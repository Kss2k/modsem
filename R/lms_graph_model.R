# Model setup for the recursive (graph) LMS backend.
#
# This backend integrates a standard-normal innovation for every latent
# variable. A single joint covariance matrix transforms those innovations, and
# latent variables are then evaluated in structural (topological) order, which
# is what lets it handle n-way interactions and ordered indicators. Indicators
# are conditionally independent given the latent variables.
#
# The E-step and the likelihood live in R/lms_graph_estep.R.

lmsGraphOrderedColumns <- function(data, model.syntax, ordered = NULL) {
  observed <- getOVs(modsemify(model.syntax))
  automatic <- names(data)[vapply(data, is.ordered, logical(1L))]
  columns <- unique(c(ordered, automatic))
  missing <- setdiff(columns, names(data))
  mod_stopif(length(missing), paste("Ordered variables not found:",
                                    paste(missing, collapse = ", ")))
  intersect(columns, observed)
}


lmsGraphThresholdSpec <- function(data, ordered, group = NULL,
                                  link = c("logit", "probit")) {
  link <- match.arg(link)
  groups <- if (is.null(group)) rep(1L, NROW(data)) else
    match(data[[group]], unique(data[[group]]))
  specs <- list()
  delta <- numeric()
  for (g in sort(unique(groups))) for (variable in ordered) {
    values <- data[[variable]][groups == g]
    code <- if (is.factor(values)) as.integer(values) else
      match(values, sort(unique(values[!is.na(values)])))
    categories <- max(code, na.rm = TRUE)
    mod_stopif(!is.finite(categories) || categories < 2L,
      sprintf("Ordered indicator `%s` has fewer than two categories in group %d.",
              variable, g))
    probability <- vapply(seq_len(categories - 1L), function(k)
      mean(code <= k, na.rm = TRUE), numeric(1L))
    probability <- pmin(pmax(probability, 1e-5), 1 - 1e-5)
    thresholds <- if (link == "logit") stats::qlogis(probability) else
      stats::qnorm(probability)
    gaps <- pmax(diff(thresholds), 1e-8)
    values.delta <- c(thresholds[1L], if (length(gaps)) log(expm1(gaps)))
    labels <- paste0(variable, "|t", seq_len(categories - 1L))
    index <- length(delta) + seq_along(values.delta)
    delta <- c(delta, values.delta)
    specs[[paste(g, variable, sep = "::")]] <- list(
      group = g, variable = variable, K = categories,
      index = index, labels = labels
    )
  }
  list(specs = specs, delta = delta)
}


lmsGraphPrepareOrdered <- function(model, data, model.syntax, ordered = NULL,
                                   group = NULL,
                                   link = c("logit", "probit")) {
  link <- match.arg(link)
  ordered <- lmsGraphOrderedColumns(data, model.syntax, ordered)
  if (!length(ordered)) return(model)
  threshold.info <- lmsGraphThresholdSpec(data, ordered, group, link)
  old.theta <- model$theta
  threshold.starts <- numeric()

  for (g in seq_along(model$models)) {
    submodel <- model$models[[g]]
    indicators <- rownames(submodel$matrices$lambdaX)
    ordered.g <- intersect(ordered, indicators)
    indicator.index <- match(ordered.g, indicators)

    # Identification for ordinal responses: zero intercept, and a response
    # scale fixed by the link rather than estimated. Off-diagonal response
    # residuals remain unsupported.
    #
    # The likelihood never reads this entry -- every `theta(column, column)`
    # in the kernels and scores sits in the continuous branch -- so it exists
    # for whatever reconstructs implied variances from the matrices, chiefly
    # the standardized solution. It therefore has to carry the variance the
    # link actually implies: a standard logistic residual has variance pi^2/3,
    # not 1. Recording 1 under `link = "logit"` left standardized loadings and
    # R-squared for ordered indicators computed against a residual roughly
    # 3.3x too small.
    response.variance <- if (identical(link, "logit")) pi^2 / 3 else 1
    submodel$matrices$tauX[indicator.index, 1L] <- 0
    submodel$matrices$thetaDelta[indicator.index, ] <- 0
    submodel$matrices$thetaDelta[, indicator.index] <- 0
    diag(submodel$matrices$thetaDelta)[indicator.index] <- response.variance
    submodel$labelMatrices$tauX[indicator.index, 1L] <- ""
    submodel$labelMatrices$thetaDelta[indicator.index, ] <- ""
    submodel$labelMatrices$thetaDelta[, indicator.index] <- ""

    specs <- Filter(function(x) identical(x$group, g), threshold.info$specs)
    maximum <- max(c(0L, vapply(specs, function(x) x$K - 1L, integer(1L))))
    threshold.delta <- matrix(NaN, length(indicators), maximum,
      dimnames = list(indicators, paste0("t", seq_len(maximum))))
    for (spec in specs)
      threshold.delta[spec$variable, seq_len(spec$K - 1L)] <-
        threshold.info$delta[spec$index]
    submodel$matrices$thresholdDelta <- threshold.delta
    submodel$matrices$thresholds <- thresholdDeltaToThresholdMatrix(threshold.delta)
    submodel$labelMatrices$thresholdDelta <- matrix("", NROW(threshold.delta),
      NCOL(threshold.delta), dimnames = dimnames(threshold.delta))
    submodel$labelMatrices$thresholds <- submodel$labelMatrices$thresholdDelta

    free <- is.finite(threshold.delta)
    starts <- threshold.delta[free]
    names(starts) <- getParamNamesMatrix(threshold.delta, "thresholdDelta")[free]
    if (g > 1L) names(starts) <- sprintf("%s.g%d", names(starts), g)
    threshold.starts <- c(threshold.starts, starts)
    submodel$matrices$thresholdDelta[free] <- NA_real_
    submodel$info$ordered <- ordered.g
    model$models[[g]] <- submodel
  }

  params <- createTheta(model, parTable.in = model$parTable)
  common <- intersect(names(old.theta), names(params$theta))
  params$theta[common] <- old.theta[common]
  params$theta[names(threshold.starts)] <- threshold.starts
  model$params[names(params)] <- params
  model$theta <- params$theta
  model$params$bounds <- getParamBounds(model)
  model$params$gradientStruct <- getGradientStruct(model, model$theta, method = "lms")
  model$info$ordered <- ordered
  model$info$lms.graph.link <- link
  model
}


lmsGraphBackend <- function(link = c("logit", "probit")) {
  link <- match.arg(link)
  list(
    name = "graph",
    link = link,
    pstep = function(...) pstepLmsGraph(..., link = link),
    complete = compLogLikLmsGraph,
    observed = obsLogLikLmsGraph,
    gradient.complete = gradientCompLogLikLmsGraph,
    gradient.observed = gradientObsLogLikLmsGraph,
    hessian.complete = hessianCompLogLikLmsGraph,
    hessian.observed = hessianObsLogLikLmsGraph,
    # The joint M-step, NOT the ECM split (#27). The split assumed the
    # structural block was separable from the measurement kernel once the
    # E-step's latent nodes were held fixed. It is not: the nodes are
    # (xi, zeta), so eta is reconstructed from them through alpha/Gamma/Omega
    # and every structural parameter moves the indicator means. Perturbing
    # structural parameters alone changed the measurement objective in 39 of 39
    # probes, and the joint objective's whole response to such a move lived in
    # the measurement half. The structural block was optimising a surrogate
    # that missed that by a median of 205%, and lowered the true Q in 6 of 39
    # M-steps here (worst -16.6) -- see test_em_monotonicity.R.
    #
    # `mstepLmsGraphEcm` is kept, unwired, because the split becomes exact
    # under (xi, eta) nodes, where p(y|xi,eta), p(eta|xi) and p(xi) separate
    # cleanly and the latter two have closed forms. That reparameterisation is
    # the open work; until then the joint step is the correct one.
    mstep = mstepLms
  )
}


lmsGraphValidate <- function(submodel, tolerance = 1e-12) {
  M <- submodel$matrices
  mod_stopif(isTRUE(submodel$info$hasComposites),
             "The experimental LMS graph backend does not yet support composites.")
  residual <- M$thetaDelta
  if (length(residual)) {
    offdiag <- residual
    diag(offdiag) <- 0
    mod_stopif(any(abs(offdiag) > tolerance, na.rm = TRUE),
      paste0("The experimental LMS graph backend requires conditionally independent ",
             "indicators; residual covariances between indicators are not supported."))
    mod_stopif(any(!is.finite(diag(residual))) || any(diag(residual) <= 0),
               "All indicator residual variances must be finite and positive.")
  }

  invisible(TRUE)
}


lmsGraphLatentCovariance <- function(M) {
  exogenous <- if (length(M$A)) M$A %*% t(M$A) else matrix(numeric(), 0L, 0L)
  cross <- M$covZetaXi
  endogenous <- M$psi

  if (!NROW(exogenous)) return(endogenous)
  if (!NROW(endogenous)) return(exogenous)

  rbind(cbind(exogenous, t(cross)),
        cbind(cross, endogenous))
}


# The graph backend integrates over every latent variable, unlike the legacy
# backend which only integrates the nonlinear dimensions.
lmsGraphDimension <- function(submodel)
  submodel$info$numXis + submodel$info$numEtas


# `thresholds` is absent whenever a model has no ordered indicators, but the
# compiled code always expects the matrix to exist.
lmsGraphMatrices <- function(submodel) {
  M <- submodel$matrices
  if (is.null(M$thresholds))
    M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
  M
}


lmsGraphOrderedIndex <- function(submodel, M = lmsGraphMatrices(submodel)) {
  index <- match(submodel$info$ordered %||% character(), rownames(M$lambdaX)) - 1L
  index[!is.na(index)]
}


lmsGraphSamplingWeights <- function(submodel) {
  weights <- submodel$data$weights
  if (is.null(weights)) rep(1, submodel$data$n) else weights
}


lmsGraphStates <- function(submodel, nodes) {
  lmsGraphValidate(submodel)
  states <- lmsGraphStatesCpp(
    submodel$matrices, nodes, submodel$info$numXis, submodel$info$numEtas
  )
  colnames(states) <- c(submodel$info$xis, submodel$info$etas)
  states
}


lmsGraphWorkspace <- function(submodel) {
  lmsGraphWorkspaceCpp(submodel$data$data.split, submodel$data$colidx0,
                       lmsGraphOrderedIndex(submodel))
}


lmsGraphLogKernel <- function(submodel, nodes, link = c("logit", "probit"),
                              workspace = NULL, shared = FALSE) {
  link <- match.arg(link)
  lmsGraphValidate(submodel)
  if (is.null(workspace)) workspace <- lmsGraphWorkspace(submodel)
  lmsGraphLogKernelWorkspaceCpp(
    lmsGraphMatrices(submodel), nodes, submodel$info$numXis,
    submodel$info$numEtas, workspace, identical(link, "logit"),
    ThreadEnv$n.threads %||% 1L, shared)
}
