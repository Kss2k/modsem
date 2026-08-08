# Experimental recursive LMS backend ---------------------------------------
#
# This backend integrates a standard-normal innovation for every latent
# variable.  A single joint covariance matrix transforms those innovations;
# latent variables are subsequently evaluated in structural (topological)
# order.  It intentionally starts with conditionally independent continuous
# indicators so that the numerical and structural parts can be validated
# independently of the legacy conditional-MVN implementation.

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

    # Identification for ordinal responses: zero intercept and unit response
    # scale. Off-diagonal response residuals remain unsupported.
    submodel$matrices$tauX[indicator.index, 1L] <- 0
    submodel$matrices$thetaDelta[indicator.index, ] <- 0
    submodel$matrices$thetaDelta[, indicator.index] <- 0
    diag(submodel$matrices$thetaDelta)[indicator.index] <- 1
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
    hessian.observed = hessianObsLogLikLmsGraph
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


lmsGraphBaseRule <- function(submodel, max.nodes = 1e6) {
  dimension <- submodel$info$numXis + submodel$info$numEtas
  order <- submodel$quad$m
  nodes.per.observation <- order ^ dimension
  mod_stopif(!is.finite(nodes.per.observation) ||
               nodes.per.observation > max.nodes,
    sprintf(paste0("The LMS graph product rule requires %.0f nodes in %d latent ",
                   "dimensions. Reduce `nodes` or use the legacy LMS backend."),
            nodes.per.observation, dimension))
  base <- quadrature(m = order, k = dimension, adaptive = FALSE,
                     quad.range = submodel$quad$quad.range)
  observations <- submodel$data$n
  index <- rep(seq_len(NROW(base$n)), times = observations)
  base$n <- base$n[index, , drop = FALSE]
  base$w <- rep(base$w, times = observations)
  base$Q <- nodes.per.observation
  base$N <- observations
  base$packed <- TRUE
  base
}


lmsGraphObservationRows <- function(submodel) {
  rows <- vector("list", submodel$data$n)
  offset <- 0L
  for (pattern in seq_along(submodel$data$data.split)) {
    values <- submodel$data$data.split[[pattern]]
    columns <- submodel$data$colidx0[[pattern]]
    for (i in seq_len(NROW(values))) {
      offset <- offset + 1L
      rows[[offset]] <- list(
        values = matrix(values[i, ], nrow = 1L),
        columns = columns
      )
    }
  }
  rows
}


lmsGraphRowLogPosterior <- function(z, matrices, row, num.xis, num.etas,
                                    ordered, logistic, ncores = 1L) {
  kernel <- lmsGraphLogKernelCpp(
    matrices, matrix(z, nrow = 1L), num.xis, num.etas,
    list(row$values), list(row$columns), ordered, logistic, ncores
  )
  kernel[1L, 1L] - .5 * sum(z * z)
}


lmsGraphAdaptiveRuleR <- function(submodel, previous = NULL,
                                  link = c("logit", "probit"),
                                  tolerance = 1e-8,
                                  mode.max.iter = 25L,
                                  curvature.floor = 1e-6,
                                  max.nodes = 1e6) {
  link <- match.arg(link)
  dimension <- submodel$info$numXis + submodel$info$numEtas
  if (!dimension) return(lmsGraphBaseRule(submodel, max.nodes))

  base <- quadrature(
    m = submodel$quad$m, k = dimension, adaptive = FALSE,
    quad.range = submodel$quad$quad.range
  )
  mod_stopif(NROW(base$n) > max.nodes,
             "The adaptive LMS graph rule exceeds the node limit.")

  M <- submodel$matrices
  if (is.null(M$thresholds))
    M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
  ordered <- match(submodel$info$ordered %||% character(),
                   rownames(M$lambdaX)) - 1L
  ordered <- ordered[!is.na(ordered)]
  observations <- lmsGraphObservationRows(submodel)
  N <- length(observations)
  Q <- NROW(base$n)
  starts <- if (!is.null(previous$modes) &&
                identical(dim(previous$modes), c(N, dimension)))
    previous$modes else matrix(0, N, dimension)

  modes <- matrix(0, N, dimension)
  packed.nodes <- matrix(0, N * Q, dimension)
  packed.log.weights <- numeric(N * Q)
  convergence <- integer(N)
  curvature.adjusted <- logical(N)
  log.base.prior <- -.5 * rowSums(base$n^2) -
    .5 * dimension * log(2 * pi)

  for (i in seq_len(N)) {
    objective <- function(z) -lmsGraphRowLogPosterior(
      z, M, observations[[i]], submodel$info$numXis,
      submodel$info$numEtas, ordered, identical(link, "logit")
    )
    fit <- tryCatch(stats::nlminb(
      start = starts[i, ], objective = objective,
      control = list(iter.max = mode.max.iter,
                     eval.max = max(100L, 5L * mode.max.iter),
                     rel.tol = tolerance, x.tol = sqrt(tolerance))
    ), error = identity)
    if (inherits(fit, "error") || !all(is.finite(fit$par))) {
      mode <- starts[i, ]
      convergence[i] <- 99L
    } else {
      mode <- fit$par
      convergence[i] <- fit$convergence
    }
    modes[i, ] <- mode

    H <- tryCatch(stats::optimHess(mode, objective), error = function(e)
      diag(dimension))
    H <- .5 * (H + t(H))
    decomposition <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
    if (is.null(decomposition) || any(!is.finite(decomposition$values))) {
      decomposition <- list(values = rep(1, dimension),
                            vectors = diag(dimension))
      curvature.adjusted[i] <- TRUE
    }
    adjusted <- pmax(decomposition$values, curvature.floor)
    curvature.adjusted[i] <- curvature.adjusted[i] ||
      any(adjusted != decomposition$values)
    transform <- decomposition$vectors %*% diag(1 / sqrt(adjusted), dimension)
    adaptive.nodes <- sweep(base$n %*% t(transform), 2L, mode, "+")
    use <- (i - 1L) * Q + seq_len(Q)
    packed.nodes[use, ] <- adaptive.nodes
    log.target.prior <- -.5 * rowSums(adaptive.nodes^2) -
      .5 * dimension * log(2 * pi)
    log.jacobian <- -.5 * sum(log(adjusted))
    packed.log.weights[use] <- log(base$w) + log.target.prior +
      log.jacobian - log.base.prior
  }

  list(
    n = packed.nodes,
    w = exp(packed.log.weights),
    log.w = packed.log.weights,
    modes = modes,
    convergence = convergence,
    curvature.adjusted = curvature.adjusted,
    Q = Q, N = N, k = dimension, m = submodel$quad$m,
    quad.range = submodel$quad$quad.range,
    adaptive = TRUE,
    adaptive.frequency = submodel$quad$adaptive.frequency,
    packed = TRUE
  )
}


lmsGraphAdaptiveRule <- function(submodel, previous = NULL,
                                 link = c("logit", "probit"),
                                 tolerance = 1e-8,
                                 mode.max.iter = 25L,
                                 curvature.floor = 1e-6,
                                 derivative.step = 1e-4,
                                 max.nodes = 1e6,
                                 workspace = NULL) {
  link <- match.arg(link)
  dimension <- submodel$info$numXis + submodel$info$numEtas
  if (!dimension) return(lmsGraphBaseRule(submodel, max.nodes))
  base <- quadrature(
    m = submodel$quad$m, k = dimension, adaptive = FALSE,
    quad.range = submodel$quad$quad.range
  )
  mod_stopif(NROW(base$n) > max.nodes,
             "The adaptive LMS graph rule exceeds the node limit.")
  M <- submodel$matrices
  if (is.null(M$thresholds))
    M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
  ordered <- match(submodel$info$ordered %||% character(),
                   rownames(M$lambdaX)) - 1L
  ordered <- ordered[!is.na(ordered)]
  N <- submodel$data$n
  Q <- NROW(base$n)
  starts <- if (!is.null(previous$modes) &&
                identical(dim(previous$modes), c(N, dimension)))
    previous$modes else matrix(0, N, dimension)
  compiled <- lmsGraphAdaptiveRuleCpp(
    M, base$n, base$w, starts,
    submodel$info$numXis, submodel$info$numEtas,
    submodel$data$data.split, submodel$data$colidx0, ordered,
    identical(link, "logit"), mode.max.iter, tolerance,
    curvature.floor, derivative.step, ThreadEnv$n.threads %||% 1L,
    workspace
  )
  failed <- sum(compiled$convergence != 0L)
  adjusted <- sum(compiled$curvatureAdjusted != 0L)
  mod_warnif(failed > 0L,
    sprintf("AGHQ mode optimization did not fully converge for %d of %d observations.",
            failed, N))
  mod_warnif(adjusted > 0L,
    sprintf("AGHQ posterior curvature was regularized for %d of %d observations.",
            adjusted, N))
  list(
    n = compiled$nodes,
    w = exp(compiled$logWeights),
    log.w = compiled$logWeights,
    modes = compiled$modes,
    convergence = as.integer(compiled$convergence),
    curvature.adjusted = as.logical(compiled$curvatureAdjusted),
    Q = Q, N = N, k = dimension, m = submodel$quad$m,
    quad.range = submodel$quad$quad.range,
    adaptive = TRUE,
    adaptive.frequency = submodel$quad$adaptive.frequency,
    packed = TRUE
  )
}


lmsGraphRule <- function(submodel, previous = NULL,
                         link = c("logit", "probit"),
                         tolerance = 1e-8, max.nodes = 1e6,
                         workspace = NULL) {
  link <- match.arg(link)
  if (isTRUE(submodel$quad$adaptive))
    lmsGraphAdaptiveRule(submodel, previous, link, tolerance,
                         max.nodes = max.nodes, workspace = workspace)
  else
    lmsGraphBaseRule(submodel, max.nodes)
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
  M <- submodel$matrices
  ordered <- match(submodel$info$ordered %||% character(),
                   rownames(M$lambdaX)) - 1L
  ordered <- ordered[!is.na(ordered)]
  lmsGraphWorkspaceCpp(submodel$data$data.split, submodel$data$colidx0,
                       ordered)
}


lmsGraphLogKernel <- function(submodel, nodes, link = c("logit", "probit"),
                              workspace = NULL) {
  link <- match.arg(link)
  lmsGraphValidate(submodel)
  M <- submodel$matrices
  data <- submodel$data
  if (is.null(M$thresholds))
    M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
  ordered <- match(submodel$info$ordered %||% character(), rownames(M$lambdaX)) - 1L
  ordered <- ordered[!is.na(ordered)]
  if (is.null(workspace)) workspace <- lmsGraphWorkspace(submodel)
  lmsGraphLogKernelWorkspaceCpp(
    M, nodes, submodel$info$numXis, submodel$info$numEtas, workspace,
    identical(link, "logit"), ThreadEnv$n.threads %||% 1L)
}


lmsGraphLogSumExp <- function(x) {
  maxima <- apply(x, 1L, max)
  maxima + log(rowSums(exp(x - maxima)))
}


pstepLmsGraph <- function(model, theta, lastQuad = NULL, recalcQuad = FALSE,
                          link = c("logit", "probit"),
                          adaptive.quad.tol = 1e-8, ...) {
  link <- match.arg(link)
  filled <- fillModel(model = model, theta = theta, method = "lms")
  groups <- vector("list", model$info$n.groups)

  for (g in seq_len(model$info$n.groups)) {
    submodel <- filled$models[[g]]
    submodel$quad <- model$models[[g]]$quad
    previous <- if (!is.null(lastQuad)) lastQuad[[g]] else NULL
    workspace <- previous$workspace %||% lmsGraphWorkspace(submodel)
    rule <- if (!is.null(previous) && !recalcQuad) previous else
      lmsGraphRule(submodel, previous = previous, link = link,
                   tolerance = adaptive.quad.tol, workspace = workspace)
    rule$workspace <- workspace
    sampling.weights <- submodel$data$weights
    if (is.null(sampling.weights)) sampling.weights <- rep(1, submodel$data$n)
    M <- submodel$matrices
    if (is.null(M$thresholds))
      M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
    aggregate <- lmsGraphPstepWorkspaceCpp(
      M, rule$n, submodel$info$numXis, submodel$info$numEtas, workspace,
      rule$w, sampling.weights, identical(link, "logit"),
      ThreadEnv$n.threads %||% 1L
    )
    log.kernel <- aggregate$logKernel
    posterior <- aggregate$posterior

    groups[[g]] <- list(
      P = posterior * sampling.weights,
      posterior = posterior,
      log.kernel = log.kernel,
      V = rule$n, w = rule$w, quad = rule,
      obsLL = aggregate$logLik,
      sampling.weights = sampling.weights, link = link,
      workspace = workspace
    )
  }

  list(P_GROUPS = groups, quad = lapply(groups, `[[`, "quad"),
       obsLL = sum(vapply(groups, `[[`, numeric(1L), "obsLL")),
       lms.backend = "graph")
}


lmsGraphObjective <- function(theta, model, P, observed = FALSE, sign = -1) {
  filled <- fillModel(model = model, theta = theta, method = "lms")
  ll <- 0
  for (g in seq_len(model$info$n.groups)) {
    Pg <- P$P_GROUPS[[g]]
    if (observed) {
      log.kernel <- lmsGraphLogKernel(
        filled$models[[g]], Pg$V, Pg$link %||% "logit", Pg$workspace
      )
      ll <- ll + lmsGraphAggregateCpp(log.kernel, Pg$w, Pg$sampling.weights)$logLik
    } else {
      submodel <- filled$models[[g]]
      M <- submodel$matrices
      if (is.null(M$thresholds))
        M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
      workspace <- Pg$workspace %||% lmsGraphWorkspace(submodel)
      ll <- ll + lmsGraphCompleteWorkspaceCpp(
        M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
        workspace, Pg$P, identical(Pg$link %||% "logit", "logit"),
        ThreadEnv$n.threads %||% 1L
      )
    }
  }
  sign * ll
}


compLogLikLmsGraph <- function(theta, model, P, sign = -1, ...)
  tryCatch(lmsGraphObjective(theta, model, P, FALSE, sign), error = function(e) NA_real_)

obsLogLikLmsGraph <- function(theta, model, P, sign = -1, ...)
  tryCatch(lmsGraphObjective(theta, model, P, TRUE, sign), error = function(e) NA_real_)


lmsGraphGradient <- function(theta, objective, epsilon = 1e-6) {
  step <- epsilon * pmax(1, abs(theta))
  out <- numeric(length(theta))
  for (j in seq_along(theta)) {
    plus <- minus <- theta
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    out[j] <- (objective(plus) - objective(minus)) / (2 * step[j])
  }
  names(out) <- names(theta)
  out
}


lmsGraphHessian <- function(theta, objective,
                            epsilon = .Machine$double.eps ^ (1 / 4)) {
  step <- epsilon * pmax(1, abs(theta))
  p <- length(theta)
  H <- matrix(0, p, p, dimnames = list(names(theta), names(theta)))
  f0 <- objective(theta)
  for (j in seq_len(p)) {
    plus <- minus <- theta
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    H[j, j] <- (objective(plus) - 2 * f0 + objective(minus)) / step[j]^2
    if (j < p) for (k in (j + 1L):p) {
      pp <- pm <- mp <- mm <- theta
      pp[c(j, k)] <- pp[c(j, k)] + step[c(j, k)]
      pm[j] <- pm[j] + step[j]; pm[k] <- pm[k] - step[k]
      mp[j] <- mp[j] - step[j]; mp[k] <- mp[k] + step[k]
      mm[c(j, k)] <- mm[c(j, k)] - step[c(j, k)]
      H[j, k] <- H[k, j] <- (objective(pp) - objective(pm) -
                              objective(mp) + objective(mm)) /
                             (4 * step[j] * step[k])
    }
  }
  H
}


lmsGraphHessianFromGradient <- function(theta, gradient,
                                        epsilon = .Machine$double.eps ^ (1 / 4)) {
  step <- epsilon * pmax(1, abs(theta))
  H <- matrix(0, length(theta), length(theta),
              dimnames = list(names(theta), names(theta)))
  for (j in seq_along(theta)) {
    plus <- minus <- theta
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    H[, j] <- (gradient(plus) - gradient(minus)) / (2 * step[j])
  }
  0.5 * (H + t(H))
}


lmsGraphAnalyticalGradient <- function(theta, model, P, observed, sign) {
  filled <- fillModel(model = model, theta = theta, method = "lms")
  locations <- model$params$gradientStruct$locations
  raw <- stats::setNames(numeric(NROW(locations)), locations$param)
  for (g in seq_len(model$info$n.groups)) {
    loc <- locations[locations$group == g, , drop = FALSE]
    if (!NROW(loc)) next
    submodel <- filled$models[[g]]
    M <- submodel$matrices
    if (is.null(M$thresholds)) {
      M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
      M$thresholdDelta <- M$thresholds
    }
    ordered <- match(submodel$info$ordered %||% character(),
                     rownames(M$lambdaX)) - 1L
    ordered <- ordered[!is.na(ordered)]
    Pg <- P$P_GROUPS[[g]]
    workspace <- Pg$workspace %||% lmsGraphWorkspace(submodel)
    score <- lmsGraphReverseScoreCpp(
      M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
      submodel$data$data.split, submodel$data$colidx0, ordered,
      Pg$w, Pg$sampling.weights,
      if (observed) matrix(numeric(), 0L, 0L) else Pg$P,
      observed, loc$block, loc$row, loc$col, loc$symmetric,
      identical(Pg$link %||% "logit", "logit"),
      ThreadEnv$n.threads %||% 1L, workspace
    )
    raw[loc$param] <- score
  }
  J <- lmsFirstDerivativeJacobian(theta, model)
  drop(sign * J %*% raw)
}


gradientCompLogLikLmsGraph <- function(theta, model, P, sign = -1,
                                       epsilon = 1e-6, ...) {
  lmsGraphAnalyticalGradient(theta, model, P, observed = FALSE, sign = sign)
}

gradientObsLogLikLmsGraph <- function(theta, model, P, sign = -1,
                                      epsilon = 1e-6, ...) {
  lmsGraphAnalyticalGradient(theta, model, P, observed = TRUE, sign = sign)
}

hessianCompLogLikLmsGraph <- function(theta, model, P, sign = -1,
                                      .relStep = .Machine$double.eps ^ (1/4), ...) {
  lmsGraphHessianFromGradient(theta, function(x)
    gradientCompLogLikLmsGraph(x, model, P, sign), .relStep)
}

hessianObsLogLikLmsGraph <- function(theta, model, P, sign = -1,
                                     .relStep = .Machine$double.eps ^ (1/4), ...) {
  lmsGraphHessianFromGradient(theta, function(x)
    gradientObsLogLikLmsGraph(x, model, P, sign), .relStep)
}
