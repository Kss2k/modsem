# E-step, likelihood and derivatives for the recursive (graph) LMS backend.
#
# Two evaluation paths exist, and the distinction that matters is whether the
# quadrature rule is shared by every observation:
#
#   common  one rule for the whole group (adaptive "none" or "quasi"). The
#           posterior is reduced to sufficient statistics as it is formed, so
#           the M-step never materialises an N by Q matrix and costs O(Q).
#   packed  one rule per observation (adaptive "full"). Nodes are stacked into
#           an N * Q matrix and the M-step is O(N * Q).
#
# Model setup lives in R/lms_graph_model.R.


# Density on a common set of nodes, in the N by Q layout the quasi-adaptive
# node pruner expects.
lmsGraphDensity <- function(nodes, submodel, workspace, link) {
  result <- lmsGraphCommonPstepCpp(
    lmsGraphMatrices(submodel), nodes, submodel$info$numXis,
    submodel$info$numEtas, workspace, rep(1, NROW(nodes)),
    rep(1, submodel$data$n), identical(link, "logit"),
    ThreadEnv$n.threads %||% 1L
  )
  exp(result$logKernel)
}


# Transform a common base rule into one rule per observation, centred on each
# observation's posterior mode and scaled by its Hessian. This is the only
# thing `adaptive = "full"` adds, and it is deliberately independent of which
# integration rule produced `base`.
lmsGraphAdaptPerObservation <- function(submodel, base, link, workspace,
                                        previous = NULL,
                                        tolerance = 1e-8,
                                        mode.max.iter = 25L,
                                        curvature.floor = 1e-6,
                                        derivative.step = 1e-4) {
  dimension <- lmsGraphDimension(submodel)
  N <- submodel$data$n
  Q <- NROW(base$n)

  # The base rule's size was already checked against the node budget. What is
  # left is the N * Q packed rule, whose cost is linear in the sample rather
  # than exponential in the dimension, so it warrants a warning rather than the
  # hard error that guards product-rule explosion.
  mod_warnif(as.numeric(N) * Q > 5e6, sprintf(paste0(
    "`adaptive = \"full\"` builds %.0f nodes (%d observations x %d points). ",
    "Expect slow iterations; `adaptive = \"quasi\"` shares one rule instead."),
    as.numeric(N) * Q, N, Q))

  starts <- if (!is.null(previous$modes) &&
                identical(dim(previous$modes), c(N, dimension)))
    previous$modes else matrix(0, N, dimension)

  M <- lmsGraphMatrices(submodel)
  compiled <- lmsGraphAdaptiveRuleCpp(
    M, base$n, base$w, starts, submodel$info$numXis, submodel$info$numEtas,
    submodel$data$data.split, submodel$data$colidx0,
    lmsGraphOrderedIndex(submodel, M), identical(link, "logit"),
    mode.max.iter, tolerance, curvature.floor, derivative.step,
    ThreadEnv$n.threads %||% 1L, workspace
  )

  failed <- sum(compiled$convergence != 0L)
  adjusted <- sum(compiled$curvatureAdjusted != 0L)
  mod_warnif(failed > 0L, sprintf(
    "Adaptive mode optimization did not fully converge for %d of %d observations.",
    failed, N))
  mod_warnif(adjusted > 0L, sprintf(
    "Adaptive posterior curvature was regularized for %d of %d observations.",
    adjusted, N))

  base$n <- compiled$nodes
  base$w <- exp(compiled$logWeights)
  base$log.w <- compiled$logWeights
  base$modes <- compiled$modes
  base$convergence <- as.integer(compiled$convergence)
  base$curvature.adjusted <- as.logical(compiled$curvatureAdjusted)
  base$Q <- Q
  base$N <- N
  base$common <- FALSE
  base$packed <- TRUE
  base
}


lmsGraphBuildRule <- function(submodel, workspace, link, previous = NULL) {
  spec <- submodel$quad
  spec$k <- lmsGraphDimension(submodel)

  rule <- buildQuadRule(
    spec, density = lmsGraphDensity, previous = previous,
    submodel = submodel, workspace = workspace, link = link
  )
  rule$N <- submodel$data$n

  # The pruner's densities are only needed while selecting nodes. The E-step
  # below recomputes them, so keeping them would duplicate N * Q storage in the
  # rule that is cached across every later iteration.
  rule$F <- NULL

  if (!identical(spec$adaptive, "full") || spec$k < 1L) return(rule)
  lmsGraphAdaptPerObservation(submodel, rule, link, workspace, previous)
}


pstepLmsGraph <- function(model, theta, lastQuad = NULL, recalcQuad = FALSE,
                          link = c("logit", "probit"), ...) {
  link <- match.arg(link)
  filled <- fillModel(model = model, theta = theta, method = "lms")
  groups <- vector("list", model$info$n.groups)

  for (g in seq_len(model$info$n.groups)) {
    submodel <- filled$models[[g]]
    submodel$quad <- model$models[[g]]$quad
    previous <- if (!is.null(lastQuad)) lastQuad[[g]] else NULL

    workspace <- previous$workspace %||% lmsGraphWorkspace(submodel)
    rule <- if (!is.null(previous) && !recalcQuad) previous else
      lmsGraphBuildRule(submodel, workspace, link, previous)
    rule$workspace <- workspace

    M <- lmsGraphMatrices(submodel)
    sampling.weights <- lmsGraphSamplingWeights(submodel)

    aggregate <- if (isTRUE(rule$common))
      lmsGraphCommonEstepCpp(
        M, rule$n, submodel$info$numXis, submodel$info$numEtas, workspace,
        rule$w, sampling.weights, identical(link, "logit"),
        256L, ThreadEnv$n.threads %||% 1L
      ) else lmsGraphPstepWorkspaceCpp(
        M, rule$n, submodel$info$numXis, submodel$info$numEtas, workspace,
        rule$w, sampling.weights, identical(link, "logit"),
        ThreadEnv$n.threads %||% 1L
      )

    posterior <- aggregate$posterior %||% NULL
    groups[[g]] <- list(
      P = if (is.null(posterior)) NULL else posterior * sampling.weights,
      posterior = posterior,
      log.kernel = aggregate$logKernel %||% NULL,
      V = rule$n, w = rule$w, quad = rule,
      obsLL = aggregate$logLik,
      sampling.weights = sampling.weights, link = link,
      workspace = workspace,
      sufficient = if (isTRUE(rule$common)) aggregate$statistics else NULL,
      common = isTRUE(rule$common)
    )
  }

  quad.err <- vapply(groups, function(group) {
    err <- group$quad$error.abs
    if (!length(err) || !is.finite(err[[1L]])) 0.0 else abs(err[[1L]])
  }, numeric(1L))

  list(P_GROUPS = groups, quad = lapply(groups, `[[`, "quad"),
       quad.err = sum(quad.err),
       obsLL = sum(vapply(groups, `[[`, numeric(1L), "obsLL")),
       lms.backend = "graph")
}


lmsGraphObjective <- function(theta, model, P, observed = FALSE, sign = -1) {
  filled <- fillModel(model = model, theta = theta, method = "lms")
  ll <- 0
  for (g in seq_len(model$info$n.groups)) {
    Pg <- P$P_GROUPS[[g]]
    submodel <- filled$models[[g]]
    M <- lmsGraphMatrices(submodel)
    logistic <- identical(Pg$link %||% "logit", "logit")
    workspace <- Pg$workspace %||% lmsGraphWorkspace(submodel)

    ll <- ll + if (observed && Pg$common) {
      lmsGraphCommonEstepCpp(
        M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
        workspace, Pg$w, Pg$sampling.weights, logistic,
        256L, ThreadEnv$n.threads %||% 1L
      )$logLik

    } else if (observed) {
      lmsGraphAggregateCpp(
        lmsGraphLogKernel(submodel, Pg$V, Pg$link %||% "logit", workspace),
        Pg$w, Pg$sampling.weights
      )$logLik

    } else if (Pg$common) {
      lmsGraphSufficientCompleteCpp(
        M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
        workspace, Pg$sufficient, logistic, ThreadEnv$n.threads %||% 1L
      )

    } else {
      lmsGraphCompleteWorkspaceCpp(
        M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
        workspace, Pg$P, logistic, ThreadEnv$n.threads %||% 1L
      )
    }
  }
  sign * ll
}


compLogLikLmsGraph <- function(theta, model, P, sign = -1, ...)
  tryCatch(lmsGraphObjective(theta, model, P, FALSE, sign),
           error = function(e) NA_real_)

obsLogLikLmsGraph <- function(theta, model, P, sign = -1, ...)
  tryCatch(lmsGraphObjective(theta, model, P, TRUE, sign),
           error = function(e) NA_real_)


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
    M <- lmsGraphMatrices(submodel)
    if (!NCOL(M$thresholds)) M$thresholdDelta <- M$thresholds
    ordered <- lmsGraphOrderedIndex(submodel, M)
    Pg <- P$P_GROUPS[[g]]
    logistic <- identical(Pg$link %||% "logit", "logit")
    workspace <- Pg$workspace %||% lmsGraphWorkspace(submodel)

    score <- if (isTRUE(Pg$common)) {
      statistics <- Pg$sufficient
      if (observed) statistics <- lmsGraphCommonEstepCpp(
        M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
        workspace, Pg$w, Pg$sampling.weights, logistic,
        256L, ThreadEnv$n.threads %||% 1L
      )$statistics
      lmsGraphSufficientScoreCpp(
        M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
        statistics, ordered, loc$block, loc$row, loc$col, loc$symmetric,
        logistic
      )
    } else lmsGraphReverseScoreCpp(
      M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
      submodel$data$data.split, submodel$data$colidx0, ordered,
      Pg$w, Pg$sampling.weights,
      if (observed) matrix(numeric(), 0L, 0L) else Pg$P,
      observed, loc$block, loc$row, loc$col, loc$symmetric,
      logistic, ThreadEnv$n.threads %||% 1L, workspace
    )
    raw[loc$param] <- score
  }
  J <- lmsFirstDerivativeJacobian(theta, model)
  drop(sign * J %*% raw)
}


# The Cholesky-parameterised covariance blocks have no closed-form score in the
# recursive formulation, so those directions fall back to finite differences.
lmsGraphReplaceCovarianceGradient <- function(theta, gradient, objective,
                                              epsilon = 1e-6) {
  covariance <- grepl("^(A|psi|covZetaXi)[0-9]+$", names(theta))
  if (!any(covariance)) return(gradient)
  baseline <- objective(theta)
  step <- epsilon * pmax(1, abs(theta))
  for (j in which(covariance)) {
    plus <- theta
    plus[j] <- plus[j] + step[j]
    gradient[j] <- (objective(plus) - baseline) / step[j]
  }
  gradient
}


gradientCompLogLikLmsGraph <- function(theta, model, P, sign = -1,
                                       epsilon = 1e-6, ...) {
  gradient <- lmsGraphAnalyticalGradient(
    theta, model, P, observed = FALSE, sign = sign
  )
  lmsGraphReplaceCovarianceGradient(
    theta, gradient,
    function(value) compLogLikLmsGraph(value, model, P, sign), epsilon
  )
}


gradientObsLogLikLmsGraph <- function(theta, model, P, sign = -1,
                                      epsilon = 1e-6, ...) {
  gradient <- lmsGraphAnalyticalGradient(
    theta, model, P, observed = TRUE, sign = sign
  )
  lmsGraphReplaceCovarianceGradient(
    theta, gradient,
    function(value) obsLogLikLmsGraph(value, model, P, sign), epsilon
  )
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
