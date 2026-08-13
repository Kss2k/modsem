# E-step, likelihood and derivatives for the recursive (graph) LMS backend.
#
# There is one evaluation path: a quadrature rule per observation. A rule that
# happens to be shared by every observation (adaptive "none" or "quasi") is not
# replicated into an N * Q matrix -- it is stored once, as Q nodes, and read by
# every observation through a zero stride (`shared = TRUE` below). That keeps
# non-adaptive rules affordable in memory at large N without a second
# implementation to keep in step.
#
# The sufficient-statistic M-step that shared rules used to get, which reduced
# the posterior as it was formed and cost O(Q) rather than O(N * Q), has been
# removed along with the kernel that went with it. Shared rules therefore pay
# the same per-iteration cost as per-observation ones; `adaptive = "quasi"`
# remains worthwhile for accuracy per node under a shared rule, not for speed.
#
# Model setup lives in R/lms_graph_model.R.


# Density on a common set of nodes, in the N by Q layout the quasi-adaptive
# node pruner expects.
lmsGraphDensity <- function(nodes, submodel, workspace, link) {
  result <- lmsGraphPstepWorkspaceCpp(
    lmsGraphMatrices(submodel), nodes, submodel$info$numXis,
    submodel$info$numEtas, workspace, rep(1, NROW(nodes)),
    rep(1, submodel$data$n), identical(link, "logit"),
    ThreadEnv$n.threads %||% 1L, TRUE
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

    aggregate <- lmsGraphPstepWorkspaceCpp(
      M, rule$n, submodel$info$numXis, submodel$info$numEtas, workspace,
      rule$w, sampling.weights, identical(link, "logit"),
      ThreadEnv$n.threads %||% 1L, isTRUE(rule$common)
    )

    posterior <- aggregate$posterior %||% NULL
    P.g <- if (is.null(posterior)) NULL else posterior * sampling.weights

    groups[[g]] <- list(
      P = P.g,
      posterior = posterior,
      log.kernel = aggregate$logKernel %||% NULL,
      V = rule$n, w = rule$w, quad = rule,
      obsLL = aggregate$logLik,
      sampling.weights = sampling.weights, link = link,
      workspace = workspace,
      common = isTRUE(rule$common)
    )
  }

  quad.err <- vapply(groups, function(group) {
    err <- group$quad$error.abs
    if (!length(err) || !is.finite(err[[1L]])) 0.0 else abs(err[[1L]])
  }, numeric(1L))

  # Identifies THIS E-step, for `lmsGraphCompletePass`'s memo. The posterior and
  # the quadrature rule both change here while theta need not, so a value keyed
  # on theta alone would be read back against the wrong E-step.
  LmsGraphCache$stamp <- LmsGraphCache$stamp + 1L

  list(P_GROUPS = groups, quad = lapply(groups, `[[`, "quad"),
       quad.err = sum(quad.err),
       obsLL = sum(vapply(groups, `[[`, numeric(1L), "obsLL")),
       stamp = LmsGraphCache$stamp,
       lms.backend = "graph")
}


lmsGraphObjective <- function(theta, model, P, observed = FALSE, sign = -1) {
  if (!observed)
    return(sign * lmsGraphCompletePass(theta, model, P)$objective)

  filled <- fillModel(model = model, theta = theta, method = "lms")
  ll <- 0
  for (g in seq_len(model$info$n.groups)) {
    Pg <- P$P_GROUPS[[g]]
    submodel <- filled$models[[g]]
    ll <- ll + lmsGraphAggregateCpp(
      lmsGraphLogKernel(submodel, Pg$V, Pg$link %||% "logit",
                        Pg$workspace %||% lmsGraphWorkspace(submodel),
                        shared = isTRUE(Pg$common)),
      Pg$w, Pg$sampling.weights
    )$logLik
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
  if (!observed) {
    pass <- lmsGraphCompletePass(theta, model, P)
    return(drop(sign * lmsFirstDerivativeJacobian(theta, model) %*% pass$raw))
  }

  filled <- fillModel(model = model, theta = theta, method = "lms")
  locations <- model$params$gradientStruct$locations
  raw <- stats::setNames(numeric(NROW(locations)), locations$param)
  for (g in seq_len(model$info$n.groups)) {
    loc <- locations[locations$group == g, , drop = FALSE]
    if (!NROW(loc)) next
    submodel <- filled$models[[g]]
    M <- lmsGraphMatrices(submodel)
    if (!NCOL(M$thresholds)) M$thresholdDelta <- M$thresholds
    Pg <- P$P_GROUPS[[g]]
    raw[loc$param] <- lmsGraphReverseScoreCpp(
      M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
      submodel$data$data.split, submodel$data$colidx0,
      lmsGraphOrderedIndex(submodel, M),
      Pg$w, Pg$sampling.weights, matrix(numeric(), 0L, 0L),
      TRUE, loc$block, loc$row, loc$col, loc$symmetric,
      identical(Pg$link %||% "logit", "logit"),
      ThreadEnv$n.threads %||% 1L,
      Pg$workspace %||% lmsGraphWorkspace(submodel),
      isTRUE(Pg$common)
    )$score
  }
  drop(sign * lmsFirstDerivativeJacobian(theta, model) %*% raw)
}


# One slot, holding the last complete-data pass. See lmsGraphCompletePass.
LmsGraphCache <- new.env(parent = emptyenv())
LmsGraphCache$stamp <- 0L
LmsGraphCache$last <- NULL


# The complete-data objective AND its score, from a single traversal, memoised
# on (E-step, theta).
#
# The M-step's optimiser asks for the objective and then the gradient at the very
# same theta -- measured at 100% of gradient calls -- and then frequently asks
# for that objective again. The observed call pattern over one M-step is
#     o1 g1 o2 g2 o2 o3 g3 o4 g4 o4 o5 g5 o6 o7
# i.e. 9 objective and 5 gradient calls across 7 distinct points. Evaluated
# independently that is 14 full N by Q passes where 7 suffice; the reverse
# traversal already computes every response probability the objective needs.
#
# Note the objective always comes FIRST, so a cache fed only by the gradient
# could never serve it -- the two have to be produced together on first touch.
#
# Keyed on `P$stamp` as well as theta, because the posterior and the quadrature
# rule both change at an E-step while theta may not; keyed on theta alone a
# stale entry would be read back against the wrong E-step.
lmsGraphCompletePass <- function(theta, model, P) {
  cached <- LmsGraphCache$last
  if (!is.null(cached) && identical(cached$stamp, P$stamp) &&
      length(cached$theta) == length(theta) && all(cached$theta == theta))
    return(cached)

  filled <- fillModel(model = model, theta = theta, method = "lms")
  locations <- model$params$gradientStruct$locations
  raw <- stats::setNames(numeric(NROW(locations)), locations$param)
  objective <- 0
  for (g in seq_len(model$info$n.groups)) {
    Pg <- P$P_GROUPS[[g]]
    submodel <- filled$models[[g]]
    M <- lmsGraphMatrices(submodel)
    if (!NCOL(M$thresholds)) M$thresholdDelta <- M$thresholds
    loc <- locations[locations$group == g, , drop = FALSE]
    out <- lmsGraphReverseScoreCpp(
      M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
      submodel$data$data.split, submodel$data$colidx0,
      lmsGraphOrderedIndex(submodel, M),
      Pg$w, Pg$sampling.weights, Pg$P,
      FALSE, loc$block, loc$row, loc$col, loc$symmetric,
      identical(Pg$link %||% "logit", "logit"),
      ThreadEnv$n.threads %||% 1L,
      Pg$workspace %||% lmsGraphWorkspace(submodel),
      isTRUE(Pg$common)
    )
    objective <- objective + out$objective
    if (NROW(loc)) raw[loc$param] <- out$score
  }

  pass <- list(stamp = P$stamp, theta = theta, objective = objective, raw = raw)
  LmsGraphCache$last <- pass
  pass
}


# The Cholesky-parameterised covariance blocks (A, psi, covZetaXi) do have a
# closed-form score: `cholDirectional()` in src/lms_graph.cpp gives the
# directional derivative of the Cholesky factor, and both score routines use it
# (cases 6, 7 and 17). These directions used to be overwritten with forward
# differences, one full N-by-Q pass per covariance parameter per gradient call.
# That cost nothing on the common-rule path, where an objective evaluation is
# O(Q) via sufficient statistics, but on the per-observation path every
# evaluation traverses N-by-Q -- it was ~88% of all objective work and roughly
# two thirds of total EM time.
#
# The analytic scores were checked against central differences over the cross
# product {full, quasi} x {complete, observed} x {continuous, ordered}, plus a
# two-equation model to exercise `psi`. The covariance directions agreed to
# 1e-10..2e-8 relative, better in every configuration than the worst
# non-covariance parameter (thresholds and residual variances, 5e-9..5e-7),
# which are the entries this gradient was already trusted for.
#
# `epsilon` is retained in the signatures because callers pass it positionally.


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
