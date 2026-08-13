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

    # The frozen latent states (#28). The complete data is z = (xi, eta), NOT
    # the standardised innovations the rule holds, so z is derived once here
    # with the E-step parameters and then held fixed for the whole M-step. That
    # is what makes the complete-data objective separate:
    #   Q(theta) = sum P log f(x | z; measurement) + sum P log p(z; structural)
    # With z instead re-derived from the current parameters at every evaluation
    # -- which is what `lmsGraphStatesCpp` does, and what the objective used to
    # do -- the second term is absent, the structural parameters act through the
    # first, and the split silently optimises the wrong function (#27).
    #
    # This costs a second N by Q by d matrix alongside `rule$n`. That is the
    # same order as the nodes already stored, and only on the per-observation
    # adaptive path; a shared rule keeps both at Q by d.
    # One traversal yields both: the states to freeze, and the product values
    # the structural statistics need. Deriving them separately ran
    # `evaluateGraphStates` twice over the same N by Q nodes with the same
    # matrices, which measured 10% of a continuous run.
    derived <- if (is.null(P.g)) NULL else lmsGraphStatesProductsCpp(
      M, rule$n, submodel$info$numXis, submodel$info$numEtas
    )
    states <- derived$states

    # Structural sufficient statistics for the ECM M-step. Formed here, at the
    # E-step parameters, because the whole point is that the latent nodes are
    # then held fixed: with z fixed the complete-data objective separates into
    # a measurement half that needs the N by Q kernel and a structural half
    # that needs only these moments. See lmsGraphStructuralStatsCpp.
    structural <- if (is.null(P.g)) NULL else lmsGraphStructuralStatsCpp(
      states, derived$products, P.g, isTRUE(rule$common)
    )

    # The structural half at the E-step parameters. Kept so that
    # `completeAtEstepLmsGraph` can still report Q without a pass over the data:
    # the measurement half is the log kernel already in hand, and this is the
    # rest. O(d^3), so it costs nothing to form here.
    structural.ll <- if (is.null(structural)) NULL else
      lmsGraphStructuralCompleteCpp(
        M, structural$S, structural$W, submodel$info$numXis,
        submodel$info$numEtas, structural$numProducts
      )

    groups[[g]] <- list(
      P = P.g,
      posterior = posterior,
      log.kernel = aggregate$logKernel %||% NULL,
      V = rule$n, w = rule$w, quad = rule,
      states = states,
      obsLL = aggregate$logLik,
      sampling.weights = sampling.weights, link = link,
      workspace = workspace,
      common = isTRUE(rule$common),
      structural = structural,
      structural.ll = structural.ll
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


# The OBSERVED branch re-derives the latent states from the current parameters,
# because it is the marginal likelihood: the nodes are innovations and the map
# from them to z is exactly what theta controls. The COMPLETE branch does not,
# because the complete data is z itself (#28) -- it reads the states frozen at
# the E-step and adds the latent density that the freezing leaves behind.
lmsGraphObjective <- function(theta, model, P, observed = FALSE, sign = -1) {
  filled <- fillModel(model = model, theta = theta, method = "lms")
  ll <- 0
  for (g in seq_len(model$info$n.groups)) {
    Pg <- P$P_GROUPS[[g]]
    submodel <- filled$models[[g]]
    M <- lmsGraphMatrices(submodel)
    logistic <- identical(Pg$link %||% "logit", "logit")
    workspace <- Pg$workspace %||% lmsGraphWorkspace(submodel)

    ll <- ll + if (observed) {
      lmsGraphAggregateCpp(
        lmsGraphLogKernel(submodel, Pg$V, Pg$link %||% "logit", workspace,
                          shared = isTRUE(Pg$common)),
        Pg$w, Pg$sampling.weights
      )$logLik

    } else {
      stats <- Pg$structural
      lmsGraphCompleteStatesCpp(
        M, Pg$states, workspace, Pg$P, logistic,
        ThreadEnv$n.threads %||% 1L, isTRUE(Pg$common)
      ) + if (is.null(stats)) 0 else lmsGraphStructuralCompleteCpp(
        M, stats$S, stats$W, submodel$info$numXis, submodel$info$numEtas,
        stats$numProducts
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

    # Complete data is z, frozen (#28), so this half sees the measurement
    # parameters only -- no reverse traversal of the structural graph, and no
    # Cholesky direction. The observed likelihood is the marginal one and still
    # needs the full reverse pass, because there the nodes ARE innovations.
    score <- if (observed) lmsGraphReverseScoreCpp(
      M, Pg$V, submodel$info$numXis, submodel$info$numEtas,
      submodel$data$data.split, submodel$data$colidx0, ordered,
      Pg$w, Pg$sampling.weights, matrix(numeric(), 0L, 0L),
      TRUE, loc$block, loc$row, loc$col, loc$symmetric,
      logistic, ThreadEnv$n.threads %||% 1L, workspace,
      isTRUE(Pg$common)
    ) else lmsGraphMeasurementScoreCpp(
      M, Pg$states, workspace, Pg$P,
      loc$block, loc$row, loc$col,
      logistic, ThreadEnv$n.threads %||% 1L, isTRUE(Pg$common)
    )
    raw[loc$param] <- score
  }
  J <- lmsFirstDerivativeJacobian(theta, model)
  out <- drop(sign * J %*% raw)
  if (observed) out else out + lmsGraphStructuralGradient(theta, model, P,
                                                          sign, J)
}


# Gradient of the latent-density half of the complete-data objective, by central
# differences.
#
# The term is O(d^3) with d = 1 + numXis + numEtas + numProducts and touches no
# data at all -- about 15 microseconds an evaluation -- and it is differenced
# only over the parameters that actually reach a structural location, read off
# the Jacobian's sparsity rather than off names so that label constraints are
# handled. A parameter shared with a measurement location is included here too:
# the two halves are added, so each contributes its own part of the chain rule.
#
# An analytic score would need dB/dtheta and dSigma/dtheta for every structural
# block. Worth writing if this shows up in a profile against the N by Q
# measurement pass; it has not.
lmsGraphStructuralGradient <- function(theta, model, P, sign, J = NULL,
                                       epsilon = .Machine$double.eps ^ (1 / 3)) {
  out <- stats::setNames(numeric(length(theta)), names(theta))
  locations <- model$params$gradientStruct$locations
  structural.location <- !(locations$block %in% lmsGraphMeasurementBlocks())
  if (!any(structural.location)) return(out)
  if (is.null(J)) J <- lmsFirstDerivativeJacobian(theta, model)
  index <- which(as.logical(
    abs(J[, structural.location, drop = FALSE]) %*%
      rep(1, sum(structural.location)) > 0
  ))
  if (!length(index)) return(out)

  # A step can put a covariance parameter somewhere the Cholesky fails. That is
  # a probe outside the feasible set, not an error in the objective, so fall
  # back to the one-sided difference rather than propagating it.
  value <- function(x) tryCatch(
    structuralCompLogLikLmsGraph(x, model, P, sign),
    error = function(e) NA_real_
  )
  base <- value(theta)
  step <- epsilon * pmax(1, abs(theta[index]))
  for (k in seq_along(index)) {
    j <- index[k]
    plus <- minus <- theta
    plus[j] <- plus[j] + step[k]
    minus[j] <- minus[j] - step[k]
    up <- value(plus)
    down <- value(minus)
    out[j] <- if (is.finite(up) && is.finite(down))
      (up - down) / (2 * step[k])
    else if (is.finite(up) && is.finite(base)) (up - base) / step[k]
    else if (is.finite(down) && is.finite(base)) (base - down) / step[k]
    else 0
  }
  out
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

# Blocks whose parameters enter the likelihood only through the measurement
# model, i.e. through `means = tau + states %*% t(lambda)` and the thresholds.
# Everything else -- A, psi, alpha, beta0, gammaXi, gammaEta, covZetaXi, omega
# -- enters only through the latent density once the nodes are held fixed.
#
# A function rather than a constant: `DA_BLOCKS` lives in
# R/model_parameters_da.R, which is sourced after this file, so a top-level
# reference would not resolve at load time.
lmsGraphMeasurementBlocks <- function() {
  unlist(DA_BLOCKS[c("lambdaX", "lambdaY", "tauX", "tauY", "thetaDelta",
                     "thetaEpsilon", "thresholdDelta")], use.names = FALSE)
}


# Split the free parameter vector into the two ECM blocks. A label constraint
# can tie one free parameter to several matrix entries, so the split is taken
# from the Jacobian's sparsity rather than from names: a parameter counts as
# measurement if it touches ANY measurement location. That keeps a parameter
# shared between the two halves in the measurement step, which uses the full
# objective and is therefore correct for it either way -- the structural step
# is the one that would be wrong, since it cannot see the kernel.
lmsGraphEcmPartition <- function(theta, model) {
  locations <- model$params$gradientStruct$locations
  J <- lmsFirstDerivativeJacobian(theta, model)
  measurement.location <- locations$block %in% lmsGraphMeasurementBlocks()
  touches <- function(mask) {
    if (!any(mask)) return(rep(FALSE, length(theta)))
    abs(J[, mask, drop = FALSE]) %*% rep(1, sum(mask)) > 0
  }
  measurement <- as.logical(touches(measurement.location))
  structural <- as.logical(touches(!measurement.location)) & !measurement
  list(measurement = which(measurement), structural = which(structural))
}


# The complete-data objective at the E-step parameters, for free.
#
# The measurement half of `lmsGraphObjective()` reduces at the E-step point to
#   sum_i sum_q P_iq * sum_j log f(x_ij | z_iq)
# and the inner sum over indicators is exactly the log kernel the E-step
# already formed. So the value at the point the M-step starts from needs no
# pass over the data at all -- it is a dot product of two matrices in hand.
# Line searches need that starting value, and paying an N by Q pass for it
# would eat most of what a cheap step saves.
#
# The latent-density half (#28) is not in the kernel, so it is added from the
# value the E-step recorded. Dropping it here would make this disagree with
# `compLogLikLmsGraph` at the same theta by exactly that term, which is the
# quantity a line search is trying to improve.
completeAtEstepLmsGraph <- function(P, sign = -1) {
  ll <- 0
  for (Pg in P$P_GROUPS) {
    if (is.null(Pg$P) || is.null(Pg$log.kernel)) return(NA_real_)
    ll <- ll + sum(Pg$P * Pg$log.kernel) + (Pg$structural.ll %||% 0)
  }
  sign * ll
}


# Matrices `lmsGraphStructuralCompleteCpp` actually reads: the latent covariance
# blocks (A, psi, covZetaXi) and the structural map (gammaXi, gammaEta, omega,
# beta0, alpha). `productDesign` is a static design and is never filled, so it
# survives untouched from the template.
LMS_GRAPH_STRUCTURAL_MATRICES <- c("A", "psi", "covZetaXi", "gammaXi",
                                   "gammaEta", "omega", "beta0", "alpha")


# `fillModel` restricted to the structural half.
#
# The structural objective is O(d^3) and touches no data -- 15 microseconds in
# the compiled routine -- but the ECM structural block evaluates it thousands of
# times per M-step, because nlminb finite-differences it. A full `fillModel`
# costs ~445 of those microseconds, thirty times the callee's own work, spent
# filling matrices it never reads: lambdaX, tauX, thetaDelta, and on ordered data
# the whole thresholdDelta -> thresholds cumulative softplus, which is an R loop
# over indicators. That was 12.6% of an ordered run.
#
# This mirrors `fillMainModel` for the eight matrices above and nothing else. It
# must keep mirroring it: a change to how those are filled has to land here too,
# which is what the equivalence test in test_lms_graph.R checks.
fillStructuralLmsGraph <- function(model, theta) {
  params <- model$params
  if (is.null(names(theta))) names(theta) <- names(params$theta)

  thetaLabel <- NULL
  if (length(params$SELECT_THETA_LAB))
    thetaLabel <- suppressWarnings(calcThetaLabel(
      theta[params$SELECT_THETA_LAB[[1L]]], params$constrExprs
    ))

  lapply(seq_len(model$info$n.groups), function(g) {
    submodel <- model$models[[g]]
    M <- submodel$matrices
    keep <- intersect(LMS_GRAPH_STRUCTURAL_MATRICES, names(M))
    M[keep] <- fillMatricesLabels(M[keep], submodel$labelMatrices[keep],
                                  thetaLabel)
    thetaMain <- theta[params$SELECT_THETA_MAIN[[g]]]

    if (!is.null(submodel$covModel$matrices)) {
      thetaCov <- if (length(params$SELECT_THETA_COV[[g]]))
        theta[params$SELECT_THETA_COV[[g]]] else NULL
      M$A <- expectedCovModel(
        fillCovModel(submodel$covModel, thetaCov, thetaLabel),
        method = "lms", sortedXis = submodel$info$xis
      )
    } else {
      M$A <- fillNA_Matrix(M$A, theta = thetaMain, pattern = "^A([0-9]*)")
    }

    M$covZetaXi <- fillNA_Matrix(M$covZetaXi, theta = thetaMain,
                                 pattern = "^covZetaXi")
    M$psi       <- fillSymmetric(M$psi, fetch(thetaMain, "^psi"))
    M$alpha     <- fillNA_Matrix(M$alpha, theta = thetaMain, pattern = "^alpha")
    M$beta0     <- fillNA_Matrix(M$beta0, theta = thetaMain, pattern = "^beta0")
    M$gammaEta  <- fillNA_Matrix(M$gammaEta, theta = thetaMain,
                                 pattern = "^gammaEta")
    M$gammaXi   <- fillNA_Matrix(M$gammaXi, theta = thetaMain,
                                 pattern = "^gammaXi")
    if (!is.null(M$omega))
      M$omega <- fillNA_Matrix(M$omega, theta = thetaMain,
                               pattern = "^omega[0-9]*$")
    M
  })
}


# Structural half of the complete-data objective, from the statistics the
# E-step formed. Costs O(d^3) with d = 1 + numXis + numEtas + numProducts --
# no pass over the data at all.
structuralCompLogLikLmsGraph <- function(theta, model, P, sign = -1, ...) {
  matrices <- fillStructuralLmsGraph(model, theta)
  ll <- 0
  for (g in seq_len(model$info$n.groups)) {
    stats <- P$P_GROUPS[[g]]$structural
    if (is.null(stats)) next
    info <- model$models[[g]]$info
    ll <- ll + lmsGraphStructuralCompleteCpp(
      matrices[[g]], stats$S, stats$W, info$numXis, info$numEtas,
      stats$numProducts
    )
  }
  sign * ll
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


# Blocks the measurement Newton step understands. `lambdaY`, `tauY` and
# `thetaEpsilon` are absent by construction on the graph backend, which routes
# every indicator through the X matrices.
LMS_GRAPH_NEWTON_BLOCKS <- c(lambdaX = 0, tauX = 2, thetaDelta = 4,
                             thresholdDelta = 18)


# Latent columns whose loading on each indicator can move. The traversal
# carries one Hessian slot per entry, so keeping this to the loadings that are
# actually nonzero or free is what keeps the blocks small under simple
# structure -- one latent per indicator rather than all of them.
lmsGraphActiveLatents <- function(lambda, locations) {
  loadings <- locations[locations$block ==
                        LMS_GRAPH_NEWTON_BLOCKS[["lambdaX"]], , drop = FALSE]
  lapply(seq_len(NROW(lambda)), function(indicator) {
    free <- loadings$col[loadings$row == indicator - 1L]
    sort(unique(c(which(lambda[indicator, ] != 0) - 1L, free)))
  })
}


# Express each measurement free parameter as a direction in its indicator's
# block coordinates ([intercept, loadings, variance | thresholds]).
#
# Returns NULL whenever the block split is not visibly true, and the caller
# then keeps the nlminb path -- always correct, merely slower. It is refused
# for a free parameter outside the four matrices above, an off-diagonal
# residual covariance, a parameter spanning two indicators (a shared label or
# equality constraint), or a nonlinear parameter map, whose second-order term
# this does not carry.
lmsGraphNewtonPlan <- function(theta, model, measurement, submodel, M,
                               activeLatents, group = 1L) {
  gradientStruct <- model$params$gradientStruct
  if (length(gradientStruct$nlinDerivs)) return(NULL)

  locations <- gradientStruct$locations
  J <- lmsFirstDerivativeJacobian(theta, model)
  if (is.null(locations) || !length(measurement)) return(NULL)

  inGroup <- locations$group == group
  thresholdDelta <- M$thresholdDelta
  indicators <- NROW(M$lambdaX)
  # `lmsGraphOrderedIndex()` returns 0-based indicator indices, not a mask.
  ordered <- logical(indicators)
  ordered[lmsGraphOrderedIndex(submodel, M) + 1L] <- TRUE

  # slots per indicator, matching the layout the traversal returns
  slots <- vapply(seq_len(indicators), FUN.VALUE = integer(1L),
                  function(l) {
    extra <- if (ordered[l]) sum(is.finite(thresholdDelta[l, ])) else 1L
    1L + length(activeLatents[[l]]) + as.integer(extra)
  })

  direction <- vector("list", indicators) # slots x |theta| per indicator
  # Second-order terms from the threshold map, one entry per `thresholdDelta`
  # location: they are only known once the block gradient is in hand.
  deltaTerms <- vector("list", indicators)
  owner <- integer(length(measurement))

  for (i in seq_along(measurement)) {
    t <- measurement[i]
    touched <- which(J[t, ] != 0 & inGroup)
    if (!length(touched)) return(NULL)
    if (!all(locations$block[touched] %in% LMS_GRAPH_NEWTON_BLOCKS)) return(NULL)

    rows <- unique(locations$row[touched])
    if (length(rows) != 1L) return(NULL) # spans two indicators
    indicator <- rows + 1L
    owner[i] <- indicator

    if (is.null(direction[[indicator]])) {
      direction[[indicator]] <- matrix(0, slots[indicator], length(measurement))
      deltaTerms[[indicator]] <- list()
    }

    active <- activeLatents[[indicator]]
    meanSlots <- 1L + length(active)
    use <- if (ordered[indicator]) which(is.finite(thresholdDelta[indicator, ]))
           else integer(0)

    for (loc in touched) {
      weight <- J[t, loc]
      block <- locations$block[loc]
      col <- locations$col[loc]

      if (block == LMS_GRAPH_NEWTON_BLOCKS[["tauX"]]) {
        direction[[indicator]][1L, i] <- direction[[indicator]][1L, i] + weight

      } else if (block == LMS_GRAPH_NEWTON_BLOCKS[["lambdaX"]]) {
        slot <- match(col, active)
        if (is.na(slot)) return(NULL)
        direction[[indicator]][1L + slot, i] <-
          direction[[indicator]][1L + slot, i] + weight

      } else if (block == LMS_GRAPH_NEWTON_BLOCKS[["thetaDelta"]]) {
        if (locations$row[loc] != col) return(NULL) # off-diagonal residual
        direction[[indicator]][slots[indicator], i] <-
          direction[[indicator]][slots[indicator], i] + weight

      } else { # thresholdDelta: thresholds are a cumulative softplus of delta
        delta <- thresholdDelta[indicator, col + 1L]
        if (!is.finite(delta)) return(NULL)
        slope <- if (col == 0L) 1 else 1 / (1 + exp(-delta))
        reached <- which(use - 1L >= col)
        if (!length(reached)) next
        direction[[indicator]][meanSlots + reached, i] <-
          direction[[indicator]][meanSlots + reached, i] + weight * slope
        # The curvature of that map contributes a term that J' H J cannot
        # see, sum_k (dQ/dtau_k) * d2 tau_k / d delta^2, with
        # d2 tau / d delta^2 = slope * (1 - slope). It needs the gradient, so
        # only its shape is recorded here.
        if (col > 0L) {
          key <- paste(loc)
          term <- deltaTerms[[indicator]][[key]]
          if (is.null(term))
            term <- list(v = numeric(length(measurement)),
                         slots = meanSlots + reached,
                         factor = slope * (1 - slope))
          term$v[i] <- term$v[i] + weight
          deltaTerms[[indicator]][[key]] <- term
        }
      }
    }
  }

  keep <- which(!vapply(direction, is.null, logical(1L)))
  if (!length(keep)) return(NULL)
  list(theta = measurement, owner = owner, indicators = keep,
       direction = direction, deltaTerms = deltaTerms, slots = slots)
}


# Newton system for one indicator: its free parameters, the gradient and the
# Hessian of the complete-data objective with respect to them.
lmsGraphNewtonSystem <- function(plan, indicator, gradient, hessian) {
  B <- plan$direction[[indicator]]
  keep <- which(colSums(abs(B)) > 0)
  if (!length(keep)) return(NULL)
  B <- B[, keep, drop = FALSE]

  H <- crossprod(B, hessian %*% B)
  for (term in plan$deltaTerms[[indicator]]) {
    v <- term$v[keep]
    H <- H + (sum(gradient[term$slots]) * term$factor) * tcrossprod(v)
  }

  list(theta = plan$theta[keep], gradient = drop(crossprod(B, gradient)),
       hessian = H)
}


# One Newton step on the measurement block.
#
# With the latent nodes held fixed the measurement objective separates across
# indicators, so its Hessian is block diagonal and a single traversal yields
# every block (see lmsGraphMeasurementNewtonCpp).
#
# Mplus does the same thing rather than anything cleverer. Measured against the
# demo binary by varying MUITERATIONS at fixed MITERATIONS: one of its M-step
# iterations costs 0.110 s at N=500 and 0.465 s at N=2000, so N x 4.0 gives
# cost x 4.23 -- its M-step iterations re-traverse the data exactly as this
# does, and there is no N-independent sufficient statistic to exploit (#25).
#
# ONE step is also the right count, not Newton to convergence. A default-
# settings TECH5 trace prints exactly one measurement iteration in every
# section -- 499 of them under EM, 289 under EMA, counted mechanically over the
# whole file (#26). Iterating this block to 3 or 5 steps was tried and moved
# the EMA endgame by 15 logLik against a 110 gap, so the single step is not
# what holds `newton` back; see mstepLmsGraphEcm for the full table.
#
# Returns NULL when the block split does not hold, and the caller then keeps
# the nlminb path.
lmsGraphNewtonMeasurement <- function(theta, model, P, measurement, lower,
                                      upper, objective, link = "logit") {
  if (model$info$n.groups > 1L) return(NULL) # a shared label would span groups

  filled <- fillModel(model = model, theta = theta, method = "lms")
  submodel <- filled$models[[1L]]
  M <- lmsGraphMatrices(submodel)
  if (!NCOL(M$thresholds)) M$thresholdDelta <- M$thresholds

  locations <- model$params$gradientStruct$locations
  active <- lmsGraphActiveLatents(
    M$lambdaX, locations[locations$group == 1L, , drop = FALSE]
  )
  plan <- lmsGraphNewtonPlan(theta, model, measurement, submodel, M, active)
  if (is.null(plan)) return(NULL)

  Pg <- P$P_GROUPS[[1L]]
  fused <- lmsGraphMeasurementNewtonCpp(
    M, Pg$states,
    submodel$data$data.split, submodel$data$colidx0,
    lmsGraphOrderedIndex(submodel, M), Pg$P, active,
    identical(link, "logit"), ThreadEnv$n.threads %||% 1L,
    Pg$workspace %||% lmsGraphWorkspace(submodel), isTRUE(Pg$common)
  )

  direction <- numeric(length(theta))
  for (indicator in plan$indicators) {
    system <- lmsGraphNewtonSystem(plan, indicator, fused$gradient[[indicator]],
                                   fused$hessian[[indicator]])
    if (is.null(system)) next

    # The objective is maximised here, so its Hessian must be negative
    # definite for the Newton direction to ascend. `chol()` errors otherwise,
    # and rather than patch a single block we hand the whole M-step back to
    # nlminb -- indefiniteness means the quadratic model is not trustworthy,
    # and this is rare near the optimum.
    step <- tryCatch({
      R <- chol(-system$hessian)
      backsolve(R, backsolve(R, system$gradient, transpose = TRUE))
    }, error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) return(NULL)

    direction[system$theta] <- step
  }
  if (!any(direction != 0)) return(NULL)

  # The measurement block runs first in the ECM, so the starting objective is
  # the one at the E-step parameters -- available without a pass.
  start <- completeAtEstepLmsGraph(P, sign = -1)
  if (!is.finite(start)) start <- objective(theta)

  # Backtrack on the whole direction. A full Newton step is right when the
  # quadratic model holds; halving twice covers the early iterations where it
  # does not, and is still far cheaper than the line search it replaces.
  for (scale in c(1, 0.5, 0.25)) {
    candidate <- theta
    candidate[measurement] <- pmin(pmax(
      theta[measurement] + scale * direction[measurement],
      lower[measurement]), upper[measurement])
    value <- objective(candidate)
    if (is.finite(value) && value <= start)
      return(list(theta = candidate, objective = value, iterations = 1L,
                  convergence = 0L, message = "newton"))
  }
  NULL
}
