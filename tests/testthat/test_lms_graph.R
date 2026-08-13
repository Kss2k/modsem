testthat::test_that("LMS graph backend transforms one unified latent covariance", {
  M <- list(
    A = matrix(c(1.2, 0.3, 0, 0.8), 2L, byrow = TRUE),
    covZetaXi = matrix(c(.15, -.1), 1L),
    psi = matrix(.7, 1L)
  )
  expected <- rbind(
    cbind(M$A %*% t(M$A), t(M$covZetaXi)),
    cbind(M$covZetaXi, M$psi)
  )
  testthat::expect_equal(lmsGraphLatentCovariance(M), expected)
  testthat::expect_equal(t(lmsGraphLatentCovariance(M)), expected)
})


# The covariance directions (blocks 6, 7 and 17) used to be overwritten with
# forward differences, and this is where that hybrid was tested. The analytic
# scores turned out to be correct -- they are checked against central
# differences in "analytical LMS graph scores match central differences" below,
# whose `blocks` vector already covers 6, 7 and 17 -- so the override was pure
# cost on the per-observation path, where every objective evaluation is a full
# N-by-Q pass. Removed; the coverage lives in that test.


testthat::test_that("LMS graph states follow recursive structural order", {
  M <- list(
    A = matrix(1, 1L), covZetaXi = matrix(0, 2L, 1L), psi = diag(2L),
    beta0 = matrix(.2), alpha = matrix(c(.1, -.3), 2L),
    gammaXi = matrix(c(.5, -.2), 2L),
    gammaEta = matrix(c(0, 0, .4, 0), 2L, byrow = TRUE),
    productDesign = matrix(c(2L, 0L, 0L,
                            1L, 1L, 0L), 2L, 3L, byrow = TRUE),
    omega = matrix(c(.25, 0, 0, .3), 2L, 2L),
    thetaDelta = diag(3L), lambdaX = diag(3L), tauX = matrix(0, 3L, 1L)
  )
  submodel <- list(
    matrices = M,
    info = list(numXis = 1L, numEtas = 2L, xis = "x", etas = c("e1", "e2"),
                hasComposites = FALSE),
    quad = list(adaptive = FALSE)
  )
  nodes <- rbind(c(0, 0, 0), c(1, .5, -.25))
  states <- lmsGraphStates(submodel, nodes)

  x <- nodes[, 1] + .2
  e1 <- .1 + .5 * x + .25 * x^2 + nodes[, 2]
  e2 <- -.3 - .2 * x + .4 * e1 + .3 * x * e1 + nodes[, 3]
  testthat::expect_equal(states, cbind(x = x, e1 = e1, e2 = e2),
                         tolerance = 1e-12)
})


testthat::test_that("LMS graph rejects indicator residual covariances", {
  submodel <- list(
    matrices = list(thetaDelta = matrix(c(1, .1, .1, 1), 2L)),
    info = list(hasComposites = FALSE), quad = list(adaptive = FALSE)
  )
  testthat::expect_error(lmsGraphValidate(submodel),
                         "residual covariances between indicators")
})


testthat::test_that("compiled LMS graph kernel matches mixed logit and probit likelihoods", {
  M <- list(
    A = matrix(1, 1L), covZetaXi = matrix(numeric(), 0L, 1L),
    psi = matrix(numeric(), 0L, 0L), beta0 = matrix(0),
    alpha = matrix(numeric(), 0L, 1L), gammaXi = matrix(numeric(), 0L, 1L),
    gammaEta = matrix(numeric(), 0L, 0L),
    productDesign = matrix(0L, 0L, 1L), omega = matrix(numeric(), 0L, 0L),
    lambdaX = matrix(c(1, .5), 2L), tauX = matrix(c(0, .2), 2L),
    thetaDelta = diag(c(1, .7)),
    thresholds = matrix(c(-.4, .8, NaN, NaN), 2L, byrow = TRUE)
  )
  base.nodes <- matrix(c(-1, .25, 1.2), ncol = 1L)
  values <- matrix(c(1, .1, 2, -.3, 3, .8), 3L, byrow = TRUE)
  nodes <- base.nodes[rep(seq_len(NROW(base.nodes)), times = NROW(values)),
                      , drop = FALSE]

  direct <- function(link) {
    eta <- base.nodes[, 1L]
    out <- matrix(0, 3L, 3L)
    F <- if (link == "logit") plogis else pnorm
    for (i in seq_len(3L)) for (q in seq_len(3L)) {
      code <- values[i, 1L]
      lower <- if (code == 1L) -Inf else c(-.4, .8)[code - 1L]
      upper <- if (code == 3L) Inf else c(-.4, .8)[code]
      ordinal <- log(F(upper - eta[q]) - F(lower - eta[q]))
      continuous <- dnorm(values[i, 2L], .2 + .5 * eta[q], sqrt(.7), log = TRUE)
      out[i, q] <- ordinal + continuous
    }
    out
  }

  workspace <- lmsGraphWorkspaceCpp(list(values), list(0:1), 0L)
  for (link in c("logit", "probit")) {
    compiled <- lmsGraphLogKernelCpp(
      M, nodes, 1L, 0L, list(values), list(0:1), 0L, link == "logit"
    )
    testthat::expect_equal(compiled, direct(link), tolerance = 1e-12)
    cached <- lmsGraphLogKernelWorkspaceCpp(
      M, nodes, 1L, 0L, workspace, link == "logit"
    )
    testthat::expect_equal(cached, compiled, tolerance = 1e-12)
    complete.weights <- matrix(seq_len(length(cached)), NROW(cached))
    testthat::expect_equal(
      lmsGraphCompleteStatesCpp(
        M, lmsGraphStatesCpp(M, nodes, 1L, 0L), workspace, complete.weights,
        link == "logit"
      ),
      sum(compiled * complete.weights), tolerance = 1e-12
    )
  }
})


testthat::test_that("packed graph rules support observation-specific nodes and weights", {
  M <- list(
    A = matrix(1), covZetaXi = matrix(numeric(), 0L, 1L),
    psi = matrix(numeric(), 0L, 0L), beta0 = matrix(0),
    alpha = matrix(numeric(), 0L, 1L), gammaXi = matrix(numeric(), 0L, 1L),
    gammaEta = matrix(numeric(), 0L, 0L),
    productDesign = matrix(0L, 0L, 1L), omega = matrix(numeric(), 0L, 0L),
    lambdaX = matrix(1), tauX = matrix(0), thetaDelta = matrix(1),
    thresholds = matrix(NaN, 1L, 0L)
  )
  values <- matrix(c(.2, -.4), 2L)
  nodes <- matrix(c(-1, 1, -2, 2), ncol = 1L) # two nodes per observation
  kernel <- lmsGraphLogKernelCpp(
    M, nodes, 1L, 0L, list(values), list(0L), integer(), TRUE
  )
  expected <- rbind(
    dnorm(values[1L], c(-1, 1), log = TRUE),
    dnorm(values[2L], c(-2, 2), log = TRUE)
  )
  testthat::expect_equal(kernel, expected, tolerance = 1e-12)

  packed.weights <- c(.25, .75, .6, .4)
  aggregate <- lmsGraphAggregateCpp(kernel, packed.weights, c(1, 1))
  expected.density <- c(
    log(sum(exp(expected[1L, ]) * packed.weights[1:2])),
    log(sum(exp(expected[2L, ]) * packed.weights[3:4]))
  )
  testthat::expect_equal(as.numeric(aggregate$logDensity), expected.density,
                         tolerance = 1e-12)
})


# The two tests that lived here ("common graph sufficient statistics reproduce
# packed M-step" and "common graph statistics respect missingness patterns")
# checked that the shared-rule sufficient-statistic M-step agreed with the
# packed one. That second implementation is gone -- a shared rule is now read
# from the same Q nodes through a zero stride rather than through its own
# kernel and M-step.
#
# What replaces them is the equivalence the stride actually claims: reading Q
# shared nodes with stride 0 must give exactly what replicating those nodes N
# times and reading them with stride Q gives. The replicated side is the
# ordinary packed path, which the score and likelihood tests below already
# pin down, so this anchors the shared path to something independently tested.
testthat::test_that("a shared rule read with zero stride equals the replicated rule", {
  delta <- matrix(c(-.4, log(expm1(1.2)), NaN, NaN), 2L, byrow = TRUE)
  M <- list(
    A = matrix(1.1), covZetaXi = matrix(.12), psi = matrix(.8),
    beta0 = matrix(.15), alpha = matrix(-.1),
    gammaXi = matrix(.35), gammaEta = matrix(0),
    productDesign = matrix(c(2L, 0L), 1L, 2L, byrow = TRUE),
    omega = matrix(.18, 1L, 1L),
    lambdaX = matrix(c(1, .7, .45, 1.1), 2L, byrow = TRUE),
    tauX = matrix(c(0, .2)), thetaDelta = diag(c(1, .65)),
    thresholdDelta = delta, thresholds = thresholdDeltaToThresholdMatrix(delta)
  )
  shared.nodes <- as.matrix(expand.grid(c(-1, .2, 1.1), c(-.8, .4, 1.3)))
  Q <- NROW(shared.nodes)
  values <- matrix(c(1, -.2, 2, .4, 3, 1.0), 3L, byrow = TRUE)
  N <- NROW(values)
  row.weights <- c(1, 1.5, .7)
  shared.weights <- rep(1 / Q, Q)

  # Same rule, the two ways of storing it.
  packed.nodes <- shared.nodes[rep(seq_len(Q), times = N), , drop = FALSE]
  packed.weights <- rep(shared.weights, times = N)

  workspace <- lmsGraphWorkspaceCpp(list(values), list(0:1), 0L)

  shared.kernel <- lmsGraphLogKernelWorkspaceCpp(
    M, shared.nodes, 1L, 1L, workspace, TRUE, 1L, TRUE)
  packed.kernel <- lmsGraphLogKernelWorkspaceCpp(
    M, packed.nodes, 1L, 1L, workspace, TRUE, 1L, FALSE)
  testthat::expect_equal(shared.kernel, packed.kernel, tolerance = 1e-12)

  shared.step <- lmsGraphPstepWorkspaceCpp(
    M, shared.nodes, 1L, 1L, workspace, shared.weights, row.weights,
    TRUE, 1L, TRUE)
  packed.step <- lmsGraphPstepWorkspaceCpp(
    M, packed.nodes, 1L, 1L, workspace, packed.weights, row.weights,
    TRUE, 1L, FALSE)
  testthat::expect_equal(shared.step$logLik, packed.step$logLik,
                         tolerance = 1e-12)
  testthat::expect_equal(shared.step$posterior, packed.step$posterior,
                         tolerance = 1e-12)

  posterior <- packed.step$posterior * row.weights
  testthat::expect_equal(
    lmsGraphCompleteStatesCpp(M, lmsGraphStatesCpp(M, shared.nodes, 1L, 1L),
                              workspace, posterior, TRUE, 1L, TRUE),
    lmsGraphCompleteStatesCpp(M, lmsGraphStatesCpp(M, packed.nodes, 1L, 1L),
                              workspace, posterior, TRUE, 1L, FALSE),
    tolerance = 1e-12)

  # The score is the path with genuinely new code: with a shared rule every
  # observation accumulates into the same Q rows of the response adjoint, so
  # it is reduced per thread rather than written in place.
  blocks <- c(0L, 2L, 4L, 6L, 7L, 8L, 9L, 10L, 11L, 19L, 17L, 18L, 18L)
  rows <- c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)
  cols <- c(1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L)
  symmetric <- blocks %in% c(4L, 7L)
  for (observed in c(TRUE, FALSE)) {
    complete.weights <- if (observed) matrix(numeric(), 0L, 0L) else posterior
    for (threads in c(1L, 2L)) {
      shared.score <- lmsGraphReverseScoreCpp(
        M, shared.nodes, 1L, 1L, list(values), list(0:1), 0L,
        shared.weights, row.weights, complete.weights, observed,
        blocks, rows, cols, symmetric, TRUE, threads, workspace, TRUE)
      packed.score <- lmsGraphReverseScoreCpp(
        M, packed.nodes, 1L, 1L, list(values), list(0:1), 0L,
        packed.weights, row.weights, complete.weights, observed,
        blocks, rows, cols, symmetric, TRUE, threads, workspace, FALSE)
      testthat::expect_equal(as.numeric(shared.score),
                             as.numeric(packed.score), tolerance = 1e-10)
    }
  }
})

testthat::test_that("per-observation adaptive graph quadrature is exact for Gaussian rows", {
  values <- matrix(c(-1, .4, 1.5), ncol = 1L)
  M <- list(
    A = matrix(1), covZetaXi = matrix(numeric(), 0L, 1L),
    psi = matrix(numeric(), 0L, 0L), beta0 = matrix(0),
    alpha = matrix(numeric(), 0L, 1L), gammaXi = matrix(numeric(), 0L, 1L),
    gammaEta = matrix(numeric(), 0L, 0L),
    productDesign = matrix(0L, 0L, 1L), omega = matrix(numeric(), 0L, 0L),
    lambdaX = matrix(1, dimnames = list("y", "x")),
    tauX = matrix(0), thetaDelta = matrix(1),
    thresholds = matrix(NaN, 1L, 0L)
  )
  submodel <- list(
    matrices = M,
    info = list(numXis = 1L, numEtas = 0L, ordered = character()),
    data = list(n = NROW(values), data.split = list(values),
                colidx0 = list(0L)),
    quad = quadSpec(integration = "gh", adaptive = "full", nodes = 3L, k = 1L)
  )
  workspace <- lmsGraphWorkspace(submodel)
  base <- buildQuadRule(setAdaptiveQuad(submodel$quad, "none"))
  rule <- lmsGraphAdaptPerObservation(submodel, base, "logit", workspace,
                                      tolerance = 1e-10)
  kernel <- lmsGraphLogKernelCpp(
    M, rule$n, 1L, 0L, list(values), list(0L), integer(), TRUE
  )
  aggregate <- lmsGraphAggregateCpp(kernel, rule$w, rep(1, NROW(values)))

  testthat::expect_false(rule$common)
  testthat::expect_true(rule$packed)
  # For a Gaussian row the posterior mode and the marginal density are known
  # in closed form, so these are exact checks rather than agreement checks.
  testthat::expect_equal(rule$modes[, 1L], values[, 1L] / 2,
                         tolerance = 1e-5)
  testthat::expect_equal(as.numeric(aggregate$logDensity),
                         dnorm(values[, 1L], 0, sqrt(2), log = TRUE),
                         tolerance = 1e-7)
})


testthat::test_that("adaptive ordered likelihood agrees with a dense fixed rule", {
  values <- matrix(1:3, ncol = 1L)
  M <- list(
    A = matrix(1), covZetaXi = matrix(numeric(), 0L, 1L),
    psi = matrix(numeric(), 0L, 0L), beta0 = matrix(0),
    alpha = matrix(numeric(), 0L, 1L), gammaXi = matrix(numeric(), 0L, 1L),
    gammaEta = matrix(numeric(), 0L, 0L),
    productDesign = matrix(0L, 0L, 1L), omega = matrix(numeric(), 0L, 0L),
    lambdaX = matrix(.9, dimnames = list("y", "x")), tauX = matrix(0),
    thetaDelta = matrix(1), thresholds = matrix(c(-.6, .7), 1L),
    thresholdDelta = matrix(c(-.6, log(expm1(1.3))), 1L)
  )
  submodel <- list(
    matrices = M,
    info = list(numXis = 1L, numEtas = 0L, ordered = "y"),
    data = list(n = 3L, data.split = list(values), colidx0 = list(0L)),
    quad = quadSpec(integration = "gh", adaptive = "full", nodes = 5L, k = 1L)
  )
  workspace <- lmsGraphWorkspace(submodel)
  adaptive <- lmsGraphBuildRule(submodel, workspace, "logit")
  adaptive.kernel <- lmsGraphLogKernelCpp(
    M, adaptive$n, 1L, 0L, list(values), list(0L), 0L, TRUE
  )
  adaptive.ll <- lmsGraphAggregateCpp(
    adaptive.kernel, adaptive$w, rep(1, 3L)
  )$logDensity

  dense <- buildQuadRule(quadSpec(integration = "gh", adaptive = "none",
                                  nodes = 41L, k = 1L))
  dense.nodes <- dense$n[rep(seq_len(NROW(dense$n)), times = 3L),,
                             drop = FALSE]
  dense.weights <- rep(dense$w, times = 3L)
  dense.kernel <- lmsGraphLogKernelCpp(
    M, dense.nodes, 1L, 0L, list(values), list(0L), 0L, TRUE
  )
  dense.ll <- lmsGraphAggregateCpp(
    dense.kernel, dense.weights, rep(1, 3L)
  )$logDensity
  testthat::expect_equal(as.numeric(adaptive.ll), as.numeric(dense.ll),
                         tolerance = 2e-5)
})


testthat::test_that("analytical LMS graph scores match central differences", {
  delta <- matrix(c(-.4, log(expm1(1.2)), NaN, NaN), 2L, byrow = TRUE)
  M <- list(
    A = matrix(1.1), covZetaXi = matrix(.12), psi = matrix(.8),
    beta0 = matrix(.15), alpha = matrix(-.1),
    gammaXi = matrix(.35), gammaEta = matrix(0),
    productDesign = matrix(c(2L, 0L), 1L, 2L, byrow = TRUE),
    omega = matrix(.18, 1L, 1L),
    lambdaX = matrix(c(1, .7, .45, 1.1), 2L, byrow = TRUE),
    tauX = matrix(c(0, .2)), thetaDelta = diag(c(1, .65)),
    thresholdDelta = delta, thresholds = thresholdDeltaToThresholdMatrix(delta)
  )
  base.nodes <- as.matrix(expand.grid(c(-1, .2, 1.1), c(-.8, .4, 1.3)))
  row.weights <- c(1, 1.5, .7)
  values <- matrix(c(1, -.2, 2, .4, 3, 1.0), 3L, byrow = TRUE)
  nodes <- base.nodes[rep(seq_len(NROW(base.nodes)), times = NROW(values)),
                      , drop = FALSE]
  nodes <- nodes + rep(c(0, .15, -.2), each = NROW(base.nodes))
  weight.blocks <- rbind(
    rep(1, NROW(base.nodes)),
    seq_len(NROW(base.nodes)),
    rev(seq_len(NROW(base.nodes)))
  )
  weight.blocks <- weight.blocks / rowSums(weight.blocks)
  weights <- as.numeric(t(weight.blocks))
  blocks <- c(0L, 2L, 4L, 6L, 7L, 8L, 9L, 10L, 11L, 19L, 17L, 18L, 18L)
  rows <- c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)
  cols <- c(1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L)
  symmetric <- blocks %in% c(4L, 7L)

  objective <- function(matrices) {
    kernel <- lmsGraphLogKernelCpp(
      matrices, nodes, 1L, 1L, list(values), list(0:1), 0L, TRUE
    )
    lmsGraphAggregateCpp(kernel, weights, row.weights)$logLik
  }
  perturb <- function(matrices, index, amount) {
    # Indexed by block + 1; blocks 12/13 (the old Kronecker omegas) are gone.
    names <- c("lambdaX", NA, "tauX", NA, "thetaDelta", NA, "A", "psi",
               "alpha", "beta0", "gammaXi", "gammaEta", NA, NA, NA, NA, NA,
               "covZetaXi", "thresholdDelta", "omega")
    name <- names[blocks[index] + 1L]
    matrices[[name]][rows[index] + 1L, cols[index] + 1L] <-
      matrices[[name]][rows[index] + 1L, cols[index] + 1L] + amount
    if (name == "thresholdDelta")
      matrices$thresholds <- thresholdDeltaToThresholdMatrix(matrices$thresholdDelta)
    matrices
  }
  analytical <- lmsGraphScoreCpp(
    M, nodes, 1L, 1L, list(values), list(0:1), 0L,
    weights, row.weights, matrix(numeric(), 0L, 0L), TRUE,
    blocks, rows, cols, symmetric, TRUE
  )
  reverse <- lmsGraphReverseScoreCpp(
    M, nodes, 1L, 1L, list(values), list(0:1), 0L,
    weights, row.weights, matrix(numeric(), 0L, 0L), TRUE,
    blocks, rows, cols, symmetric, TRUE
  )
  testthat::expect_equal(as.numeric(reverse), as.numeric(analytical),
                         tolerance = 1e-10)
  epsilon <- 1e-6
  numerical <- vapply(seq_along(blocks), function(index) {
    (objective(perturb(M, index, epsilon)) -
       objective(perturb(M, index, -epsilon))) / (2 * epsilon)
  }, numeric(1L))
  # The direct measurement and acyclic structural blocks have exact scores.
  # Covariance-Cholesky directions are tested separately, while self-referential
  # eta terms are outside the recursive graph model.
  exact <- c(1L, 2L, 3L, 6L, 8L, 10L, 12L, 13L)
  testthat::expect_equal(as.numeric(analytical)[exact], numerical[exact],
                         tolerance = 2e-5, scale = 1)

  kernel <- lmsGraphLogKernelCpp(M, nodes, 1L, 1L, list(values),
                                 list(0:1), 0L, TRUE)
  posterior <- lmsGraphAggregateCpp(kernel, weights, row.weights)$posterior
  posterior <- posterior * row.weights
  complete <- lmsGraphScoreCpp(
    M, nodes, 1L, 1L, list(values), list(0:1), 0L,
    weights, row.weights, posterior, FALSE,
    blocks, rows, cols, symmetric, TRUE
  )
  complete.reverse <- lmsGraphReverseScoreCpp(
    M, nodes, 1L, 1L, list(values), list(0:1), 0L,
    weights, row.weights, posterior, FALSE,
    blocks, rows, cols, symmetric, TRUE
  )
  testthat::expect_equal(as.numeric(complete.reverse), as.numeric(complete),
                         tolerance = 1e-10)
  complete.objective <- function(matrices) lmsGraphWeightedKernelCpp(
    lmsGraphLogKernelCpp(matrices, nodes, 1L, 1L, list(values),
                         list(0:1), 0L, TRUE), posterior
  )
  complete.numerical <- vapply(seq_along(blocks), function(index) {
    (complete.objective(perturb(M, index, epsilon)) -
       complete.objective(perturb(M, index, -epsilon))) / (2 * epsilon)
  }, numeric(1L))
  testthat::expect_equal(as.numeric(complete)[exact], complete.numerical[exact],
                         tolerance = 2e-5, scale = 1)
})


testthat::test_that("reverse graph score ignores structurally absent thresholds", {
  delta <- matrix(c(-.4, .2, NaN, NaN), 2L, byrow = TRUE)
  M <- list(
    A = matrix(1), covZetaXi = matrix(numeric(), 0L, 1L),
    psi = matrix(numeric(), 0L, 0L), beta0 = matrix(0),
    alpha = matrix(numeric(), 0L, 1L), gammaXi = matrix(numeric(), 0L, 1L),
    gammaEta = matrix(numeric(), 0L, 0L),
    productDesign = matrix(0L, 0L, 1L), omega = matrix(numeric(), 0L, 0L),
    lambdaX = matrix(c(1, .5), 2L), tauX = matrix(c(0, 0), 2L),
    thetaDelta = diag(2L), thresholdDelta = delta,
    thresholds = thresholdDeltaToThresholdMatrix(delta)
  )
  score <- lmsGraphReverseScoreCpp(
    M, matrix(0, 1L, 1L), 1L, 0L,
    list(matrix(c(1, 0), 1L)), list(0:1), 0L,
    1, 1, matrix(numeric(), 0L, 0L), TRUE,
    18L, 1L, 1L, FALSE, TRUE
  )
  testthat::expect_equal(as.numeric(score), 0)
})
