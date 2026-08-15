# PML: what it is checked against, and why.
#
# k = 0. With no interaction there is nothing to condition on, the quadrature
# rule collapses to a single node, and the model is an ordinary linear ordinal
# SEM -- which lavaan fits with `estimator = "PML"`. That makes lavaan an
# INDEPENDENT reference for everything else is built on: the map from model
# matrices to bivariate x* moments, the threshold parameterisation, and the
# rectangle probabilities.
#
# k > 0. No other implementation of an interaction PML exists, so the reference
# is the generating process itself: simulate from the model the matrices
# describe and check the quadrature reproduces the cell frequencies.
#
# The bivariate normal CDF is checked separately against pbivnorm, which is the
# only part with a well-known independent implementation.
#
# lavaan matching notes:
#   * `parameterization = "theta"` fixes residual variances at 1, which is the
#     probit identification PML uses.
#   * lavaan's reported `fx` is a discrepancy against the SATURATED pairwise
#     log-likelihood, not -pl/N. Undo that before comparing objectives, or the
#     two are on different scales by a constant.

skip_if_not_installed("lavaan")
skip_if_not_installed("pbivnorm")
skip_if(Sys.getenv("MODSEM_SKIP_SLOW_TESTS") == "true")

ordinalFixture <- function(categories = 4L, seed = 42L) {
  vars <- c("x1", "x2", "x3", "z1", "z2", "z3", "y1", "y2", "y3")
  set.seed(seed)
  data <- oneInt[vars]
  for (v in vars) {
    z <- (data[[v]] - mean(data[[v]])) / stats::sd(data[[v]])
    cuts <- stats::quantile(z, probs = seq(0, 1, length.out = categories + 1L))
    cuts[1L] <- -Inf; cuts[length(cuts)] <- Inf
    data[[v]] <- cut(z, breaks = cuts, ordered_result = TRUE)
  }
  list(data = data, vars = vars)
}


# The kernel works off a whole submodel, because the moments come from LMSModel.
pmlSubmodel <- function(m, numXis, numEtas, k = 0L) {
  none <- matrix(0, 0L, 0L)
  default <- function(name, value) if (is.null(m[[name]])) value else m[[name]]
  matrices <- list(
    A = m$A, lambdaX = m$lambdaX, tauX = m$tauX, thetaDelta = m$thetaDelta,
    gammaXi = m$gammaXi, gammaEta = m$gammaEta, alpha = m$alpha,
    beta0 = m$beta0, psi = m$psi, Ieta = diag(numEtas),
    omegaXiXi  = default("omegaXiXi",  matrix(0, numEtas * numXis, numXis)),
    omegaEtaXi = default("omegaEtaXi", matrix(0, numEtas * numXis, numEtas)),
    covZetaXi  = default("covZetaXi",  matrix(0, numEtas, numXis)),
    lambdaY = none, tauY = none, thetaEpsilon = none, W = none, T = none
  )
  list(matrices = matrices,
       info = list(numXis = numXis, numEtas = numEtas, hasComposites = FALSE),
       quad = list(k = k))
}


# Cell probabilities for every pair, either through the node loop (`hoist =
# FALSE`) or from the unconditional moments (`hoist = TRUE`). Both routes go
# through the kernel the fit itself runs.
pmlProbs <- function(sub, thresholds, pairs, m = 30L, hoist = FALSE) {
  rule <- pmlQuadRule(if (hoist) 0L else sub$quad$k, m = m)
  all <- seq_len(NROW(pairs)) - 1L
  pmlProbabilitiesCpp(
    sub, nodes = rule$n, weights = rule$w, thresholdList = thresholds,
    rows = seq_along(thresholds) - 1L, pairs = pairs - 1L,
    hoisted    = if (hoist) all else integer(0L),
    integrated = if (hoist) integer(0L) else all)
}


# A k = 1 model: `Y ~ X + Z + X:Z`, two indicators each, conditioning on xi1.
# Indicators 1-4 measure the xis and are untouched by the interaction;
# indicators 5-6 measure eta and are not.
interactionFixture <- function(omega = 0.4) {
  Phi <- matrix(c(1.0, 0.3, 0.3, 1.2), 2L)
  lambda <- matrix(0, 6L, 3L)
  lambda[1:2, 1L] <- c(1, 0.8)
  lambda[3:4, 2L] <- c(1, 0.9)
  lambda[5:6, 3L] <- c(1, 0.7)
  Oxx <- matrix(0, 2L, 2L)
  Oxx[1L, 2L] <- omega
  pmlSubmodel(list(
    lambdaX = lambda, tauX = matrix(0, 6L, 1L), thetaDelta = diag(6L),
    A = t(chol(Phi)), psi = matrix(0.8, 1L),
    gammaXi = matrix(c(0.5, 0.4), 1L, 2L), gammaEta = matrix(0, 1L, 1L),
    alpha = matrix(0, 1L), beta0 = matrix(0, 2L, 1L),
    omegaXiXi = Oxx), numXis = 2L, numEtas = 1L, k = 1L)
}


testthat::test_that("the bivariate normal CDF matches pbivnorm", {
  grid <- expand.grid(a = seq(-4, 4, by = 0.5), b = seq(-4, 4, by = 0.5))
  # Both sides of the 0.925 switch between Genz's quadrature and his
  # near-singular expansion, and both sides of the 0.3/0.75 rule changes.
  for (rho in c(-0.999, -0.95, -0.8, -0.5, -0.29, 0, 0.29, 0.74, 0.8,
                0.924, 0.93, 0.99, 0.999)) {
    # ABSOLUTE agreement: cells out in the tail are ~1e-8, so a relative
    # tolerance would be testing floating-point noise rather than the routine.
    testthat::expect_lt(
      max(abs(pmlBivariateCpp(grid$a, grid$b, rho) -
                pbivnorm::pbivnorm(grid$a, grid$b, rho = rho))),
      1e-12, label = sprintf("max |diff| at rho = %.3f", rho))
  }
})


testthat::test_that("the k = 0 quadrature rule is a single unit-weight node", {
  rule <- pmlQuadRule(0L)
  testthat::expect_equal(NROW(rule$n), 1L)
  testthat::expect_equal(rule$w, 1)
  # and a k = 1 rule still integrates the standard normal exactly
  r1 <- pmlQuadRule(1L, m = 20L)
  testthat::expect_equal(sum(r1$w), 1, tolerance = 1e-12)
  testthat::expect_equal(sum(r1$w * r1$n[, 1L]), 0, tolerance = 1e-12)
  testthat::expect_equal(sum(r1$w * r1$n[, 1L]^2), 1, tolerance = 1e-12)
})


testthat::test_that("affectedness is downstream of an interaction, not exogeneity", {
  M <- interactionFixture()$matrices
  testthat::expect_equal(pmlAffectedEtas(M, 2L, 1L), TRUE)
  testthat::expect_equal(pmlCleanIndicators(M, 2L, 1L),
                         c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))

  # No omega at all: everything is clean, and the node loop should never run.
  linear <- interactionFixture(omega = 0)$matrices
  testthat::expect_equal(pmlAffectedEtas(linear, 2L, 1L), FALSE)
  testthat::expect_true(all(pmlCleanIndicators(linear, 2L, 1L)))

  # A FREE omega is NA in the unfilled matrices and must count as nonzero --
  # otherwise the split would be built at a zero starting value and be wrong for
  # every iteration after the first.
  unfilled <- linear
  unfilled$omegaXiXi[1L, 2L] <- NA_real_
  testthat::expect_equal(pmlAffectedEtas(unfilled, 2L, 1L), TRUE)

  # Two etas, omega only in eta1's block, eta2 ~ eta1. eta2 inherits the
  # nonlinearity through gammaEta even though its own block is empty.
  Oxx <- matrix(0, 4L, 2L); Oxx[1L, 2L] <- 0.4
  chained <- list(gammaEta = matrix(c(0, 0.5, 0, 0), 2L, 2L),
                  omegaXiXi = Oxx, omegaEtaXi = matrix(0, 4L, 2L))
  testthat::expect_equal(pmlAffectedEtas(chained, 2L, 2L), c(TRUE, TRUE))

  chained$gammaEta[] <- 0
  testthat::expect_equal(pmlAffectedEtas(chained, 2L, 2L), c(TRUE, FALSE))
})


testthat::test_that("pair probabilities are a proper distribution over cells", {
  sub <- interactionFixture()
  pairs <- t(utils::combn(6L, 2L))
  thresholds <- rep(list(c(-0.8, 0, 0.9)), 6L)
  for (P in pmlProbs(sub, thresholds, pairs)) {
    testthat::expect_true(all(P > 0))
    testthat::expect_equal(sum(P), 1, tolerance = 1e-9)
  }
})


testthat::test_that("hoisting a clean pair is exact, not an approximation", {
  sub <- interactionFixture()
  pairs <- t(utils::combn(6L, 2L))
  thresholds <- rep(list(c(-0.7, 0.1, 0.8)), 6L)
  index <- function(j, k) which(pairs[, 1L] == j & pairs[, 2L] == k)

  fine <- pmlProbs(sub, thresholds, pairs, m = 60L)
  coarse <- pmlProbs(sub, thresholds, pairs, m = 5L)
  hoisted <- pmlProbs(sub, thresholds, pairs, hoist = TRUE)

  # A clean pair: the quadrature is CONVERGING ON the closed form, so the
  # hoisted value is the exact answer the node loop can only approach.
  p <- index(1L, 2L)
  testthat::expect_lt(max(abs(fine[[p]] - hoisted[[p]])), 1e-9)
  testthat::expect_gt(max(abs(coarse[[p]] - hoisted[[p]])), 1e-6)

  # A pair on the interaction's outcome is genuinely non-normal, so the
  # unconditional moments are NOT the answer and the loop is required.
  p <- index(5L, 6L)
  testthat::expect_gt(max(abs(fine[[p]] - hoisted[[p]])), 1e-4)
})


testthat::test_that("the k > 0 conditioning reproduces the generating process", {
  sub <- interactionFixture()
  M <- sub$matrices
  tau <- c(-0.7, 0.1, 0.8)
  thresholds <- rep(list(tau), 6L)
  omega <- M$omegaXiXi[1L, 2L]
  pairs <- t(utils::combn(6L, 2L))

  set.seed(1L)
  chunks <- 4L; n <- 5e5L
  N <- chunks * n
  counts <- array(0, dim = c(4L, 4L, NROW(pairs)))
  for (ch in seq_len(chunks)) {
    xi <- matrix(stats::rnorm(2L * n), n, 2L) %*% t(M$A)
    eta <- as.vector(xi %*% t(M$gammaXi)) + omega * xi[, 1L] * xi[, 2L] +
      stats::rnorm(n, sd = sqrt(M$psi[[1L]]))
    xstar <- cbind(xi, eta) %*% t(M$lambdaX) + matrix(stats::rnorm(6L * n), n, 6L)
    code <- apply(xstar, 2L, function(v) findInterval(v, tau) + 1L)
    for (p in seq_len(NROW(pairs)))
      counts[, , p] <- counts[, , p] +
        table(factor(code[, pairs[p, 1L]], levels = 1:4),
              factor(code[, pairs[p, 2L]], levels = 1:4))
  }

  fine <- pmlProbs(sub, thresholds, pairs, m = 60L)
  hoisted <- pmlProbs(sub, thresholds, pairs, hoist = TRUE)
  clean <- pmlCleanIndicators(M, 2L, 1L)

  # Per-cell Monte Carlo SE is ~1.7e-4 here. Measured across all 15 pairs the
  # quadrature sits at 2.0-4.3e-4 (1-3 SE), the hoist on a clean pair is
  # identical to it, and the hoist on a dirty pair is 1.0-2.1e-2. The bounds
  # below are set from those, with the two regimes an order of magnitude apart.
  # At 2e7 draws the quadrature comparison tightens to 1.0-2.2e-4, SE 5.4e-5.
  for (p in seq_len(NROW(pairs))) {
    j <- pairs[p, 1L]; k <- pairs[p, 2L]
    observed <- counts[, , p] / N
    testthat::expect_lt(max(abs(fine[[p]] - observed)), 1e-3,
                        label = sprintf("quadrature error on pair (%d,%d)", j, k))

    # And the split is not over-eager: dropping the integral on a pair marked
    # `integrated` is wrong by more than an order of magnitude above MC noise,
    # while on a `hoisted` pair it changes nothing.
    gap <- max(abs(hoisted[[p]] - observed))
    if (clean[[j]] && clean[[k]]) testthat::expect_lt(gap, 1e-3)
    else                          testthat::expect_gt(gap, 5e-3)
  }
})


testthat::test_that("PML matches lavaan on a linear ordinal model (k = 0)", {
  fx <- ordinalFixture()
  syntax <- '
    X =~ x1 + x2 + x3
    Z =~ z1 + z2 + z3
    Y =~ y1 + y2 + y3
    Y ~ X + Z
  '
  lav <- lavaan::sem(syntax, data = fx$data, ordered = fx$vars,
                     estimator = "PML", parameterization = "theta")
  testthat::expect_true(lavaan::lavInspect(lav, "converged"))

  # Build the same model's matrices from lavaan's estimates, then check the
  # objective agrees with lavaan's own -- which is the saturated pairwise
  # log-likelihood minus fx * N.
  categories <- vapply(fx$data[fx$vars],
                       function(x) length(levels(x)), integer(1L))
  tables <- pmlPairTables(fx$data, fx$vars, categories)
  saturated <- sum(vapply(tables$tables, function(tb) {
    n <- as.vector(tb); p <- n / sum(n); sum(n[n > 0] * log(p[n > 0]))
  }, numeric(1L)))

  pe <- lavaan::parameterEstimates(lav)
  pick <- function(l, o, r) {
    v <- pe$est[pe$lhs == l & pe$op == o & pe$rhs == r]
    if (length(v)) v[[1L]] else 0
  }
  lambda <- matrix(0, 9L, 3L)
  lambda[1:3, 1L] <- c(1, pick("X", "=~", "x2"), pick("X", "=~", "x3"))
  lambda[4:6, 2L] <- c(1, pick("Z", "=~", "z2"), pick("Z", "=~", "z3"))
  lambda[7:9, 3L] <- c(1, pick("Y", "=~", "y2"), pick("Y", "=~", "y3"))
  Phi <- matrix(c(pick("X", "~~", "X"), pick("X", "~~", "Z"),
                  pick("X", "~~", "Z"), pick("Z", "~~", "Z")), 2L)
  sub <- pmlSubmodel(list(
    lambdaX = lambda, tauX = matrix(0, 9L, 1L), thetaDelta = diag(9L),
    A = t(chol(Phi)), psi = matrix(pick("Y", "~~", "Y"), 1L),
    gammaXi = matrix(c(pick("Y", "~", "X"), pick("Y", "~", "Z")), 1L, 2L),
    gammaEta = matrix(0, 1L, 1L), alpha = matrix(0, 1L),
    beta0 = matrix(0, 2L, 1L)), numXis = 2L, numEtas = 1L, k = 0L)
  thresholds <- split(pe$est[pe$op == "|"],
                      rep(seq_len(9L), each = categories[[1L]] - 1L))

  rows <- seq_len(9L)
  partition <- pmlPartition(sub$matrices, 2L, 1L, tables$pairs, rows)
  testthat::expect_length(partition$integrated, 0L)

  rule <- pmlQuadRule(0L)
  objective <- pmlObjectiveCpp(
    sub, nodes = rule$n, weights = rule$w, thresholdList = thresholds,
    rows = rows - 1L, pairs = tables$pairs - 1L, countList = tables$tables,
    hoisted = partition$hoisted - 1L, integrated = partition$integrated - 1L)
  expected <- saturated - lavaan::lavInspect(lav, "optim")$fx * NROW(fx$data)
  testthat::expect_equal(objective, expected, tolerance = 1e-6,
                         info = sprintf("PML %.6f vs lavaan-implied %.6f",
                                        objective, expected))
})
