# PML against lavaan, at k = 0.
#
# With no interaction there is nothing to condition on, the quadrature rule
# collapses to a single node, and the model is an ordinary linear ordinal SEM --
# which lavaan fits with `estimator = "PML"`. That makes lavaan an INDEPENDENT
# reference for the piece everything else is built on: the map from model
# matrices to bivariate x* moments, the threshold parameterisation, and the
# rectangle probabilities.
#
# It is deliberately the k = 0 case. The conditioning path is written for k > 0
# but is NOT verified against an external reference here; that needs a reference
# which does not exist yet, and is the next stage.
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


# The moments come from LMSModel, so a fixture has to be a whole submodel, not
# just the handful of matrices the arithmetic touches.
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
  implied <- pmlSelectMoments(pmlUnconditionalMomentsCpp(sub), seq_len(6L))
  thresholds <- rep(list(c(-0.8, 0, 0.9)), 6L)
  P <- pmlPairProbabilities(implied, thresholds, 1L, 3L)
  testthat::expect_true(all(P > 0))
  testthat::expect_equal(sum(P), 1, tolerance = 1e-9)
})


testthat::test_that("hoisting a clean pair is exact, not an approximation", {
  sub <- interactionFixture()
  rows <- seq_len(6L)
  thresholds <- rep(list(c(-0.7, 0.1, 0.8)), 6L)

  mixture <- function(j, k, m) {
    rule <- pmlQuadRule(1L, m = m)
    moments <- pmlMomentsCpp(sub, rule$n)
    out <- 0
    for (q in seq_along(moments))
      out <- out + rule$w[[q]] * pmlPairProbabilities(
        pmlSelectMoments(moments[[q]], rows), thresholds, j, k)
    out
  }
  hoisted <- function(j, k)
    pmlPairProbabilities(pmlSelectMoments(pmlUnconditionalMomentsCpp(sub), rows),
                         thresholds, j, k)

  # A clean pair: the quadrature is CONVERGING ON the closed form, so the
  # hoisted value is the exact answer the node loop can only approach.
  testthat::expect_lt(max(abs(mixture(1L, 2L, 60L) - hoisted(1L, 2L))), 1e-9)
  testthat::expect_gt(max(abs(mixture(1L, 2L, 5L) - hoisted(1L, 2L))), 1e-6)

  # A pair on the interaction's outcome is genuinely non-normal, so the
  # unconditional moments are NOT the answer and the loop is required.
  testthat::expect_gt(max(abs(mixture(5L, 6L, 60L) - hoisted(5L, 6L))), 1e-4)
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

  objective <- pmlObjective(sub, thresholds = thresholds, tables = tables,
                            rule = pmlQuadRule(0L), partition = partition,
                            rows = rows, sign = -1)
  expected <- saturated - lavaan::lavInspect(lav, "optim")$fx * NROW(fx$data)
  testthat::expect_equal(-objective, expected, tolerance = 1e-6,
                         info = sprintf("PML %.6f vs lavaan-implied %.6f",
                                        -objective, expected))
})
