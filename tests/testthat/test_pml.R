# PML against lavaan, at k = 0.
#
# With no interaction there is nothing to condition on, the quadrature rule
# collapses to a single node, and the model is an ordinary linear ordinal SEM --
# which lavaan fits with `estimator = "PML"`. That makes lavaan an INDEPENDENT
# reference for the piece everything else is built on: the map from model
# matrices to bivariate x* moments, the threshold parameterisation, and the
# rectangle probabilities.
#
# It is deliberately the k = 0 case. The conditioning path in pmlLatentMoments()
# is written for k > 0 but is NOT verified here; validating it needs a reference
# that does not exist yet, and that is the next stage.
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


testthat::test_that("pair probabilities are a proper distribution over cells", {
  fx <- ordinalFixture()
  M <- list(
    lambdaX = matrix(c(1, 0.8, 0, 0, 0, 1), 3L, 2L),
    tauX = matrix(0, 3L, 1L), thetaDelta = diag(3L),
    A = matrix(c(1.2, 0, 0.3, 0.9), 2L), psi = matrix(0.7, 1L),
    gammaXi = matrix(0.5, 1L, 1L), gammaEta = matrix(0, 1L, 1L),
    alpha = matrix(0, 1L), beta0 = matrix(0, 1L, 1L),
    covZetaXi = matrix(0, 1L, 1L)
  )
  # a 2-latent measurement model, no eta: build moments directly
  latent <- list(mean = c(0, 0), cov = M$A %*% t(M$A))
  implied <- pmlImpliedMoments(M, latent)
  thresholds <- rep(list(c(-0.8, 0, 0.9)), 3L)
  P <- pmlPairProbabilities(implied, thresholds, 1L, 2L)
  testthat::expect_true(all(P > 0))
  testthat::expect_equal(sum(P), 1, tolerance = 1e-9)
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
  M <- list(
    lambdaX = lambda, tauX = matrix(0, 9L, 1L), thetaDelta = diag(9L),
    A = t(chol(Phi)), psi = matrix(pick("Y", "~~", "Y"), 1L),
    gammaXi = matrix(c(pick("Y", "~", "X"), pick("Y", "~", "Z")), 1L, 2L),
    gammaEta = matrix(0, 1L, 1L), alpha = matrix(0, 1L),
    beta0 = matrix(0, 2L, 1L), covZetaXi = matrix(0, 1L, 2L)
  )
  thresholds <- split(pe$est[pe$op == "|"],
                      rep(seq_len(9L), each = categories[[1L]] - 1L))

  objective <- pmlObjective(M, numXis = 2L, numEtas = 1L,
                            thresholds = thresholds, tables = tables,
                            rule = pmlQuadRule(0L), sign = -1)
  expected <- saturated - lavaan::lavInspect(lav, "optim")$fx * NROW(fx$data)
  testthat::expect_equal(-objective, expected, tolerance = 1e-6,
                         info = sprintf("PML %.6f vs lavaan-implied %.6f",
                                        -objective, expected))
})
