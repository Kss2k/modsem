devtools::load_all()

# Pruning nodes by their weights before evaluating the integrand should not
# change the integral, also when an outlier has its posterior mass beyond the
# weight cutoff (see `prunedQuadrature()`)

set.seed(2938)
s <- 0.2
y <- rnorm(200, mean = 0, sd = sqrt(1 + s^2))

counter <- new.env()
densityFun <- function(y) {
  function(z, ...) {
    counter$n <- counter$n + NROW(z)
    outer(y, z[, 1], FUN = \(yy, zz) dnorm(yy, mean = zz, sd = s))
  }
}

integral1D <- function(fun, m, wtol) {
  q <- quadrature(m = m, k = 1)
  counter$n <- 0
  pruned <- prunedQuadrature(fun = fun, quadw = q$w, quadn = q$n,
                             a = -Inf, b = Inf, tol = 1e-12, wtol = wtol)
  list(I = pruned$I.cur, evaluated = counter$n, total = length(q$w))
}

wtol <- QUAD_WEIGHT_TOL

for (m in c(32, 100, 300)) {
  # No outlier: tails are trimmed well inside the cutoff, so the dropped
  # nodes should never be evaluated
  fun  <- densityFun(y)
  full <- integral1D(fun, m = m, wtol = NULL)
  wcut <- integral1D(fun, m = m, wtol = wtol)

  testthat::expect_equal(wcut$I, full$I, tolerance = 1e-10)
  testthat::expect_lt(wcut$evaluated, wcut$total)

  # Outlier beyond the weight cutoff (posterior mode at ~11.5, cutoff at
  # ~9-10). GH(32) is a coarse rule, where neighbouring weights in the tails
  # differ by more than `edge.ratio`
  fun  <- densityFun(c(y, 12))
  full <- integral1D(fun, m = m, wtol = NULL)
  wcut <- integral1D(fun, m = m, wtol = wtol)

  testthat::expect_equal(wcut$I, full$I, tolerance = 1e-10)
  testthat::expect_equal(wcut$evaluated, wcut$total)
}

# Adaptive quadrature, single dimension
for (y.out in c(10, 12)) {
  fun  <- densityFun(c(y, y.out))
  full <- suppressWarnings(adaptiveGaussQuadrature(fun, a = -Inf, b = Inf, m = 32, k = 1,
                                                   m.ceil = NULL, wtol = NULL))
  wcut <- suppressWarnings(adaptiveGaussQuadrature(fun, a = -Inf, b = Inf, m = 32, k = 1,
                                                   m.ceil = NULL, wtol = wtol))

  testthat::expect_equal(wcut$integral, full$integral, tolerance = 1e-10)
  testthat::expect_equal(wcut$m.ceil, full$m.ceil)
}

# Adaptive quadrature, product grid (two dimensions)
s2 <- 0.3
Y  <- rbind(matrix(rnorm(400, mean = 0, sd = sqrt(1 + s2^2)), ncol = 2), c(11, 0.5))
fun2D <- function(z, ...) {
  vapply(seq_len(NROW(z)), FUN.VALUE = numeric(NROW(Y)), FUN = \(j)
         dnorm(Y[, 1], mean = z[j, 1], sd = s2) * dnorm(Y[, 2], mean = z[j, 2], sd = s2))
}

full <- suppressWarnings(adaptiveGaussQuadrature(fun2D, a = -Inf, b = Inf, m = 16, k = 2,
                                                 m.ceil = NULL, wtol = NULL))
wcut <- suppressWarnings(adaptiveGaussQuadrature(fun2D, a = -Inf, b = Inf, m = 16, k = 2,
                                                 m.ceil = NULL, wtol = wtol))

testthat::expect_equal(wcut$integral, full$integral, tolerance = 1e-10)
testthat::expect_equal(NROW(wcut$n), NROW(full$n))
