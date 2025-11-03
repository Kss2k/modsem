devtools::load_all()

m0 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z
'

set.seed(2934829)
oneInt2 <- oneInt
oneInt2$weights <- rnorm(NROW(oneInt2), mean = 1, sd = 0.2)
oneInt2$group   <- sample(2, NROW(oneInt2), replace = TRUE)

testthat::expect_no_condition({
  fit.lms.h0 <- modsem(m0, oneInt2, method = "lms", sampling.weights = "weights")
  fit.qml.h0 <- modsem(m0, oneInt2, method = "qml", sampling.weights = "weights")
  fit.lms.h0.mg1 <- modsem(m0, oneInt2, method = "lms", group = "group", sampling.weights = "weights", sampling.weights.normalization = "total")
  fit.lms.h0.mg2 <- modsem(m0, oneInt2, method = "lms", group = "group", sampling.weights = "weights", sampling.weights.normalization = "group")
})

testthat::expect_true(fit.lms.h0$iterations == 2L)
testthat::expect_true(fit.qml.h0$iterations == 1L)
testthat::expect_true(fit.lms.h0.mg1$iterations == 2L)
testthat::expect_true(fit.lms.h0.mg2$iterations == 2L)

w.g1.1 <- fit.lms.h0.mg1$model$models[[1L]]$data$weights
w.g2.1 <- fit.lms.h0.mg1$model$models[[2L]]$data$weights
w.g1.2 <- fit.lms.h0.mg2$model$models[[1L]]$data$weights
w.g2.2 <- fit.lms.h0.mg2$model$models[[2L]]$data$weights

testthat::expect_true(sum(w.g1.1) + sum(w.g2.1) == length(w.g1.1) + length(w.g2.1))
testthat::expect_true(sum(w.g1.1) != length(w.g1.1))
testthat::expect_true(sum(w.g2.1) != length(w.g2.1))

testthat::expect_true(sum(w.g1.2) + sum(w.g2.2) == length(w.g1.2) + length(w.g2.2))
testthat::expect_true(sum(w.g1.2) == length(w.g1.2))
testthat::expect_true(sum(w.g2.2) == length(w.g2.2))

m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

testthat::expect_no_condition({
  fit.lms.h1 <- modsem(m1, oneInt2, method = "lms", sampling.weights = "weights")
  fit.qml.h1 <- modsem(m1, oneInt2, method = "qml", sampling.weights = "weights")
})


testthat::expect_no_condition({
  bootstrap_modsem(fit.lms.h0, R = 5L, optimize = TRUE, type = "nonparametric")
  bootstrap_modsem(fit.qml.h0, R = 5L, optimize = TRUE, type = "nonparametric")
  bootstrap_modsem(fit.lms.h0, R = 5L, optimize = TRUE, type = "parametric")
  bootstrap_modsem(fit.qml.h0, R = 5L, optimize = TRUE, type = "parametric")
  bootstrap_modsem(fit.lms.h0.mg1, R = 5L, optimize = TRUE, type = "nonparametric")
  bootstrap_modsem(fit.lms.h0.mg1, R = 5L, optimize = TRUE, type = "parametric")
})


# Check missing sampling weights
oneInt3 <- oneInt2
oneInt3$weights[c(2, 3, 5, 18)] <- NA

testthat::expect_error(
  modsem(m1, oneInt3, method = "lms", sampling.weights = "weights"),
  regexp = ".*sampling.weights.*cannot have missing.*"
)

# Check negative sampling weights
oneInt3 <- oneInt2
oneInt3$weights[c(2, 3, 5, 18)] <- c(-1.2, -0.001, -0.2, -3.2)

testthat::expect_error(
  modsem(m1, oneInt3, method = "lms", sampling.weights = "weights"),
  regexp = ".*sampling.weights.*cannot have negative values.*"
)
