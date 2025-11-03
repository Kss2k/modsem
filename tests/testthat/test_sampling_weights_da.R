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

testthat::expect_no_condition({
  fit.lms.h0 <- modsem(m0, oneInt2, method = "lms", sampling.weights = "weights")
  fit.qml.h0 <- modsem(m0, oneInt2, method = "qml", sampling.weights = "weights")
})

testthat::expect_true(fit.lms.h0$iterations == 2L)
testthat::expect_true(fit.qml.h0$iterations == 1L)


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
})
