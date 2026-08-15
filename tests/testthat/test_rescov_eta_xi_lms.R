devtools::load_all()

testthat::test_that("Residual covariances between xis and etas in the LMS approach", {
  set.seed(4201)
  n <- 2000L
  sigma <- matrix(c(1.00, 0.35, 0.20,
                    0.35, 0.50, 0.40,
                    0.20, 0.40, 0.60), nrow = 3L)
  latent <- matrix(stats::rnorm(n * 3L), ncol = 3L) %*% chol(sigma)
  x <- latent[, 1L]
  y <- x + latent[, 2L]
  z <- y + latent[, 3L]

  make_indicator <- function(v, loading) {
    loading * v + stats::rnorm(n, sd = sqrt(0.25))
  }

  data <- data.frame(
    x1 = make_indicator(x, 1.00),
    x2 = make_indicator(x, 0.85),
    x3 = make_indicator(x, 0.90),
    y1 = make_indicator(y, 1.00),
    y2 = make_indicator(y, 0.85),
    y3 = make_indicator(y, 0.90),
    z1 = make_indicator(z, 1.00),
    z2 = make_indicator(z, 0.85),
    z3 = make_indicator(z, 0.90)
  )

  syntax <- "
    X =~ x1 + x2 + x3
    Y =~ y1 + y2 + y3

    Y ~ 1*X
    X ~~ Y
  "

  fit.lav <- lavaan::sem(syntax, data, meanstructure = TRUE)
  fit.mod <- modsem(syntax, data, method = "lms", optimize = FALSE, convergence.abs = 1e-6)

  ebm.mod <- c(modsem_predict(fit.mod, method = "EBM"))
  ebm.lav <- c(lavaan::lavPredict(fit.lav, method = "EBM"))
  testthat::expect_equal(ebm.mod, ebm.lav, tol = 1e-4)

  ml.mod  <- c(modsem_predict(fit.mod, method = "ML"))
  ml.lav  <- c(lavaan::lavPredict(fit.lav, method = "ML"))
  testthat::expect_equal(ml.mod, ml.lav, tol = 1e-4)

  est.mod <- c(coef(fit.mod))
  est.lav <- lavaan::coef(fit.lav)[names(est.mod)]

  testthat::expect_equal(est.mod, est.lav, tol = 1e-4)

  syntax.int <- "
    X =~ x1 + x2 + x3
    Y =~ y1 + y2 + y3

    Y ~ 1*X + X:X
    X ~~ Y
  "

  y.int <- x + change + 0.2*x^2

  data.int <- data.frame(
    x1 = make_indicator(x, 1.00),
    x2 = make_indicator(x, 0.85),
    x3 = make_indicator(x, 0.90),
    y1 = make_indicator(y.int, 1.00),
    y2 = make_indicator(y.int, 0.85),
    y3 = make_indicator(y.int, 0.90)
  )

  fit.mod.int <- modsem(syntax.int, data, method = "lms")
  est.mod.int <- coef(fit.mod.int)
  testthat::expect_equal(
    est.mod.int[["X~~Y"]],
    est.mod[["X~~Y"]],
    tol = 1e-3
  )

  # Check that residual covariances work with cov.syntax models
  syntax.z.main <- "
    X =~ x1 + x2 + x3
    Y =~ y1 + y2 + y3
    Z =~ z1 + z2 + z3

    Z ~ 1*Y

    X ~~ Y
    Y ~~ Z
  "
  
  syntax.z.cov <- "
    Y ~ 1*X
  "

  fit.lav <- lavaan::sem(paste0(syntax.z.main, syntax.z.cov), data, meanstructure = TRUE)
  fit.mod <- modsem(syntax.z.main, data, method = "lms", cov.syntax = syntax.z.cov, optimize = FALSE, convergence.abs = 1e-8)

  est.mod <- c(coef(fit.mod))
  est.lav <- lavaan::coef(fit.lav)[names(est.mod)]

  testthat::expect_equal(est.mod, est.lav, tol = 1e-4)
})
