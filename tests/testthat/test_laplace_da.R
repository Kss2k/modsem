devtools::load_all()

testthat::test_that("Laplace DA estimator returns finite continuous-model estimates", {
  syntax <- "
    X =~ x1 + x2 + x3
    Y =~ y1 + y2 + y3
    Z =~ z1 + z2 + z3
    Y ~ X + Z + X:Z
  "

  fit <- modsem(
    syntax,
    data = oneInt[1:200, ],
    method = "laplace",
    calc.se = FALSE
  )

  testthat::expect_s3_class(fit, "modsem_da")
  testthat::expect_identical(fit$method, "laplace")
  testthat::expect_true(is.finite(fit$logLik))
  testthat::expect_true(all(is.finite(coef(fit, type = "free"))))
})


testthat::test_that("Laplace DA cold starts retain posterior variance", {
  syntax <- "
    X =~ x1 + x2 + x3
    Y =~ y1 + y2 + y3
    Z =~ z1 + z2 + z3
    Y ~ X + Z + X:Z
  "

  fit <- modsem(
    syntax,
    data = oneInt[1:80, ],
    method = "laplace",
    calc.se = FALSE
  )

  latent.vars <- coef(fit)[c("X~~X", "Y~~Y", "Z~~Z")]
  testthat::expect_true(is.finite(fit$logLik))
  testthat::expect_true(all(is.finite(latent.vars)))
  testthat::expect_true(all(latent.vars > 0.3 & latent.vars < 2.5))
})


testthat::test_that("Laplace DA handles ordered indicators with probit theta identification", {
  syntax <- "
    X =~ x1 + x2 + x3
    Y =~ y1 + y2 + y3
    Z =~ z1 + z2 + z3
    Y ~ X + Z + X:Z
  "
  ordered.data <- oneInt[1:80, ]
  ordered.vars <- colnames(ordered.data)
  ordered.data[ordered.vars] <- lapply(ordered.data[ordered.vars], function(x) {
    ordered(cut(x, quantile(x, probs = seq(0, 1, length.out = 5)),
                include.lowest = TRUE))
  })

  fit <- modsem(
    syntax,
    data = ordered.data,
    method = "laplace",
    ordered = ordered.vars,
    calc.se = FALSE
  )

  matrices <- fit$start.model$models[[1L]]$matrices
  testthat::expect_true(is.finite(fit$logLik))
  testthat::expect_true(all(diag(matrices$thetaDelta) == 1))
  testthat::expect_true(all(drop(matrices$tauX) == 0))
  testthat::expect_equal(sum(fit$parTable$op == "|"), 27L)
})


testthat::test_that("Laplace DA combines continuous and ordered indicator likelihoods", {
  syntax <- "
    X =~ x1 + x2
    Y =~ y1 + y2
    Y ~ X
  "
  mixed.data <- oneInt[1:60, c("x1", "x2", "y1", "y2")]
  ordered.vars <- c("x1", "y1")
  mixed.data[ordered.vars] <- lapply(mixed.data[ordered.vars], function(x) {
    ordered(cut(x, quantile(x, probs = seq(0, 1, length.out = 4)),
                include.lowest = TRUE))
  })

  fit <- modsem(
    syntax,
    data = mixed.data,
    method = "laplace",
    ordered = ordered.vars,
    calc.se = FALSE
  )

  testthat::expect_true(is.finite(fit$logLik))
  testthat::expect_true(all(is.finite(coef(fit, type = "free"))))
})


testthat::test_that("Laplace DA rejects residual covariances involving ordered indicators", {
  syntax <- "
    X =~ x1 + x2
    Y =~ y1 + y2
    Y ~ X
    x1 ~~ x2
  "
  data <- oneInt[1:30, c("x1", "x2", "y1", "y2")]
  data$x1 <- ordered(cut(
    data$x1,
    quantile(data$x1, probs = seq(0, 1, length.out = 4)),
    include.lowest = TRUE
  ))

  testthat::expect_error(
    modsem(syntax, data, method = "laplace", ordered = "x1", calc.se = FALSE),
    regexp = "Residual covariances involving ordered indicators"
  )
})
