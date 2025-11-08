library(testthat)
devtools::load_all()

grad_inputs <- lms_model_inputs


compare_gradients <- function(block_id, row, col, symmetric = 0L, tol = 2e-2) {
  inputs <- grad_inputs()

  gn <- \() modsem:::gradLogLikLmsCpp_cd(
    inputs$modelR, inputs$P,
    block = as.integer(block_id), row = as.integer(row),
    col = as.integer(col), symmetric = as.integer(symmetric),
    colidxR = inputs$colidx, n = inputs$n, d = inputs$d,
    npatterns = 1L, eps = 1e-7, ncores = 1L
  )

  ga <- \() modsem:::gradLogLikLmsCpp(
    inputs$modelR, inputs$P,
    block = as.integer(block_id), row = as.integer(row),
    col = as.integer(col), symmetric = as.integer(symmetric),
    colidxR = inputs$colidx, n = inputs$n, d = inputs$d,
    npatterns = 1L, eps = 1e-7, ncores = 1L
  )

  grad_numeric <- gn()
  grad_analytic <- ga()

  expect_equal(as.numeric(grad_analytic),
               as.numeric(grad_numeric),
               tolerance = tol)
}


test_that("analytic alpha gradient matches numeric finite differences", {
  compare_gradients(block_id = 8, row = 0, col = 0)
})


test_that("analytic beta0 gradient matches numeric finite differences", {
  compare_gradients(block_id = 9, row = 0, col = 0)
})


test_that("analytic tauX gradient matches numeric finite differences", {
  compare_gradients(block_id = 2, row = 0, col = 0)
})


test_that("analytic lambdaX gradient matches numeric finite differences", {
  compare_gradients(block_id = 0, row = 0, col = 0)
})


test_that("analytic A gradient matches numeric finite differences", {
  compare_gradients(block_id = 6, row = 0, col = 0)
})

test_that("analytic A off-diagonal gradient matches numeric finite differences", {
  compare_gradients(block_id = 6, row = 0, col = 1)
})


test_that("analytic psi gradient matches numeric finite differences", {
  compare_gradients(block_id = 7, row = 0, col = 0, symmetric = 1L)
})

test_that("analytic psi off-diagonal gradient matches numeric finite differences", {
  compare_gradients(block_id = 7, row = 0, col = 1, symmetric = 1L)
})


test_that("analytic gammaXi gradient matches numeric finite differences", {
  compare_gradients(block_id = 10, row = 0, col = 0)
})

test_that("analytic gammaEta gradient matches numeric finite differences", {
  compare_gradients(block_id = 11, row = 0, col = 0)
})

test_that("analytic omegaXiXi gradient matches numeric finite differences", {
  compare_gradients(block_id = 12, row = 0, col = 0)
})

test_that("analytic omegaEtaXi gradient matches numeric finite differences", {
  compare_gradients(block_id = 13, row = 0, col = 0)
})


test_that("analytic thetaDelta gradient matches numeric finite differences", {
  compare_gradients(block_id = 4, row = 0, col = 0, symmetric = 1L)
})

test_that("analytic thetaDelta off-diagonal gradient matches numeric finite differences", {
  compare_gradients(block_id = 4, row = 0, col = 1, symmetric = 1L)
})
