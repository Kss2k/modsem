library(testthat)
devtools::load_all()

tpb_model <- "
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC + PBC:SN

  BEH ~ 1
  ATT ~ 1
  INT ~ 1
  PBC ~ 1
  SN  ~ 1
"

grad_inputs <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) return(cache)

    data("TPB", package = "modsem", envir = environment())

    fit <- suppressWarnings(
      suppressMessages(
        modsem(
          tpb_model,
          data = TPB,
          method = "lms",
          nodes = 8,
          optimize = FALSE,
          standardize.data = TRUE,
          mean.observed = FALSE,
          verbose = FALSE
        )
      )
    )

    submodel <- fit$model$models[[1L]]
    P_group  <- fit$estep$P_GROUPS[[1L]]
    dataR    <- submodel$data

    cache <<- list(
      modelR = submodel,
      P      = P_group,
      colidx = dataR$colidx0,
      n      = as.integer(dataR$n.pattern),
      d      = as.integer(dataR$d.pattern)
    )

    cache
  }
})


compare_gradients <- function(block_id, row, col, tol = 2e-2) {
  inputs <- grad_inputs()

  grad_numeric <- modsem:::gradLogLikLmsCpp_cd(
    inputs$modelR, inputs$P,
    block = as.integer(block_id), row = as.integer(row),
    col = as.integer(col), symmetric = as.integer(0),
    colidxR = inputs$colidx, n = inputs$n, d = inputs$d,
    npatterns = 1L, eps = 1e-7, ncores = 1L
  )

  grad_analytic <- modsem:::gradLogLikLmsCpp(
    inputs$modelR, inputs$P,
    block = as.integer(block_id), row = as.integer(row),
    col = as.integer(col), symmetric = as.integer(0),
    colidxR = inputs$colidx, n = inputs$n, d = inputs$d,
    npatterns = 1L, eps = 1e-7, ncores = 1L
  )

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


test_that("analytic psi gradient matches numeric finite differences", {
  compare_gradients(block_id = 7, row = 0, col = 0)
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
  compare_gradients(block_id = 4, row = 0, col = 0)
})
