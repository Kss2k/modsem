library(testthat)
devtools::load_all()

compare_block_hessian <- function(block_id, rows, cols, symmetric,
                                  tol_grad = 2e-2, tol_hess = 2e-2,
                                  rel_step = 1e-4, min_abs = 1e-2) {
  inputs <- lms_model_inputs()
  block <- rep(as.integer(block_id), length(rows))

  analytic <- modsem:::hessCompLogLikLmsCpp(
    inputs$modelR, inputs$P,
    block = block, row = rows, col = cols, symmetric = symmetric,
    colidxR = inputs$colidx, n = inputs$n, d = inputs$d,
    npatterns = 1L, relStep = rel_step, minAbs = min_abs, ncores = 1L
  )

  numeric <- modsem:::hessCompLogLikLmsCpp_fd(
    inputs$modelR, inputs$P,
    block = block, row = rows, col = cols, symmetric = symmetric,
    colidxR = inputs$colidx, n = inputs$n, d = inputs$d,
    npatterns = 1L, relStep = rel_step, minAbs = min_abs, ncores = 1L
  )

  expect_equal(as.numeric(analytic$gradient),
               as.numeric(numeric$gradient),
               tolerance = tol_grad)
  expect_equal(as.numeric(analytic$Hessian),
               as.numeric(numeric$Hessian),
               tolerance = tol_hess)
}

test_that("analytic alpha Hessian matches numeric finite differences", {
  inputs <- lms_model_inputs()
  alpha_mat <- inputs$modelR$matrices$alpha
  dims <- dim(alpha_mat)
  n_r <- if (length(dims)) dims[1] else length(alpha_mat)
  n_c <- if (length(dims) > 1) dims[2] else 1L
  rows <- as.integer(rep(seq_len(n_r) - 1L, times = n_c))
  cols <- as.integer(rep(seq_len(n_c) - 1L, each = n_r))
  symmetric <- rep(0L, length(rows))
  compare_block_hessian(block_id = 8L, rows = rows, cols = cols,
                        symmetric = symmetric)
})

test_that("analytic tauX Hessian matches numeric finite differences", {
  rows <- as.integer(c(0L, 1L))
  cols <- as.integer(c(0L, 0L))
  symmetric <- rep(0L, length(rows))
  compare_block_hessian(block_id = 2L, rows = rows, cols = cols,
                        symmetric = symmetric)
})
