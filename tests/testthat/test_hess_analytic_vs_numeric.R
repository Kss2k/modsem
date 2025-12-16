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

  numeric <- modsem:::hessCompLogLikLmsCpp_central(
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

test_that("analytic lambdaX Hessian matches numeric finite differences", {
  inputs <- lms_model_inputs()
  lambda <- inputs$modelR$matrices$lambdaX
  dims <- dim(lambda)
  n_r <- if (length(dims)) dims[1] else length(lambda)
  n_c <- if (length(dims) > 1) dims[2] else 1L
  n_test_r <- min(2L, n_r)
  n_test_c <- min(2L, n_c)
  rows <- as.integer(rep(seq_len(n_test_r) - 1L, times = n_test_c))
  cols <- as.integer(rep(seq_len(n_test_c) - 1L, each = n_test_r))
  symmetric <- rep(0L, length(rows))
  compare_block_hessian(block_id = 0L, rows = rows, cols = cols,
                        symmetric = symmetric)
})

test_that("analytic thetaDelta Hessian matches numeric finite differences", {
  inputs <- lms_model_inputs()
  theta <- inputs$modelR$matrices$thetaDelta
  dims <- dim(theta)
  n_r <- if (length(dims)) dims[1] else length(theta)
  n_c <- if (length(dims) > 1) dims[2] else 1L
  skip_if(n_r < 1L || n_c < 1L, "thetaDelta has insufficient size")
  row_vals <- c(0L, min(1L, n_r - 1L))
  col_vals <- row_vals
  rows <- as.integer(row_vals)
  cols <- as.integer(col_vals)
  symmetric <- rep(1L, length(rows))
  compare_block_hessian(block_id = 4L, rows = rows, cols = cols,
                        symmetric = symmetric)
})

test_that("analytic gammaXi Hessian matches numeric finite differences", {
  inputs <- lms_model_inputs()
  Gx <- inputs$modelR$matrices$gammaXi
  dims <- dim(Gx)
  n_r <- if (length(dims)) dims[1] else length(Gx)
  n_c <- if (length(dims) > 1) dims[2] else 1L
  skip_if(n_r == 0L || n_c == 0L, "gammaXi matrix is empty")
  rows <- as.integer(c(0L, min(1L, n_r - 1L)))
  cols <- as.integer(c(0L, min(1L, n_c - 1L)))
  symmetric <- rep(0L, length(rows))
  compare_block_hessian(block_id = 10L, rows = rows, cols = cols,
                        symmetric = symmetric)
})

test_that("analytic gammaEta Hessian matches numeric finite differences", {
  inputs <- lms_model_inputs()
  Ge <- inputs$modelR$matrices$gammaEta
  dims <- dim(Ge)
  n_r <- if (length(dims)) dims[1] else length(Ge)
  n_c <- if (length(dims) > 1) dims[2] else 1L
  skip_if(n_r == 0L || n_c == 0L, "gammaEta matrix is empty")
  rows <- as.integer(c(0L, min(1L, n_r - 1L)))
  cols <- as.integer(c(0L, min(1L, n_c - 1L)))
  symmetric <- rep(0L, length(rows))
  compare_block_hessian(block_id = 11L, rows = rows, cols = cols,
                        symmetric = symmetric)
})

test_that("analytic beta0 Hessian matches numeric finite differences", {
  inputs <- lms_model_inputs()
  beta0 <- inputs$modelR$matrices$beta0
  dims <- dim(beta0)
  n_r <- if (length(dims)) dims[1] else length(beta0)
  n_c <- if (length(dims) > 1) dims[2] else 1L
  skip_if(n_r == 0L || n_c == 0L, "beta0 matrix is empty")
  n_test_r <- min(2L, n_r)
  n_test_c <- min(2L, n_c)
  rows <- as.integer(rep(seq_len(n_test_r) - 1L, times = n_test_c))
  cols <- as.integer(rep(seq_len(n_test_c) - 1L, each = n_test_r))
  symmetric <- rep(0L, length(rows))
  compare_block_hessian(block_id = 9L, rows = rows, cols = cols,
                        symmetric = symmetric)
})

test_that("numeric A Hessian matches central difference baseline", {
  inputs <- lms_model_inputs()
  A <- inputs$modelR$matrices$A
  dims <- dim(A)
  n_r <- if (length(dims)) dims[1] else length(A)
  n_c <- if (length(dims) > 1) dims[2] else 1L
  skip_if(n_r == 0L || n_c == 0L, "A matrix is empty")
  rows <- as.integer(rep(seq_len(min(2L, n_r)) - 1L, times = min(2L, n_c)))
  cols <- as.integer(rep(seq_len(min(2L, n_c)) - 1L, each = min(2L, n_r)))
  symmetric <- rep(0L, length(rows))
  compare_block_hessian(block_id = 6L, rows = rows, cols = cols,
                        symmetric = symmetric)
})

test_that("numeric Psi Hessian matches central difference baseline", {
  inputs <- lms_model_inputs()
  Psi <- inputs$modelR$matrices$psi
  dims <- dim(Psi)
  n_r <- if (length(dims)) dims[1] else length(Psi)
  skip_if(n_r == 0L, "Psi matrix is empty")
  rows <- as.integer(c(0L, min(1L, n_r - 1L)))
  cols <- rows
  symmetric <- rep(1L, length(rows))
  compare_block_hessian(block_id = 7L, rows = rows, cols = cols,
                        symmetric = symmetric)
})
