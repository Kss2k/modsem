devtools::load_all()


m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'


rthreshold <- \(k, offset = runif(1, min = -1, max = 1), sigma = 0.35) {
  t <- seq_len(k) - mean(seq_len(k)) + offset
  t <- t + runif(k, min = -sigma, max = sigma)
  c(-Inf, t, Inf)
}


cut_data <- function(data, k = 5, choose = NULL) {
  if (is.null(choose))
    choose <- colnames(data)

  standardize <- \(x) (x - mean(x)) / sd(x)

  for (var in choose) {
    x <- standardize(data[[var]])
    t <- rthreshold(k)
    data[[var]] <- cut(x, breaks = t, ordered_result = TRUE)
  }

  data
}


testthat::test_that("ordered LMS fits without conditions and returns stable interaction estimate", {
  choose <- colnames(oneInt)

  set.seed(2837290)
  data <- cut_data(oneInt, choose = choose)

  fit <- NULL
  testthat::expect_no_error({
    fit <- suppressWarnings(suppressMessages(
      modsem(
        m1,
        data,
        method = "lms",
        ordered = choose,
        mean.observed = FALSE,
        ordered.fixed.seed = FALSE,
        ordered.polyak.juditsky = TRUE,
        verbose = FALSE
      )
    ))
  })

  # Standardized, not `parameter_estimates()`. These bounds are on the
  # standardized scale: they were set when ordered models went through a
  # different estimator whose latent variances sat near 1, so the two scales
  # coincided. The graph backend leaves the latents in the y* metric the link
  # implies (variances ~17), where the interaction carries sd(X)sd(Z)/sd(Y)
  # rather than sd(X)/sd(Y) -- a factor of ~2.7 on this model. The standardized
  # estimate is the scale-free quantity the test is actually about.
  par_table <- standardized_estimates(fit)
  int_est <- par_table[
    par_table$lhs == "Y" &
      par_table$op == "~" &
      par_table$rhs == "X:Z",
    "est"
  ]

  testthat::expect_length(int_est, 1L)
  testthat::expect_gte(int_est, 0.440)
  testthat::expect_lte(int_est, 0.460)
})
