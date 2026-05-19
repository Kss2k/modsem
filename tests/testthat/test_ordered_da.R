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
  choose <- c("x1", "x2", "z1", "y1")

  set.seed(2837290)
  data <- cut_data(oneInt, choose = choose)

  fit <- NULL
  testthat::expect_no_condition({
    fit <- suppressMessages(
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
    )
  })

  par_table <- parameter_estimates(fit)
  int_est <- par_table[
    par_table$lhs == "Y" &
      par_table$op == "~" &
      par_table$rhs == "X:Z",
    "est"
  ]

  testthat::expect_length(int_est, 1L)
  testthat::expect_gte(int_est, 0.440)
  testthat::expect_lte(int_est, 0.462)
})
