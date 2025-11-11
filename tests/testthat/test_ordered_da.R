devtools::load_all()


m1 <- '
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
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

  thresholds <- list()
  for (var in choose) {
    x <- standardize(data[[var]])
    t <- rthreshold(k)
    y <- cut(x, breaks = t, ordered_result = TRUE)
    y <- as.ordered(as.integer(y))

    min.x <- min(x)
    max.x <- max(x)

    data[[var]]       <- y
    thresholds[[var]] <- t[t >= min.x & t <= max.x]
  }

  list(data = data, thresholds = thresholds)
}



CHOOSE <- list(c("y1", "y2", "y3"))

for (choose in CHOOSE) {
  set.seed(2837290)
  CUTS <- cut_data(oneInt, choose = choose)
  oneInt2 <- CUTS$data
  devtools::load_all();lms1 <- modsem(m1, oneInt2, method = "lms", ordered = choose)
  thresholds <- CUTS$thresholds


  thresholds.table <- NULL
  parTable <- parameter_estimates(lms1)
  for (col in choose) {
    tau.true   <- thresholds[[col]]
    tau.true   <- tau.true[is.finite(tau.true)]
    mask       <- parTable$lhs == col & parTable$op == "|"
    tau.est    <- parTable[mask, "est"]
    tau.lower  <- parTable[mask, "ci.lower"]
    tau.upper  <- parTable[mask, "ci.upper"]
    pars <- paste0(col, "|t", seq_along(tau.true))

    rows <- data.frame(parameter = pars, true = tau.true,
                       est = tau.est, diff = tau.true - tau.est,
                       ci.lower = tau.lower, ci.upper = tau.upper,
                       ok = tau.true >= tau.lower & tau.true <= tau.upper)
    thresholds.table <- rbind(thresholds.table, rows)
  }

  print(modsemParTable(thresholds.table))
  testthat::expect_true(sum(thresholds.table$ok) / NROW(thresholds.table) >= 0.95) # 95% confidence
}
