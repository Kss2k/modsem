devtools::load_all()
library(lavaan)

m1 <- '
# Outer Model
  X =~ x1 + x2 # + x3
  Z =~ z1 + z2 # + z3
  Y =~ y1 + y2 # + y3

# Inner Model
  Y ~ X + Z # + X:Z
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

    min.x <- min(x)
    max.x <- max(x)

    data[[var]]       <- y
    thresholds[[var]] <- t[t >= min.x & t <= max.x]
  }

  list(data = data, thresholds = thresholds)
}

set.seed(2837290)
ovs <- getOVs(modsemify(m1))
choose <- intersect(colnames(oneInt), ovs)
CUTS <- cut_data(oneInt, choose = choose)
oneInt2 <- CUTS$data
m2 <- '
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z #+ X:Z

  x1 ~~ 0.128*x1
  x2 ~~ 0.209*x2
  x3 ~~ 0.187*x3                          
  z1 ~~ 0.134*z1                          
  z2 ~~ 0.191*z2                          
  z3 ~~ 0.169*z3                          
  y1 ~~ 0.049*y1                          
  y2 ~~ 0.085*y2                          
  y3 ~~ 0.090*y3                          
'
lms2 <- modsem(m1, oneInt2, method = "lms", ordered = choose, estimator = "PML",
               optimize = TRUE, n.threads = 7)
lms3 <- modsem(m1, oneInt2, method = "lms", ordered = choose, estimator = "PML",
               optimize = FALSE, start = lms2$theta, convergence.rel = 1e-17, n.threads = 5)
lms4 <- modsem(m1, oneInt2, method = "lms", ordered = choose, estimator = "PML",
               optimize = FALSE, start = lms2$theta, convergence.rel = 1e-17, n.threads = 5)

reOrder <- \(x) as.ordered(as.integer(x))
oneInt3 <- oneInt2
oneInt3[choose] <- lapply(oneInt2[choose], reOrder)
fit_lav <- sem(m1, oneInt3, estimator = "PML")
#   lms1 <- modsem(m1, oneInt2, method = "lms", ordered = choose,
#                  ordered.iter = 75, ordered.warmup = 20)


CHOOSE <- list(c("x1", "x2", "z1", "y1"),
               colnames(oneInt))

for (choose in CHOOSE) {
  set.seed(2837290)
  CUTS <- cut_data(oneInt, choose = choose)
  oneInt2 <- CUTS$data
  lms1 <- modsem(m1, oneInt2, method = "lms", ordered = choose, estimator = "PML",
                 optimize = TRUE)
#   ,
#                  ordered.iter = 75, ordered.warmup = 20)
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

# Compare internals
tau.x <- c(-1.5, -0.8)
tau.y <- c(-0.9, 0.1)


var.x <- 1.2
var.y <- 1.4
cov.xy <- 0.8
cov.m <- matrix(c(var.x, cov.xy, cov.xy, var.y), nrow=2)
cor.m <- cov2cor(cov.m)
rho <- cor.m[2, 1]
x <- 1
y <- 1

# cov.xy / (sqrt(var.x) * sqrt(var.y))
sd.x <- sqrt(var.x)
sd.y <- sqrt(var.y)
lavaan:::pbinorm(lower.x = tau.x[1] / sd.x, upper.x = tau.x[2] / sd.x,
                 lower.y = tau.y[1] / sd.y, upper.y = tau.y[2] / sd.y,
                 rho = rho)

X <- matrix(c(x, y), nrow = 1, ncol = 2)
TAU <- matrix(c(tau.x, tau.y), nrow = 2, byrow = TRUE)
exp(probPML(X, mu = c(0, 0), Sigma = cov.m, isOrderedEnum = c(1, 2), thresholds = TAU))
