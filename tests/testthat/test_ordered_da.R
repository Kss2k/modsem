devtools::load_all()
library(lavaan)

m1 <- '
# Outer Model
  X =~ x1 + x2 # + x3
  Z =~ z1 + z2 # + z3
  Y =~ y1 + y2 # + y3

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
lms2 <- modsem(m1, oneInt2, method = "lms", ordered = choose, estimator = "PML",
               optimize = TRUE, n.threads = 7, robust.se = TRUE,
               adaptive.quad = TRUE, nodes = 16, adaptive.frequency = 20)

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
TAU.X <- c(-100, -1.5, -0.8)
TAU.Y <- c(-100, -0.9, 0.1)


var.x <- 1.2
var.y <- 1.4
cov.xy <- 0.8
cov.m <- matrix(c(var.x, cov.xy, cov.xy, var.y), nrow=2)
cor.m <- cov2cor(cov.m)
rho <- cor.m[2, 1]
x <- 2
y <- 2

tau.x <- TAU.X[x:(x+1L)]
tau.y <- TAU.Y[y:(y+1L)]

# cov.xy / (sqrt(var.x) * sqrt(var.y))
sd.x <- sqrt(var.x)
sd.y <- sqrt(var.y)
lavaan:::pbinorm(lower.x = tau.x[1] / sd.x, upper.x = tau.x[2] / sd.x,
                 lower.y = tau.y[1] / sd.y, upper.y = tau.y[2] / sd.y,
                 rho = rho)

X <- matrix(c(1, 1), nrow = 1, ncol = 2)
TAU <- matrix(c(tau.x, tau.y), nrow = 2, byrow = TRUE)

c<-1.4
exp(probPML(X, mu = c(0, 0)+c, Sigma = cov.m, isOrderedEnum = c(1, 2), thresholds = TAU))
exp(probPML_Fast(X, mu = c(0, 0), Sigma = cov.m, isOrderedEnum = c(1, 2), thresholds = TAU))


