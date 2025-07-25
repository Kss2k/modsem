devtools::load_all()
library(mvtnorm)


cov_xz <- 0.8
cov_mx <- 0.7
cov_mz <- 0.4
var_x <- 1.2
var_z <- 1.8
var_m <- 0.9

sigma <- matrix(
 c( var_x, cov_xz, cov_mx,
   cov_xz,  var_z, cov_mz,
   cov_mx, cov_mz,  var_m),
  nrow = 3, ncol = 3
)

set.seed(139427)
N <- 1e7
XI <- rmvnorm(N, mean = rep(0, 3), sigma = sigma)

XI
X <- XI[, 1]
Z <- XI[, 2]
M <- XI[, 3]

XX <- X * X
XZ <- X * Z
ZZ <- Z * Z
MX <- M * X
MZ <- M * Z

pt.cov <- data.frame(
  lhs = c("X", "X", "X", "Z", "Z", "M"),
  op = "~~",
  rhs = c("X", "Z", "M", "Z", "M", "M"),
  est = c(var_x, cov_xz, cov_mx, var_z, cov_mz, var_m)
)
names <- c("X", "Z", "M", "X:X", "X:Z", "Z:Z", "M:X", "M:Z")

pt.int <- data.frame(
  lhs = "Y",
  op = "~",
  rhs = names,
  est = 0
)

pt <- rbind(pt.cov, pt.int)
pt.var <- var_interactions(pt)


INT <- data.frame(X, Z, M, XX, XZ, ZZ, MX, MZ)
colnames(INT) <- names
observed <- cov(INT)
expected <- observed
expected[TRUE] <- NA
for (x in colnames(observed)) {
  for (y in rownames(observed)) {
    match <- pt.var[pt.var$op == "~~" & pt.var$lhs == x & pt.var$rhs == y, "est"]

    if (length(match))
      expected[y, x] <- expected[x, y] <- match
  }
}

testthat::expect_equal(round(expected, 2), round(observed, 2))
