devtools::load_all()
set.seed(29723234)

library(mvtnorm)

n <- 1e5
Sigma <- matrix(c(
  1.2, 0.7, 0.8, 0.2,
  0.7, 1.8, 0.6, 0.3,
  0.8, 0.6, 1.4, 0.6,
  0.2, 0.3, 0.6, 2.0
), nrow = 4)

R <- 10L

S <- 0
center <- \(x) x - mean(x, na.rm = TRUE)

for (i in seq_len(R)) {
  printf("\riter %d/%d", i, R)
  XI <- rmvnorm(n, sigma = Sigma)
  x <- center(XI[, 1])
  z <- center(XI[, 2])
  w <- center(XI[, 3])
  m <- center(XI[, 4])

  xz <- x * z
  xw <- x * w
  zw <- z * w
  xm <- x * m
  zm <- z * m
  wm <- w * m

  xzw <- x * z * w # recurisve centering algorithm
  xzm <- x * z * m
  xwm <- x * w * m
  zwm <- z * w * m
  xzwm <- x * z * w * m

  Si <- cov(data.frame(x, z, w, xz, xw, zw, 
                       xm, zm, wm, xzw, xzm,
                       xwm, zwm, xzwm))

  S <- S + Si
}
cat("\n")

S <- S/R
round(S, 2)
  
X <- XI[, 1]
Z <- XI[, 2]
W <- XI[, 3]

Y <- 1.2 * X + 0.4 * Z + 0.7 * W + 
  0.2 * W * Z +
  0.7 * W * X +
  1.2 * X * Z +
  2.2 * X * Z * W + rnorm(n, sd = sqrt(2))

createInd <- \(x, lambda, epsilon = 0.2) lambda * x + rnorm(n, sd = sqrt(epsilon))


x1 <- createInd(X, 1)
x2 <- createInd(X, 0.8)
x3 <- createInd(X, 0.9)

z1 <- createInd(Z, 1)
z2 <- createInd(Z, 0.8)
z3 <- createInd(Z, 0.9)

w1 <- createInd(W, 1)
w2 <- createInd(W, 0.8)
w3 <- createInd(W, 0.9)

y1 <- createInd(Y, 1)
y2 <- createInd(Y, 0.8)
y3 <- createInd(Y, 0.9)

data <- data.frame(x1, x2, x3, 
                   z1, z2, z3,
                   w1, w2, w3,
                   y1, y2, y3)

model <- '
 X =~ x1 + x2 + x3
 Z =~ z1 + z2 + z3
 W =~ w1 + w2 + w3
 Y =~ y1 + y2 + y3

 Y ~ X + Z + W + X:Z + X:W + Z:W + X:Z:W
 # True values are
 #   Y ~ 1.2 *     X +
 #       0.4 *     Z +
 #       0.7 *     W +
 #       1.2 *   X:Z +
 #       0.7 *   X:W +
 #       0.2 *   Z:W +
 #       2.2 * X:Z:W +
'

est <- modsem(model, data, method = "dblcent")
summary(est, H0 = FALSE)
# n = 10,000
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   Y ~                                                 
#>     X   (true: 1.2)   1.186    0.012  101.886    0.000
#>     Z   (true: 0.4)   0.411    0.008   50.541    0.000
#>     W   (true: 0.7)   0.701    0.010   68.565    0.000
#>     XZ  (true: 1.2)   1.193    0.008  147.196    0.000
#>     XW  (true: 0.7)   0.678    0.008   90.412    0.000
#>     ZW  (true: 0.2)   0.213    0.008   26.954    0.000
#>     XZW (true: 2.2)   2.216    0.004  503.945    0.000
est <- modsem(model, data, method = "dblcent", match = TRUE)
summary(est, H0 = FALSE)
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   Y ~                                                 
#>     X   (true: 1.2)   1.192    0.012   97.960    0.000
#>     Z   (true: 0.4)   0.397    0.008   47.377    0.000
#>     W   (true: 0.7)   0.699    0.011   65.709    0.000
#>     XZ  (true: 1.2)   1.195    0.008  145.091    0.000
#>     XW  (true: 0.7)   0.688    0.008   89.628    0.000
#>     ZW  (true: 0.2)   0.178    0.008   21.835    0.000
#>     XZW (true: 2.2)   2.176    0.004  486.566    0.000


est <- modsem(model, data, method = "dblcent", rcs = TRUE)
summary(est, H0 = FALSE)


#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   Y ~                                                 
#>     X   (true: 1.2)   0.979    0.013   76.464    0.000
#>     Z   (true: 0.4)   0.069    0.009    7.280    0.000
#>     W   (true: 0.7)   0.384    0.012   33.288    0.000
#>     XZ  (true: 1.2)   1.388    0.010  142.664    0.000
#>     XW  (true: 0.7)   0.770    0.009   84.243    0.000
#>     ZW  (true: 0.2)   0.220    0.010   22.208    0.000
#>     XZW (true: 2.2)   3.032    0.006  486.156    0.000
