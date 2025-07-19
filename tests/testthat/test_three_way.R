devtools::load_all()
set.seed(29723234)

library(mvtnorm)

n <- 1e4
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

modsem:::timeExpr(
est <- modsem(model, data, method = "dblcent")
)
summary(est, H0 = FALSE)
#n = 100,000
#> Regressions:
#>                     Estimate  Std.Err  z-value  P(>|z|)
#>  Y ~                                                 
#>    X     (true: 1.2)   1.199    0.012  103.820    0.000
#>    Z     (true: 0.4)   0.383    0.008   46.859    0.000
#>    W     (true: 0.7)   0.699    0.010   68.139    0.000
#>    XZ    (true: 1.2)   1.181    0.008  145.896    0.000
#>    XW    (true: 0.7)   0.710    0.007   95.378    0.000
#>    ZW    (true: 0.2)   0.199    0.008   25.132    0.000
#>    XZW   (true: 2.2)   2.198    0.005  486.690    0.000

est <- modsem(model, data, method = "dblcent", match = TRUE)
summary(est, H0 = FALSE)
#n = 100,000
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
     
est.rcs.1 <- modsem(model, data, method = "dblcent", rcs = TRUE,
                    rcs.scale.corrected = FALSE)
summary(est.rcs.1, H0 = FALSE)
#n = 100,000
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   Y ~                                                 
#>     X   (true: 1.2)   0.987    0.013   77.095    0.000
#>     Z   (true: 0.4)   0.010    0.010    1.058    0.290
#>     W   (true: 0.7)   0.349    0.012   29.759    0.000
#>     XZ  (true: 1.2)   1.340    0.010  134.033    0.000
#>     XW  (true: 0.7)   0.819    0.009   88.521    0.000
#>     ZW  (true: 0.2)   0.226    0.010   22.749    0.000
#>     XZW (true: 2.2)   3.076    0.007  458.627    0.000

est.rcs.2 <- modsem(model, data, method = "dblcent", rcs = TRUE,
                    rcs.scale.corrected = TRUE)
summary(est.rcs.2, H0 = FALSE)
#n = 100,000
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   Y ~                                                 
#>     X   (true: 1.2)   1.135    0.011   98.696    0.000
#>     Z   (true: 0.4)   0.253    0.008   29.929    0.000
#>     W   (true: 0.7)   0.570    0.010   54.762    0.000
#>     XZ  (true: 1.2)   1.163    0.008  149.563    0.000
#>     XW  (true: 0.7)   0.702    0.007   96.105    0.000
#>     ZW  (true: 0.2)   0.211    0.008   27.286    0.000