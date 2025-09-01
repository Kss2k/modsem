# devtools::install_github("kss2k/modsem", ref = "higher-order-rcs")
# library(modsem)
devtools::load_all()
library(lavaan)
library(mvtnorm)


# Simulate data from two higher order variables X and Z
N <- 2000L
Sigma <- diag(c(1.2, 1.8))
Sigma[1, 2] <- Sigma[2, 1] <- 0.6
Mu <- rep(0, 2L)

set.seed(2394827)
XI <- rmvnorm(N, mean = Mu, sigma = Sigma)

disturbance <- \(eps) rnorm(N, mean = 0, sd = sqrt(eps))

# Second Order LVs
X <- XI[, 1]
Z <- XI[, 2]
Y <- 0.2 * X + 0.5 * Z + 0.5 * X * Z + disturbance(0.5)

# First Order LVs
X1 <- 1   * X + disturbance(0.2)
X2 <- 0.8 * X + disturbance(0.2)

Z1 <- 1   * Z + disturbance(0.2)
Z2 <- 0.8 * Z + disturbance(0.2)

# OV Indicators
x1 <- 1   * X1 + disturbance(0.2)
x2 <- 0.8 * X1 + disturbance(0.2)

x3 <- 1   * X2 + disturbance(0.2)
x4 <- 0.8 * X2 + disturbance(0.2)

z1 <- 1   * Z1 + disturbance(0.2)
z2 <- 0.8 * Z1 + disturbance(0.2)

z3 <- 1   * Z2 + disturbance(0.2)
z4 <- 0.8 * Z2 + disturbance(0.2)


y1 <- 1   * Y + disturbance(0.2)
y2 <- 0.8 * Y + disturbance(0.2)

data <- data.frame(x1, x2, x3, x4,
                   z1, z2, z3, z4,
                   y1, y2)


sem.syntax <- '
  LX1 =~ x1 + x2
  LX2 =~ x3 + x4

  LZ1 =~ z1 + z2
  LZ2 =~ z3 + z4

  LX =~ LX1 + LX2
  LZ =~ LZ1 + LZ2

  LY =~ y1 + y2

  LY ~ LX + LZ + LX:LZ
'

fit_lms_rcs <- modsem(sem.syntax, data, method = "lms", rcs = TRUE)
summary(fit_lms_rcs)
# n = 2000
#> Regressions:
#>                  Estimate  Std.Error  z.value  P(>|z|)
#>   LY ~
#>     LX              0.161      0.026    6.077    0.000
#>     LZ              0.487      0.023   21.349    0.000
#>     LX:LZ           0.491      0.018   27.764    0.000

runMplus <- tryCatch({
    MplusAutomation::detectMplus()
    TRUE
  },
  error = function(e) FALSE
)

if (runMplus) {
  fit_mplus <- modsem(sem.syntax, data = data, method = "mplus")
  summary(fit_mplus)
  #> Regressions:
  #>                  Estimate  Std.Error  z.value  P(>|z|)
  #>   LY ~
  #>     LX   (LY<-LX)   0.166      0.027    6.148    0.000
  #>     LZ   (LY<-LZ)   0.481      0.023   20.913    0.000
  #>     LX:LZ  (LY<-)   0.499      0.020   24.950    0.000
}

testthat::expect_warning({
  fit_lms <- modsem(sem.syntax, data = data, method = "lms")
}, regexp = "Unable to optimize.*")
standardized_estimates(fit_lms)

summary(fit_lms)
#> Regressions:
#>                  Estimate  Std.Error  z.value  P(>|z|)
#>   LY ~
#>     LX              0.166      0.027    6.081    0.000
#>     LZ              0.481      0.023   20.582    0.000
#>     LX:LZ           0.498      0.020   24.608    0.000

testthat::expect_error(modsem(sem.syntax, data = data, method = "qml"),
                       regexp = "Higher-order.*qml.*Try.*lms.*")
