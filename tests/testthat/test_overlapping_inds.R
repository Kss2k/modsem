devtools::load_all()

library(mvtnorm)
set.seed(29347098)

var_X      <- 1
var_Z      <- 1
cov_X_Z    <- 0.2 

gamma_Y_X  <- 0.3
gamma_Y_Z  <- 0.5
gamma_Y_XZ <- 0.3

zeta_Y     <- 0.8
lambda_1   <- 1
lambda_2   <- .7
lambda_3   <- .8

epsilon    <- 0.2
beta_1     <- 1.2
beta_2     <- 0.8
beta_3     <- 1.5
N          <- 2000

residual <- \(epsilon) rnorm(N, sd = sqrt(epsilon))
create_ind <- \(lv, beta, lambda, epsilon) beta + lambda * lv + residual(epsilon)

SXI <- matrix(c(var_X, cov_X_Z,
                cov_X_Z, var_Z), nrow = 2)
XI <- rmvnorm(N, sigma = SXI)

X <- XI[, 1]
Z <- XI[, 2]

Y <- 
  gamma_Y_X * X + 
  gamma_Y_Z * Z + 
  gamma_Y_XZ * X * Z +
  residual(zeta_Y)

x1 <- create_ind(X, beta_1, lambda_1, epsilon)
x2 <- create_ind(X, beta_2, lambda_2, epsilon)
x3 <- create_ind(X, beta_3, lambda_3, epsilon)

z1 <- create_ind(Z, beta_1, lambda_1, epsilon)
z2 <- create_ind(Z, beta_2, lambda_2, epsilon)
z3 <- create_ind(Z, beta_3, lambda_3, epsilon)

y1 <- create_ind(Y, beta_1, lambda_1, epsilon)
y2 <- beta_1 + lambda_2 * Y + lambda_3 * X + residual(epsilon)
y3 <- beta_2 + lambda_3 * Y + lambda_2 * X + residual(epsilon)

data <- data.frame(
   x1, x2, x3,
   z1, z2, z3,
   y1, y2, y3
)

model <- '
  X =~ x1 + x2 + x3 + y2 + y3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

testthat::expect_no_condition({
  modsem(model, data = data, method = "lms")
  summary(fit)
})

testthat::expect_error({
  modsem(model, data = data, method = "qml")
}, regexp = "The same indicator cannot .*")
