
# devtools::install_github("kss2k/modsem")

library(mvtnorm)
devtools::load_all()
library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cSEM)
library(lme4)

# ──────────────────────────────────────────────────────────────────────────────
# Setup Simulation
# ──────────────────────────────────────────────────────────────────────────────


model.lvl1 <- '
X =~ x1 + x2 + x3 # + x4 + x5
Z =~ z1 + z2 + z3 # + z4 + z5
Y =~ y1 + y2 + y3 # + y4 + y5

Y ~ X + Z # + X:Z
'

model.lme4 <- 'Y ~ X + Z + (1 + X + Z | cluster)'

var_X      <- 1
var_Z      <- 1
cov_X_Z    <- 0.2 

beta_Y     <- 0
gamma_Y_X  <- 0.3
gamma_Y_Z  <- 0.5
gamma_Y_XZ <- 0 # exclude for now
gamma_Y_ZZ <- 0 # exclude for now
gamma_Y_XX <- 0 # exclude for now

var_XZ     <- var_X*var_Z + cov_X_Z^2
zeta_Y     <- 1 - (gamma_Y_X^2 * var_X +
                   gamma_Y_Z^2 * var_Z +
                   2 * gamma_Y_X*gamma_Y_Z * cov_X_Z +
                   gamma_Y_XZ^2 * var_XZ)

var_gamma_Y_X  <- 0.04
var_gamma_Y_Z  <- 0.06
var_gamma_Y_XZ <- 0
var_beta_Y     <- 0.1

lambda_1   <- 1
lambda_2   <- .7
lambda_3   <- .8

epsilon    <- 0.2
beta_1     <- 1.2
beta_2     <- 0.8
beta_3     <- 1.5
n          <- 2500
k          <- 5 # number of categories in cluster

lhs <- "Y"
rhs <- c("X", "Z", "X:Z", "X:X", "Z:Z")
parTable.true <- data.frame(
  lhs = lhs,
  op = "~",
  rhs = rhs,
  par = sprintf("%s~%s", lhs, rhs),
  est = c(gamma_Y_X, gamma_Y_Z, gamma_Y_XZ, gamma_Y_XX, gamma_Y_ZZ),
  se = NA,
  pvalue = NA,
  method = "True value",
  id = 1,
  elapsed = 0
)


sim_data <- function(N = n, K = k) {
  residual <- function(epsilon) rnorm(N, sd = sqrt(epsilon))


  create_ind <- function(lv, beta, lambda, epsilon) {
    beta + lambda * lv + residual(epsilon)
  }


  SXI <- matrix(c(var_X, cov_X_Z,
                  cov_X_Z, var_Z), nrow = 2)
  XI <- rmvnorm(N, sigma = SXI)
  
  X <- XI[, 1]
  Z <- XI[, 2]
  
  Y <- rep(0, N)
    # gamma_Y_X * X + 
    # gamma_Y_Z * Z + 
    # gamma_Y_XZ * X * Z +
    # gamma_Y_XX * X * X +
    # gamma_Y_ZZ * Z * Z 
  
  cluster <- sample(K, n, replace = TRUE)

  for (i in unique(cluster)) {
    cond <- cluster == i
    dbeta <- beta_Y + rnorm(1L, mean = 0, sd = sqrt(var_beta_Y))
    dgammax <- gamma_Y_X + rnorm(1L, mean = 0, sd = sqrt(var_gamma_Y_X))
    dgammaz <- gamma_Y_Z + rnorm(1L, mean = 0, sd = sqrt(var_gamma_Y_Z))
    dgammaxz <- gamma_Y_XZ + rnorm(1L, mean = 0, sd = sqrt(var_gamma_Y_XZ))
    
    Y[cond] <- Y[cond] + dbeta + 
      dgammax * X[cond] + dgammaz * Z[cond] + 
      dgammaxz * (X * Z)[cond]
  }

  print(var(Y))
  zeta_Y <- 1 - var(Y) # empirical solution
  Y <- Y + residual(zeta_Y)

  x1 <- create_ind(X, beta_1, lambda_1, epsilon)
  x2 <- create_ind(X, beta_2, lambda_2, epsilon)
  x3 <- create_ind(X, beta_3, lambda_3, epsilon)
  
  z1 <- create_ind(Z, beta_1, lambda_1, epsilon)
  z2 <- create_ind(Z, beta_2, lambda_2, epsilon)
  z3 <- create_ind(Z, beta_3, lambda_3, epsilon)
  
  y1 <- create_ind(Y, beta_1, lambda_1, epsilon)
  y2 <- create_ind(Y, beta_2, lambda_2, epsilon)
  y3 <- create_ind(Y, beta_3, lambda_3, epsilon)
  
  data <- data.frame(
    x1, x2, x3,
    z1, z2, z3,
    y1, y2, y3,
    cluster
  )

  print(
    lmer('Y ~ X + Z + (1 + X + Z | cluster)', data = data.frame(X, Z, Y, cluster))
  )
  data
}
