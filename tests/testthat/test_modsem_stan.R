devtools::load_all()

library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE, threads_per_chain = 4)      # cache compiled models
options(mc.cores = parallel::detectCores()) 

# Here we estimate the models using STAN.
# We're essentially specifying the models as
# Bayesian models with flat priors for the model
# parameters, mimicking ML-estimation


# ------------------------------------------------------------------------------
# Two-way Interaction
# ------------------------------------------------------------------------------

m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

# First we compile the STAN model, this can be slow
# and is therefore done once, such that the compiled
# STAN code can be reused for the same model syntax later
compiled_model.2way <- compile_stan_model(m1)
cat(compiled_model$syntax)

# Fit the model based on the compiled STAN code
fit.2way <- modsem_stan(compiled_model = compiled_model.2way, 
                        data = oneInt, iter = 2000, chains = 2)
# We can get a summary
summary(fit.2way)

# We can also get the standardized (and unstandardized)
# estimates
standardized_estimates(fit.2way)
parameter_estimates(fit.2way)


# ------------------------------------------------------------------------------
# Three-way Interaction
# ------------------------------------------------------------------------------

# Simulate data

set.seed(29723234)
n <- 500
Sigma <- matrix(c(
  1.2, 0.7, 0.8,
  0.7, 1.8, 0.6,
  0.8, 0.6, 1.4
), nrow = 3)

XI <- rmvnorm(n, sigma = Sigma)

X <- XI[, 1]
Z <- XI[, 2]
W <- XI[, 3]

Y <- 1.2 * X + 0.4 * Z + 0.7 * W +
  0.2 * W * Z +
  0.7 * W * X +
  1.2 * X * Z +
  2.2 * X * Z * W + rnorm(n, sd = sqrt(2))

createInd <- \(x, lambda, tau, epsilon = 0.2) {
  tau + lambda * x + rnorm(n, sd = sqrt(epsilon))
}

x1 <- createInd(X, 1.0, 1.2)
x2 <- createInd(X, 0.8, 0.8)
x3 <- createInd(X, 0.9, 1.0)

z1 <- createInd(Z, 1.0, 1.2)
z2 <- createInd(Z, 0.8, 0.8)
z3 <- createInd(Z, 0.9, 1.0)

w1 <- createInd(W, 1.0, 1.2)
w2 <- createInd(W, 0.8, 0.8)
w3 <- createInd(W, 0.9, 1.0)

y1 <- createInd(Y, 1.0, 1.2)
y2 <- createInd(Y, 0.8, 0.8)
y3 <- createInd(Y, 0.9, 1.0)

data.3way <- data.frame(x1, x2, x3,
                        z1, z2, z3,
                        w1, w2, w3,
                        y1, y2, y3)
m.3way <- '
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

# First we compile the STAN model, this can be slow
# and is therefore done once, such that the compiled
# STAN code can be reused for the same model syntax later
compiled_model.3way <- compile_stan_model(m.3way)

# Fit the model based on the compiled STAN code
fit.3way <- modsem_stan(
  compiled_model = compiled_model_3way,
  data   = data.3way,
  chains = 2,
  iter   = 4000 # More iterations should yield more stable estimates
)
#> Regressions:
#>                  Estimate  Std.Error  z.value  P(>|z|)
#>   Y ~           
#>     X               1.040      0.159    6.535    0.000
#>     Z               0.423      0.128    3.296    0.001
#>     W               0.696      0.126    5.524    0.000
#>     X:Z             1.385      0.173    8.023    0.000
#>     X:W             0.733      0.175    4.193    0.000
#>     Z:W             0.190      0.178    1.070    0.284
#>     X:Z:W           2.227      0.107   20.731    0.000

summary(fit.3way)
standardized_estimates(fit.3way)
parameter_estimates(fit.3way)
