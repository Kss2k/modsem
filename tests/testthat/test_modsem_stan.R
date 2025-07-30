devtools::load_all()

library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)      # cache compiled models
options(mc.cores = parallel::detectCores()) 

m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

compiled_model <- compile_stan_model(m1)

fit <- modsem_stan(compiled_model = compiled_model, 
                   data = oneInt, iter = 4000)
summary(fit)
standardized_estimates(fit)
set.seed(29723234)


n <- 2000
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

# Compile a STAN model based on the lavaan syntax to save 
# time when re-estimating the model,
# I.e., we don't want to compile each iteraction of the simulation
compiled_model_3way <- compile_stan_model(m.3way)
fit.3way <- modsem_stan(
  compiled_model = compiled_model_3way,
  data   = data.3way,
  chains = 2,
  iter   = 4000
)

standardized_estimates(fit.3way)
parameter_estimates(fit.3way)
