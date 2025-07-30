devtools::load_all()
library(rstan)
rstan_options(auto_write = TRUE)      # cache compiled models
options(mc.cores = parallel::detectCores()) 


stan.syntax <- '
// sem_latent_interaction.stan (updated)
// SEM with three latent factors (X, Z, Y) and latent interaction X×Z → Y
// Identification: first loading fixed to 1; *latent means constrained to 0*.
// All nine indicator intercepts τ are now free.

functions {
  // empirical covariance helper (for generated quantities)
  real cov_vector(vector a, vector b) {
    return (dot_self(a - mean(a)))/(num_elements(a)-1);
  }
}

data {
  int<lower=1> N;                       // sample size
  matrix[N, 9] Y;                       // observed indicators (x1–x3, z1–z3, y1–y3)
}

parameters {
  //--------------------------------------
  // Measurement part
  //--------------------------------------
  vector[2] lambda_x;                   // loadings: x2, x3
  vector[2] lambda_z;                   // loadings: z2, z3
  vector[2] lambda_y;                   // loadings: y2, y3

  vector[9] tau;                        // ALL indicator intercepts free
  vector<lower=0>[9] sigma_e;           // residual SDs

  //--------------------------------------
  // Latent part
  //--------------------------------------
  real<lower=0> sigma_X;
  real<lower=0> sigma_Z;
  real<lower=0> sigma_Y;

  real beta_X;                          // structural slopes
  real beta_Z;
  real beta_XZ;

  vector[N] X;                          // person‑level latent scores
  vector[N] Z;
  vector[N] Y_lat;                      // endogenous latent Y
}

transformed parameters {
  vector[N] XZ = X .* Z;                // latent interaction term
}

model {
  //--------------------------------------
  // Priors
  //--------------------------------------
  // lambda_x ~ normal(1, 1);
  // lambda_z ~ normal(1, 1);
  // lambda_y ~ normal(1, 1);

  // tau      ~ normal(0, 2);

  // sigma_X  ~ exponential(1);
  // sigma_Z  ~ exponential(1);
  // sigma_Y  ~ exponential(1);
  // sigma_e  ~ exponential(1);

  // beta_X   ~ normal(0, 1);
  // beta_Z   ~ normal(0, 1);
  // beta_XZ  ~ normal(0, 1);

  //--------------------------------------
  // Latent variable distributions (means fixed to 0)
  //--------------------------------------
  X     ~ normal(0, sigma_X);
  Z     ~ normal(0, sigma_Z);
  Y_lat ~ normal(beta_X * X + beta_Z * Z + beta_XZ * XZ, sigma_Y);

  //--------------------------------------
  // Measurement model
  //--------------------------------------
  for (n in 1:N) {
    // X indicators
    Y[n,1] ~ normal(tau[1] + X[n],                    sigma_e[1]);   // x1 (λ fixed = 1)
    Y[n,2] ~ normal(tau[2] + lambda_x[1] * X[n],      sigma_e[2]);   // x2
    Y[n,3] ~ normal(tau[3] + lambda_x[2] * X[n],      sigma_e[3]);   // x3

    // Z indicators
    Y[n,4] ~ normal(tau[4] + Z[n],                    sigma_e[4]);   // z1 (λ fixed = 1)
    Y[n,5] ~ normal(tau[5] + lambda_z[1] * Z[n],      sigma_e[5]);   // z2
    Y[n,6] ~ normal(tau[6] + lambda_z[2] * Z[n],      sigma_e[6]);   // z3

    // Y indicators
    Y[n,7] ~ normal(tau[7] + Y_lat[n],                sigma_e[7]);   // y1 (λ fixed = 1)
    Y[n,8] ~ normal(tau[8] + lambda_y[1] * Y_lat[n],  sigma_e[8]);   // y2
    Y[n,9] ~ normal(tau[9] + lambda_y[2] * Y_lat[n],  sigma_e[9]);   // y3
  }
}

generated quantities {
  real corr_XZ = cov_vector(X, Z) / (sd(X) * sd(Z));
  real corr_XY = cov_vector(X, Y_lat) / (sd(X) * sd(Y_lat));
  real corr_ZY = cov_vector(Z, Y_lat) / (sd(Z) * sd(Y_lat));
}
'

devtools::load_all()
library(rstan)
rstan_options(auto_write = TRUE)      # cache compiled models
options(mc.cores = parallel::detectCores()) 

model <- rstan::stan_model(model_code = stan.syntax)
Y <- oneInt[c("x1", "x2", "x3",
              "z1", "z2", "z3",
              "y1", "y2", "y3")]
stan_data <- list(
  N = nrow(Y),
  Y = Y 
)

fit <- sampling(
  object = model,
  data   = stan_data,
  chains = 2,
  iter   = 2000,
  warmup = 1000
)

summary(fit, c("beta_X", "beta_Z", "beta_XZ"))
stan_rhat(fit)
stan_trace(fit, "beta_X")
stan_trace(fit, "beta_XZ")
stan_trace(fit, "beta_Z")


stan.syntax <- "
// sem_latent_interaction.stan (updated)
// SEM with three latent factors (X, Z, Y) and latent interaction X×Z → Y
// Identification: first loading fixed to 1; *latent means constrained to 0*.
// All nine indicator intercepts τ are now free.

functions {
  vector getIthProduct(int i, int K, int N, N_INT, matrix OMEGA, matrix ETA) {

    vector[N] product;
    product[:] = 1; // fill with ones

    for (j in 1:N_INT) {
      if (OMEGA[i, j]) {
        product = product .* ETA[:, i];
      }
    }

    return product; 
  }
}


data {
  int<lower=1>                 N; // Sample size
  int<lower=1>                 K; // Number of indicators
  int<lower=1>                 N_XIS; // Number of xis
  int<lower=1>                 N_ETAS; // Number of etas
  int<lower=2>                 N_LVS; // Number of lVs
  int<lower=0>                 N_INT; // Number of interaction terms

  matrix[K, N_LVS]             LAMBDA; // Structure of measurement model
  matrix[N_ETAS, N_LVS + N_INT] GAMMA; // Structure of structural model
  matrix[N_ETAS, N_INT]         OMEGA; // Structure of interaction-coefficients
  matrix[N_INT, N_LVS]       PRODUCTS; // Structure of interaction terms

  matrix[N_LVS, N_LVS]         PSI; // Structure of latent covariances
  matrix[K, K]                 THETA; // Structure of indicator covariances
  vector[K]                    TAU; // Structure of indicator intercepts
  vector[N_LVS]                ALPHA; // Structure of LV intercepts

  int<lower=0> N_FREE_LAMBDA;     // sum(abs(LAMBDA))
  int<lower=0> N_FREE_GAMMA;      // sum(abs(GAMMA))
  int<lower=0> N_FREE_OMEGA; 

  int<lower=0> N_FREE_DIAG_THETA;  // sum(abs(diag(THETA)))
  int<lower=0> N_FREE_LOWER_THETA; // sum(abs(THETA[is.lower(THETA)]))

  int<lower=0> N_FREE_DIAG_PSI;  // sum(abs(diag(PSI)))
  int<lower=0> N_FREE_LOWER_PSI; // sum(abs(GAMMA[is.lower(GAMMA)]))

  int<lower=0> N_FREE_TAU;         // sum(abs(TAU))
  int<lower=0> N_FREE_ALPHA;       // sum(abs(ALPHA))

  // Observed Data
  matrix[N, K] Y;
}


parameters {
  // Measurement model
  vector[N_FREE_LAMBDA] lambda;
  vector[N_FREE_TAU] tau;
  vector<lower=0>[N_FREE_DIAG_THETA] theta_d;
  vector[N_FREE_LOWER_THETA] theta_l;
 
  // Structural model
  vector[N_FREE_GAMMA] gamma;
  vector[N_FREE_ALPHA] ALPHA;
  vector<lower=0>[N_FREE_DIAG_PSI] psi_d;
  vector[N_FREE_LOWER_PSI] psi_l;

  // LVs
  matrix[N, N_LVS] XI;

  // Indicator disturbances
  // matrix[N, K] EPSILON;
}

transformed_parameters {
  // Fill matrices
  matrix[K, N_LVS]              Lambda; // Structure of measurement model
  matrix[N_ETAS, N_LVS + N_INT] Gamma; // Structure of structural model

  matrix[N_LVS, N_LVS]         Psi;   // Structure of latent covariances
  matrix[K, K]                 Theta; // Structure of indicator covariances
  vector[K]                    Tau;   // Structure of indicator intercepts
  vector[N_LVS]                Alpha; // Structure of LV intercepts

  // Fill Lambda
  int k = 1;
  for (i in 1:K) {
    int filledFirst = 0;

    for (j in 1:N_LVS) {
      int fill = LAMBDA[i, j];

      if (fill && !filledFirst) {
        Lambda[i, j] = 1;
      } else if (fill) {
        Lambda[i, j] = lambda[k];
        k = k + 1;
      }
    }
  }

  // Fill Gamma
  k = 1;
  for (i in 1:N_ETAS) {
    for (j in 1:N_LVS) {
      int fill = GAMMA[i, j];

      if (fill) {
        Gamma[i, j] = gamma[k];
        k = k + 1;
      }
    }
  }
  
  // Fill OMEGA
  k = 1;
  for (i in 1:N_ETAS) {
    for (j in 1:N_INT) {
      int fill = OMEGA[i, j];

      if (fill) {
        Omega[i, j] = omega[k];
        k = k + 1;
      }
    }
  }

  // Fill Diagonal Theta
  k = 1;
  for (i in 1:K) {
    int fill = THETA[i, i];

    if (fill) {
      Theta[i, i] = theta_d[k];
      k = k + 1;
    }
  }
  
  // Fill Off-Diagonal Theta
  k = 1;
  for (i in 1:K) {

    for (j in 1:i) {
      if (j >= i) {
        continue;
      }
    
      int fill = THETA[i, j];
    
      if (fill) {
      Theta[i, j] = theta_l[k];
      Theta[j, i] = theta_l[k];
      k = k + 1;
    }
  }

  // Fill Diagonal Psi
  k = 1;
  for (i in 1:N_LVS) {
    int fill = PSI[i, i];

    if (fill) {
      Psi[i, i] = psi_d[k];
      k = k + 1;
    }
  }
  
  // Fill Off-Diagonal Theta
  k = 1;
  for (i in 1:N_LVS) {

    for (j in 1:i) {
      if (j >= i) {
        continue;
      }
    
      int fill = PSI[i, j];
    
      if (fill) {
      Psi[i, j] = psi_l[k];
      Psi[j, i] = psi_l[k];
      k = k + 1;
    }
  }

  // Fill Alpha
  k = 1;
  for (i in 1:N_LVS) {
    int fill = ALPHA[i];
    if (fill) {
      Alpha[i] = alpha[k];
      k = k + 1;
    }
  }
  
  // Fill Tau
  k = 1;
  for (i in 1:K) {
    int fill = TAU[i];
    if (fill) {
      Tau[i] = tau[k];
      k = k + 1;
    }
  }

  // Calculate LV values
  matrix[N, N_LVS + N_INT] ETA;

  for (i in 1:N_LVS) {
    ETA[, i] = XI[, i]; // For xis the actual values are filled
                          // For etas we start with the disturbances
  }

  for (i in 1:N_ETAS) {
    for (j in 1:N_LVS) {
      if (GAMMA[i, j]) {
        ETA[, i] = ETA[, i] + GAMMA[i, j] * ETA[, j];
      }
    }
  
    for (j in 1:N_INT) {
      if (OMEGA[i, j]) {
        ETA[, i] = 
          ETA[, i] + OMEGA[i, j] * getIthProduct(j, K, N, N_INT, OMEGA, ETA);
      }
    }
  }

  // Calculate Expected Indicator values
  matrix[N, K] X = ETA * Lambda';
  for (i in 1:K) {
    X[, i] = Tau[i] + X[, i];
  }
}


model {
  // Priors
  // No priors (yet)

  // Latent variable distributions
  row_vector[N_LVS] alpha_row = to_row_vector(Alpha);   // convert mean
  for (n in 1:N) {
    to_row_vector(XI[n]) ~ multi_normal(alpha_row, Psi);
  }
  matrix[N, K] mu = X;                 // expected scores (built in TP block)
  for (n in 1:N) {
    to_row_vector(Y[n]) ~ multi_normal(to_row_vector(mu[n]), Theta);
  }
}

"
devtools::load_all()
library(rstan)
rstan_options(auto_write = TRUE)      # cache compiled models
options(mc.cores = parallel::detectCores()) 

model <- rstan::stan_model(model_code = stan.syntax)
Y <- oneInt[c("x1", "x2", "x3",
              "z1", "z2", "z3",
              "y1", "y2", "y3")]
stan_data <- list(
  N = nrow(Y),
  Y = Y 
)

fit <- sampling(
  object = model,
  data   = stan_data,
  chains = 2,
  iter   = 2000,
  warmup = 1000
)

summary(fit, c("beta_X", "beta_Z", "beta_XZ"))
stan_rhat(fit)
stan_trace(fit, "beta_X")
stan_trace(fit, "beta_XZ")
stan_trace(fit, "beta_Z")



stan.syntax.3way <- '
// sem_latent_interaction.stan (updated)
// SEM with three latent factors (X, Z, Y) and latent interaction X×Z → Y
// Identification: first loading fixed to 1; *latent means constrained to 0*.
// All nine indicator intercepts τ are now free.

functions {
  // empirical covariance helper (for generated quantities)
  real cov_vector(vector a, vector b) {
    return (dot_self(a - mean(a)))/(num_elements(a)-1);
  }
}

data {
  int<lower=1> N;                       // sample size
  matrix[N, 12] Y;                       // observed indicators (x1–x3, z1–z3, y1–y3)
}

parameters {
  //--------------------------------------
  // Measurement part
  //--------------------------------------
  vector[2] lambda_x;                   // loadings: x2, x3
  vector[2] lambda_z;                   // loadings: z2, z3
  vector[2] lambda_w;                   // loadings: w2, w3
  vector[2] lambda_y;                   // loadings: y2, y3

  vector[12] tau;                        // ALL indicator intercepts free
  vector<lower=0>[12] sigma_e;           // residual SDs

  //--------------------------------------
  // Latent part
  //--------------------------------------
  real<lower=0> sigma_X;
  real<lower=0> sigma_Z;
  real<lower=0> sigma_W;
  real<lower=0> sigma_Y;

  real beta_X;                          // structural slopes
  real beta_Z;
  real beta_W;
  real beta_XZ;
  real beta_XW;
  real beta_ZW;
  real beta_XZW;

  vector[N] X;                          // person‑level latent scores
  vector[N] Z;
  vector[N] W;
  vector[N] Y_lat;                      // endogenous latent Y
}

transformed parameters {
  vector[N] XZ = X .* Z;                // latent interaction term
  vector[N] XW = X .* W;                // latent interaction term
  vector[N] ZW = Z .* W;                // latent interaction term
  vector[N] XZW = X .* ZW;                // latent interaction term
}

model {
  //--------------------------------------
  // Priors
  //--------------------------------------
  // lambda_x ~ normal(1, 1);
  // lambda_z ~ normal(1, 1);
  // lambda_y ~ normal(1, 1);

  // tau      ~ normal(0, 2);

  // sigma_X  ~ exponential(1);
  // sigma_Z  ~ exponential(1);
  // sigma_Y  ~ exponential(1);
  // sigma_e  ~ exponential(1);

  // beta_X   ~ normal(0, 1);
  // beta_Z   ~ normal(0, 1);
  // beta_XZ  ~ normal(0, 1);

  //--------------------------------------
  // Latent variable distributions (means fixed to 0)
  //--------------------------------------
  X     ~ normal(0, sigma_X);
  Z     ~ normal(0, sigma_Z);
  W     ~ normal(0, sigma_W);
  Y_lat ~ normal(beta_X * X + 
                 beta_Z * Z +
                 beta_W * W +
                 beta_XZ * XZ +
                 beta_XW * XW +
                 beta_ZW * ZW +
                 beta_XZW * XZW, 
                 sigma_Y);

  //--------------------------------------
  // Measurement model
  //--------------------------------------
    // X indicators
  Y[,1] ~ normal(tau[1] + X,                    sigma_e[1]);   // x1 (λ fixed = 1)
  Y[,2] ~ normal(tau[2] + lambda_x[1] * X,      sigma_e[2]);   // x2
  Y[,3] ~ normal(tau[3] + lambda_x[2] * X,      sigma_e[3]);   // x3

  // Z indicators
  Y[,4] ~ normal(tau[4] + Z,                    sigma_e[4]);   // z1 (λ fixed = 1)
  Y[,5] ~ normal(tau[5] + lambda_z[1] * Z,      sigma_e[5]);   // z2
  Y[,6] ~ normal(tau[6] + lambda_z[2] * Z,      sigma_e[6]);   // z3

  // W indicators
  Y[,7] ~ normal(tau[7] + W,                sigma_e[7]);   // y1 (λ fixed = 1)
  Y[,8] ~ normal(tau[8] + lambda_w[1] * W,  sigma_e[8]);   // y2
  Y[,9] ~ normal(tau[9] + lambda_w[2] * W,  sigma_e[9]);   // y3

  // Y indicators
  Y[,10] ~ normal(tau[10] + Y_lat,                sigma_e[10]);   // y1 (λ fixed = 1)
  Y[,11] ~ normal(tau[11] + lambda_y[1] * Y_lat[n],  sigma_e[11]);   // y2
  Y[,12] ~ normal(tau[12] + lambda_y[2] * Y_lat[n],  sigma_e[12]);   // y3
}

generated quantities {
  real cov_XZ = cov_vector(X, Z);
  real cov_XW = cov_vector(X, W);
  real cov_ZW = cov_vector(Z, W);
}
'


model <- rstan::stan_model(model_code = stan.syntax)
Y <- oneInt[c("x1", "x2", "x3",
              "z1", "z2", "z3",
              "y1", "y2", "y3")]
stan_data <- list(
  N = nrow(Y),
  Y = Y 
)

fit <- sampling(
  object = model,
  data   = stan_data,
  chains = 2,
  iter   = 2000,
  warmup = 1000
)

summary(fit, c("beta_X", "beta_Z", "beta_XZ"))
stan_rhat(fit)
stan_trace(fit, "beta_X")
stan_trace(fit, "beta_XZ")
stan_trace(fit, "beta_Z")


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


model <- rstan::stan_model(model_code = stan.syntax.3way)
Y <- data.3way[c("x1", "x2", "x3",
                 "z1", "z2", "z3",
                 "w1", "w2", "w3",
                 "y1", "y2", "y3")]
stan_data <- list(
  N = nrow(Y),
  Y = Y 
)

fit <- sampling(
  object = model,
  data   = stan_data,
  chains = 2,
  iter   = 5000,
  warmup = 3000
)

beta <- c("beta_X", "beta_Z", "beta_W",
          "beta_XZ", "beta_XW", "beta_ZW",
          "beta_XZW")
lambda <- c("lambda_x", "lambda_z",
            "lambda_w", "lambda_y")
sigma_e <- c("sigma_e")
sigma <- c("sigma_X", "sigma_Z", "sigma_W", "sigma_Y")
tau  <- c("tau")
pars <- c(lambda, beta, sigma, sigma_e, tau)
summary(fit, pars)
stan_rhat(fit)
stan_trace(fit, pars)
stan_trace(fit, "beta_XZ")
stan_trace(fit, "beta_Z")

