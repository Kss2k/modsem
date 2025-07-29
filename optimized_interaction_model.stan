// ──────────────────────────────────────────────────────────────────────────────
//    sem_latent_interaction_fast.stan
//    (flat / undefined priors – add yours where desired)
// ──────────────────────────────────────────────────────────────────────────────
functions {
  // nothing needed – the product is vectorised now
}

data {
  // ── general sizes ───────────────────────────────────────────────────────────
  int<lower=1> N;                       // sample size
  int<lower=1> K;                       // # indicators
  int<lower=1> N_XIS;                   // # ξ's (exogenous LVs)
  int<lower=1> N_LVS;                   // total # latent variables (= ξ + η)
  int<lower=1> N_ETAS;                  // # η's  (endogenous LVs)
  int<lower=0> N_INT;                   // # latent interaction terms

  // ── index arrays (pre‑built in R) ───────────────────────────────────────────
  int<lower=1> N_FREE_LAMBDA;
  int<lower=1> lambda_row[N_FREE_LAMBDA];
  int<lower=1> lambda_col[N_FREE_LAMBDA];

  int<lower=1> N_FREE_GAMMA;
  int<lower=1> gamma_row[N_FREE_GAMMA];
  int<lower=1> gamma_col[N_FREE_GAMMA];

  int<lower=1> N_FREE_OMEGA;
  int<lower=1> omega_row[N_FREE_OMEGA];
  int<lower=1> omega_col[N_FREE_OMEGA];

  // first loading index (for identification, fixed to 1)
  int<lower=1> first_loading_row[N_LVS];

  // latent‑interaction column pairs:   PROD_q = eta[, left_q] .* eta[, right_q]
  int<lower=1> prod_left[N_INT];
  int<lower=1> prod_right[N_INT];

  // ── priors that are *structural zeros/ones*, NOT parameters ────────────────
  // (same as before but now just boolean arrays if you still need them)

  // ── observed data ───────────────────────────────────────────────────────────
  matrix[N, K] Y;
}

parameters {
  // measurement model
  vector[N_FREE_LAMBDA] lambda;

  // structural model
  vector[N_FREE_GAMMA]  gamma;
  vector[N_FREE_OMEGA]  omega;

  // LV means
  vector[N_LVS] alpha;          // you may set some of these to zero in data

  // covariance factors (Cholesky) – free, flat
  cholesky_factor_cov[K]  L_theta;
  cholesky_factor_cov[N_LVS] L_psi;

  // latent scores
  matrix[N, N_LVS] XI;
}

transformed parameters {
  // ── measurement matrices ───────────────────────────────────────────────────
  matrix[K, N_LVS] Lambda = rep_matrix(0, K, N_LVS);
  for (k in 1:N_FREE_LAMBDA)
    Lambda[ lambda_row[k] , lambda_col[k] ] = lambda[k];
  for (j in 1:N_LVS)                       // set first loading = 1
    Lambda[ first_loading_row[j] , j ] = 1;

  // ── structural matrices ────────────────────────────────────────────────────
  matrix[N_ETAS, N_LVS] Gamma = rep_matrix(0, N_ETAS, N_LVS);
  for (k in 1:N_FREE_GAMMA)
    Gamma[ gamma_row[k] , gamma_col[k] ] = gamma[k];

  matrix[N_ETAS, N_INT] Omega = rep_matrix(0, N_ETAS, N_INT);
  for (k in 1:N_FREE_OMEGA)
    Omega[ omega_row[k] , omega_col[k] ] = omega[k];

  // ── latent variables (row‑wise for vectorisation) ──────────────────────────
  matrix[N, N_LVS] eta = XI;                 // start with ξ and disturbances

  // build N×N_INT matrix of element‑wise products in one shot
  matrix[N, N_INT] PROD;
  for (q in 1:N_INT)
    PROD[, q] = eta[, prod_left[q]] .* eta[, prod_right[q]];

  // propagate structural part (matrix products, no loops)
  eta[, (N_XIS + 1):N_LVS] =
      eta * Gamma'                 // linear effects
    + PROD * Omega';               // latent interactions

  // ── expected indicators ────────────────────────────────────────────────────
  matrix[N, K] mu = eta * Lambda';           // N × K
}

model {
  // flat / undefined priors – **add your own here if desired**

  // latent variables
  {
    row_vector[N_LVS] alpha_row = to_row_vector(alpha);
    for (n in 1:N)
      XI[n] ~ multi_normal_cholesky(alpha_row, L_psi);
  }

  // observed indicators
  for (n in 1:N)
    Y[n] ~ multi_normal_cholesky(mu[n]', L_theta);
}

/*–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  generated quantities – optional: cov matrices, replicated data, …
––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––*/
