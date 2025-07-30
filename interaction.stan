// sem_latent_interaction_fast.stan
//   Three latent factors (X, Z, Y) plus a latent interaction X×Z → Y
//   Identification: first loading per factor fixed to 1
//   All indicator intercepts τ are free
//   Fast version: non-centred Ξ, Cholesky likelihoods, vectorised products
//--------------------------------------------------------------------------

functions {
  /**
   * Multiply the N-length latent vectors whose indices are listed in
   * idx_row (zeros mean “not used”) and return the N-vector product.
   */
  vector get_product(int N_LVS,
                     int N,
                     int idx_row[N_LVS],   // <- fixed-size int array
                     matrix ETA) {
    vector[N] out = rep_vector(1, N);
    for (j in 1:N_LVS)
      if (idx_row[j] > 0)
        out = out .* ETA[, idx_row[j]];
    return out;
  }
}

//////////////////////////////////////////////////////////
// DATA
//////////////////////////////////////////////////////////
data {
  int<lower=1>  N;                   // sample size
  int<lower=1>  K;                   // number of indicators
  int<lower=1>  N_XIS;               // number of ξ (exogenous) factors
  int<lower=1>  N_ETAS;              // number of η  (endogenous) factors
  int<lower=2>  N_LVS;               // total latent factors = N_XIS + N_ETAS
  int<lower=0>  N_INT;               // number of latent interactions

  //--------------------------------------------------------------------
  //  mask matrices (0 = fixed to 0, 1 = free, 2 = reference loading 1)
  //--------------------------------------------------------------------
  matrix[K,         N_LVS]  LAMBDA_mask;     // measurement model
  matrix[N_ETAS,    N_LVS]  GAMMA_mask;      // structural paths
  matrix[N_ETAS,    N_INT]  OMEGA_mask;      // interaction paths
  matrix[N_INT,     N_LVS]  PRODUCTS_mask;   // which LVs in each product
  matrix[N_LVS,     N_LVS]  PSI_mask;        // latent (co)variances
  matrix[K,         K]      THETA_mask;      // indicator (co)variances
  vector[K]                   TAU_mask;      // indicator intercepts
  vector[N_LVS]               ALPHA_mask;    // latent means

  //--------------------------------------------------------------------
  //  helper counts = number of free parameters in each block
  //--------------------------------------------------------------------
  int N_FREE_LAMBDA;
  int N_FREE_GAMMA;
  int N_FREE_OMEGA;
  int N_FREE_PSI_D;
  int N_FREE_PSI_L;
  int N_FREE_THETA_D;
  int N_FREE_THETA_L;
  int N_FREE_TAU;
  int N_FREE_ALPHA;

  //--------------------------------------------------------------------
  //  index array for products (built once in R/Python would be better)
  //--------------------------------------------------------------------
  array[N_INT, N_LVS] int PROD_IDX;     // 0 where not used

  // observed data
  matrix[N, K] Y;
}

//////////////////////////////////////////////////////////
// PARAMETERS
//////////////////////////////////////////////////////////
parameters {
  //---------------- measurement ----------------
  vector[N_FREE_LAMBDA]      lambda;
  vector[N_FREE_TAU]         tau;
  vector<lower=0>[N_FREE_THETA_D] theta_d;
  vector[N_FREE_THETA_L]          theta_l;

  //---------------- structural -----------------
  vector[N_FREE_GAMMA]       gamma;
  vector[N_FREE_OMEGA]       omega;
  vector[N_FREE_ALPHA]       alpha;
  vector<lower=0>[N_FREE_PSI_D] psi_d;
  vector[N_FREE_PSI_L]            psi_l;

  //---------------- latent variables -----------
  matrix[N, N_LVS] XI_raw;                     // non-centred
  cholesky_factor_corr[N_LVS] Lcorr_Psi;
  vector<lower=0>[N_LVS]      sigma_Psi;
}

//////////////////////////////////////////////////////////
// TRANSFORMED PARAMETERS
//////////////////////////////////////////////////////////
transformed parameters {
  //-------------------------------------------------
  // 1.  free-parameter matrices / vectors
  //-------------------------------------------------
  matrix[K,          N_LVS] Lambda   = rep_matrix(0, K, N_LVS);
  matrix[N_ETAS,     N_LVS] Gamma    = rep_matrix(0, N_ETAS, N_LVS);
  matrix[N_ETAS,     N_INT] Omega    = rep_matrix(0, N_ETAS, N_INT);
  matrix[N_LVS,      N_LVS] Psi      = rep_matrix(0, N_LVS, N_LVS);
  matrix[K,          K]     Theta    = rep_matrix(0, K, K);
  vector[K]                 Tau      = rep_vector(0, K);
  vector[N_LVS]             Alpha    = rep_vector(0, N_LVS);

  { //--------------------------------------------
    int k;

    // Lambda (loadings)
    k = 1;
    for (r in 1:K)
      for (c in 1:N_LVS)
        if (LAMBDA_mask[r, c] == 2)         // reference loading
          Lambda[r, c] = 1;
        else if (LAMBDA_mask[r, c] == 1) {
          Lambda[r, c] = lambda[k];
          ++k;
        }

    //--- Gamma -----------------------------------------------------------------
    k = 1;
    for (r in 1:N_ETAS) {
      for (c in 1:N_LVS) {
        if (GAMMA_mask[r, c]) {
          Gamma[r, c] = gamma[k];
          ++k;
        }
      }
    }

    //--- Omega -----------------------------------------------------------------
    k = 1;
    for (r in 1:N_ETAS) {
      for (c in 1:N_INT) {
        if (OMEGA_mask[r, c]) {
          Omega[r, c] = omega[k];
          ++k;
        }
      }
    }

    //--- Theta lower triangle ---------------------------------------------------
    k = 1;
    for (i in 2:K) {
      for (j in 1:(i - 1)) {
        if (THETA_mask[i, j]) {
          Theta[i, j] = theta_l[k];
          Theta[j, i] = theta_l[k];
          ++k;
        }
      }
    }

    //--- Psi lower triangle -----------------------------------------------------
    k = 1;
    for (i in 2:N_LVS) {
      for (j in 1:(i - 1)) {
        if (PSI_mask[i, j]) {
          Psi[i, j] = psi_l[k];
          Psi[j, i] = psi_l[k];
          ++k;
        }
      }
    }
    // Alpha
    k = 1;
    for (i in 1:N_LVS)
      if (ALPHA_mask[i])
        Alpha[i] = alpha[k++];

    // Tau
    k = 1;
    for (i in 1:K)
      if (TAU_mask[i])
        Tau[i] = tau[k++];
  }

  //-------------------------------------------------
  // 2.  latent variables (non-centred)
  //-------------------------------------------------
  matrix[N_LVS, N_LVS] L_Psi   = diag_pre_multiply(sigma_Psi, Lcorr_Psi);
  matrix[N,       N_LVS] XI    = rep_matrix(Alpha', N) + XI_raw * L_Psi';

  //-------------------------------------------------
  // 3.  build η factors including interactions
  //-------------------------------------------------
  matrix[N, N_LVS] ETA = XI;                       // start with ξ & η residuals

  for (e in 1:N_ETAS) {
    int idx = N_XIS + e;                           // column in ETA for ηₑ

    // structural main-effects
    for (j in 1:N_LVS)
      if (GAMMA_mask[e, j])
        ETA[, idx] += Gamma[e, j] * ETA[, j];

    // latent interactions
    for (p in 1:N_INT)
      if (OMEGA_mask[e, p])
        ETA[, idx] += Omega[e, p] *
          get_product(N_LVS, N, PROD_IDX[p], ETA);
  }

  //-------------------------------------------------
  // 4.  expected indicators
  //-------------------------------------------------
  matrix[N_LVS, K] Lambda_T = Lambda';             // cache transpose
  matrix[N, K]     Mu       = ETA * Lambda_T;
  for (k in 1:K)
    Mu[, k] += Tau[k];
}

//////////////////////////////////////////////////////////
// MODEL
//////////////////////////////////////////////////////////
model {
  //-------------------- priors (weakly informative) ---------
  lambda    ~ normal(0, 0.5);
  gamma     ~ normal(0, 0.5);
  omega     ~ normal(0, 0.5);
  alpha     ~ normal(0, 1);

  sigma_Psi ~ lognormal(0, 0.5);
  Lcorr_Psi ~ lkj_corr_cholesky(2);

  theta_d   ~ lognormal(-1, 0.5);
  theta_l   ~ normal(0, 0.3);
  psi_d     ~ lognormal(-1, 0.5);
  psi_l     ~ normal(0, 0.3);

  //-------------------- latent variables (non-centred) -----
  to_vector(XI_raw) ~ normal(0, 1);        // standard normal i.i.d.

  //-------------------- likelihood --------------------------
  {
    matrix[K, K] L_Theta = cholesky_decompose(Theta);
    for (n in 1:N)
      Y[n] ~ multi_normal_cholesky(Mu[n], L_Theta);
  }
}

//////////////////////////////////////////////////////////
//  Optionally add a generated-quantities block to monitor
//  the interaction scores, predicted Y, etc.
//////////////////////////////////////////////////////////
