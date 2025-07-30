devtools::load_all()

library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)      # cache compiled models
options(mc.cores = parallel::detectCores()) 

stan.syntax <- "
// sem_latent_interaction.stan (updated)
// SEM with three latent factors (X, Z, Y) and latent interaction X×Z → Y
// Identification: first loading fixed to 1; *latent means constrained to 0*.
// All nine indicator intercepts τ are now free.

functions {
  vector getIthProduct(int i, int N_LVS, int N, matrix PRODUCTS, matrix ETA) {

    vector[N] product = rep_vector(0, N);

    int firstFound = 0;
    for (j in 1:N_LVS) {

      if (PRODUCTS[i, j]) {

        if (!firstFound) {
          product = ETA[, j];
          firstFound = 1;

        } else {
          product = product .* ETA[, j];
        }
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
  matrix[N_ETAS, N_LVS]         GAMMA; // Structure of structural model
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
  vector[N_FREE_ALPHA] alpha;
  vector[N_FREE_OMEGA] omega;
  vector<lower=0>[N_FREE_DIAG_PSI] psi_d;
  vector[N_FREE_LOWER_PSI] psi_l;

  // LVs
  matrix[N, N_LVS] XI;

  // Indicator disturbances
  // matrix[N, K] EPSILON;
}


transformed parameters {
  // Declare matrices
  matrix[K, N_LVS]     Lambda; // Structure of measurement model
  matrix[N_ETAS, N_LVS] Gamma; // Structure of structural model
  matrix[N_ETAS, N_INT] Omega;

  matrix[N_LVS, N_LVS]         Psi;   // Structure of latent covariances
  matrix[K, K]                 Theta; // Structure of indicator covariances
  vector[K]                    Tau;   // Structure of indicator intercepts
  vector[N_LVS]                Alpha; // Structure of LV intercepts
  
  // Fill Matrices
  Lambda = rep_matrix(0, K, N_LVS);
  Gamma  = rep_matrix(0, N_ETAS, N_LVS);
  Omega  = rep_matrix(0, N_ETAS, N_INT);
  Theta  = rep_matrix(0, K, K);
  Psi    = rep_matrix(0, N_LVS, N_LVS);
  Tau    = rep_vector(0, K);
  Alpha  = rep_vector(0, N_LVS);

  {
    // Fill Lambda
    int k = 1;
    for (j in 1:N_LVS) {
      real filledFirst = 0;

      for (i in 1:K) {
        real fill = LAMBDA[i, j];

        if (fill && !filledFirst) {
          Lambda[i, j] = 1;
          filledFirst = 1;

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
        real fill = GAMMA[i, j];

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
        real fill = OMEGA[i, j];

        if (fill) {
          Omega[i, j] = omega[k];
          k = k + 1;
        }
      }
    }

    // Fill Diagonal Theta
    k = 1;
    for (i in 1:K) {
      real fill = THETA[i, i];

      if (fill) {
        Theta[i, i] = theta_d[k];
        k = k + 1;
      }
    }
    
    // Fill Off-Diagonal Theta
    k = 1;
    for (i in 1:K) {

      for (j in 1:(i-1)) {
        real fill = THETA[i, j];
      
        if (fill) {
          Theta[i, j] = theta_l[k];
          Theta[j, i] = theta_l[k];
          k = k + 1;
        }
      }
    }

    // Fill Diagonal Psi
    k = 1;
    for (i in 1:N_LVS) {
      real fill = PSI[i, i];

      if (fill) {
        Psi[i, i] = psi_d[k];
        k = k + 1;
      }
    }
    
    // Fill Off-Diagonal Psi
    k = 1;
    for (i in 1:N_LVS) {

      for (j in 1:(i-1)) {
        real fill = PSI[i, j];
      
        if (fill) {
          Psi[i, j] = psi_l[k];
          Psi[j, i] = psi_l[k];
          k = k + 1;
        }
      }
    }

    // Fill Alpha
    k = 1;
    for (i in 1:N_LVS) {
      real fill = ALPHA[i];
      if (fill) {
        Alpha[i] = alpha[k];
        k = k + 1;
      }
    }
    
    // Fill Tau
    k = 1;
    for (i in 1:K) {
      real fill = TAU[i];

      if (fill) {
        Tau[i] = tau[k];
        k = k + 1;
      }
    }
  }

  // Calculate LV values
  matrix[N, N_LVS] ETA;

  for (i in 1:N_LVS) {
    ETA[, i] = XI[, i]; // For xis the actual values are filled
                          // For etas we start with the disturbances
  }

  for (i in 1:N_ETAS) {
    int idx = N_XIS + i;

    for (j in 1:N_LVS) {
      if (GAMMA[i, j]) {
        ETA[, idx] = ETA[, idx] + Gamma[i, j] * ETA[, j];
      }
    }
  
    for (j in 1:N_INT) {
      if (OMEGA[i, j]) {
        ETA[, idx] = 
          ETA[, idx] + Omega[i, j] * getIthProduct(j, N_LVS, N, PRODUCTS, ETA);
      }
    }
  }

  // Calculate Expected Indicator values
  matrix[N, K] X = ETA * Lambda';
  for (i in 1:K) {
    X[, i] = Tau[i] + X[, i];
  }

  // DEBUGGING
  // print(\"Omega\");
  // print(Omega);
  // print(\"Gamma\");
  // print(Gamma);
  // print(\"Lambda\");
  // print(Lambda);
  // print(\"PRODUCTS\");
  // print(PRODUCTS);
  // print(\"Theta\");
  // print(Theta);
  // print(\"Psi\");
  // print(Psi);

  // print(\"head(ETA)\");
  // print(ETA[1:5, ]);
  // 
  // print(\"head(XI)\");
  // print(XI[1:5, ]);
  // 
  // print(\"head(X)\");
  // print(X[1:5, ]);

  // print(\"head(XZ)\");
  // print(getIthProduct(1, N_LVS, N, PRODUCTS, ETA)[1:5]);
}

model {
  // Priors
  // No priors (yet)

  // Latent variable distributions
  {
    row_vector[N_LVS] alpha_row = to_row_vector(Alpha);   // convert mean
    for (n in 1:N) {
      XI[n] ~ multi_normal(alpha_row, Psi);
    }
  } 

  {
    matrix[N, K] mu = X;                 // expected scores (built in TP block)
    for (n in 1:N) {
      Y[n] ~ multi_normal(to_row_vector(mu[n]), Theta);
    }
  }
}

"

# Functions
specifyModelSTAN <- function(syntax = NULL,
                             data = NULL,
                             auto.fix.first = TRUE,
                             auto.fix.single = TRUE,
                             orthogonal.y = FALSE,
                             orthogonal.x = FALSE,
                             mean.observed = TRUE) {
  if (!is.null(syntax)) parTable <- modsemify(syntax)
  stopif(is.null(parTable), "No parTable found")

  # endogenous variables (etas)model
  etas    <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  etas    <- etas[length(etas):1] # reverse
  numEtas <- length(etas)

  indsEtas    <- getIndsLVs(parTable, etas)
  numIndsEtas <- vapply(indsEtas, FUN.VALUE = vector("integer", 1L),
                        FUN = length)
  allIndsEtas    <- unique(unlist(indsEtas))
  numAllIndsEtas <- length(allIndsEtas)

  # exogenouts variables (xis) and interaction terms
  intTerms      <- unique(parTable[grepl(":", parTable$rhs), "rhs"])
  numInts       <- length(intTerms)
  varsInts      <- stringr::str_split(intTerms, pattern = ":")
  varsInts      <- stats::setNames(varsInts, nm = intTerms)
  allVarsInInts <- unique(unlist(varsInts))

  xis           <- getXis(parTable, checkAny = TRUE)
  numXis        <- length(xis)

  indsXis    <- getIndsLVs(parTable, xis)
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis    <- unique(unlist(indsXis))
  numAllIndsXis <- length(allIndsXis)

  lVs <- c(xis, etas)
  numLVs <- length(lVs)
  indsLVs <- getIndsLVs(parTable, lVs)
  allInds <- c(allIndsXis, allIndsEtas)
  numAllInds <- length(allInds)

  # clean data
  data <- cleanAndSortData(data, allIndsXis, allIndsEtas, impute.na = impute.na)


  # measurement model x
  listLambda <- constructLambda(lVs, indsLVs, parTable = parTable,
                                auto.fix.first = auto.fix.first)
  LAMBDA <- listLambda$numeric
  LAMBDA[is.na(LAMBDA)] <- 1

  listTau <- constructTau(lVs, indsLVs, parTable = parTable)
  TAU <- listTau$numeric
  TAU[is.na(TAU)] <- 1

  listAlpha <- constructAlpha(lVs, parTable = parTable,
                              mean.observed = mean.observed)
  ALPHA <- listAlpha$numeric
  ALPHA[is.na(ALPHA)] <- 1

  listTheta <- constructTheta(lVs, indsLVs, parTable = parTable,
                              auto.fix.single = auto.fix.single)
  THETA <- listTheta$numeric
  THETA[is.na(THETA)] <- 1

  listGamma <- constructGamma(etas, lVs, parTable = parTable)
  GAMMA <- listGamma$numeric
  GAMMA[is.na(GAMMA)] <- 1

  listOmega <- constructGamma(etas, intTerms, parTable = parTable)
  OMEGA <- listOmega$numeric
  OMEGA[is.na(OMEGA)] <- 1
  
  listPsi  <- constructPsi(etas, parTable = parTable, orthogonal.y = orthogonal.y)
  psi      <- listPsi$numeric

  listPhi <- constructPhi(xis, method = "qml", cov.syntax = NULL,
                          parTable = parTable, orthogonal.x = orthogonal.x)
  phi <- listPhi$numeric

  PSI <- diagPartitionedMat(phi, psi)
  PSI[is.na(PSI)] <- 1

  PRODUCTS <- matrix(0, nrow = numInts, ncol = numLVs, dimnames = list(intTerms, lVs))
  for (intTerm in intTerms) {
    elems <- varsInts[[intTerm]]

    PRODUCTS[intTerm, elems] <- 1
  }

  lowerTri <- \(X) X[lower.tri(X, diag = FALSE)]

  list(N      = nrow(data),
       K      = numAllInds,
       N_XIS  = numXis,
       N_ETAS = numEtas,
       N_LVS  = numLVs,
       N_INT  = numInts,
  
       LAMBDA   = LAMBDA,
       GAMMA    = GAMMA,
       OMEGA    = OMEGA,
       PRODUCTS = PRODUCTS,

       PSI   = PSI,
       THETA = THETA,

       ALPHA = as.vector(ALPHA),
       TAU   = as.vector(TAU),

       N_FREE_LAMBDA = sum(is.na(listLambda$numeric)),
       N_FREE_GAMMA = sum(GAMMA),
       N_FREE_OMEGA  = sum(OMEGA),

       N_FREE_DIAG_THETA  = sum(diag(THETA)),
       N_FREE_LOWER_THETA = sum(lowerTri(THETA)),

       N_FREE_DIAG_PSI  = sum(diag(PSI)),
       N_FREE_LOWER_PSI = sum(lowerTri(PSI)),

       N_FREE_ALPHA = sum(ALPHA),
       N_FREE_TAU   = sum(TAU),

       Y = as.matrix(data)
  )
}

m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

stan_data <- specifyModelSTAN(m1, data = oneInt)
stan_model <- stan_model(model_code = stan.syntax)

fit <- sampling(
  object = stan_model,
  data   = stan_data,
  chains = 2,
  iter   = 2000,
  warmup = 1000
)

summary(fit, c("gamma"))
summary(fit, c("omega"))
stan_rhat(fit)
stan_trace(fit, "omega")
stan_trace(fit, "gamma")


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

stan_data.3way <- specifyModelSTAN(m.3way, data = data.3way)

fit.3way <- sampling(
  object = stan_model,
  data   = stan_data.3way,
  chains = 2,
  iter   = 2000,
  warmup = 1000
)

summary(fit.3way, c("gamma"))
summary(fit.3way, c("omega"))
stan_rhat(fit.3way)
stan_trace(fit.3way, "omega")
stan_trace(fit.3way, "gamma")
