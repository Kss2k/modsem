STAN_LAVAAN_MODELS <- rlang::env(
  syntaxes = NULL,
  compiled = NULL
)


STAN_OPERATOR_LABELS <- c(
  "__INTERCEPT" = "~1",
  "__REGRESSION__" = "~",
  "__COVARIANCE__" = "~~",
  "__MEASUREMENT__" = "=~",
  "__XWITH__" = ":"
)


STAN_SYNTAX_BLOCKS <- "
functions {
%s

}

data {
%s

}

parameters {
%s

}

transformed parameters {
%s

}

model {
%s

}

generated quantities {
%s

}
"


#' Compile \code{STAN} model based on a \code{lavaan} model
#'
#' @param model.syntax \code{lavaan} syntax.
#' @param compile Should compilation be performed? If \code{FALSE} only the \code{STAN}
#'   is generated, and not compiled.
#' @param force Should compilation of previously compiled models be forced?
#' @export
compile_stan_model <- function(model.syntax, compile = TRUE, force = FALSE,
                               ordered = NULL) {
  if (is.null(ordered))
    ordered <- character(0)

  if (length(ordered)) {
    model.syntax <- paste(
      model.syntax,
      sprintf("#ORDERED %s", paste0(ordered, collapse = " ")),
      sep = "\n"
    )
  }

  parTable <- modsemify(model.syntax)

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

  collapse <- function(..., sep = "\n") {
    args <- list(...)
    do.call(
      "paste",
      args = c(lapply(args, FUN = \(arg) paste0(arg, collapse = sep)),
               list(sep = sep))
    )
  }

  FUNCTIONS <-  "
    // empirical covariance helper (for generated quantities)
    real cov_vector(vector a, vector b) {
      int N = num_elements(a);

      // Center the vectors
      vector[N] a_centered = a - mean(a);
      vector[N] b_centered = b - mean(b);

      // dot_product is faster than manual loop and clearer than dot_self
      return dot_product(a_centered, b_centered) / (N - 1);
    }
  "

  DATA <- "int<lower=1> N;  // sample size"
  PARAMETERS <- NULL
  TRANSFORMED_PARAMETERS <- NULL
  MODEL <- NULL
  GENERATED_QUANTITIES <- NULL
  EXCLUDE.PARS <- NULL


  STAN_INDS_LV <- function(lV) {
    inds <- indsLVs[[lV]]
    k <- length(inds)

    # --- data {}: keep the matrix for legacy pipelines (cols unused if ordinal)
    dataInds <- sprintf("matrix[N, %d] INDICATORS_%s;", k, lV)

    # We will accumulate code per-indicator to allow mixing ordinal/continuous.
    par_lines   <- character()
    model_lines <- character()
    data_lines  <- dataInds

    # Names for parameters that exist irrespective of ordinal/continuous
    # Intercepts:
    labTau <- sprintf("%s__INTERCEPT", inds)
    par_lines <- c(par_lines, sprintf("real %s;", labTau))

    # Loadings: first is fixed to 1 by construction; free for the rest
    # (we keep your original labels X__MEASUREMENT__x2 etc.)
    free_inds  <- inds[-1]
    labLambda  <- if (length(free_inds)) sprintf("%s__MEASUREMENT__%s", lV, free_inds) else character(0)
    if (length(labLambda)) {
      par_lines <- c(par_lines, sprintf("real %s;", labLambda))
    }

    # Residual SDs (continuous-only) and ordinal cutpoints (ordinal-only) will be added in the loop

    for (j in seq_along(inds)) {
      ind <- inds[j]
      # loading term: first indicator fixed to 1, others use their free parameter
      loading_term <- if (j == 1L) "1" else sprintf("%s__MEASUREMENT__%s", lV, ind)

      if (ind %in% ordered) {
        # ---- ORDINAL INDICATOR ----
        # data: category count + integer responses
        data_lines <- paste0(
          data_lines, "\n",
          sprintf("int<lower=2> K_%s;", ind), "\n",
          sprintf("int<lower=1, upper=K_%s> INDICATORS_%s[N];", ind, ind)
        )

        # parameters: NO residual SD; add cutpoints
        par_lines <- c(
          par_lines,
          sprintf("ordered[K_%s - 1] %s__CUTPOINTS;", ind, ind)
        )

        # model: cumulative logit with linear predictor
        model_lines <- c(
          model_lines,
          sprintf("{"),
          sprintf("  vector[N] eta_%s = %s__INTERCEPT + %s * %s;", ind, ind, loading_term, lV),
          sprintf("  INDICATORS_%s ~ ordered_logistic(eta_%s, %s__CUTPOINTS);", ind, ind, ind),
          sprintf("}")
        )

      } else {
        # ---- CONTINUOUS INDICATOR ----
        # parameters: residual SD present
        par_lines <- c(
          par_lines,
          sprintf("real<lower=0> %s__COVARIANCE__%s;", ind, ind)
        )

        # model: Gaussian measurement eq.
        model_lines <- c(
          model_lines,
          sprintf("INDICATORS_%s[,%d] ~ normal(%s__INTERCEPT + %s * %s, %s__COVARIANCE__%s);",
                  lV, j, ind, loading_term, lV, ind, ind)
        )
      }
    }

    list(
      parameters = collapse(par_lines),
      model      = collapse(model_lines),
      data       = collapse(data_lines)
    )
  }

  STAN_PAR_XIS <- function(xis) {
    k <- length(xis)

    # parameters {}
    parLOmega   <- sprintf("cholesky_factor_corr[%d] L_Omega;", k)
    parSqrtDPhi <- sprintf("vector<lower=0>[%d] sqrtD_Phi;", k)
    # parZXiMat   <- sprintf("matrix[N, %d] Z_XI_Matrix;", k)
    parXiMat   <- sprintf("matrix[N, %d] XI_Matrix;", k)

    # transformed parameters {}
    tparLSigma    <- sprintf("matrix[%d, %d] L_Sigma = diag_pre_multiply(sqrtD_Phi, L_Omega);", k, k)
    tparMuXi      <- sprintf("row_vector[%d] MU_XI = rep_row_vector(0, %d);", k, k)
    tparXiArr     <- sprintf("array[N] vector[%d] XI_Array;", k)
    tparXiArrFill <- "for (i in 1:N) {XI_Array[i] = (XI_Matrix[i, ])';}"
    tparXi       <- sprintf("matrix[N, %d] XI_Matrix = Z_XI_Matrix * L_Sigma';", k)

    xiVectors <- NULL
    for (i in seq_along(xis)) {
      name     <- xis[[i]]
      xiVector <- sprintf("vector[N] %s = col(XI_Matrix, %d);", name, i)
      xiVectors <- c(xiVectors, xiVector)
    }

    tparXiVectors <- collapse(xiVectors)

    # model {}
    modLOmega <- sprintf("L_Omega ~ lkj_corr_cholesky(1);")
    # modZXiMat <- "to_vector(Z_XI_Matrix) ~ normal(0, 1);"
    # modXiMat <- "for (n in 1:N) {XI_Matrix[n,] ~ multi_normal_cholesky(MU_XI, L_Sigma);}"
    modXiArr <- "XI_Array ~ multi_normal_cholesky(MU_XI, L_Sigma);"

    parameters  <- collapse(parLOmega, parSqrtDPhi, parXiMat)
    tparameters <- collapse(tparLSigma, tparXiVectors, tparMuXi, tparXiArr, tparXiArrFill)
    model       <- collapse(modLOmega, modXiArr)

    EXCLUDE.PARS <<- c(EXCLUDE.PARS, "XI_Array", "XI_Matrix", xis)

    list(parameters = parameters, transformed_parameters = tparameters,
         model = model)
  }


  STAN_PAR_ETA <- function(eta) {
    indeps <- unique(parTable[parTable$lhs == eta & parTable$op == "~", "rhs"])
    indeps <- stringr::str_replace_all(indeps, pattern = ":",
                                       replacement = "__XWITH__")

    # parameters {}
    labBeta <- sprintf("%s__REGRESSION__%s", eta, indeps)
    labSD <- sprintf("%s__COVARIANCE__%s", eta, eta)

    parBeta   <- sprintf("real %s;", labBeta)
    parValues <- sprintf("vector[N] %s;", eta)
    parSD  <- sprintf("real<lower=0> %s;", labSD)

    # model {}
    projEta <- collapse(sprintf("%s * %s", labBeta, indeps), sep = " + ")
    modEta <- sprintf("%s ~ normal(%s, %s);", eta, projEta, labSD)

    parameters <- collapse(parValues, parSD, parBeta)

    EXCLUDE.PARS <<- c(EXCLUDE.PARS, eta)

    list(parameters = parameters, model = modEta)
  }


  STAN_COMPUTED_PRODUCTS <- function(intTerms) {
    transformed_parameters <- NULL

    for (intTerm in intTerms) {
      elems <- stringr::str_split(intTerm, pattern = ":")[[1L]]
      prodEq <- collapse(elems, sep = " .* ")
      nameIntTerm <- collapse(elems, sep = "__XWITH__")

      transformed_parameters <- c(
        transformed_parameters,
        sprintf("vector[N] %s = %s;", nameIntTerm, prodEq)
      )

      EXCLUDE.PARS <<- c(EXCLUDE.PARS, nameIntTerm)
    }

    list(transformed_parameters = collapse(transformed_parameters))
  }


  STAN_COMPUTED_COVARIANCES <- function(vars) {
    vars   <- stringr::str_replace_all(vars, pattern = ":", replacement = "__XWITH__")
    combos <- getUniqueCombos(vars)
    generated_quantities <- NULL

    for (i in seq_len(NROW(combos))) {
      lhs <- combos[[1L]][[i]]
      rhs <- combos[[2L]][[i]]

      generated_quantities <- c(
        generated_quantities,
        sprintf("real %s__COVARIANCE__%s = cov_vector(%s, %s);",
                lhs, rhs, lhs, rhs)
      )
    }

    list(generated_quantities = collapse(generated_quantities))
  }


  STAN_COMPUTED_VARIANCES <- function(vars) {
    vars   <- stringr::str_replace_all(vars, pattern = ":", replacement = "__XWITH__")
    generated_quantities <- NULL

    for (var in vars) {
      generated_quantities <- c(
        generated_quantities,
        sprintf("real %s__COVARIANCE__%s = cov_vector(%s, %s);",
                var, var, var, var)
      )
    }

    list(generated_quantities = collapse(generated_quantities))
  }


  add2block <- function(FUN, ...) {
    blocks <- FUN(...)

    for (name in names(blocks)) {
      block <- blocks[[name]]
      switch(name,
             functions = {
               FUNCTIONS <<- collapse(FUNCTIONS, block)
             },

             data = {
               DATA <<- collapse(DATA, block)
             },

             parameters = {
               PARAMETERS <<- collapse(PARAMETERS, block)
             },

             transformed_parameters = {
               TRANSFORMED_PARAMETERS <<- collapse(TRANSFORMED_PARAMETERS, block)
             },

             model = {
               MODEL <<- collapse(MODEL, block)
             },

             generated_quantities = {
               GENERATED_QUANTITIES <<- collapse(GENERATED_QUANTITIES, block)
             }
      )
    }
  }

  for (lV in lVs)   add2block(STAN_INDS_LV, lV = lV)
  for (eta in etas) add2block(STAN_PAR_ETA, eta = eta)

  add2block(STAN_PAR_XIS, xis = xis)
  add2block(STAN_COMPUTED_PRODUCTS, intTerms = intTerms)
  add2block(STAN_COMPUTED_COVARIANCES, vars = c(xis, intTerms))
  add2block(STAN_COMPUTED_VARIANCES, vars = c(xis, intTerms))

  stanModelSyntax <- sprintf(STAN_SYNTAX_BLOCKS,
                             FUNCTIONS, DATA, PARAMETERS,
                             TRANSFORMED_PARAMETERS,
                             MODEL, GENERATED_QUANTITIES)

  SYNTAXES <- STAN_LAVAAN_MODELS$syntaxes
  COMPILED <- STAN_LAVAAN_MODELS$compiled
  match    <- STAN_LAVAAN_MODELS$syntaxes == model.syntax

  if (compile && any(match) && !force) {
    message("Reusing compiled Stan model...")
    stanModel <- last(COMPILED[match]) # if a duplicate somehow appears, pick last/newest match

  } else if (compile) {
    message("Compiling Stan model...")

    stanModel <- rstan::stan_model(model_code = stanModelSyntax)

    SYNTAXES <- c(SYNTAXES, model.syntax)
    COMPILED <- c(COMPILED, stanModel)

    STAN_LAVAAN_MODELS$syntaxes <- SYNTAXES
    STAN_LAVAAN_MODELS$compiled <- COMPILED

  } else stanModel <- NULL


  list(syntax = stanModelSyntax,
       stan_model = stanModel,
       info = list(lVs = lVs,
                   xis = xis,
                   etas = etas,
                   indsLVs = indsLVs,
                   allIndsXis = allIndsXis,
                   allIndsEtas = allIndsEtas,
                   exclude.pars = unique(EXCLUDE.PARS)))
}


getStanData <- function(compiled_model, data, missing = "listwise", ordered = NULL) {
  if (is.null(ordered)) ordered <- character(0)

  lVs         <- compiled_model$info$lVs
  indsLVs     <- compiled_model$info$indsLVs
  allIndsXis  <- compiled_model$info$allIndsXis
  allIndsEtas <- compiled_model$info$allIndsEtas

  # 1) Pre-coerce requested ordinal columns in the raw data
  #    (safe even if columns are already numeric; ensures stable ordering)
  for (col in ordered) {
    if (!col %in% names(data)) {
      stop("`ordered` indicator '", col, "' not found in data.")
    }
    data[[col]] <- as.integer(as.ordered(data[[col]]))
  }

  # 2) Run your existing missing-data preparation (listwise or otherwise)
  INDICATORS <- prepDataModsemDA(data, allIndsXis, allIndsEtas,
                                 missing = missing)$data.full

  stan_data <- list(N = nrow(INDICATORS))

  # 3) Emit the latent-specific indicator matrices (unchanged)
  for (lV in lVs) {
    name <- sprintf("INDICATORS_%s", lV)
    inds <- indsLVs[[lV]]
    if (!all(inds %in% colnames(INDICATORS))) {
      missing_cols <- paste(setdiff(inds, colnames(INDICATORS)), collapse = ", ")
      stop("Indicators missing from prepared data for latent '", lV, "': ", missing_cols)
    }
    stan_data[[name]] <- as.matrix(INDICATORS[, inds, drop = FALSE])
  }

  # 4) For each ordered indicator, add integer vector 1..K and K
  remap_to_consecutive <- function(x) {
    # x should be atomic, no NAs expected after listwise handling
    u <- sort(unique(x))
    # Create a 1..K mapping even if labels were not consecutive (e.g., 0/2/5)
    map <- setNames(seq_along(u), as.character(u))
    as.integer(unname(map[as.character(x)]))
  }

  for (ind in ordered) {
    if (!ind %in% colnames(INDICATORS)) {
      stop("`ordered` indicator '", ind, "' not found after preprocessing.")
    }
    x_raw <- INDICATORS[, ind]
    # Ensure integer 1..K coding regardless of original labels
    x_int <- remap_to_consecutive(x_raw)
    K     <- as.integer(max(x_int))

    stan_data[[sprintf("INDICATORS_%s", ind)]] <- x_int
    stan_data[[sprintf("K_%s", ind)]]          <- K
  }

  stan_data
}


# General Stan Model for estimating SEMs with interaction terms
STAN_MODEL_GENERAL <- "
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

  array[N] vector[N_LVS] marginalXI;
  row_vector[N_LVS] marginalMeanXI = rep_row_vector(0, N_LVS);

  for (n in 1:N) {
    marginalXI[n] = (XI[n, ])' - Alpha;
  }

  array[N] vector[K] marginalY;
  row_vector[K] marginalMeanY = rep_row_vector(0, K);

  for (n in 1:N) {
    marginalY[n] = (Y[n, ] - X[n, ])';
  }
}

model {
  // Priors
  // No priors (yet)

  marginalXI ~ multi_normal(marginalMeanXI, Psi);
  marginalY  ~ multi_normal(marginalMeanY, Theta);
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

## Example
# model.syntax <- '
#   X =~ x1 + x2 + x3
#   Z =~ z1 + z2 + z3
#   Y =~ y1 + y2 + y3
#
#   Y ~ X + Z + X:Z
# '
# stan_data <- specifyModelSTAN(model.syntax, data = oneInt)
# stan_model <- stan_model(model_code = STAN_MODEL_GENERAL)
#
# fit <- sampling(
#   object = stan_model,
#   data   = stan_data,
#   chains = 2,
#   iter   = 2000,
#   warmup = 1000
# )
#
# summary(fit, c("gamma"))
# summary(fit, c("omega"))
# summary(fit, c("Psi"))
# stan_rhat(fit)
# stan_trace(fit, "omega")
# stan_trace(fit, "gamma")
