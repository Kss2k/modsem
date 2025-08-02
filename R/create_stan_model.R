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
compile_stan_model <- function(model.syntax, compile = TRUE) {
  message("Compiling STAN model...")

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


  STAN_INDS_LV <- function(lV) {
    inds <- indsLVs[[lV]]
    k <- length(inds)

    # data {}
    dataInds <- sprintf("matrix[N, %d] INDICATORS_%s;", k, lV)

    # parameters {}
    labTau <- sprintf("%s__INTERCEPT", inds)
    labResSD <- sprintf("%s__COVARIANCE__%s", inds, inds)
    labLambda <- sprintf("%s__MEASUREMENT__%s", lV, inds[-1]) # first loading is constrained

    parTau    <- sprintf("real %s;", labTau)
    parResSD  <- sprintf("real<lower=0> %s;", labResSD)
    parLambda <- sprintf("real %s;", labLambda)

    # model {}
    idx <- seq_along(inds)
    modInd <- sprintf("INDICATORS_%s[,%d] ~ normal(%s + %s * %s, %s);", 
                      lV, idx, labTau, c("1", labLambda), lV, labResSD)

    parameters <- collapse(parTau, parResSD, parLambda)
    model      <- collapse(modInd)
    data       <- collapse(dataInds)

    list(parameters = parameters, model = model, data = data)
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

  if (compile) stanModel <- rstan::stan_model(model_code = stanModelSyntax)
  else         stanModel <- NULL

  list(syntax = stanModelSyntax,
       stan_model = stanModel,
       info = list(lVs = lVs,
                   xis = xis,
                   etas = etas,
                   indsLVs = indsLVs,
                   allIndsXis = allIndsXis,
                   allIndsEtas = allIndsEtas))
}


get_stan_data <- function(compiled_model, data, impute.na = FALSE) {
  lVs         <- compiled_model$info$lVs
  indsLVs     <- compiled_model$info$indsLVs
  allIndsXis  <- compiled_model$info$allIndsXis
  allIndsEtas <- compiled_model$info$allIndsEtas

  # clean data
  INDICATORS <- cleanAndSortData(data, allIndsXis, allIndsEtas, impute.na = impute.na)

  stan_data <- list(N = nrow(data))

  for (lV in lVs) {
    name <- sprintf("INDICATORS_%s", lV)
    inds <- indsLVs[[lV]]
    stan_data[[name]] = INDICATORS[, inds, drop = FALSE]
  }

  stan_data
}


EXAMPLE_STAN_SYNTAX_3WAY <- '
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
  Y[,11] ~ normal(tau[11] + lambda_y[1] * Y_lat,  sigma_e[11]);   // y2
  Y[,12] ~ normal(tau[12] + lambda_y[2] * Y_lat,  sigma_e[12]);   // y3
}

generated quantities {
  real cov_XZ = cov_vector(X, Z);
  real cov_XW = cov_vector(X, W);
  real cov_ZW = cov_vector(Z, W);
}
'
