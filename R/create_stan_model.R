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
compile_stan_model <- function(model.syntax, compile = TRUE, force = FALSE) {
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


getStanData <- function(compiled_model, data, impute.na = FALSE) {
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
