STAN_LAVAAN_MODELS <- rlang::env(
  syntaxes = NULL,
  compiled = NULL
)


STAN_OPERATOR_LABELS <- c(
  "__INTERCEPT" = "~1",
  "__REGRESSION__" = "~",
  "__COVARIANCE__" = "~~",
  # "__STDDEV__"     = "~~~",
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
                               ordered = NULL, ordered.link = c("logit", "probit")) {
  ordered.link <- tolower(ordered.link)
  ordered.link <- match.arg(ordered.link)

  if (is.null(ordered))
    ordered <- character(0)

  if (length(ordered)) {
    model.syntax <- paste(
      model.syntax,
      sprintf("#ORDERED %s", paste0(ordered, collapse = " ")),
      sep = "\n"
    )
  }

  if (ordered.link == "probit") {
    resVarFormulaInd <- "1"
    resSDFormulaInd  <- "1"
    orderedLinkFun   <- "ordered_probit"
  } else if (ordered.link == "logit") {
    resVarFormulaInd <- "(pi()^2)/3"
    resSDFormulaInd  <- sprintf("sqrt(%s)", resVarFormulaInd)
    orderedLinkFun   <- "ordered_logistic"
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

  # NEW: normalize and precompute ordinal sets
  if (is.null(ordered)) ordered <- character(0)
  ord_set <- unique(ordered)
  is_allOrdinal_lv <- vapply(lVs, function(lv) {
    inds <- indsLVs[[lv]]
    length(inds) > 0 && all(inds %in% ord_set)
  }, logical(1))

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

  getMod <- \(lhs, op, rhs) getModParTable(lhs, op, rhs, parTable = parTable)

  DATA <- "int<lower=1> N;  // sample size"
  PARAMETERS <- NULL
  TRANSFORMED_PARAMETERS <- NULL
  MODEL <- NULL
  GENERATED_QUANTITIES <- NULL
  EXCLUDE.PARS <- NULL

  STAN_INDS_LV <- function(lV) {
    inds <- indsLVs[[lV]]
    k <- length(inds)

    parLines   <- character()
    modelLines <- character()
    quantLines <- character()
    dataInds   <- sprintf("matrix[N, %d] INDICATORS_%s;", length(inds), lV)
    dataLines  <- dataInds
    tparLines  <- character()

    # Intercepts only for continuous indicators
    contInds <- inds[!inds %in% ord_set]
    ordInds  <- inds[inds %in% ord_set]

    fixContVar <- length(contInds) + length(ordInds) <= 1L

    for (i in seq_along(contInds)) {
      ind             <- contInds[[i]]
      parInterceptInd <- sprintf("real %s__INTERCEPT;", ind)
      parLambdaInd    <- sprintf("real %s__MEASUREMENT__%s;", ind, ind)

      modTauInd    <- getMod(lhs = ind, op = "~1", rhs = "")
      modLambdaInd <- getMod(lhs = lV,  op = "=~", rhs = ind)
      modVarInd    <- getMod(lhs = ind, op = "~~", rhs = ind)

      fixTauInd     <- !is.na(modTauInd)
      fixLambdaInd  <- !is.na(modLambdaInd) || i == 1L
      fixContVarInd <- !is.na(modVarInd) || fixContVar

      if (fixTauInd) {
        parTauInd  <- NULL
        tparTauInd <- sprintf("real %s__INTERCEPT = %s;", ind, modTauInd)
      } else {
        parTauInd  <- sprintf("real %s__INTERCEPT;", ind)
        tparTauInd <- NULL
      }

      if (fixLambdaInd) {
        mod           <- if (!is.na(modLambdaInd)) modLambdaInd else "1"
        parLambdaInd  <- NULL
        tparLambdaInd <- sprintf("real %s__MEASUREMENT__%s = %s;", lV, ind, mod)
      } else {
        parLambdaInd  <- sprintf("real %s__MEASUREMENT__%s;", lV, ind)
        tparLambdaInd <- NULL
      }

      if (fixContVarInd) {
        mod        <- if (!is.na(modVarInd)) modVarInd else "0"
        parResSDInd  <- NULL
        tparResSDInd <- sprintf("real %s__STDDEV__%s = sqrt(%s);", ind, ind, mod)
      } else {
        parResSDInd  <- sprintf("real<lower=0> %s__STDDEV__%s;", ind, ind)
        tparResSDInd <- NULL
      }

      tparResVarInd <- sprintf(
        "real %s__COVARIANCE__%s = %s__STDDEV__%s^2;",
        ind, ind, ind, ind
      )

      modMeasrInd <- sprintf(
        "INDICATORS_%s[,%d] ~ normal(%s__INTERCEPT + %s__MEASUREMENT__%s * %s, %s__STDDEV__%s);",
        lV, i, ind, lV, ind, lV, ind, ind
      )

      parLines <- c(
        parLines,
        parTauInd,
        parLambdaInd,
        parResSDInd
      )

      tparLines <- c(
        tparLines,
        tparTauInd,
        tparLambdaInd,
        tparResSDInd,
        tparResVarInd
      )

      modelLines <- c(
        modelLines,
        modMeasrInd
      )
    }

    for (j in seq_along(ordInds)) {
      ind <- ordInds[[j]]

      if (TRUE)    parLambdaInd  <- sprintf("real<lower=0> %s__UMEASUREMENT__%s;", lV, ind)
      else         parLambdaInd  <- sprintf("real %s__UMEASUREMENT__%s;", lV, ind)

      dataK_Ind     <- sprintf("int<lower=2> K_%s;", ind)
      dataInd       <- sprintf("int<lower=1, upper=K_%s> ORD_INDICATOR_%s[N];", ind, ind)
      parCutInd     <- sprintf("ordered[K_%s - 1] %s__UCUTPOINTS;", ind, ind)
      responseInd   <- sprintf("vector[N] LV_RESPONSE_%s_ = %s__UMEASUREMENT__%s * %s;", ind, lV, ind, lV)
      likelihoodInd <- sprintf("ORD_INDICATOR_%s ~ %s(LV_RESPONSE_%s_, %s__UCUTPOINTS);",
                               ind, orderedLinkFun, ind, ind)
      modMeasrInd <- sprintf("{\n  %s\n  %s\n}", responseInd, likelihoodInd)

      quantTotalVarInd <- sprintf(
        "real TOTAL_VAR_%s_ = (%s__UMEASUREMENT__%s ^ 2) * cov_vector(%s, %s) + %s;",
        ind, lV, ind, lV, lV, resVarFormulaInd
      )

      quantStdCutPoints <- sprintf(
        "ordered[K_%s - 1] %s__CUTPOINTS = %s__UCUTPOINTS / sqrt(TOTAL_VAR_%s_);",
        ind, ind, ind, ind
      )

      quantLambdaInd <- sprintf(
        "real %s__MEASUREMENT__%s = %s__UMEASUREMENT__%s / sqrt(TOTAL_VAR_%s_);",
        lV, ind, lV, ind, ind
      )

      quantResVarInd <- sprintf(
        "real %s__COVARIANCE__%s = %s / TOTAL_VAR_%s_;",
        ind, ind, resVarFormulaInd, ind
      )

      dataLines  <- c(dataLines, dataK_Ind, dataInd)
      parLines   <- c(parLines, parLambdaInd, parCutInd)
      quantLines <- c(quantLines, quantTotalVarInd, quantStdCutPoints,
                      quantResVarInd, quantLambdaInd)
      modelLines <- c(modelLines, modMeasrInd)
    }

    list(
         parameters = collapse(parLines),
         model      = collapse(modelLines),
         data       = collapse(dataLines),
         transformed_parameters = collapse(tparLines),
         generated_quantities = collapse(quantLines)
    )
  }

  # Centered parameterization
  STAN_PAR_XIS <- function(xis) {
    k <- length(xis)

    # Identify which xi's are all-ordinal
    fix_idx <- which(vapply(xis, function(lv) is_allOrdinal_lv[lv], logical(1)))
    free_idx <- setdiff(seq_len(k), fix_idx)

    parLOmega   <- sprintf("cholesky_factor_corr[%d] L_Omega;", k)

    if (length(free_idx) > 0L) {
      parSqrtDPhi <- sprintf("vector<lower=0>[%d] sqrtD_Phi_free;", length(free_idx))
    } else {
      parSqrtDPhi <- NULL
    }

    parXiMat   <- sprintf("matrix[N, %d] XI_Matrix;", k)

    # transformed parameters: rebuild sqrtD_Phi with fixed 1's at fix_idx
    tparBuildSqrt <- c(
      sprintf("vector[%d] sqrtD_Phi;", k),
      " {",
      if (length(free_idx) > 0L) "  int c = 1;" else NULL,
      if (length(seq_len(k)) > 0L) {
        paste0(
          vapply(seq_len(k), function(i) {
            if (i %in% fix_idx) {
              sprintf("  sqrtD_Phi[%d] = 1;", i)
            } else {
              sprintf("  sqrtD_Phi[%d] = sqrtD_Phi_free[c]; c += 1;", i)
            }
          }, character(1L)),
          collapse = "\n"
        )
      } else NULL,
      " }"
    )

    tparLSigma    <- sprintf("matrix[%d, %d] L_Sigma = diag_pre_multiply(sqrtD_Phi, L_Omega);", k, k)
    tparMuXi      <- sprintf("row_vector[%d] MU_XI = rep_row_vector(0, %d);", k, k)
    tparXiArr     <- sprintf("array[N] vector[%d] XI_Array;", k)
    tparXiArrFill <- "for (i in 1:N) {XI_Array[i] = (XI_Matrix[i, ])';}"

    xiVectors <- NULL
    for (i in seq_along(xis)) {
      name     <- xis[[i]]
      xiVector <- sprintf("vector[N] %s = col(XI_Matrix, %d);", name, i)
      xiVectors <- c(xiVectors, xiVector)
    }
    tparXiVectors <- collapse(xiVectors)

    modLOmega <- "L_Omega ~ lkj_corr_cholesky(1);"
    modXiArr  <- "XI_Array ~ multi_normal_cholesky(MU_XI, L_Sigma);"

    parameters  <- collapse(c(parLOmega, parSqrtDPhi, parXiMat))
    tparameters <- collapse(c(tparBuildSqrt, tparLSigma, tparXiVectors, tparMuXi, tparXiArr, tparXiArrFill))
    model       <- collapse(c(modLOmega, modXiArr))

    EXCLUDE.PARS <<- c(EXCLUDE.PARS, "XI_Array", "XI_Matrix", xis)

    list(parameters = parameters, transformed_parameters = tparameters, model = model)
  }

  # Non-centered parameterization
  # STAN_PAR_XIS <- function(xis) {
  #   k <- length(xis)
  #
  #   # Identify which xi's are all-ordinal (variance fixed to 1)
  #   fix_idx  <- which(vapply(xis, function(lv) is_allOrdinal_lv[lv], logical(1)))
  #   free_idx <- setdiff(seq_len(k), fix_idx)
  #
  #   # ---------------- parameters ----------------
  #   parLOmega <- sprintf("cholesky_factor_corr[%d] L_Omega;", k)
  #
  #   if (length(free_idx) > 0L) {
  #     parSqrtDPhi <- sprintf("vector<lower=0>[%d] sqrtD_Phi_free;", length(free_idx))
  #   } else {
  #     parSqrtDPhi <- NULL
  #   }
  #
  #   # Non-centered: matrix of standard normals (N x k)
  #   parZ <- sprintf("matrix[N, %d] z_XI;", k)
  #
  #   # ---------------- transformed parameters ----------------
  #   # Rebuild sqrtD_Phi with fixed 1's at fix_idx
  #   tparBuildSqrt <- c(
  #     sprintf("vector[%d] sqrtD_Phi;", k),
  #     "{",
  #     if (length(free_idx) > 0L) "  int c = 1;" else NULL,
  #     if (k > 0L) {
  #       paste0(
  #         vapply(seq_len(k), function(i) {
  #           if (i %in% fix_idx) {
  #             sprintf("  sqrtD_Phi[%d] = 1;", i)
  #           } else {
  #             sprintf("  sqrtD_Phi[%d] = sqrtD_Phi_free[c]; c += 1;", i)
  #           }
  #         }, character(1L)),
  #         collapse = "\n"
  #       )
  #     } else NULL,
  #     "}"
  #   )
  #
  #   tparLSigma <- sprintf("matrix[%d, %d] L_Sigma = diag_pre_multiply(sqrtD_Phi, L_Omega);", k, k)
  #
  #   # Row-vector mean (0 by default). Change if you want non-zero LV means.
  #   tparMuXi <- sprintf("row_vector[%d] MU_XI = rep_row_vector(0, %d);", k, k)
  #
  #   # Build XI_Matrix in one vectorised statement:
  #   # each row: MU_XI + z_i * L_Sigma'
  #   tparXiMatDecl <- sprintf("matrix[N, %d] XI_Matrix;", k)
  #   tparXiMatFill <- "XI_Matrix = rep_matrix(MU_XI, N) + z_XI * L_Sigma';"
  #
  #   # Also expose named column vectors (X, Z, ...) for downstream code
  #   xiVectors <- NULL
  #   for (i in seq_along(xis)) {
  #     xiVectors <- c(xiVectors, sprintf("vector[N] %s = col(XI_Matrix, %d);", xis[[i]], i))
  #   }
  #   tparXiVectors <- collapse(xiVectors)
  #
  #   # ---------------- model ----------------
  #   modLOmega <- "L_Omega ~ lkj_corr_cholesky(1);"
  #   # Non-centered prior: standard normals (vectorised via to_vector)
  #   modZ      <- "to_vector(z_XI) ~ std_normal();"
  #   # (Optionally add a weak prior for sqrtD_Phi_free if desired.)
  #
  #   parameters  <- collapse(c(parLOmega, parSqrtDPhi, parZ))
  #   transformed <- collapse(c(
  #     tparBuildSqrt, tparLSigma, tparMuXi,
  #     tparXiMatDecl, tparXiMatFill,
  #     tparXiVectors
  #   ))
  #   model <- collapse(c(modLOmega, modZ))
  #
  #   # Hide internals if you use an exclusion mechanism
  #   EXCLUDE.PARS <<- c(EXCLUDE.PARS, "z_XI", "XI_Matrix", xis)
  #
  #   list(
  #     parameters = parameters,
  #     transformed_parameters = transformed,
  #     model = model
  #   )
  # }

  STAN_PAR_ETA <- function(eta) {
    indeps <- unique(parTable[parTable$lhs == eta & parTable$op == "~", "rhs"])
    indeps <- stringr::str_replace_all(indeps, ":", "__XWITH__")

    labBeta <- sprintf("%s__REGRESSION__%s", eta, indeps)
    allOrd <- is_allOrdinal_lv[eta]

    parBeta   <- if (length(labBeta)) sprintf("real %s;", labBeta) else NULL
    parValues <- sprintf("vector[N] %s;", eta)

    tparLines <- NULL

    if (allOrd) {
      # FIXED disturbance SD = 1
      parSD <- NULL
      projEta <- if (length(indeps)) collapse(sprintf("%s * %s", labBeta, indeps), sep = " + ") else "0"
      modEta <- sprintf("%s ~ normal(%s, 1);", eta, projEta)

      # Expose the fixed residual variance under the usual name
      tparLines <- sprintf("real %s__COVARIANCE__%s = 1;", eta, eta)

    } else {
      # Free disturbance SD (parameter remains as before)
      labSD <- sprintf("%s__STDDEV__%s", eta, eta)
      parSD <- sprintf("real<lower=0> %s;", labSD)
      tparLines <- sprintf("real %s__COVARIANCE__%s = %s^2;", eta, eta, labSD)
      projEta <- if (length(indeps)) collapse(sprintf("%s * %s", labBeta, indeps), sep = " + ") else "0"
      modEta <- sprintf("%s ~ normal(%s, %s);", eta, projEta, labSD)
    }

    parameters  <- collapse(c(parValues, parSD, parBeta))
    EXCLUDE.PARS <<- c(EXCLUDE.PARS, eta)

    list(
         parameters = parameters,
         model = modEta,
         transformed_parameters = collapse(tparLines)
    )
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
                   exclude.pars = unique(EXCLUDE.PARS),
                   parTable = parTable))
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

    stan_data[[sprintf("ORD_INDICATOR_%s", ind)]] <- x_int
    stan_data[[sprintf("K_%s", ind)]]             <- K
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


getModParTable <- function(lhs, op, rhs, parTable, .default = NA) {
  parTable[parTable$mod == "", "mod"] <- .default

  if (op == "~~") {
    out <- parTable[((parTable$lhs == lhs & parTable$rhs == rhs) |
                     (parTable$lhs == rhs & parTable$rhs == lhs)) &
                    parTable$op == op, "mod"]
  } else {
    out <- parTable[parTable$lhs == lhs & parTable$rhs == rhs &
                    parTable$op == op, "mod"]
  }

  if (!length(out)) .default else out
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
