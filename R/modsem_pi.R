#' Interaction between latent variables using product indicators
#'
#' @param model.syntax \code{lavaan} syntax
#' 
#' @param data dataframe
#' 
#' @param method method to use:
#' \code{"rca"} = residual centering approach (passed to \code{lavaan}),
#' \code{"uca"} = unconstrained approach (passed to \code{lavaan}),
#' \code{"dblcent"} = double centering approach (passed to \code{lavaan}),
#' \code{"pind"} = prod ind approach, with no constraints or centering (passed to \code{lavaan}),
#' \code{"custom"} = use parameters specified in the function call (passed to \code{lavaan})
#' 
#' @param match should the product indicators be created by using the match-strategy
#' 
#' @param standardize.data should data be scaled before fitting model
#' 
#' @param first.loading.fixed Should the first factor loading in the latent product be fixed to one?
#' 
#' @param center.before should indicators in products be centered before computing products (overwritten by \code{method}, if \code{method != NULL})
#' 
#' @param center.after should indicator products be centered after they have been computed?
#' 
#' @param residuals.prods should indicator products be centered using residuals (overwritten by \code{method}, if \code{method != NULL})
#' 
#' @param residual.cov.syntax should syntax for residual covariances be produced (overwritten by \code{method}, if \code{method != NULL})
#' 
#' @param constrained.prod.mean should syntax for product mean be produced (overwritten by \code{method}, if \code{method != NULL})
#' 
#' @param center.data should data be centered before fitting model
#' 
#' @param constrained.loadings should syntax for constrained loadings be produced (overwritten by \code{method}, if \code{method != NULL})
#' 
#' @param constrained.var should syntax for constrained variances be produced (overwritten by \code{method}, if \code{method != NULL})
#' 
#' @param constrained.res.cov.method method for constraining residual covariances
#' 
#' @param auto.scale methods which should be scaled automatically (usually not useful)
#' 
#' @param auto.center methods which should be centered automatically (usually not useful)
#' 
#' @param estimator estimator to use in \code{lavaan}
#' 
#' @param group group variable for multigroup analysis
#' 
#' @param run should the model be run via \code{lavaan}, if \code{FALSE} only modified syntax and data is returned
#'
#' @param suppress.warnings.lavaan should warnings from \code{lavaan} be suppressed?
#' 
#' @param ... arguments passed to other functions, e.g., \code{lavaan}
#' 
#' @return \code{modsem} object
#' @export 
#' @description
#' \code{modsem_pi()} is a function for estimating interaction effects between latent variables, 
#' in structural equation models (SEMs), using product indicators.
#' Methods for estimating interaction effects in SEMs can basically be split into 
#' two frameworks: 
#' 1. Product Indicator based approaches (\code{"dblcent"}, \code{"rca"}, \code{"uca"}, 
#' \code{"ca"}, \code{"pind"}), and 
#' 2. Distributionally based approaches (\code{"lms"}, \code{"qml"}).
#' \code{modsem_pi()} is essentially a fancy wrapper for \code{lavaan::sem()} which generates the 
#' necessary syntax and variables for the estimation of models with latent product indicators.
#' Use \code{default_settings_pi()} to get the default settings for the different methods.
#' 
#' @examples
#' library(modsem)
#' # For more examples, check README and/or GitHub.
#' # One interaction
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 +x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'   
#'   # Inner model
#'   Y ~ X + Z + X:Z 
#' '
#' 
#' # Double centering approach
#' est1 <- modsem_pi(m1, oneInt)
#' summary(est1)
#' 
#' \dontrun{
#' # The Constrained Approach 
#' est1Constrained <- modsem_pi(m1, oneInt, method = "ca")
#' summary(est1Constrained)
#' }
#' 
#' # Theory Of Planned Behavior
#' tpb <- ' 
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   # Covariances
#'   ATT ~~ SN + PBC
#'   PBC ~~ SN 
#'   # Causal Relationships
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC 
#'   BEH ~ INT:PBC  
#' '
#' 
#' # Double centering approach
#' estTpb <- modsem_pi(tpb, data = TPB)
#' summary(estTpb)
#'
#' \dontrun{
#' # The Constrained Approach 
#' estTpbConstrained <- modsem_pi(tpb, data = TPB, method = "ca")
#' summary(estTpbConstrained)
#' }
modsem_pi <- function(model.syntax = NULL,
                      data = NULL,
                      method = "dblcent",
                      match = NULL,
                      standardize.data = FALSE,
                      center.data = FALSE,
                      first.loading.fixed = TRUE,
                      center.before = NULL,
                      center.after = NULL,
                      residuals.prods = NULL,
                      residual.cov.syntax = NULL,
                      constrained.prod.mean = NULL,
                      constrained.loadings = NULL,
                      constrained.var = NULL,
                      constrained.res.cov.method = NULL,
                      auto.scale = "none",
                      auto.center = "none",
                      estimator = "ML", 
                      group = NULL,
                      run = TRUE, 
                      suppress.warnings.lavaan = FALSE,
                      ...) {
  stopif(is.null(model.syntax), "No model syntax provided in modsem")
  stopif(is.null(data), "No data provided in modsem")
  if (!is.data.frame(data)) data <- as.data.frame(data)

  methodSettings <-
    getMethodSettingsPI(method, args =
                        list(center.before = center.before,
                             center.after = center.after,
                             residuals.prods = residuals.prods,
                             residual.cov.syntax = residual.cov.syntax,
                             constrained.prod.mean = constrained.prod.mean,
                             constrained.loadings = constrained.loadings,
                             constrained.var = constrained.var,
                             constrained.res.cov.method = constrained.res.cov.method,
                             first.loading.fixed = first.loading.fixed,
                             match = match))

  # Get the specifications of the model 
  modelSpec <- parseLavaan(model.syntax, colnames(data), 
                           match = methodSettings$match)

  # Data Processing  -----------------------------------------------------------
  data          <- data[c(modelSpec$oVs, group)]
  completeCases <- stats::complete.cases(data)
  if (any(!completeCases)) {
    warning2("Removing missing values case-wise.")
    data <- data[completeCases, ]
  }

  if (standardize.data || method %in% auto.scale) {
    data <- lapplyDf(data, FUN = scaleIfNumeric, scaleFactor = FALSE)
  }

  if (center.data || method %in% auto.center) {
    data <- lapplyDf(data, FUN = function(x) x - mean(x))
  }

  prodInds <-
    createProdInds(modelSpec,
                   data = data,
                   center.before = methodSettings$center.before,
                   center.after = methodSettings$center.after,
                   residuals.prods = methodSettings$residuals.prods)
  mergedProdInds <- combineListDf(prodInds)

  # using list_cbind so that mergedProdInds can be NULL
  newData <- purrr::list_cbind(list(data, mergedProdInds))

  # Genereating a new syntax with constraints and measurmentmodel
  parTable <- addSpecsParTable(modelSpec,
                               residual.cov.syntax = methodSettings$residual.cov.syntax,
                               constrained.res.cov.method = methodSettings$constrained.res.cov.method,
                               constrained.prod.mean = methodSettings$constrained.prod.mean,
                               constrained.loadings = methodSettings$constrained.loadings,
                               constrained.var = methodSettings$constrained.var,
                               firstFixed = first.loading.fixed,
                               ...)

  newSyntax <- parTableToSyntax(parTable, removeColon = TRUE)

  modelSpec$prodInds <- prodInds
  modelSpec$syntax   <- newSyntax
  modelSpec$data     <- newData
  modelSpec$parTable <- parTable

  if (run) {
    lavWrapper <- getWarningWrapper(silent = suppress.warnings.lavaan)
    lavEst <- tryCatch(lavaan::sem(newSyntax, newData, estimator = estimator, 
                                   group = group, ...) |> lavWrapper(),
                       error = function(cnd) {
                         warning2(capturePrint(cnd))
                         NULL
                       })
    coefParTable <- tryCatch(lavaan::parameterEstimates(lavEst),
                             error = function(cnd) NULL)
    modelSpec$lavaan       <- lavEst
    modelSpec$coefParTable <- coefParTable
  }

  structure(modelSpec, class = c("modsem_pi", "modsem"), method = method)
}


createProdInds <- function(modelSpec,
                           data,
                           center.before = FALSE,
                           center.after = FALSE,
                           residuals.prods = FALSE) {
  indProds <- purrr::map2(.x = modelSpec$relDfs,
                          .y = modelSpec$indsInLatentProds,
                          .f = createIndProds,
                          data = data,
                          centered = center.before)
  if (residuals.prods) {
    indProds <- purrr::map2(.x = indProds, .y = modelSpec$indsInLatentProds,
                            .f = calculateResidualsDf, data = data)

  } else if (!is.logical(residuals.prods)) {
    stop2("residualProds was neither FALSE nor TRUE in createProdInds")
  }

  if (center.after) {
    indProds <- lapply(indProds, FUN = function(df)
                       lapplyDf(df, FUN = function(x) x - mean(x)))

  }

  indProds
}


createIndProds <- function(relDf, indNames, data, centered = FALSE) {
  # Getting the indProd names
  varnames <- unname(colnames(relDf))
  # Selecting the inds from the dataset
  inds <- data[indNames]
  # Check if inds are numeric
  isNumeric <- sapply(inds, is.numeric)

  stopif(any(!isNumeric), "Expected inds to be numeric when creating prods")

  # Centering them
  if (centered) inds <- lapplyDf(inds, FUN = function(x) x - mean(x))

  prods <- lapplyNamed(varnames,
                FUN = function(varname, data, relDf)
                  multiplyIndicatorsCpp(data[relDf[[varname]]]),
                data = inds, relDf = relDf, names = varnames)

  # return as data.frame()
  structure(prods, row.names = seq_len(nrow(data)),
            class = "data.frame")
}


calculateResidualsDf <- function(dependentDf, independentNames, data) {
  # Using purrr::list_cbind() is more efficient than cbind()
  combinedData <- purrr::list_cbind(list(dependentDf, data))

  # Getting the names of the dependent variables
  dependentNames <- colnames(dependentDf)
  # Getting formula
  formula <- getResidualsFormula(dependentNames, independentNames)

  if (length(dependentNames <= 1)) {
    res <- as.data.frame(stats::residuals(stats::lm(formula = formula, 
                                                    combinedData)))
    colnames(res) <- dependentNames
    return(res)
  }

  stats::residuals(stats::lm(formula = formula, combinedData))
}


getResidualsFormula <- function(dependendtNames, indepNames) {
  formulaDep   <- paste0("cbind(", stringr::str_c(dependendtNames,
                                collapse = ", "), ")")
  formulaIndep <- stringr::str_c(indepNames, collapse = " + ")
  paste0(formulaDep, " ~ ", formulaIndep)
}


addSpecsParTable <- function(modelSpec,
                             residual.cov.syntax = FALSE,
                             constrained.res.cov.method = "equality",
                             constrained.prod.mean = FALSE,
                             constrained.loadings = FALSE,
                             constrained.var = FALSE,
                             firstFixed = TRUE,
                             ...) {
  relDfs       <- modelSpec$relDfs
  latentProds  <- modelSpec$latentProds
  indProdNames <- modelSpec$indProdNames
  parTable     <- modelSpec$parTable

  if (is.null(relDfs) || length(relDfs) < 1) return(parTable)

  measureParTable <- purrr::map2(.x = latentProds, .y = indProdNames,
                                 .f = getParTableMeasure, operator = "=~", 
                                 firstFixed = firstFixed) |>
    purrr::list_rbind()
  parTable <- rbindParTable(parTable, measureParTable)

  if (constrained.var || constrained.loadings || constrained.prod.mean) {
    parTable <- addVariances(parTable) |> 
      addCovariances() |>
      labelParameters() |>
      labelFactorLoadings() 
  }

  if (!is.logical(residual.cov.syntax)) {
    stop2("residual.cov.syntax is not FALSE or TRUE in generateSyntax")

  } else if (residual.cov.syntax) {
    residualCovariances <- purrr::map(.x = relDfs, .f = getParTableResCov,
                                      method = constrained.res.cov.method,
                                      pt = parTable, ...)  |>
      purrr::list_rbind()
    parTable <- rbindParTable(parTable, residualCovariances)
  }

  if (constrained.var) parTable <- specifyVarCov(parTable, relDfs)
  if (constrained.loadings) parTable <- specifyFactorLoadings(parTable, relDfs)

  if (constrained.prod.mean) {
    restrictedMeans <- purrr::map2(modelSpec$prodNames,
                                   modelSpec$elementsInProdNames,
                                   getParTableRestrictedMean,
                                   createLabels = !constrained.var,
                                   pt = parTable) |>
      purrr::list_rbind()
    parTable <- rbindParTable(parTable, restrictedMeans)
  }

  modEnv$parTable <- parTable
  parTable
}


# this function assumes a prod of only two latent variables no more
getParTableRestrictedMean <- function(prodName, elementsInProdName, 
                                      createLabels = TRUE, pt) {
  stopif(length(elementsInProdName) > 2,
         "The mean of a latent prod should not be constrained when there",
         " are more than two variables in the prod term. Please use a",
         " different method \n")

  meanLabel     <- createLabelMean(prodName)
  meanStructure <- createParTableRow(vecLhsRhs = c(prodName, ""),
                                     op = "~1", mod = meanLabel)
  covEquation   <- trace_path(pt, elementsInProdName[[1]], elementsInProdName[[2]])
  meanFormula   <- createParTableRow(vecLhsRhs = c(meanLabel, covEquation), 
                                     op = "==") 
  rbind(meanStructure, meanFormula)
}


multiplyInds <- function(df) {
  if (is.null(df)) return(NULL)
  if (ncol(df) <= 1) return(df[[1]])

  y <- cbind.data.frame(df[[1]] * df[[2]],
                        df[,-(1:2),drop = FALSE])

  multiplyInds(y)
}


getParTableMeasure <- function(dependentName,
                               predictorNames,
                               operator = "=~",
                               firstFixed = FALSE) {
  stopif(length(dependentName) > 1, "Expected dependentName ",
         "to be a single string in getParTableMeasure")

  if (length(predictorNames) == 1 && dependentName == predictorNames) {
    # In this case it should be seen as an observed variable, 
    # and should not have a measurement model
    return(NULL)
  }

  nRows    <- length(predictorNames)
  parTable <- data.frame(lhs = rep(dependentName, nRows),
                         op = rep(operator, nRows),
                         rhs = predictorNames,
                         mod = vector("character", nRows))

  if (firstFixed) parTable[["mod"]][[1]] <- "1"
  parTable
}


createParTableRow <- function(vecLhsRhs, op, mod = "") {
  data.frame(lhs = vecLhsRhs[[1]], op = op, rhs = vecLhsRhs[[2]], mod = mod)
}


#' Get \code{lavaan} syntax for product indicator approaches
#'
#' @param model.syntax \code{lavaan} syntax
#' 
#' @param method method to use:
#' \code{"rca"} = residual centering approach,
#' \code{"uca"} = unconstrained approach,
#' \code{"dblcent"} = double centering approach,
#' \code{"pind"} = prod ind approach, with no constraints or centering,
#' \code{"custom"} = use parameters specified in the function call
#' 
#' @param match should the product indicators be created by using the match-strategy
#' 
#' @param ... arguments passed to other functions (e.g., \link{modsem_pi})
#' 
#' @return \code{character} vector
#' @export 
#' @description
#' \code{get_pi_syntax()} is a function for creating the \code{lavaan} syntax used for estimating 
#' latent interaction models using one of the product indicators in \code{lavaan}.
#' 
#' @examples
#' library(modsem)
#' library(lavaan)
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 + x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'   
#'   # Inner model
#'   Y ~ X + Z + X:Z 
#' '
#' syntax <- get_pi_syntax(m1)
#' data <- get_pi_data(m1, oneInt)
#' est <- sem(syntax, data)
get_pi_syntax <- function(model.syntax,
                          method = "dblcent",
                          match = FALSE,
                          ...) {
  oVs       <- getOVs(model.syntax = model.syntax)
  emptyData <- as.data.frame(matrix(0, nrow = 1, ncol = length(oVs), 
                                    dimnames = list(NULL, oVs)))
  modsem_pi(model.syntax, method = method, match = match, 
            data = emptyData, run = FALSE, ...)$syntax
}


#' Get data with product indicators for different approaches
#'
#' @param model.syntax \code{lavaan} syntax
#' @param data data to create product indicators from 
#' @param method method to use:
#' \code{"rca"} = residual centering approach,
#' \code{"uca"} = unconstrained approach,
#' \code{"dblcent"} = double centering approach,
#' \code{"pind"} = prod ind approach, with no constraints or centering,
#' \code{"custom"} = use parameters specified in the function call
#' @param match should the product indicators be created by using the match-strategy
#' 
#' @param ... arguments passed to other functions (e.g., \link{modsem_pi})
#' 
#' @return \code{data.frame}
#' @export 
#' @description
#' \code{get_pi_syntax()} is a function for creating the \code{lavaan} syntax used for estimating 
#' latent interaction models using one of the product indicators in \code{lavaan}.
#' 
#' @examples
#' library(modsem)
#' library(lavaan)
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 +x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'   
#'   # Inner model
#'   Y ~ X + Z + X:Z 
#' '
#' syntax <- get_pi_syntax(m1)
#' data <- get_pi_data(m1, oneInt)
#' est <- sem(syntax, data)
get_pi_data <- function(model.syntax, data, method = "dblcent",
                        match = FALSE, ...) {
  modsem_pi(model.syntax, data = data, method = method, match = match, 
            run = FALSE, ...)$data
}



#' extract lavaan object from modsem object estimated using product indicators
#'
#' @param object modsem object
#' 
#' @return lavaan object
#' @export 
#' @examples
#' library(modsem)
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 + x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'   
#'   # Inner model
#'   Y ~ X + Z + X:Z 
#' '
#' est <- modsem_pi(m1, oneInt)
#' lav_est <- extract_lavaan(est) 
extract_lavaan <- function(object) {
  if (!inherits(object, "modsem_pi")) {
    stop2("object is not of class modsem_pi")
  }
  object$lavaan
}
