#' Interaction between latent variables using product indicators
#'
#' @param modelSyntax lavaan syntax
#' @param data dataframe
#' @param method method to use:
#' "rca" = residual centering approach (passed to lavaan),
#' "uca" = unconstrained approach (passed to lavaan),
#' "dblcent" = double centering approach (passed to lavaan),
#' "pind" = prod ind approach, with no constraints or centering (passed to lavaan),
#' "custom" = use parameters specified in the function call (passed to lavaan)
#' @param match should the product indicators be created by using the match-strategy
#' @param standardizeData should data be scaled before fitting model
#' @param firstLoadingFixed Sould the first factorloading in the latent prod be fixed to one?
#' @param centerBefore should inds in prods be centered before computing prods (overwritten by method, if method != NULL)
#' @param centerAfter should ind prods be centered after they have been computed?
#' @param residualsProds should ind prods be centered using residuals (overwritten by method, if method != NULL)
#' @param residualCovSyntax should syntax for residual covariances be produced (overwritten by method, if method != NULL)
#' @param constrainedProdMean should syntax prod mean be produced (overwritten by method, if method != NULL)
#' @param centerData should data be centered before fitting model
#' @param constrainedLoadings should syntax for constrained loadings be produced (overwritten by method, if method != NULL)
#' @param constrainedVar should syntax for constrained variances be produced (overwritten by method, if method != NULL)
#' @param constrainedResCovMethod method for constraining residual covariances
#' @param auto.scale methods which should be scaled automatically (usually not useful)
#' @param auto.center methods which should be centered automatically (usually not useful)
#' @param estimator estimator to use in lavaan
#' @param group group variable for multigroup analysis
#' @param run should the model be run via lavaan, if FALSE only modified syntax and data is returned
#' @param ... arguments passed to other functions, e.g,. lavaan
#' @return ModSEM object
#' @export 
#' @description
#' modsem_pi is a function for estimating interaction effects between latent variables, 
#' in structural equation models (SEMs), using product indicators.
#' Methods for estimating interaction effects in SEM's can basically be split into 
#' two frameworks: 1. Product Indicator based approaches ("dblcent", "rca", "uca", 
#' "ca", "pind"), and 2. Distributionally based approaches ("lms", "qml").
#' modsem_pi() is essentially a just 
#' a fancy wrapper for lavaan::sem()  which generates the 
#' necessary syntax, and variables for the estimation of models with latent product indicators.
#' @examples
#' library(modsem)
#' # For more examples check README and/or GitHub.
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
#'   LATT =~ att1 + att2 + att3 + att4 + att5
#'   LSN =~ sn1 + sn2
#'   LPBC =~ pbc1 + pbc2 + pbc3
#'   LINT =~ int1 + int2 + int3
#'   LBEH =~ b1 + b2
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   # Covariances
#'   LATT ~~ LSN + LPBC
#'   LPBC ~~ LSN 
#'   # Causal Relationsships
#'   LINT ~ LATT + LSN + LPBC
#'   LBEH ~ LINT + LPBC 
#'   LBEH ~ LINT:LPBC  
#' '
#' 
#' # double centering approach
#' estTpb <- modsem_pi(tpb, data = TPB)
#' summary(estTpb)
#'
#' \dontrun{
#' # The Constrained Approach 
#' estTpbConstrained <- modsem_pi(tpb, data = TPB, method = "ca")
#' summary(estTpbConstrained)
#' }
modsem_pi <- function(modelSyntax = NULL,
                   data = NULL,
                   method = "dblcent",
                   match = FALSE,
                   standardizeData = FALSE,
                   centerData = FALSE,
                   firstLoadingFixed = TRUE,
                   centerBefore = NULL,
                   centerAfter = NULL,
                   residualsProds = NULL,
                   residualCovSyntax = NULL,
                   constrainedProdMean = NULL,
                   constrainedLoadings = NULL,
                   constrainedVar = NULL,
                   constrainedResCovMethod = NULL,
                   auto.scale = "none",
                   auto.center = "none",
                   estimator = "ML", 
                   group = NULL,
                   run = TRUE, 
                   ...) {
  if (is.null(modelSyntax)) stop("No model syntax provided in modsem")
  if (is.null(data)) stop("No data provided in modsem")
  if (!is.data.frame(data)) data <- as.data.frame(data)
  
  # Get the specifications of the model 
  modelSpec <- parseLavaan(modelSyntax, colnames(data), match = match)

  # Data Processing  -----------------------------------------------------------
  data <- data[c(modelSpec$oVs, group)]
  completeCases <- stats::complete.cases(data)
  if (any(!completeCases)) {
    warning("Removing missing values case-wise.")
    data <- data[completeCases, ]
  }

  ## Standardizing data
  if (standardizeData || method %in% auto.scale) {
    data <- lapplyDf(data, FUN = scaleIfNumeric, scaleFactor = FALSE)
  }
  ## Centering Data (should not be paired with standardize data)
  if (centerData || method %in% auto.center) {
    data <- lapplyDf(data, FUN = function(x) x - mean(x))
  }

  methodSettings <-
    getMethodSettings(method,
                      args =
                        list(centerBefore = centerBefore,
                             centerAfter = centerAfter,
                             residualsProds = residualsProds,
                             residualCovSyntax = residualCovSyntax,
                             constrainedProdMean = constrainedProdMean,
                             constrainedLoadings = constrainedLoadings,
                             constrainedVar = constrainedVar,
                             constrainedResCovMethod = constrainedResCovMethod,
                             firstLoadingFixed = firstLoadingFixed))

  # ModSEM-algorithm for prod ind based approaches -----------------------------
  prodInds <-
    createProdInds(modelSpec,
                   data = data,
                   centerBefore = methodSettings$centerBefore,
                   centerAfter = methodSettings$centerAfter,
                   residualsProds = methodSettings$residualsProds)
  mergedProdInds <- combineListDf(prodInds)

  # using list_cbind so that mergedProdInds can be NULL
  newData <- purrr::list_cbind(list(data, mergedProdInds))
  # Genereating a new syntax with constraints and measurmentmodel --------------
  parTable <- addSpecsParTable(modelSpec,
                               residualCovSyntax = methodSettings$residualCovSyntax,
                               constrainedResCovMethod = methodSettings$constrainedResCovMethod,
                               constrainedProdMean = methodSettings$constrainedProdMean,
                               constrainedLoadings = methodSettings$constrainedLoadings,
                               constrainedVar = methodSettings$constrainedVar,
                               firstFixed = firstLoadingFixed,
                               ...)

  newSyntax <- parTableToSyntax(parTable, removeColon = TRUE)

  modelSpec$prodInds <- prodInds
  modelSpec$syntax <- newSyntax
  modelSpec$data <- newData
  modelSpec$parTable <- parTable

  # Estimating the model via lavaan::sem() 
  if (run) {
    lavEst <- tryCatch(lavaan::sem(newSyntax, newData, estimator = estimator, 
                                   group = group, ...),
                       error = function(cnd) {
                         warning("Error in Lavaan: \n")
                         warning(capturePrint(cnd))
                         NULL
                       })
    coefParTable <- tryCatch(lavaan::parameterEstimates(lavEst),
                             error = function(cnd) NULL)
    modelSpec$lavaan <- lavEst
    modelSpec$coefParTable <- coefParTable
  }
  structure(modelSpec, class = "modsem_pi", method = method)
}


createProdInds <- function(modelSpec,
                           data,
                           centerBefore = FALSE,
                           centerAfter = FALSE,
                           residualsProds = FALSE) {
  
  indProds <- purrr::map2(.x = modelSpec$relDfs,
                          .y = modelSpec$indsInLatentProds,
                          .f = createIndProds,
                          data = data,
                          centered = centerBefore)
  if (residualsProds) {
    indProds <-
      purrr::map2(.x = indProds,
                  .y = modelSpec$indsInLatentProds,
                  .f = calculateResidualsDf,
                  data = data)

  } else if (!is.logical(residualsProds)) {
    stop("residualProds was neither FALSE nor TRUE in createProdInds")
  }
  if (centerAfter) {
    indProds <- lapply(indProds, FUN = function(df)
                       lapplyDf(df, FUN = function(x) x - mean(x)))

  }

  indProds
}


createIndProds <- function(relDf,
                           indNames,
                           data,
                           centered = FALSE) {

  # Getting the indProd names
  varnames <- unname(colnames(relDf))

  # Selecting the inds from the dataset
  inds <- data[indNames]
  # Check if inds are numeric
  isNumeric <- sapply(inds, is.numeric)

  if (sum(as.integer(!isNumeric)) > 0) {
    stop("Expected inds to be numeric when creating prods")
  }

  # Centering them, if center == TRUE
  if (centered == TRUE) {
    inds <- lapplyDf(inds,
                     FUN = function(x) x - mean(x))

  }

  prods <-
    lapplyNamed(varnames,
                FUN = function(varname, data, relDf)
                  multiplyIndicatorsCpp(data[relDf[[varname]]]),
                data = inds,
                relDf = relDf,
                names = varnames)

  # return as data.frame()
  structure(prods,
            row.names = seq_len(nrow(data)),
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
  formulaDep <- paste0("cbind(",
                       stringr::str_c(dependendtNames,
                                      collapse = ", "),
                       ")")
  formulaIndep <- stringr::str_c(indepNames, collapse = " + ")
  paste0(formulaDep, " ~ ", formulaIndep)
}


addSpecsParTable <- function(modelSpec,
                             residualCovSyntax = FALSE,
                             constrainedResCovMethod = "equality",
                             constrainedProdMean = FALSE,
                             constrainedLoadings = FALSE,
                             constrainedVar = FALSE,
                             firstFixed = TRUE,
                             ...) {
  relDfs <- modelSpec$relDfs
  latentProds <- modelSpec$latentProds
  indProdNames <- modelSpec$indProdNames
  parTable <- modelSpec$parTable

  if (is.null(relDfs) || length(relDfs) < 1) {
    return(parTable)
  }
  # Measure model latent prods ------------------------------------------
  measureParTable <-
    purrr::map2(.x = latentProds,
                .y = indProdNames,
                .f = getParTableMeasure,
                operator = "=~",
                firstFixed = firstFixed) |>
    purrr::list_rbind()
  parTable <- rbindParTable(parTable, measureParTable)

  # label parameters and add if necessary --------------------------------------
  if (constrainedVar || constrainedLoadings || constrainedProdMean) {
    parTable <- addVariances(parTable) |> 
      addCovariances() |>
      labelParameters() |>
      labelFactorLoadings() 
  }

  # Residual covariances -------------------------------------------------------
  if (!is.logical(residualCovSyntax)) {
    stop("residualCovSyntax is not FALSE or TRUE in generateSyntax")

  } else if (residualCovSyntax) {
    residualCovariances <- purrr::map(.x = relDfs,
                                      .f = getParTableResCov,
                                      method = constrainedResCovMethod,
                                      pt = parTable, # for caluclating formulas using path tracer
                                      ...)  |>
      purrr::list_rbind()
    parTable <- rbindParTable(parTable, residualCovariances)
  }

  # Constrained Vars and Covs --------------------------------------------------
  if (constrainedVar) {
    parTable <- specifyVarCov(parTable, relDfs)
  }

  # Constrained Factor loadings ------------------------------------------------
  if (constrainedLoadings) {
    parTable <- specifyFactorLoadings(parTable, relDfs)
  }

  # Constrained prod mean syntax -----------------------------------------------
  if (constrainedProdMean) {
    restrictedMeans <- purrr::map2(modelSpec$prodNames,
                                   modelSpec$elementsInProdNames,
                                   getParTableRestrictedMean,
                                   createLabels = !constrainedVar,
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
  if (length(elementsInProdName) > 2) {
    stop("The mean of a latent prod should not be constrained when there",
         " are more than two variables in the prod term. Please use a",
         " different method \n")
  }
  #label <- createLabelCov(elementsInProdName[[1]], elementsInProdName[[2]])
  meanLabel <- createLabelMean(prodName)

  meanStructure <- createParTableRow(vecLhsRhs = c(prodName, "1"),
                                     op = "~", mod = meanLabel)
  covEquation <- tracePath(pt, elementsInProdName[[1]], elementsInProdName[[2]])
  meanFormula <- createParTableRow(vecLhsRhs = c(meanLabel, covEquation), 
                                   op = "==") 
  rbind(meanStructure, meanFormula)
}


multiplyInds <- function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  if (ncol(df) <= 1){
    return(df[[1]])
  }

  y <- cbind.data.frame(df[[1]] * df[[2]],
                        df[,-(1:2),drop = FALSE])

  multiplyInds(y)
}


#  function for getting unique combinations of two values in x
getUniqueCombinations <- function(x) {
  # Base case, x is 1 length long and there are no unique combos
  if (length(x) <= 1) {
    return(NULL)
  }
  rest <- getUniqueCombinations(x[-1])
  combos <- data.frame(V1 = rep(x[[1]], length(x) - 1),
                       V2 = x[-1])
  rbind(combos, rest)
}


getParTableMeasure <- function(dependentName,
                               predictorNames,
                               operator = "=~",
                               firstFixed = FALSE) {
  if (length(dependentName) > 1) {
    stop("Expected dependentName to be a single string in getParTableMeasure")
  }
  if (length(predictorNames) == 1 && dependentName == predictorNames) {
    # In this case it should be seen as an observed variable, 
    # and should not have a measurement model
    return(NULL)
  }
  nRows <- length(predictorNames)
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
