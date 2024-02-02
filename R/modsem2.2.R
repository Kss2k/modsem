
#' Interaction between latent variables
#'
#' @param modelSyntax lavaan syntax
#' @param data dataframe
#' @param method method to use:
#' "rca" = residual centering approach (passed to lavaan),
#' "uca" = unconstrained approach (passed to lavaan),
#' "dblcent" = double centering approach (passed to lavaan),
#' "pind" = prod ind approach (passed to lavaan),
#' "lms" = laten model structural equations (passed to nlsem),
#' "custom" = use parameters specified in the function call (passed to lavaan)
#' @param standardizeData should data be scaled before fitting model
#' @param isMeasureSpecified have you specified the measure model for the latent prod
#' @param firstLoadingFixed Sould the first factorloading in the latent prod be fixed to one?
#' @param centerBefore should inds in prods be centered before computing prods (overwritten by method, if method != NULL)
#' @param centerAfter should ind prods be centered after they have been computed?
#' @param residualsProds should ind prods be centered using residuals (overwritten by method, if method != NULL)
#' @param residualCovSyntax should syntax for residual covariances be produced (overwritten by method, if method != NULL)
#' @param constrainedProdMean should syntax prod mean be produced (overwritten by method, if method != NULL)
#' @param ... arguments passed to other functions, e.g,. lavaan
#'
#' @return
#' @export
#'
#' @examples
#' library(modsem)
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 +x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'
#'   # Inner model
#'   Y ~ X + Z + X:Z
#''
#' est1 <- modsem(m1, oneInt)
#' summary(est1)

modsem <- function(modelSyntax = NULL,
                   data = NULL,
                   method = "rca",
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
                   qml = FALSE,
                   auto.scale = "none",
                   auto.center = "none",
                   estimator = "MLM", # We shold correct for non-normal Y
                   ...) {
  if (is.null(modelSyntax)) {
    stop("No model syntax provided in modsem")
  }
  if (is.null(data)) {
    stop("No data provided in modsem")
  }
  # PreSteps -------------------------------------------------------------------
  modEnv$data <- data
  ## Standardizing data
  if (standardizeData == TRUE || method %in% auto.scale) {
    modEnv$data <- lapplyDf(modEnv$data,
                     FUN = scaleIfNumeric,
                     scaleFactor = FALSE)
  }
    ## Centering Data (should not be paired with standardize data)
  if (centerData == TRUE || method %in% auto.center) {
    modEnv$data <- lapplyDf(modEnv$data,
                     FUN = function(x) x - mean(x))
  }

    # Get the specifications of the model --------------------------------------
  modelSpec <- parseLavaan(modelSyntax, colnames(data))
  
    # Setting parameters according to method -----------------------------------
  if (method == "lms") {
    # If method is LMS we pass it own to its own version of modsem()
    LMS <- modsem.LMS(modelSpec,
                      data = modEnv$data,
                      qml = qml,
                      standardizeData = standardizeData)
    return(LMS)

  } else if (method == "mplus") {
    mplus <- modsem.mplus(modelSyntax,
                          data = modEnv$data)
    return(mplus)
  }

  methodSettings <- getMethodSettings(method, args = list(centerBefore = centerBefore,
                                     centerAfter = centerAfter,
                                     residualsProds = residualsProds,
                                     residualCovSyntax = residualCovSyntax,
                                     constrainedProdMean = constrainedProdMean,
                                     constrainedLoadings = constrainedLoadings,
                                     constrainedVar = constrainedVar,
                                     constrainedResCovMethod = constrainedResCovMethod,
                                     firstLoadingFixed = firstLoadingFixed))

  # ModSEM-algorithm for prod ind based approaches --------------------

  # Calculating prodinidicators based on method specifications --------------
  prodInds <-
    createProdInds(modelSpec,
                            data = modEnv$data,
                            centerBefore = methodSettings$centerBefore,
                            centerAfter = methodSettings$centerAfter,
                            residualsProds = methodSettings$residualsProds)
  # Merging prodInds into a single dataset ----------------------------
    # Old solution (does not handle/warn about duplicates)
      # mergedProdInds <- purrr::list_cbind(unname(prodInds))
  mergedProdInds <- combineListDf(prodInds)

  # Creating a new dataset with the prodinds
    # using list_cbind so that mergedProdInds can be NULL
  newData <- purrr::list_cbind(list(modEnv$data, mergedProdInds))
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

  # Estimating the model via lavaan::sem()
  lavaanEstimation <- lavaan::sem(newSyntax, newData, estimator = estimator, ...)
  coefParTable <- lavaan::parameterEstimates(lavaanEstimation)
  # Adding a bunch of stuff to the model, before i return it
  modelSpec$prodInds <- prodInds
  modelSpec$newSyntax <- newSyntax
  modelSpec$newData <- newData
  modelSpec$lavaan <- lavaanEstimation
  modelSpec$parTable <- parTable
  modelSpec$coefParTable <- coefParTable
  # this is not pretty either
  structure(modelSpec,
            class = "ModSEM",
            method = method)

}




createProdInds <- function(modelSpec,
                                    data,
                                    centerBefore = FALSE,
                                    centerAfter = FALSE,
                                    residualsProds = FALSE) {

  indProds <- purrr::map2(.x = modelSpec$relDfs,
                                   .y = modelSpec$indsInProdTerms,
                                   .f = createIndProds,
                                   data = data,
                                   centered = centerBefore)
  if (residualsProds == TRUE) {
    indProds <-
      purrr::map2(.x = indProds,
                  .y = modelSpec$indsInProdTerms,
                  .f = calculateResidualsDf,
                  data = data)

  } else if (!is.logical(residualsProds)) {
    stop("residualProds was neither FALSE nor TRUE in createProdInds")
  }
  # new additon
  if (centerAfter == TRUE) {
    indProds <-
      lapply(indProds,
             FUN = function(df)
               lapplyDf(df,
                        FUN = function(x) x - mean(x)))

  }

  indProds
}



createIndProds <- function(relationDf,
                                    indNames,
                                    data,
                                    centered = FALSE) {

  # Getting the indProd names
  varnames <- unname(colnames(relationDf))

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
                FUN = function(varname, data, relationDf)
                  multiplyIndicatorsCpp(data[relationDf[[varname]]]),
                data = inds,
                relationDf = relationDf,
                names = varnames)

  # return as data.frame()
  structure(prods,
            row.names = 1:nrow(data),
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

      res <- as.data.frame(residuals(lm(formula = formula, combinedData)))
      colnames(res) <- dependentNames
      return(res)
  }
  residuals(lm(formula = formula, combinedData))

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
  unspecifiedLatentProds <- modelSpec$unspecifiedLatentProds
  indProdNamesNonSpec <- modelSpec$indProdNamesUnspecifiedLatents
  parTable <- modelSpec$parTable

  if (is.null(relDfs) || length(relDfs) < 1) {
    return(parTable)
  }

  # Measure model latent prods ------------------------------------------
  measureParTable <-
    purrr::map2(.x = unspecifiedLatentProds,
                .y = indProdNamesNonSpec,
                .f = getParTableMeasure,
                operator = "=~",
                firstFixed = firstFixed) |>
    purrr::list_rbind()

  parTable <- rbindParTable(parTable, measureParTable)
  # Residual covariances -------------------------------------------------------
  if (!is.logical(residualCovSyntax)) {
    stop("residualCovSyntax is not FALSE or TRUE in generateSyntax")

  } else if (residualCovSyntax) {
    residualCovariances <- purrr::map(.x = relDfs,
                                      .f = getParTableResCov,
                                      method = constrainedResCovMethod, 
                                      ...)  |>
      purrr::list_rbind()
    parTable <- rbindParTable(parTable, residualCovariances)
  }

  # Constrained Vars and Covs --------------------------------------------------
  if (constrainedVar == TRUE) {
    parTable <- labelVarCov(parTable, relDfs)
    parTable <- specifyVarCov(parTable, relDfs)
  }

  # Constrained Factor loadings ------------------------------------------------
  if (constrainedLoadings) {
    if (!constrainedVar) {
      parTable <- labelVarCov(parTable, relDfs)
    }
    parTable <- labelFactorLoadings(parTable)
    parTable <- specifyFactorLoadings(parTable, relDfs)
  }

  # Constrained prod mean syntax --------------------------------------------
  if (constrainedProdMean) {
    restrictedMeans <- purrr::map2(modelSpec$prodNames,
                modelSpec$elementsInProdNames,
                getParTableRestrictedMean,
                createLabels = !constrainedVar) |>
      purrr::list_rbind()
    parTable <- rbindParTable(parTable, restrictedMeans)
  }

  modEnv$parTable <- parTable
  parTable
}



# this function assumes a prod of only two latent variables no more
getParTableRestrictedMean <- function(prodName, elementsInProdName, createLabels = TRUE) {
  if (length(elementsInProdName) > 2) {
    stop("The mean of a latent prod should not be constrained when there",
             " are more than two variables in the prod term. Please use a",
             " different method \n")
  }
  label <- createLabelCov(elementsInProdName[[1]], elementsInProdName[[2]])

  if (createLabels) {
    covariance <- createParTableRow(vecLhsRhs = elementsInProdName[1:2],
                                    op = "~~",
                                    mod = label)
  } else if (!createLabels) {
    covariance <- NULL
  }
  meanStructure <- createParTableRow(vecLhsRhs = c(prodName, "1"),
                                     op = "~",
                                     mod = label)
  rbind(covariance, meanStructure)
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

  nRows <- length(predictorNames)
  parTable <- data.frame(lhs = rep(dependentName, nRows),
                         op = rep(operator, nRows),
                         rhs = predictorNames,
                         mod = vector("character", nRows))
  if (firstFixed) {
    parTable[["mod"]][[1]] <- "1"
  }

  parTable
}

ModSEM <- setClass("ModSEM")


createParTableRow <- function(vecLhsRhs, op, mod = "") {
  data.frame(lhs = vecLhsRhs[[1]], op = op, rhs = vecLhsRhs[[2]], mod = mod)
}


#' summary.ModSEM
#'
#' @param object modsem object
#' @rdname summary
#' @export
summary.ModSEM <- function(object, ...) {
  cat("ModSEM: \nMethod =", attributes(object)$method, "\n")
  if (attributes(object)$method == "Mplus") {
    object$coefParTable

  } else lavaan::summary(object$lavaan)
}


#' summary.ModSEM
#'
#' @param modelSyntax
#' @rdname summary modsem object
#' @export
setMethod("summary", "ModSEM", summary.ModSEM)



