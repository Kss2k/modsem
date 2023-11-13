
#' Interaction between latent variables
#'
#' @param modelSyntax lavaan syntax
#' @param data dataframe
#' @param method method to use:
#' "rca" = residual centering approach (passed to lavaan),
#' "uca" = unconstrained approach (passed to lavaan),
#' "dblcent" = double centering approach (passed to lavaan),
#' "pind" = product indicator approach (passed to lavaan),
#' "lms" = laten model structural equations (passed to nlsem),
#' "custom" = use parameters specified in the function call (passed to lavaan)
#' @param standardizeData should data be scaled before fitting model
#' @param isMeasurementSpecified have you specified the measurement model for the latent product
#' @param firstLoadingFixed Sould the first factorloading in the latent product be fixed to one?
#' @param centerBefore should indicators in products be centered before computing products (overwritten by method, if method != NULL)
#' @param centerAfter should indicator products be centered after they have been computed?
#' @param residualsProducts should indicator products be centered using residuals (overwritten by method, if method != NULL)
#' @param residualCovSyntax should syntax for residual covariances be produced (overwritten by method, if method != NULL)
#' @param constrainedProductMean should syntax product mean be produced (overwritten by method, if method != NULL)
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
                   isMeasurementSpecified = FALSE,
                   firstLoadingFixed = TRUE,
                   centerBefore = FALSE,
                   centerAfter = FALSE,
                   residualsProducts = FALSE,
                   residualCovSyntax = FALSE,
                   constrainedProductMean = FALSE,
                   qml = FALSE,
                   ...) {
  if (is.null(modelSyntax)) {
    stop2("No model syntax provided in modsem")
  }
  if (is.null(data)) {
    stop2("No data provided in modsem")
  }

  # PreSteps -------------------------------------------------------------------

  LavDataToBeModified <<- data
  ## Standardizing data
  if (standardizeData == TRUE) {
    LavDataToBeModified <- lapplyDf(LavDataToBeModified,
                     FUN = scaleIfNumeric,
                     scaleFactor = FALSE)
  }

    ## Centering Data (should not be paired with standardize data)
  if (centerData == TRUE) {
    LavDataToBeModified <- lapplyDf(LavDataToBeModified,
                     FUN = function(x) x - mean(x))
  }


  # parcelInfo <- getParcelInfo(modelSyntax)
  #
  # if (!is.null(parcelInfo)) {
  #   if (isMeasurementSpecified != TRUE) {
  #     warning2("isMeasurementSpecified should be TRUE, when using parceling \n")
  #   }
  #   modelSyntax <- fixParcelSyntax(modelSyntax, parcelInfo = parcelInfo)
  #   data <- computeParcels(parcelInfo, data = data)
  # }

    # Get the specifications of the model --------------------------------------
  modelSpecification <-
    parseLavaan(modelSyntax, isMeasurementSpecified = isMeasurementSpecified)

    # Setting parameters according to method -----------------------------------
  switch(method,
    "rca" = {
      centerBefore <- FALSE
      residualsProducts <- TRUE
      residualCovSyntax <- TRUE
      constrainedProductMean <- FALSE
    },
    "uca" = {
      centerBefore <- TRUE
      residualsProducts <- FALSE
      residualCovSyntax <- TRUE
      constrainedProductMean <- TRUE
    },
    "pind" = {
      centerBefore <- FALSE
      residualsProducts <- FALSE
      residualCovSyntax <- FALSE
      constrainedProductMean <- FALSE
    },
    "dblcent" = {
      centerBefore <- TRUE
      centerAfter <- TRUE
      residualsProducts <- FALSE
      residualCovSyntax <- TRUE
    },
    "lms" = {
      # If method is LMS we pass it own to its own version of modsem()
      if (isMeasurementSpecified == TRUE) {
        stop2("Measurement for products should not be specified when using LMS")
      }
      LMS <- modsem.LMS(modelSpecification,
                        data = data,
                        qml = qml,
                        standardizeData = standardizeData)
      return(LMS)
    },
    "custom" = {
      warning2(
        'Method was "custom", using specifications set indside the function call'
        )
    },
    stop2("Unkown method in modsem, have you made a typo?")
    )

  # ModSEM-algorithm for product indicator based approaches --------------------

  # Calculating productinidicators based on method specifications --------------
  productIndicators <-
    createProductIndicators(modelSpecification,
                            data = LavDataToBeModified,
                            centerBefore = centerBefore,
                            centerAfter = centerAfter,
                            residualsProducts = residualsProducts)

  # Merging productIndicators into a single dataset ----------------------------
    # Old solution (does not handle/warn about duplicates)
      # mergedProductIndicators <- purrr::list_cbind(unname(productIndicators))
  mergedProductIndicators <- combineListDf(productIndicators)

  # Creating a new dataset with the productindicators
  newData <- cbind.data.frame(LavDataToBeModified, mergedProductIndicators)

  # Genereating a new syntax with constraints and measurmentmodel --------------

  newSyntax <- generateSyntax(modelSpecification,
                              isMeasurementSpecified = isMeasurementSpecified,
                              residualCovSyntax = residualCovSyntax,
                              constrainedProductMean = constrainedProductMean,
                              firstFixed = firstLoadingFixed)

  # Estimating the model via lavaan::sem()
  lavaanEstimation <- lavaan::sem(newSyntax, newData, ...)

  # Adding a bunch of stuff to the model, befur i return it
  modelSpecification$productIndicators <- productIndicators
  modelSpecification$newSyntax <- newSyntax
  modelSpecification$newData <- newData
  modelSpecification$lavaan <- lavaanEstimation
  # this is not pretty either
  structure(modelSpecification,
            class = "ModSEM",
            method = method)

}




createProductIndicators <- function(modelSpecification,
                                    data,
                                    centerBefore = FALSE,
                                    centerAfter = FALSE,
                                    residualsProducts = FALSE) {

  indicatorProducts <- purrr::map2(.x = modelSpecification$relationalDfs,
                                   .y = modelSpecification$indicatorsInProductTerms,
                                   .f = createIndicatorProducts,
                                   data = data,
                                   centered = centerBefore)
  if (residualsProducts == TRUE) {
    indicatorProducts <-
      purrr::map2(.x = indicatorProducts,
                  .y = modelSpecification$indicatorsInProductTerms,
                  .f = calculateResidualsDf,
                  data = data)

  } else if (!is.logical(residualsProducts)) {
    stop2("residualProducts was neither FALSE nor TRUE in createProductIndicators")
  }
  # new additon
  if (centerAfter == TRUE) {
    indicatorProducts <-
      lapply(indicatorProducts,
             FUN = function(df)
               lapplyDf(df,
                        FUN = function(x) x - mean(x)))

  }

  indicatorProducts
}



createIndicatorProducts <- function(relationDf,
                                    indicatorNames,
                                    data,
                                    centered = FALSE) {

  # Getting the indicatorProduct names
  varnames <- unname(colnames(relationDf))

  # Selecting the indicators from the dataset
  indicators <- data[indicatorNames]
  # Check if indicators are numeric
  isNumeric <- sapply(indicators, is.numeric)

  if (sum(as.integer(!isNumeric)) > 0) {
    stop2("Expected indicators to be numeric when creating products")
  }

  # Centering them, if center == TRUE
  if (centered == TRUE) {
    indicators <- lapplyDf(indicators,
                           FUN = function(x) x - mean(x))

  }

  products <-
    lapplyNamed(varnames,
                FUN = function(varname, data, relationDf)
                  multiplyIndicatorsCpp(data[relationDf[[varname]]]),
                data = indicators,
                relationDf = relationDf,
                names = varnames)

  # return as data.frame()
  structure(products,
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



generateSyntax <- function(modelSpecification,
                           isMeasurementSpecified = FALSE,
                           residualCovSyntax = FALSE,
                           constrainedProductMean = FALSE,
                           firstFixed = TRUE) {

  modelSyntax <- modelSpecification$lavaanSyntax
  relationalDfs <- modelSpecification$relationalDfs
  parTable <- modelSpecification$parTable
  productNames <- names(relationalDfs)

  # Measurement model latent products ------------------------------------------
  if (!is.logical(residualCovSyntax)) {
    stop2("residualCovSyntax is not FALSE or TRUE in generateSyntax")

  } else if (isMeasurementSpecified == FALSE) {
    measurementParTable <-
      purrr::map2(.x = modelSpecification$productNamesLatents,
                      .y = modelSpecification$indicatorProductNamesLatents,
                      .f = getParTableMeasurement,
                      operator = "=~",
                      firstFixed = firstFixed) |>
      purrr::list_rbind()
    parTable <- rbind(parTable, measurementParTable)
  }

  # Residual covariances -------------------------------------------------------
  if (!is.logical(residualCovSyntax)) {
    stop2("residualCovSyntax is not FALSE or TRUE in generateSyntax")

  } else if (residualCovSyntax == TRUE) {
    residualCovariances <- purrr::map2(.x = relationalDfs,
                               .y = productNames,
                               .f = getParTableResCov) |>
      purrr::list_rbind()
    parTable <- rbind(parTable, residualCovariances)
  }

  # Constrained product mean syntax --------------------------------------------
  if (constrainedProductMean == TRUE) {
    restrictedMeans <- purrr::map2(modelSpecification$productNames,
                modelSpecification$elementsInProductNames,
                getParTableRestrictedMean) |>
      purrr::list_rbind()

    parTable <- rbind(parTable, restrictedMeans)
  }

  parTableToSyntax(parTable, removeColon = TRUE)
}



# this function assumes a product of only two latent variables no more
getParTableRestrictedMean <- function(productName, elementsInProductName) {
  if (length(elementsInProductName) > 2) {
    warning2("The mean of a latent product should not be constrained when there",
             " are more than two variables in the product term. Please use a",
             " different method \n")
  }
  label <- paste0("COV_", productName)
  covariance <- createParTableRow(vecLhsRhs = elementsInProductName[1:2],
                                  op = "~~",
                                  mod = paste0("COV_", productName))
  meanStructure <- createParTableRow(vecLhsRhs = c(productName, "1"),
                                     op = "~",
                                     mod = label)
  rbind(covariance, meanStructure)
}




multiplyIndicators <- function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  if (ncol(df) <= 1){
    return(df[[1]])
  }

  y <- cbind.data.frame(df[[1]] * df[[2]],
                        df[,-(1:2),drop = FALSE])


  multiplyIndicators(y)
}




# specify residual covariance in lavaanify parTable-format
getParTableResCov <- function(relationalDf, productName) {
  if (ncol(relationalDf) <= 1) {
    return(NULL)
  }
  productNames <- colnames(relationalDf)
  uniqueCombinations <- getUniqueCombinations(productNames)
  # Now we want to specify the covariance based on shared indicators
  isShared <- vector("logical", length = nrow(uniqueCombinations))

  for (i in 1:nrow(uniqueCombinations)) {
    indicatorsProduct1 <- unlist(relationalDf[uniqueCombinations[i, "V1"]])
    indicatorsProduct2 <- unlist(relationalDf[uniqueCombinations[i, "V2"]])
    # Compare the Indicators in product1 and product2, and convert to integer
    sharedValues <- as.integer(indicatorsProduct1 %in% indicatorsProduct2)
    # Sum the values
    numberShared <- sum(sharedValues)

    if (numberShared >= 1) {
      isShared[[i]] <- TRUE
    } else if (numberShared == 0) {
      isShared[[i]] <- FALSE
    }
  }

  # Syntax for oblique covariances
  syntaxOblique <- apply(uniqueCombinations[isShared, c("V1", "V2")],
                         MARGIN = 1,
                          FUN = createParTableRow,
                          op = "~~") |>
    purrr::list_rbind()


  syntaxOrthogonal <- apply(uniqueCombinations[!isShared, c("V1", "V2")],
                            MARGIN = 1,
                            FUN = createParTableRow,
                            op = "~~",
                            mod = "0") |>
    purrr::list_rbind()

  rbind(syntaxOrthogonal, syntaxOblique)
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



getParTableMeasurement <- function(dependentName,
                                        predictorNames,
                                        operator = "=~",
                                        firstFixed = FALSE) {

  nRows <- length(predictorNames)
  parTable <- data.frame(lhs = rep(dependentName, nRows),
                         op = rep(operator, nRows),
                         rhs = predictorNames,
                         mod = vector("character", nRows))
  if (firstFixed == TRUE) {
    parTable[["mod"]][[1]] <- "1"
  }

  parTable
}



createParTableRow <- function(vecLhsRhs, op, mod = "") {
  data.frame(lhs = vecLhsRhs[[1]], op = op, rhs = vecLhsRhs[[2]], mod = mod)
}


ModSEM <- setClass("ModSEM")



#' summary.ModSEM
#'
#' @param object modsem object
#' @rdname summary
#' @export
summary.ModSEM <- function(object, ...) {
  cat("ModSEM: \nMethod =", attributes(object)$method, "\n")
  lavaan::summary(object$lavaan)
}


#' summary.ModSEM
#'
#' @param modelSyntax
#' @rdname summary modsem object
#' @export
setMethod("summary", "ModSEM", summary.ModSEM)



