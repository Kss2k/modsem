
modsem <- function(modelSyntax = NULL,
                   data = NULL,
                   method = "rca",
                   standardizeData = FALSE,
                   isMeasurementSpecified = FALSE,
                   centeredProducts = FALSE,
                   residualsProducts = FALSE,
                   residualCovSyntax = FALSE,
                   constrainedProductMean = FALSE,
                   ...) {
  if (is.null(modelSyntax)) {
    stop2("No model syntax provided in modsem")
  }
  if (is.null(data)) {
    stop2("No data provided in modsem")
  }

  # Setting parameters according to method -------------------------------------
  if (is.null(method)) {
    warning2("Method was NULL, using specifications set indside the function call")

  } else if (method == "rca") {
    centeredProducts <- FALSE
    residualsProducts <- TRUE
    residualCovSyntax <- TRUE
    constrainedProductMean <- FALSE

  } else if (method == "unconstrained") {
    centeredProducts <- TRUE
    residualsProducts <- FALSE
    residualCovSyntax <- TRUE
    constrainedProductMean <- TRUE

  } else if (method == "regression") {
    centeredProducts <- FALSE
    residualsProducts <- FALSE
    residualCovSyntax <- FALSE
    constrainedProductMean <- FALSE

  } else {
    stop2("Unkown method in modsem, have you made a typo?")
  }

  if (standardizeData == TRUE) {
    data <- lapplyDf(data,
                    FUN = scaleIfNumeric,
                    scaleFactor = FALSE)
  }

  # Get the specifications of the model ----------------------------------------
  modelSpecification <-
    parseLavaan(modelSyntax,isMeasurementSpecified = isMeasurementSpecified)

  # Calculating productinidicators based on method specifications --------------
  productIndicators <-
    createProductIndicators(modelSpecification,
                            data = data,
                            centeredProducts = centeredProducts,
                            residualsProducts = residualsProducts)

  # Merging productIndicators into a single dataset ----------------------------
  # (NB i should add a warning if they contain the same names), and delete copies...
  mergedProductIndicators <- purrr::list_cbind(unname(productIndicators))
  # Creating a new dataset with the productindicators
  newData <- cbind.data.frame(data, mergedProductIndicators)

  # Genereating a new syntax with constraints and measurmentmodel --------------

  newSyntax <- generateSyntax(modelSpecification,
                              isMeasurementSpecified = isMeasurementSpecified,
                              residualCovSyntax = residualCovSyntax,
                              constrainedProductMean = constrainedProductMean)

  # Estimating the model via lavaan::sem()
  lavaanEstimation <- lavaan::sem(newSyntax, newData)

  # Adding a bunch of stuff to the model, befur i return it
  modelSpecification$productIndicators <- productIndicators
  modelSpecification$newSyntax <- newSyntax
  modelSpecification$newData <- newData
  modelSpecification$lavaan <- lavaanEstimation
  # this is not pretty either
  class(modelSpecification) <- "modsem"
  modelSpecification
}




createProductIndicators <- function(modelSpecification,
                                    data,
                                    centeredProducts = FALSE,
                                    residualsProducts = FALSE) {

  indicatorProducts <- purrr::map2(.x = modelSpecification$relationalDfs,
                                   .y = modelSpecification$indicatorsInProductTerms,
                                   .f = createIndicatorProducts,
                                   data = data,
                                   centered = centeredProducts)
  if (residualsProducts == TRUE) {
    purrr::map2(.x = indicatorProducts,
                .y = modelSpecification$indicatorsInProductTerms,
                .f = calculateResidualsDf,
                data = data)
  } else if (residualsProducts == FALSE) {
    indicatorProducts
  } else {
    stop2("residualProducts was neither FALSE nor TRUE in createProductIndicators")
  }

}



createIndicatorProducts <- function(relationDf,
                                    indicatorNames,
                                    data,
                                    centered = FALSE) {
  # Getting the indicatorProduct names
  varnames <- unname(colnames(relationDf))

  # Selecting the indicators from the dataset
  indicators <- data[indicatorNames]

  # Centering them, if center == TRUE
  if (centered == TRUE) {
    indicators <- apply(indicators, 2, scale, scale = FALSE) |>
      as.data.frame()
  }
  # Creating a dataframe to hold the computed indicatorproducts
  products <- data.frame(matrix(ncol = length(varnames), nrow = nrow(data)))

  # Setting the productames (e.g., var1var2 = var1*var2) and indicator names
  products <- structure(products,
                        indicatorNames = indicatorNames)
  colnames(products) <- varnames

  # Loop to create the indicatorProducts
  for (i in seq_along(varnames)) {
    varname <- varnames[[i]]
    products[[varname]] <- multiplyIndicators(indicators[relationDf[[varname]]])
  }

  products
}


# function for calculating residuals for a dataframe of productindicators
calculateResidualsDf <- function(dependentDf, independentNames, data) {

  # Do i want to explicitly coerce this?? cbind() should return a df, it
  # it's inputs are df's
  combinedData <- cbind(dependentDf, data)

  # Getting the names of the dependent variables
  dependentNames <- colnames(dependentDf)

  # calculating the residuals for each product, using the same predictors.
  residuals <- lapplyNamed(dependentNames,
                           FUN = calculateResidualsSingle,
                           independentNames = independentNames,
                           data = combinedData)


  as.data.frame(residuals)
}



calculateResidualsSingle <- function(dependentName, independentNames, data) {
  # generating a formula
  formula <- generateFormula.lm(dependentName, independentNames, operator = "~")

  residuals(lm(formula, data))
}




generateFormula.lm <- function(dependentName,
                               predictorNames,
                               operator = "~",
                               firstFixed = FALSE) {

  predictors <- stringr::str_c(predictorNames, collapse = " + ")

  if (firstFixed == TRUE) {
    predictors <- paste0("1*", predictors)
  }
  paste(dependentName, operator, predictors)

}



generateSyntax <- function(modelSpecification,
                           isMeasurementSpecified = FALSE,
                           residualCovSyntax = FALSE,
                           constrainedProductMean = FALSE) {

  modelSyntax <- modelSpecification$lavaanSyntax
  relationalDfs <- modelSpecification$relationalDfs
  productNames <- names(relationalDfs)
  addedSyntax <- c()
  # Measurement model latent products ------------------------------------------
  if (!is.logical(residualCovSyntax)) {
    stop2("residualCovSyntax is not FALSE or TRUE in generateSyntax")

  } else if (isMeasurementSpecified == FALSE) {
    measurementSyntax <-
      purrr::map2_chr(.x = modelSpecification$productNamesLatents,
                      .y = modelSpecification$indicatorProductNamesLatents,
                      .f = generateFormula.measurement,
                      operator = "=~",
                      firstFixed = TRUE)
    addedSyntax <- c(measurementSyntax, addedSyntax)
  }

  # Residual covariances -------------------------------------------------------
  if (!is.logical(residualCovSyntax)) {
    stop2("residualCovSyntax is not FALSE or TRUE in generateSyntax")

  } else if (residualCovSyntax == TRUE) {
    residualCovariances <- purrr::map2(.x = relationalDfs,
                               .y = productNames,
                               .f = getSyntaxResidualCovariances) |>
      unlist()
    addedSyntax <- c(addedSyntax, residualCovariances)
  }

  # Constrained product mean syntax --------------------------------------------
  if (constrainedProductMean == TRUE) {
    meanSyntax <-
      generateRestrictedMeanSyntax(modelSpecification$productNames,
                                   modelSpecification$elementsInProductNames)
    addedSyntax <- c(addedSyntax, meanSyntax)
  }

  syntaxElements <- c(modelSyntax, addedSyntax)
  stringr::str_c(syntaxElements, collapse = "\n")
}




# plural
generateRestrictedMeanSyntax <- function(productNames, elementsInProductNames) {

  purrr::map2(productNames, elementsInProductNames, restrictedMeanSyntaxSingle) |>
    unlist() |>
    stringr::str_c(collapse = "\n  ")

}



# this function assumes a product of only two latent variables no more
restrictedMeanSyntaxSingle <- function(productName, elementsInProductName) {
  covarianceLabel <- paste0(" ~~ cov", productName, "*")
  covariance <- stringr::str_c(elementsInProductName[1:2], collapse = covarianceLabel)
  meanStructure <- paste0(productName, " ~ cov", productName, "*1")
  paste(covariance, meanStructure, sep = "\n  ")
}




multiplyIndicators <- function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  if (ncol(df) <= 1){
    return(df)
  }
  product <- df[[1]] * df[[2]]
  y <- cbind(product, df[,-(1:2),drop = FALSE])


  unlist(multiplyIndicators(y))
}




# specify residual covariance in lavaanify parTable-format
getSyntaxResidualCovariances <- function(relationalDf, productName) {
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
                         FUN = stringr::str_c,
                         collapse = " ~~ ",
                         simplify = TRUE) |>
    stringr::str_c(collapse = "\n")


  syntaxOrthogonal <- apply(uniqueCombinations[!isShared, c("V1", "V2")],
                            MARGIN = 1,
                            FUN = stringr::str_c,
                            collapse = " ~~ 0*",
                            simplify = TRUE) |>
    stringr::str_c(collapse = "\n")

  title <- paste0("# Residual (Co)Variances: ", productName)
  paste(title,
        syntaxOblique,
        syntaxOrthogonal,
        sep = "\n")
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



generateFormula.measurement <- function(dependentName,
                                        predictorNames,
                                        operator = "=~",
                                        firstFixed = FALSE) {

  predictors <- stringr::str_c(predictorNames, collapse = " + ")

  if (firstFixed == TRUE) {
    predictors <- paste0("1*", predictors)
  }
  title <- paste("# Measurement Model:", dependentName)
  formula <- paste(dependentName, operator, predictors)
  paste(title, formula, sep = "\n")
}


summary.modsem <- function(object, ...) {
  lavaan::summary(object$lavaan)
}
