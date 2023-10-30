
modsem <- function(modelSyntax,
                   data,
                   method = "rca",
                   isMeasurementSpecified = FALSE,
                   ...) {
  # Get the specifications of the model
  modelSpecification <-
    parseLavaan(modelSyntax,isMeasurementSpecified = isMeasurementSpecified)
  modelSpecification$data <- data
  class(modelSpecification) <- method

  productIndicators <- createProductIndicators(modelSpecification)
  mergedProductIndicators <- purrr::list_cbind(unname(productIndicators))
  newSyntax <- generateSyntax(modelSpecification,
                              isMeasurementSpecified = isMeasurementSpecified)
  newData <- cbind.data.frame(data, mergedProductIndicators)
  lavaanEstimation <- lavaan::sem(newSyntax,
                                  newData)

  modelSpecification$productIndicators <- productIndicators
  modelSpecification$newSyntax <- newSyntax
  modelSpecification$newData <- newData
  modelSpecification$lavaan <- lavaanEstimation
  # this is not pretty either
  class(modelSpecification) <- "modsem"
  modelSpecification
}


summary.modsem <- function(object, ...) {
  lavaan::summary(object$lavaan)
}


# Creating product indicators --------------------------------------------------


createProductIndicators <- function(modelSpecification) {
  UseMethod("createProductIndicators")
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



generateSyntax <- function(modelSpecification, isMeasurementSpecified) {
  isMeasurementSpecified <- isMeasurementSpecified
  UseMethod("generateSyntax")
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


