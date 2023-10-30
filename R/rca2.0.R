
createProductIndicators.rca <- function(modelSpecification, ...) {
  # Compute indicatorProducts (wihtout centering the variables)

  indicatorProducts <- purrr::map2(.x = modelSpecification$relationalDfs,
                                   .y = modelSpecification$indicatorsInProductTerms,
                                   .f = createIndicatorProducts,
                                   data = modelSpecification$data,
                                   centered = FALSE)

  purrr::map2(.x = indicatorProducts,
              .y = modelSpecification$indicatorsInProductTerms,
              .f = calculateResidualsDf,
              data = modelSpecification$data)
}



# fucntion for calculating residuals for a dataframe of productindicators
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




generateSyntax.rca <- function(modelSpecification, isMeasurementSpecified = FALSE) {
  modelSyntax <- modelSpecification$lavaanSyntax
  relationalDfs <- modelSpecification$relationalDfs
  productNames <- names(relationalDfs)
  addedSyntax <- purrr::map2(.x = relationalDfs,
                             .y = productNames,
                             .f = getSyntaxResidualCovariances) |>
    unlist()

  if (isMeasurementSpecified == FALSE) {

    measurementSyntax <-
      purrr::map2_chr(.x = modelSpecification$productNamesLatents,
                      .y = modelSpecification$indicatorProductNamesLatents,
                      .f = generateFormula.measurement,
                      operator = "=~",
                      firstFixed = TRUE)

    addedSyntax <- c(measurementSyntax, addedSyntax)


  }
  syntaxElements <- c(modelSyntax, addedSyntax)
  stringr::str_c(syntaxElements, collapse = "\n")
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



