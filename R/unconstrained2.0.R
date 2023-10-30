

createProductIndicators.unconstrained <- function(modelSpecification) {
  purrr::map2(.x = modelSpecification$relationalDfs,
              .y = modelSpecification$indicatorsInProductTerms,
              .f = createIndicatorProducts,
              data = modelSpecification$data,
              centered = TRUE)

}


generateSyntax.unconstrained <- function(modelSpecification, isMeasurementSpecified = FALSE) {
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
  meanSyntax <- generateRestrictedMeanSyntax(modelSpecification$productNames,
                                             modelSpecification$elementsInProductNames)
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
