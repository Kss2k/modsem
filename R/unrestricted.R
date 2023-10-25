

createProductIndicators.unrestricted <- function(modelSpecification) {
  purrr::map2(.x = modelSpecification$relationDf,
              .y = modelSpecification$indicatorNames,
              .f = createIndicatorProducts,
              data = modelSpecification$data,
              centered = TRUE)

}


generateSyntax.unrestricted <- function(modelSpecification, ...) {
  modelSyntax <- modelSpecification$originalSyntax
  relationDf <- modelSpecification$relationDf
  # syntax for outer model and covariances
  outerAndCovs <- purrr::map2(.x = names(relationDf),
                             .y = relationDf,
                             .f = generateSyntax.outer.covs,
                             modelSyntax = modelSyntax, ...) |>
    unlist()

  meanSyntax <- restrictedMeanSyntax(modelSpecification$originalSyntax)
  syntaxElements <- c(fixLatentNames(modelSyntax), outerAndCovs, meanSyntax)
  stringr::str_c(syntaxElements, collapse = "\n")
}



# plural
restrictedMeanSyntax <- function(modelSyntax) {
  latents <- getLatentProducts(modelSyntax)
  purrr::map2(names(latents), latents, restrictedMeanSyntaxSingle) |>
    unlist() |>
    stringr::str_c(collapse = "\n  ")

}

# this function assumes a product of only two latent variables no more
restrictedMeanSyntaxSingle <- function(productName, varnames) {
  covarianceLabel <- paste0(" ~~ cov", productName, "*")
  covariance <- stringr::str_c(varnames[1:2], collapse = covarianceLabel)
  meanStructure <- paste0(productName, " ~ cov", productName, "*1")
  paste(covariance, meanStructure, sep = "\n  ")
}
