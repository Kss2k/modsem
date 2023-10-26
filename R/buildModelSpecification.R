
# create the functions assuming a single latent variable
createProductCombos <- function(modelSyntax) {
  # List with dataframes
      # Example:
      # createPairLists(example3)
  indicators <- getProductIndicators(modelSyntax)

  lapply(indicators, expand.grid)

}




# function for getting all individual indicators used in each latent product
getIndividualIndicators <- function(modelSyntax) {
  # List with vectors
  # Example: getIndividualIndicators(example3)

  indicators <- getProductIndicators(modelSyntax)

  lapply(indicators, function(x) unname(unlist(x)))

}




createIndicatorNames <- function(modelSyntax) {
  # Create names for indicatorproducts for the respective latent products
    # Example:
      # indicatorProductNames(createPairLists(example3))

  productCombo <- createProductCombos(modelSyntax)
  lapply(productCombo, function(df) apply(df, 1, stringr::str_c, collapse = ""))
}




# There must be a more elegant way of doing this
createRelationDf <- function(modelSyntax) {
  # This creates a list with dataframes where the rownames correspond to the latent variables,
  # the columns the name of the product
  # The values represent the varibles in the combos
  productCombos <- createProductCombos(modelSyntax)
  indicatorProductNames <- createIndicatorNames(modelSyntax)

  # first transpose, then convert back to dataframe
  relationDf <- lapply(productCombos, function(x) as.data.frame(t(x)))

  # Setting names
  for (i in seq_along(relationDf)) {
    names(relationDf[[i]]) <- indicatorProductNames[[i]]
    # Creating a attribute-vector with all the indicators used to create the products in the df
    attr(relationDf[[i]], "indicatorNames") <- uniqueValuesDf(relationDf[[i]])
  }


  relationDf
}





createModelSpecification <- function(modelSyntax, data, method = "rca", returnLavaan = FALSE,
                                      ...) {
  modelSpecification <- list(
    originalSyntax = modelSyntax,

    # I might not need these:, since createRelationDf contains the same info
      #productCombos = createProductCombos(modelSyntax),
      #indicatorProductNames = createIndicatorNames(modelSyntax),
    indicatorNames = getIndividualIndicators(modelSyntax),

    # This is a combination of the two above:
      # I definetely need a shorter word for this
    relationDf = createRelationDf(modelSyntax),
    data = data
  )


  # setting class to method used
  class(modelSpecification) <- method

  # computing the productindicators based on class/method
  modelSpecification$productIndicators <- createProductIndicators(modelSpecification)

  # combining the product indicators into a single dataframe (i.e., not a list of dataframes)
    # Using unname() avoids getting names like df1.var1, ... df2.var1 etc.. but assumes unique names

  combinedProductIndicators <- purrr::list_cbind(unname(modelSpecification$productIndicators))

  # combining the product indicators with the dataset, into one coherent one
  modelSpecification$modifiedData <- cbind(data, combinedProductIndicators)

  # Generating new syntax based on method
  modelSpecification$modifiedSyntax <- generateSyntax(modelSpecification)

  if (returnLavaan == TRUE) {
    return(lavaan::sem(modelSpecification$modifiedSyntax,
                      modelSpecification$modifiedData))
  }
  modelSpecification$lavaan <- lavaan::sem(modelSpecification$modifiedSyntax,
                                           modelSpecification$modifiedData, ...)
  # returning object
  modelSpecification
}


uniqueValuesDf <- function(df) {
    df |>
    unlist() |>
    unname() |>
    unique()
}
