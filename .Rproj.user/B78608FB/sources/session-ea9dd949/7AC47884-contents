
modsem <- function(modelSyntax,
                   data,
                   method = "rca",
                   specifyCovs = TRUE,
                   firstFixed = TRUE,
                   returnLavaan = FALSE,
                   ...) {
  # Get the specifications of the model
  modelSpecification <- createModelSpecification(modelSyntax, data, method = method, returnLavaan = returnLavaan, ...)


  modelSpecification
}


readModsem <- function(modelSyntax, method = "rca", generateSyntax = TRUE) {
  modelSpecification <- list(
    originalSyntax = modelSyntax,
    indicatorNames = getIndividualIndicators(modelSyntax),
    relationDf = createRelationDf(modelSyntax)
  )
  class(modelSpecification) <- method
  if (generateSyntax == FALSE) {
    return(modelSpecification)
  }
  modelSpecification$modifiedSyntax <- generateSyntax(modelSpecification)
  modelSpecification
}

# ProductIndicators ####


createProductIndicators <- function(modelSpecification) {
  UseMethod("createProductIndicators")
}




multiplyIndicators <- function(df) {
  if (is.null(df)) return(NULL)
  if (ncol(df) <= 1) return(df)
  product <- df[[1]] * df[[2]]
  y <- cbind(product, df[,-(1:2),drop = FALSE])

  # old solution: multiplyIndicators(y)$product
  unlist(multiplyIndicators(y))
}





createIndicatorProducts <- function(relationDf, indicatorNames, data, centered = FALSE) {
  # Getting the indicatorProduct names
  varnames <- colnames(relationDf)

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



#### Syntax related ####


generateSyntax <- function(modelSpecification) {
  UseMethod("generateSyntax")
}




# function to generate outermodel syntaxes and covariances
generateSyntax.outer.covs <- function(latentName,
                                      relationDf,
                                      modelSyntax,
                                      specifyCovs = TRUE,
                                      firstFixed = TRUE) {
  productIndicators <- colnames(relationDf)


  # structural model
  structuralModel <- generateFormula(latentName, productIndicators, operator = "=~", firstFixed)

  # residual covariances
  if (specifyCovs == TRUE) {
    headerCovs <- paste0("# ModSEM: Residual (Co)Variances ", latentName)
    residualCovariances <- paste(headerCovs,
                                 syntaxResidualCovariance(relationDf), "\n",
                                 sep = "\n  ")
  } else if (specifyCovs == FALSE) {
    residualCovariances <- ""
  }
  # Adding them all together
  headerOuter <- paste0("  # ModSEM: Measurement Model ", latentName)
  newSyntax <- paste(headerOuter,
                     structuralModel,
                     "\n  ",
                     residualCovariances,
                     sep = "\n  ")

  newSyntax
}


# Used to generate formula for lm and lavaan
generateFormula <- function(dependentName,
                            predictorNames,
                            operator = "~",
                            firstFixed = FALSE) {

  predictors <- stringr::str_c(predictorNames, collapse = " + ")

  if (firstFixed == TRUE) {
    predictors <- paste0("1*", predictors)
  }
  paste(dependentName, operator, predictors)
}





# THIS IS just placholder, there is way to much nesting going on
syntaxResidualCovariance<- function(df) {
  n <- ncol(df)
  relationalMatrix <- matrix(FALSE, nrow = n, ncol = n)
  variables <- colnames(df)

  for (i in 1:n) {
    for (j in 1:i) {
      sharedValuesLogical <- sum(df[[variables[i]]] %in% df[[variables[j]]]) >= 1
      relationalMatrix[i, j] <- sharedValuesLogical
    }
  }
  # is this necessary?
  colnames(relationalMatrix) <- variables
  rownames(relationalMatrix) <- variables


  replaceLogical(relationalMatrix)

}




replaceLogical <- function(relationalMatrix) {
  colNames <- colnames(relationalMatrix)


  out <- matrix("", nrow = nrow(relationalMatrix),
                ncol = ncol(relationalMatrix))

  orthogonal <- residualCovariance(colNames,oblique = FALSE)
  oblique <- residualCovariance(colNames, oblique = TRUE)


  out[relationalMatrix == FALSE] <- orthogonal[relationalMatrix == FALSE]
  out[relationalMatrix == TRUE ] <- oblique[relationalMatrix == TRUE]
  out[lower.tri(out)] |>
    stringr::str_flatten(collapse = "\n  ")
}




residualCovariance <- function(varnames, oblique = TRUE) {
  if (oblique == TRUE) {
    operator <- ""
  } else if (oblique == FALSE) {
    operator <- "0*"
  }

  replacement <- outer(varnames, varnames, function(x,y) paste0(x, " ~~ ", operator, y)) |>
    as.matrix()
  colnames(replacement) <- varnames
  rownames(replacement) <- varnames
  replacement
}


fixLatentNames <- function(modelSyntax) {
  stringr::str_replace_all(modelSyntax, ":", "")
}

