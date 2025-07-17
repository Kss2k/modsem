parseLavaan <- function(model.syntax = NULL, 
                        variableNames = NULL, 
                        match = FALSE,
                        suppress.warnings.match = FALSE,
                        match.recycle = FALSE) {
  # Check if a model.syntax is provided, if not we should return an error
  if (is.null(model.syntax))
    stop2("No model.syntax provided")
  else if (!is.character(model.syntax))
    stop2("The provided model syntax is not a string!")
  else if (length(model.syntax) > 1)
    stop2("The provided model syntax is not of length 1")

  parTable <- modsemify(model.syntax)
  structuralExprs <- parTable[parTable$op == "~",]
  measureExprs    <- parTable[parTable$op %in% c("=~", "<~"), ]

  oVs  <- getOVs(parTable)
  lVs  <- getLVs(parTable)
  vars <- getVarsPI(parTable)
  prodNames <- getProdNames(parTable)
  prodNamesCleaned <- stringr::str_remove_all(prodNames, ":")

  # Get all the indicators in the model
  inds <- getIndicators(parTable, observed=TRUE)
  stopif(!all(inds %in% variableNames),
         "Unable to find observed variables in data: ",
         capturePrint(inds[!inds %in% variableNames]))

  # Are prods latent?
  elementsInProds <- lapplyNamed(prodNames,
                                 FUN = splitProdName,
                                 pattern = ":",
                                 names = prodNamesCleaned)
  checkElementsInProds(elementsInProds, lVs=lVs, oVs=oVs)
  checkHigherOrderInteractions(elementsInProds, parTable=parTable)

  # Inds belonging to latent variables which are specified in the syntax
  indsLatents <- structureLavExprs(measureExprs)

  if (length(prodNamesCleaned) > 0) {
    # Get inds belonging to latent variables, or if observed, just get the
    # observed variable in prod terms
    indsInLatentProds <-
      lapplyNamed(elementsInProds[prodNamesCleaned],
                  FUN = getIndsMultipleVariables,
                  indsLatents = indsLatents,
                  names = prodNamesCleaned)

    # Creating a relDF for prodTerms
    relDfs <- lapply(indsInLatentProds, FUN = createRelDf, 
                     match = match, match.recycle = match.recycle,
                     suppress.warnings.match = suppress.warnings.match)

    # Get a list with all the inds in each interactionterm
    indsInLatentProds <- lapplyNamed(indsInLatentProds,
                                   FUN = function(x) unname(unlist(x)),
                                   names = prodNamesCleaned)
    # create the names for the indProds
    indProdNames <- lapplyNamed(relDfs,
                                FUN = colnames,
                                names = names(relDfs))
  } else { # in the case where ther is no interaction effects
    indsInLatentProds <- NULL
    relDfs <- NULL
    indProdNames <- NULL
  }

  # Return modelSpec
  modelSpec <- list(model.syntax = model.syntax,
                    parTable     = parTable,

                    oVs = oVs,
                    lVs = lVs,
                    prodNames           = prodNamesCleaned,
                    elementsInProdNames = elementsInProds,
                    relDfs = relDfs,

                    indsInLatentProds = indsInLatentProds,
                    latentProds       = names(indsInLatentProds),
                    indProdNames      = indProdNames)
  modelSpec
}


# Function for structuring exprs in a parTable
structureLavExprs <- function(parTable = NULL) {
  # If empty return NULL
  if (is.null(parTable)) return(NULL)
  # Same if nrow <= 0
  else if (nrow(parTable) <= 0) return(NULL)

  # Get Dependents
  names <- unique(parTable$lhs)

  # see utils for definition of selectValuesByCol() and lapplyNamed()
  lapplyNamed(names, selectValuesByCol, parTable, "rhs" ,"lhs")
}


getIndsVariable <- function(varName = NULL,  indsLatents = NULL) {
  stopif(is.null(varName), "Error in getIndsVariable(), varName is NULL")
  # Get the names of the latent variables in the model
  lVs <- names(indsLatents)

  # Check if our varName is a latent- or an observed variable
  if (!(varName %in% lVs)) {
    # If it is not a latent variable, we should just return the observed variable
    return(varName)
  } else if (varName %in% lVs) {
    # If it is a latent variable, we should return its inds
    return(indsLatents[[varName]])
  }

  # Error if it is neither a latent- nor an observerd variable
  stop2("Something went wrong in getIndsVariable(), ",
        "varName neither observed not latent")
}


# Function for getting inds for mutliple variables
getIndsMultipleVariables <- function(varNames = NULL, indsLatents = NULL) {
  # varNames = vector with variableNames
  # output should be a list
  # Check that arguements are supplied/not NULL
  stopif(is.null(varNames), "Error in getIndsMultipleVariables(), varNames is NULL")
  # use getIndsVariable for each element in varNames, and return a named list
  lapplyNamed(varNames,
              FUN = getIndsVariable,
              indsLatents = indsLatents)
}


splitProdName <- function(prodName, pattern) {
  stopif(is.null(prodName) ,"prodNames in splitProdName was NULL")
  stopif(is.null(pattern) ,"pattern in splitProdName was NULL")
  stringr::str_split_1(prodName, pattern)
}


splitProdNamesVec <- function(prodNames, pattern) {
  stopif(is.null(prodNames) ,"prodNames in splitProdNamesVec was NULL")
  stopif(is.null(pattern) ,"pattern in splitProdNamesVec was NULL")
  unlist(stringr::str_split(prodNames, pattern))
}


fixProdNames <- function(prodName, pattern = NULL) {
  stopif(is.null(pattern) ,"pattern in fixProdNames was NULL")
  stringr::str_remove_all(prodName, pattern)
}


fixLatentNamesSyntax <- function(model.syntax, pattern) {
  stringr::str_replace_all(model.syntax, pattern, "")
}


createRelDf <- function(indsProdTerm, 
                        match = FALSE, 
                        suppress.warnings.match = FALSE,
                        match.recycle = TRUE) {
  if (match) {
    lengths <- vapply(indsProdTerm, FUN.VALUE = integer(1L),
                      FUN = length)

    if ((shortest <- min(lengths)) != (longest <- max(lengths))) {
      warnif(!suppress.warnings.match, "Unequal number of indicators for latent variables ",
             "in product term, some indicators will be recycled!")

      if (match.recycle)
        selector <- \(x) x[rep(seq_along(x), length.out = longest)]
      else
        selector <- \(x) x[seq_len(shortest)]

      indsProdTerm <- lapply(X = indsProdTerm, FUN = selector)
    }

    relDf <- t(as.data.frame(indsProdTerm))

  } else if (!match) {
    allCombos <- t(expand.grid(indsProdTerm))
    relDf <- NULL
    
    for (i in seq_len(ncol(allCombos))) {
      if (!greplRowDf(allCombos[, i], relDf)) 
        relDf <- cbind(relDf, allCombos[, i])
    }
  }

  names <- apply(relDf, MARGIN = 2, FUN = stringr::str_c, collapse = "")
  structure(as.data.frame(relDf),
            names = names,
            row.names = names(indsProdTerm))
}


# This should return just an observed variable, if it does not belong in a latent one
getIndsLatents <- function(names, latents) {
  lapplyNamed(names,
              FUN = function(name, latents) {
                if (name %in% names(latents)) {
                  latents[[name]]
                } else {
                  name
                }
              },
              latents = latents,
              names = names)
}
