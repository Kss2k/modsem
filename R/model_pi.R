parseLavaan <- function(model_syntax = NULL, variableNames = NULL, match = FALSE) {
  # Checking prerequisites -----------------------------------------------------
  # Check if a model_syntax is provided, if not we should return an error
  if (is.null(model_syntax)) 
    stop2("No model_syntax provided")
  else if (!is.character(model_syntax)) 
    stop2("The provided model syntax is not a string!")
  else if (length(model_syntax) > 1) 
    stop2("The provided model syntax is not of length 1")

  # Convert to partable --------------------------------------------------------
  parTable <- modsemify(model_syntax)
  structuralExprs <- parTable[parTable$op == "~",]
  measureExprs <- parTable[parTable$op %in% c("=~", "<~"), ]
  lVs <- unique(parTable$lhs[parTable$op == "=~"])
  vars <- unique(c(parTable$rhs[parTable$op %in% c("~", "=~") & 
                   parTable$rhs != "1"],
                   parTable$lhs[parTable$op == "~"])) |>
    stringr::str_split(pattern = ":", simplify = FALSE) |>
    unlist() |> unique()
  oVs <- vars[!vars %in% lVs]
  prodNames <- parTable$rhs[grepl(":", parTable$rhs) & parTable$op == "~"] |>
    unique()
  prodNamesCleaned <- stringr::str_remove_all(prodNames, ":")

  # Get all the indicators in the model
  inds <- unique(measureExprs$rhs[!grepl(":", measureExprs$rhs)])
  if (!all(inds %in% variableNames)) {
    stop2("Unable to find observed variables in data: ",
         capturePrint(inds[!inds %in% variableNames]))
  }

  # Are prods latent?
  elementsInProds <- lapplyNamed(prodNames,
                                 FUN = splitProdName,
                                 pattern = ":",
                                 names = prodNamesCleaned)

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
    relDfs <- lapply(indsInLatentProds,
                     FUN = createRelDf,
                     match = match)

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

  # Return modelSpec -----------------------------------------------------------
  modelSpec <- list(model_syntax = model_syntax,
                    parTable = parTable,
                    oVs = oVs,
                    lVs = lVs,
                    prodNames = prodNamesCleaned,
                    elementsInProdNames = elementsInProds,
                    relDfs = relDfs,
                    indsInLatentProds = indsInLatentProds,
                    latentProds = names(indsInLatentProds),
                    indProdNames = indProdNames)
  modelSpec
}


# Function for structuring exprs in a parTable ---------------------------
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
  if (is.null(varName)) {
    return(
           stop2("Error in getIndsVariable(), varName is NULL")
    )
  }
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


# Function for getting inds for mutliple variables -----------------------
getIndsMultipleVariables <- function(varNames = NULL, indsLatents = NULL) {
  # varNames = vector with variableNames
  # output should be a list
  # Check that arguements are supplied/not NULL
  if (is.null(varNames)) {
    stop2("Error in getIndsMultipleVariables(), varNames is NULL")
  }
  # use getIndsVariable for each element in varNames, and return a named list
  lapplyNamed(varNames,
              FUN = getIndsVariable,
              indsLatents = indsLatents)
}


splitProdName <- function(prodName, pattern) {
  if (is.null(prodName)) {
    stop2("prodNames in splitProdName was NULL")
  }
  if (is.null(pattern)) {
    stop2("pattern not supplied in splitProdName")
  }
  stringr::str_split_1(prodName, pattern)
}


splitProdNamesVec <- function(prodNames, pattern) {
  if (is.null(prodNames)) {
    stop2("prodNames in splitProdNamesVec was NULL")
  }
  if (is.null(pattern)) {
    stop2("pattern not supplied in splitProdName")
  }
  unlist(stringr::str_split(prodNames, pattern))
}


fixProdNames <- function(prodName, pattern = NULL) {
  if (is.null(pattern)) {
    stop2("pattern not supplied in fixProdNames")
  }
  stringr::str_remove_all(prodName, pattern)
}


fixLatentNamesSyntax <- function(model_syntax, pattern) {
  stringr::str_replace_all(model_syntax, pattern, "")
}


createRelDf <- function(indsProdTerm, match = FALSE) {
  if (match) {
    lengths <- vapply(indsProdTerm, FUN.VALUE = integer(1L), 
                      FUN = length) 
    if ((shortest <- min(lengths)) != (longest <- max(lengths))) {
      warning2("Unequal number of indicators for latent variables ",
              "in product term, not all indicators will be used")
      indsProdTerm <- lapply(indsProdTerm,
                             FUN = function(x, shortest) x[seq_len(shortest)],
                             shortest = shortest)
    }
    relDf <- t(as.data.frame(indsProdTerm))
  } else if (!match) {
    allCombos <- t(expand.grid(indsProdTerm))
    relDf <- NULL
    for (i in seq_len(ncol(allCombos))) {
      if (!greplRowDf(allCombos[, i], relDf)) {
        relDf <- cbind(relDf, allCombos[, i])
      }
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
