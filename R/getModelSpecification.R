parseLavaan <- function(modelSyntax = NULL, variableNames = NULL, match = FALSE) {
  # Checking prerequisites -----------------------------------------------------
  # Check if a modelSyntax is provided, if not we should return an error
  if (is.null(modelSyntax)) {
    stop("No modelSyntax provided")
    # Check if .model is string
  } else if (!is.character(modelSyntax)) {
    stop("The provided model syntax is not a string!")
  }
  # Convert to partable --------------------------------------------------------
  parTable <- modsemify(modelSyntax)

  # Extracting some general information ----------------------------------------
  # Structural
  structuralExprs <- parTable[parTable$op == "~",]

  # Interactions/Prod exprs
  interactionExprs <- parTable[grepl(":", parTable$rhs) & parTable$op == "~", ]
  prodNames <- unique(interactionExprs$rhs)
  prodNamesCleaned <- stringr::str_remove_all(prodNames, ":")

  # Logical vector indicating whether a row is a measureExpr
  isMeasure <- parTable$op %in% c("=~", "<~")
  measureExprs <- parTable[isMeasure, ]

  prodMeasureExprs <- parTable[isMeasure & grepl(":", parTable$lhs), ]
  specifiedProds <- unique(stringr::str_remove_all(prodMeasureExprs$lhs,":"))
  unspecifiedProds <- prodNamesCleaned[!(prodNamesCleaned %in% specifiedProds)]

  # Get all the inds in the model
  inds <- unique(measureExprs$rhs[!grepl(":", measureExprs$rhs)])
  if (!all(inds %in% variableNames)) {
    stop("Unable to find observed variables in data: ",
         capturePrint(inds[!inds %in% variableNames]))
  }
  # Get latent variable names (assumes first order)
  latentVariables <- unique(measureExprs$lhs)

  # Structuring the information ------------------------------------------------
  # Are prods latent?
  elementsInProds <- lapplyNamed(prodNames,
                                 FUN = splitProdName,
                                 pattern = ":",
                                 names = prodNamesCleaned)
  areElementsLatent <-
    lapply(elementsInProds,
           FUN = function(x, inds)  (x %in% latentVariables),
           inds)

  numberLatentElements <- lapply(areElementsLatent,
                                 FUN = sum)

  # A prod is considered latent if it has at least one latent variable
  isProdLatent <- lapply(numberLatentElements, FUN = as.logical)
  latentProds <- unlist(prodNamesCleaned)[unlist(isProdLatent)]

  # Inds belonging to latent variables which are specified in the syntax
  indsLatents <- structureLavExprs(measureExprs)

  relDfs <- list()
  indsInProdTerms <- list()
  indsInUnspecifedLatentProds <- list()
  indProdNamesUnspecifiedLatents <- list()
  # Computation for non-specified latent prods ---------------------------------
  if (length(unspecifiedProds) > 0) {
    # nSpec = unspecied
    # Get inds belonging to latent variables, or if observed, just get the
    # observed variable in prod terms
    nSpecIndsInProdTerms <-
      lapplyNamed(elementsInProds[unspecifiedProds],
                  FUN = getIndsMultipleVariables,
                  indsLatents = indsLatents,
                  names = unspecifiedProds)

    # Creating a relDF for prodTerms
    nSpecRelDfs <- lapply(nSpecIndsInProdTerms,
                          FUN = createRelDf,
                          match = match)

    # Get a list with all the inds in each interactionterm
    nSpecIndsInProdTerms <- lapplyNamed(nSpecIndsInProdTerms,
                                        FUN = function(x) unname(unlist(x)),
                                        names = unspecifiedProds)
    # create the names for the indProds
    nSpecIndNamesProdTerms <- lapplyNamed(nSpecRelDfs,
                                          FUN = colnames,
                                          names = names(nSpecRelDfs))
    relDfs <- c(relDfs, nSpecRelDfs)

    indProdNamesUnspecifiedLatents <-
      nSpecIndNamesProdTerms[unlist(isProdLatent[unspecifiedProds])]
    indsInProdTerms <- c(indsInProdTerms,
                         nSpecIndsInProdTerms)
    indsInUnspecifedLatentProds <-
      nSpecIndsInProdTerms[unlist(isProdLatent[unspecifiedProds])]
  }

  # Specified latent prods -----------------------------------------------------
  if (length(specifiedProds) > 0) {
    #fix names in lhs-column
    prodMeasureExprs[["lhs"]] <-
      fixProdNames(prodMeasureExprs[["lhs"]], pattern = ":")
    # Split the exprs based on what prod term they belong to
    exprsLatentProds <-
      lapplyNamed(specifiedProds, FUN = function(name, df) df[df$lhs == name, ],
                  df = prodMeasureExprs)
    # create a list with ind names in latent prod, and fix names
    specIndProdNamesLatents <- lapply(exprsLatentProds, 
                                      FUN = function(df) df[["rhs"]])

    # Get all inds in the prodTerms
    specRelDfs <- createRelDfSpecifiedProds(specIndProdNamesLatents,
                                            elementsInProds[specifiedProds]) 
    specIndsInProdTerms <-
      lapplyNamed(specRelDfs,
                  FUN = function(df) {
                    out <- vector("list", nrow(df))
                    names(out) <- rownames(df)
                    for(i in names(out)) {
                      out[[i]] <- as.vector(unlist(df[i, ]))
                    }
                    out
                  }, names = names(specRelDfs))
    specIndProdNamesLatents <-
      lapplyNamed(specIndProdNamesLatents, FUN = fixProdNames,
                  pattern = ":", names = names(specIndProdNamesLatents))
    specIndsInProdTerms <- lapplyNamed(specIndsInProdTerms,
                                       FUN = function(x) unname(unlist(x)),
                                       names = unspecifiedProds)
    # Create relational df with all the possible combos
    indsInProdTerms <- c(indsInProdTerms,
                         specIndsInProdTerms)
    relDfs <- c(relDfs, specRelDfs) 
  }
  # Return modelSpec -----------------------------------------------------------
  modelSpec <- list(modelSyntax = modelSyntax,
                    parTable = parTable,
                    prodNames = prodNamesCleaned,
                    elementsInProdNames = elementsInProds,
                    relDfs = relDfs,
                    indsInProdTerms = indsInProdTerms,
                    indsInUnspecifedLatentProds = indsInUnspecifedLatentProds,
                    unspecifiedLatentProds = names(indsInUnspecifedLatentProds),
                    indProdNamesUnspecifiedLatents = indProdNamesUnspecifiedLatents)
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
           stop("Error in getIndsVariable(), varName is NULL")
    )
  }
  # Get the names of the latent variables in the model
  latentVariables <- names(indsLatents)

  # Check if our varName is a latent- or an observed variable
  if (!(varName %in% latentVariables)) {
    # If it is not a latent variable, we should just return the observed variable
    return(varName)
  } else if (varName %in% latentVariables) {
    # If it is a latent variable, we should return its inds
    return(indsLatents[[varName]])
  }

  # Error if it is neither a latent- nor an observerd variable
  stop("Something went wrong in getIndsVariable(), ",
        "varName neither observed not latent")
}


# Function for getting inds for mutliple variables -----------------------
getIndsMultipleVariables <- function(varNames = NULL, indsLatents = NULL) {
  # varNames = vector with variableNames
  # output should be a list
  # Check that arguements are supplied/not NULL
  if (is.null(varNames)) {
    stop("Error in getIndsMultipleVariables(), varNames is NULL")
  }
  # use getIndsVariable for each element in varNames, and return a named list
  lapplyNamed(varNames,
              FUN = getIndsVariable,
              indsLatents = indsLatents)
}


splitProdName <- function(prodName, pattern) {
  if (is.null(prodName)) {
    stop("prodNames in splitProdName was NULL")
  }
  if (is.null(pattern)) {
    stop("pattern not supplied in splitProdName")
  }
  stringr::str_split_1(prodName, pattern)
}


splitProdNamesVec <- function(prodNames, pattern) {
  if (is.null(prodNames)) {
    stop("prodNames in splitProdNamesVec was NULL")
  }
  if (is.null(pattern)) {
    stop("pattern not supplied in splitProdName")
  }
  unlist(stringr::str_split(prodNames, pattern))
}


fixProdNames <- function(prodName, pattern = NULL) {
  if (is.null(pattern)) {
    stop("pattern not supplied in fixProdNames")
  }
  stringr::str_remove_all(prodName, pattern)
}


fixLatentNamesSyntax <- function(modelSyntax, pattern) {
  stringr::str_replace_all(modelSyntax, pattern, "")
}


createRelDf <- function(indsProdTerm, match = FALSE) {
  if (match) {
    lengths <- vapply(indsProdTerm, FUN.VALUE = integer(1L), 
                      FUN = length) 
    if ((shortest <- min(lengths)) != (longest <- max(lengths))) {
      warning("Unequal number of indicators for latent variables ",
              "in product term, not all indicators will be used")
      indsProdTerm <- lapply(indsProdTerm,
                             FUN = function(x, shortest) x[seq_len(shortest)],
                             shortest = shortest)
    }
    relDf <- t(as.data.frame(indsProdTerm))
  }
  else if (!match) 
    relDf <- t(expand.grid(indsProdTerm))
  names <- apply(relDf, MARGIN = 2, FUN = stringr::str_c, collapse = "")
  structure(as.data.frame(relDf),
            names = names)
}


createRelDfSpecifiedProds <- function(indsProdTerms, latents) {
  specRelDfs <- lapply(indsProdTerms, FUN = stringr::str_split,
                       pattern = ":") |>
    lapply(FUN = function(latent)
            lapply(latent, 
              FUN = function(ind)
                as.data.frame(matrix(ind, nrow = 1))) |> 
       purrr::list_rbind())
  specRelDfs <- 
    purrr::map2(.x = specRelDfs, .y = latents,
                .f = function(.x, .y) {
                  colnames(.x) <- .y
                  rownames(.x) <- apply(.x, MARGIN = 1, 
                                        stringr::str_c,
                                        collapse = "", 
                                        simplify = TRUE)
                  .x
                })
  lapply(specRelDfs, FUN = function(df) as.data.frame(t(df))) 
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
