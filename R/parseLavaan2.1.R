


#' Parse lavaan model
#'
#' @param modelSyntax lavaan syntax
#' @param isMeasureSpecified have you specified the measure model for the latent prod
#'
#' @return
#' @export
#'
#' @examples
parseLavaan <- function(
    modelSyntax   = NULL) {

  # Checking prerequisites -----------------------------------------------------
    # Check if a modelSyntax is provided, if not we should return an error
  if (is.null(modelSyntax)) {
    stop2("No modelSyntax provided")

    # Check if .model is string
  } else if (!is.character(modelSyntax)) {
    stop2("The provided model syntax is not a string!")
  }

  # Convert to lavaan partable -------------------------------------------------
  parTable <- modsemify(modelSyntax)

  # Extracting some general information ----------------------------------------
    # Structural
  structuralExprs <- parTable[parTable$op == "~",]

    # Interactions/Prod exprs
  interactionExprs <-
    parTable[grepl(":", parTable$rhs) & parTable$op == "~", ]
  prodNames <- unique(interactionExprs$rhs)
        # Get names for prodTerms without __COLON__
  prodNamesCleaned <- stringr::str_remove_all(prodNames,
                                                 ":")

  # Logical vector indicating whether a row is a measureExpr
  isMeasure <- parTable$op %in% c("=~", "<~")
  measureExprs <- parTable[isMeasure, ]

  prodMeasureExprs <-
    parTable[isMeasure & grepl(":", parTable$lhs), ]
  specifiedProds <-
    unique(stringr::str_remove_all(prodMeasureExprs$lhs,":"))
  unspecifiedProds <-
    prodNamesCleaned[!(prodNamesCleaned %in% specifiedProds)]

  # Get all the inds in the model
  inds <- unique(measureExprs$rhs)
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
  isProdLatent <- lapply(numberLatentElements,
                            FUN = as.logical)
  latentProds <- unlist(prodNamesCleaned)[unlist(isProdLatent)]

    # Inds belonging to latent variables which are specified in the syntax
  indsLatents <- structureLavExprs(measureExprs)


  relDfs <- list()
  indsInProdTerms <- list()
  indsInUnspecifedLatentProds <- list()
  indProdNamesUnspecifiedLatents <- list()
  # Computation for non-specified latent prods -----------------------------
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
                            FUN = createRelDf)

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

  # Specified latent prods --------------------------------------------------
  if (length(specifiedProds) > 0) {

    #fix names in lhs-column
    prodMeasureExprs[["lhs"]] <-
      fixProdNames(prodMeasureExprs[["lhs"]],
                      pattern = ":")

    # Split the exprs based on what prod term they belong to
    exprsLatentProds <-
      lapplyNamed(specifiedProds,
                  FUN = function(name, df) df[df$lhs == name, ],
                  df = prodMeasureExprs)
    # create a list with ind names in latent prod, and fix names
    specIndProdNamesLatents <-
      lapply(exprsLatentProds,
             FUN = function(df) df[["rhs"]])

    # Get all inds in the prodTerms

    specIndsInProdTerms <- lapply(elementsInProds[specifiedProds],
                                  FUN = getIndsLatents,
                                  latents = indsLatents)

    # Create relational df with all the possible combos
    allCombosRelDfs <- lapplyNamed(
      specIndsInProdTerms,
      FUN = createRelDf)

    # then prune these combos
    specRelDfs <- purrr::map2(
      .x = specIndProdNamesLatents,
      .y = allCombosRelDfs,
      .f = function(name, df) df[stringr::str_remove_all(name, pattern = ":")])

    # Flatten the list, so that it collectively referst to XY, not XY$X and XY$Y
    specIndsInProdTerms <- lapplyNamed(specIndsInProdTerms,
                                       FUN = function(x) unname(unlist(x)),
                                       names = specifiedProds)
    specIndProdNamesLatents <-
      lapplyNamed(specIndProdNamesLatents,
                  FUN = fixProdNames,
                  pattern = ":",
                  names = names(specIndProdNamesLatents))

    relDfs  <- c(relDfs, specRelDfs)
    indsInProdTerms <- c(indsInProdTerms,
                                  specIndsInProdTerms)
  }

  # Info nlsem -----------------------------------------------------------------
  etaNames <- unique(structuralExprs$lhs)
  allPredictors <- unique(structuralExprs$rhs)
  simplePredictors <- allPredictors[!grepl(":", allPredictors)]
  # I do this in the counter intuitive way, to keep the same order as the
  # measure exprs in the model, which is the way the nlsem model reads
  # the syntax
  latentSimplePredictors <-
    latentVariables[latentVariables %in% simplePredictors]

  nlsemInfo <- list(etaNames = etaNames,
                    indsEta = indsLatents[etaNames],
                    xiNames = latentSimplePredictors,
                    indsXi = indsLatents[latentSimplePredictors],
                    modelSyntax = modelSyntax)

  # Return modelSpec --------------------------------------------------
  modelSpec <- list(
    nlsem = nlsemInfo,
    parTable = parTable,
    prodNames = prodNamesCleaned,
    elementsInProdNames = elementsInProds,
    relDfs = relDfs,
    indsInProdTerms = indsInProdTerms,

    #unpecifiedLatentProds = latentProds[!(latentProds %in% specifiedProds)],
    indsInUnspecifedLatentProds = indsInUnspecifedLatentProds,
    unspecifiedLatentProds = names(indsInUnspecifedLatentProds),
    indProdNamesUnspecifiedLatents = indProdNamesUnspecifiedLatents
    )

  modelSpec
}



# Function for structuring exprs in a parTable ---------------------------
structureLavExprs <- function(parTable = NULL) {
  # If empty return NULL
  if (is.null(parTable)) return(NULL)
  # Same if nrow <= 0
  else if (nrow(parTable) <= 0) return(NULL)

  ### Get Dependents
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
  stop2(
    "Something went wrong in getIndsVariable(), varName neither observed not latent"
  )
}



# Function for getting inds for mutliple variables -----------------------

getIndsMultipleVariables <- function(varNames = NULL,
                                           indsLatents = NULL) {
  # varNames = vector with variableNames
  # output should be a list

  # Check that arguements are supplied/not NULL
  if (is.null(varNames)) {
    return(
      stop2("Error in getIndsMultipleVariables(), varNames is NULL")
    )
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



fixLatentNamesSyntax <- function(modelSyntax, pattern) {
  stringr::str_replace_all(modelSyntax, pattern, "")
}



createRelDf <- function(indsProdTerm) {
  relDf <- t(expand.grid(indsProdTerm))
  names <- apply(relDf, MARGIN = 2, FUN = stringr::str_c, collapse = "")

  structure(as.data.frame(relDf),
            names = names)

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




