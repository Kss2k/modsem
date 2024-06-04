getParTableResCov <- function(relDf, method, ...) {
  switch(method, 
         "simple" = getParTableResCov.simple(relDf),
         "ca" = getParTableResCov.ca(relDf, ...),
         "equality" = getParTableResCov.equality(relDf, ...))
}



# Simple -----------------------------------------------------------------------
getParTableResCov.simple <- function(relDf) {
  if (ncol(relDf) <= 1) {
    return(NULL)
  }
  prodNames <- sort(colnames(relDf))
  uniqueCombinations <- getUniqueCombinations(prodNames)
  # Now we want to specify the covariance based on shared inds
  isShared <- vector("logical", length = nrow(uniqueCombinations))

  for (i in seq_len(nrow(uniqueCombinations))) {
    indsProd1 <- unlist(relDf[uniqueCombinations[i, "V1"]])
    indsProd2 <- unlist(relDf[uniqueCombinations[i, "V2"]])
    # Compare the Inds in prod1 and prod2, and convert to integer
    sharedValues <- as.integer(indsProd1 %in% indsProd2)
    # Sum the values
    numberShared <- sum(sharedValues)
    if (numberShared >= 1) {
      isShared[[i]] <- TRUE
    } else if (numberShared == 0) {
      isShared[[i]] <- FALSE
    }
  }
  # Syntax for oblique covariances
  prodsSharingInds <- uniqueCombinations[isShared, c("V1", "V2")]
  if (nrow(prodsSharingInds) > 0) {
    syntaxOblique <- apply(prodsSharingInds,
                           MARGIN = 1,
                           FUN = createParTableRow,
                           op = "~~") |>
      purrr::list_rbind()
  } else {
    syntaxOblique <- NULL
  }
  prodsNotSharingInds <- uniqueCombinations[!isShared, c("V1", "V2")]
  if (nrow(prodsNotSharingInds) > 0) {
    syntaxOrthogonal <- apply(prodsNotSharingInds,
                              MARGIN = 1,
                              FUN = createParTableRow,
                              op = "~~",
                              mod = "0") |>
      purrr::list_rbind()
  } else {
    syntaxOrthogonal <- NULL
  }
  rbind(syntaxOrthogonal, syntaxOblique)
}



# Rescovs with same indicator constrained to equality --------------------------
getParTableResCov.equality <- function(relDf, setToZero = FALSE) {
  if (ncol(relDf) <= 1) {
    return(NULL)
  }
  prodNames <- sort(colnames(relDf))
  sharedMatrix <- matrix("", nrow = length(prodNames), ncol = length(prodNames),
                         dimnames = list(prodNames, prodNames))
  # Now we want to specify the covariance based on shared inds
  for (i in prodNames) {
    for (j in prodNames) {
      sharedIndicators <- relDf[[i]][relDf[[i]] %in% relDf[[j]]]
      sharedMatrix[i, j] <- stringr::str_c(sharedIndicators,
                                           collapse = "_")
    }
  }
  labelMatrix <- sharedMatrix
  labelMatrix <- ifelse(labelMatrix == "", "0", paste0("share_", labelMatrix))
  labelMatrix[upper.tri(labelMatrix, diag = TRUE)] <- ""

  uniqueCombos <- getUniqueCombinations(prodNames)
  uniqueCombos[["labels"]] <- vector("character", length = nrow(uniqueCombos))
  for (i in seq_len(nrow(uniqueCombos))) {
    uniqueCombos[["labels"]][[i]] <- labelMatrix[uniqueCombos[i, "V2"],
                                                 uniqueCombos[i, "V1"]]
  }

  parTable <- apply(uniqueCombos, MARGIN = 1,
                    FUN = function(x)
                      createParTableRow(x[c("V1", "V2")], op = "~~", mod = x[["labels"]])
                    ) |>
    purrr::list_rbind()
  if (!setToZero) parTable <- parTable[parTable$mod != 0, ]
  parTable
}



# Constrained Approach ---------------------------------------------------------
getParTableResCov.ca <- function(relDf, pt) {
  if (nrow(relDf) > 2) {
    stop2("Constrained approach for constraining residual covariances should ",
         "not be used with latent products with more than two components")
  }
  if (ncol(relDf) <= 1) {
    return(NULL)
  }
  prodNames <- colnames(relDf)
  labelMatrix <- matrix("", nrow = length(prodNames), ncol = length(prodNames),
                        dimnames = list(prodNames, prodNames))
  # Now we want to specify the covariance based on shared inds
  parTable <- NULL
  for (i in 2:nrow(labelMatrix)) {
    rhs <- rownames(labelMatrix)[[i]]

    for (j in 1:(i-1)) {
      lhs <- colnames(relDf)[[j]]
      sharedIndicators <- relDf[[i]] %in% relDf[[j]]

      if (sum(sharedIndicators) >= 1) {
        labelMatrix[i, j] <- createLabelCov(rhs, lhs)

      } else if (sum(sharedIndicators) < 1) {
        labelMatrix[i, j] <- "0"
      }

      parTable <- rbind(parTable,
                        createParTableRow(c(rhs, lhs),
                                          op = "~~",
                                          mod = labelMatrix[i,j]))
    }

  }

  #apply eq constraints to those which are not set to zero
  if (length(parTable$mod[parTable$mod != "0"]) > 0) {
    eqConstraints <- apply(parTable[parTable$mod != "0", c("lhs", "rhs")],
                           MARGIN = 1,
                           FUN = function(vars, relDf)
                             getFormulaResCovProdInd(vars[["lhs"]], 
                                                     vars[["rhs"]], 
                                                     relDf, pt),
                           relDf = relDf) |>
      purrr::list_rbind()
  } else {
    eqConstraints <- NULL
  }
  rbind(parTable, eqConstraints)
}



getFormulaResCovProdInd <- function(indProd1, indProd2, relDf, pt) {
  if (is.null(indProd1) || is.null(indProd2)) {
    return(NULL)
  }
  cols <- c(indProd1, indProd2)
  rowShared <- relDf[relDf[[indProd1]] %in_paired% relDf[[indProd2]], cols]
  forceRowNames(rowShared) <- rownames(relDf)[relDf[[indProd1]] %in_paired% relDf[[indProd2]]]
  rowNotShared <- relDf[!(relDf[[indProd1]] %in_paired% relDf[[indProd2]]), cols]
  forceRowNames(rowNotShared) <- rownames(relDf)[!relDf[[indProd1]] %in_paired% relDf[[indProd2]]]

  latentNotShared <- rownames(rowNotShared)
  indShared <- rowShared[1, 1]
  indsNotShared <- unlist(rowNotShared[1, 1:2])
  lambdaShared <- createLabelLambda(indsNotShared, latentNotShared)
  varLatentNotShared <- tracePath(pt, latentNotShared, latentNotShared)
  varIndShared <- createLabelVar(indShared)

  rhs <- paste(lambdaShared[[1]], lambdaShared[[2]],
               varLatentNotShared, varIndShared, sep = " * ")
  lhs <- createLabelCov(indProd1, indProd2)
  createParTableRow(c(lhs, rhs), op = "==")
}
