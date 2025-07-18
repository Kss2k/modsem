getParTableResCov <- function(relDf, 
                              method, 
                              pt = NULL,
                              explicit.zero = FALSE,
                              include.single.inds = FALSE,
                              len.diff.ignore = 1L,
                              setToZero = FALSE) {
  simple <- \() 
    getParTableResCov.simple(relDf, 
                             explicit.zero = explicit.zero,
                             include.single.inds = include.single.inds,
                             len.diff.ignore = len.diff.ignore)
  switch(method,
    simple = simple(),
    simple.no.warn = suppressWarnings(simple()),
    ca = getParTableResCov.ca(relDf, pt = pt),
    equality = getParTableResCov.equality(relDf, setToZero = setToZero)
  )
}


getParTableResCov.simple <- function(relDf, explicit.zero = FALSE, include.single.inds = FALSE,
                                     len.diff.ignore = 1L) {
  EMPTY <- data.frame(lhs = NULL, op = NULL, rhs = NULL, mod = NULL)
  attr(EMPTY, "OK") <- TRUE

  if (length(relDf) <= 1) 
    return(EMPTY)

  OK <- TRUE

  if (include.single.inds) {
    allInds <- unique(unlist(relDf))
    relDf <- as.list(relDf)
    relDf <- c(relDf, as.list(setNames(allInds, nm = allInds)))
  }

  prodNames <- sort(names(relDf))
  uniqueCombinations <- getUniqueCombos(prodNames)
  # Now we want to specify the covariance based on shared inds
  isShared <- vector("logical", length = nrow(uniqueCombinations))

  for (i in seq_len(nrow(uniqueCombinations))) {
    indsProd1 <- unlist(relDf[uniqueCombinations[i, "V1"]])
    indsProd2 <- unlist(relDf[uniqueCombinations[i, "V2"]])
    # Compare the Inds in prod1 and prod2, and convert to integer
    sharedValues <- as.integer(indsProd1 %in% indsProd2)
    
    # Sum the values
    numberShared <- sum(sharedValues)
                                     
    # if there is a difference in number of elems in `len.diff.ignore` they should be ignored...
    # See the simulation results in `tests/testthat/test_three_way.R`
    #> cov(x, z, w, xz, xw, zx, xzw)
    #>        x    z    w     xz    xw     zw    xzw
    #> x   1.20 0.70 0.80  0.000 0.000  0.000  1.280
    #> z   0.70 1.80 0.60  0.000 0.000  0.000  1.860
    #> w   0.80 0.60 1.40  0.000 0.000  0.000  0.960
    #> xz  0.00 0.00 0.00  2.651 1.280  1.860 -0.001
    #> xw  0.00 0.00 0.00  1.280 2.321  1.460  0.000
    #> zw  0.00 0.00 0.00  1.860 1.460  2.880 -0.002
    #> xzw 1.28 1.86 0.96 -0.001 0.000 -0.002  8.225
    len.diff <- abs(length(indsProd2) - length(indsProd1))
    ignore <- len.diff %in% len.diff.ignore

    if      (numberShared >= 1 && !ignore) isShared[[i]] <- TRUE
    else if (numberShared == 0 && !ignore) isShared[[i]] <- FALSE
  }

  if (all(isShared)) {
    warning2("All residual covariances between product indicators were freed!\n",
             "The model will likely not be identifiable! Please try passing:\n",
             "  `res.cov.method = \"none\"` or `res.cov.method = \"equality\"`",
             immediate. = FALSE)
    OK <- FALSE
  }

  prodsSharingInds    <- uniqueCombinations[isShared, c("V1", "V2")]
  prodsNotSharingInds <- uniqueCombinations[!isShared, c("V1", "V2")]
  
  # Syntax for oblique covariances
  if (nrow(prodsSharingInds) > 0) {
    syntaxOblique <- apply(prodsSharingInds,
                           MARGIN = 1,
                           FUN = createParTableRow,
                           op = "~~") |>
      purrr::list_rbind()

  } else syntaxOblique <- NULL

  if (explicit.zero && nrow(prodsNotSharingInds) > 0) {
    syntaxOrthogonal <- apply(prodsNotSharingInds,
                              MARGIN = 1,
                              FUN = createParTableRow,
                              op = "~~",
                              mod = "0") |>
      purrr::list_rbind()

  } else syntaxOrthogonal <- NULL

  out <- rbind(syntaxOrthogonal, syntaxOblique)
  
  if (is.null(out))
    return(EMPTY)

  attr(out, "OK") <- OK
  out
}



# Rescovs with same indicator constrained to equality
getParTableResCov.equality <- function(relDf, setToZero = FALSE) {
  EMPTY <- data.frame(lhs = NULL, op = NULL, rhs = NULL, mod = NULL)
  attr(EMPTY, "OK") <- TRUE

  if (length(relDf) <= 1) 
    return(EMPTY)

  OK <- TRUE

  prodNames <- sort(names(relDf))
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

  uniqueCombos <- getUniqueCombos(prodNames)
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
  
  attr(parTable, "OK") <- OK
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
  varLatentNotShared <- trace_path(pt, latentNotShared, latentNotShared)
  varIndShared <- createLabelVar(indShared)

  rhs <- paste(lambdaShared[[1]], lambdaShared[[2]],
               varLatentNotShared, varIndShared, sep = " * ")
  lhs <- createLabelCov(indProd1, indProd2)
  createParTableRow(c(lhs, rhs), op = "==")
}
