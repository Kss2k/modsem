getParTableResCov <- function(relDf, 
                              method, 
                              pt = NULL,
                              explicit.zero = FALSE,
                              include.single.inds = FALSE,
                              setToZero = FALSE) {
  simple <- \() 
    getParTableResCov.simple(relDf, 
                             explicit.zero = explicit.zero,
                             include.single.inds = include.single.inds)
  switch(method,
    simple = simple(),
    simple.no.warn = suppressWarnings(simple()),
    ca = getParTableResCov.ca(relDf, pt = pt),
    equality = getParTableResCov.equality(relDf, setToZero = setToZero)
  )
}


getParTableResCov.simple <- function(relDf, explicit.zero = FALSE, include.single.inds = FALSE) {
  EMPTY <- data.frame(lhs = NULL, op = NULL, rhs = NULL, mod = NULL)
  attr(EMPTY, "OK") <- TRUE

  if (length(relDf) <= 1) 
    return(EMPTY)

  OK <- TRUE

  if (include.single.inds) {
    allInds <- unique(unlist(relDf))
    relDf <- as.list(relDf)
    relDf <- c(relDf, as.list(stats::setNames(allInds, nm = allInds)))
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
    #> cov(x, z, w, xz, xw, zw, xm, zm, wm, xzw, xzm, xwm, zwm, xzwm)
    #>         x    z    w   xz   xw   zw   xm   zm   wm   xzw  xzm  xwm  zwm  xzwm
    #> x    1.20 0.70 0.80 0.00 0.00 0.00 0.00 0.00 0.00  1.84 0.64 1.04 0.78  0.00
    #> z    0.70 1.80 0.60 0.00 0.00 0.00 0.00 0.00 0.00  2.28 0.78 0.78 1.44  0.00
    #> w    0.80 0.60 1.40 0.00 0.00 0.00 0.00 0.00 0.00  1.94 0.78 1.24 1.14  0.00
    #> xz   0.00 0.00 0.00 2.65 1.28 1.86 0.50 0.57 0.36  0.00 0.00 0.00 0.00  3.36
    #> xw   0.00 0.00 0.00 1.28 2.32 1.46 0.88 0.54 0.76  0.00 0.00 0.00 0.00  3.25
    #> zw   0.00 0.00 0.00 1.86 1.46 2.88 0.66 1.26 0.78  0.00 0.00 0.00 0.00  4.09
    #> xm   0.00 0.00 0.00 0.50 0.88 0.66 2.44 1.46 1.72  0.00 0.00 0.00 0.00  4.54
    #> zm   0.00 0.00 0.00 0.57 0.54 1.26 1.46 3.69 1.38  0.00 0.00 0.00 0.00  5.56
    #> wm   0.00 0.00 0.00 0.36 0.76 0.78 1.72 1.38 3.16  0.00 0.00 0.00 0.00  4.96
    #> xzw  1.84 2.28 1.94 0.00 0.00 0.00 0.00 0.00 0.00 10.25 3.91 3.88 4.56  0.01
    #> xzm  0.64 0.78 0.78 0.00 0.00 0.00 0.00 0.00 0.00  3.91 6.98 4.69 5.79  0.01
    #> xwm  1.04 0.78 1.24 0.00 0.00 0.00 0.00 0.00 0.00  3.88 4.69 7.67 5.42  0.00
    #> zwm  0.78 1.44 1.14 0.00 0.00 0.00 0.00 0.00 0.00  4.56 5.79 5.42 8.90  0.00
    #> xzwm 0.00 0.00 0.00 3.36 3.25 4.09 4.54 5.56 4.96  0.01 0.01 0.00 0.00 28.76
    ignore <- length(indsProd1) %% 2 != length(indsProd2) %% 2

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
