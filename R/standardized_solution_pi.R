correctStdSolutionPI <- function(object, parTable.std, grouping = NULL) {
  # Function for correcting the standardized solution for interaction terms
  # lavaan standardizes the interaction term (X:Z) to have a variance of 1,
  # which is not correct if cov(X, Z) != 0.
  parTable <- parameter_estimates(object)
  intTerms <- object$elementsInProdNames
  cols <- c("est", "se", "ci.lower", "ci.upper")
  cols <- intersect(cols, names(parTable.std))

  vcov <- modsem_inspect(object, what = "cov.all")
  vcor <- modsem_inspect(object, what = "cor.all")

  if (!is.null(grouping)) {
    parTable     <- subsetByGrouping(parTable, grouping = grouping)
    parTable.std <- subsetByGrouping(parTable.std, grouping = grouping)
    if (!NROW(parTable) || !NROW(parTable.std)) return(NULL)

    stopif(!"block" %in% names(grouping), "Missing information about `block`!\n")

    vcov <- vcov[[grouping[["block"]]]]
    vcor <- vcor[[grouping[["block"]]]]
  } 

  stopif(!"est" %in% cols, 
         "The parTable must contain the 'est' column for the standardized solution.")

  for (xz in names(intTerms)) {
    elems <- intTerms[[xz]]
    vars  <- diag(vcov[elems, elems])
    sds   <- sqrt(vars)

    rowsXZ <- parTable[parTable$rhs == xz & parTable$op == "~", , drop = FALSE]
    y <- rowsXZ$lhs[[1]]

    if (!length(y)) {
      warning2("No endogenous variable found for interaction term '", xz, "'.",
               immediate. = FALSE)
      next
    }

    # Find the relevant exogenous variables
    struct <- parTable[parTable$lhs == y & parTable$op == "~", , drop = FALSE]
    xis    <- struct$rhs

    # Unstandardized terms
    vary  <- vcov[y, y]
    sdy   <- sqrt(vary)
    B3    <- getCoefs(y = y, x = xz, parTable = parTable)

    # Correlations
    combosXis <- getUniqueCombos(xis)
    combosXis <- combosXis[combosXis[[1]] == xz |
                           combosXis[[2]] == xz, , drop = FALSE]

    lxis <- combosXis[[1]]
    rxis <- combosXis[[2]]

    corrs <- vapply(seq_along(lxis), FUN.VALUE = numeric(1L), 
                    FUN = \(i) vcor[lxis[i], rxis[i]])

    # Incorrectly standardized terms
    lcoefsIncorrect <- getCoefs(x = lxis, y = y, parTable = parTable.std)
    rcoefsIncorrect <- getCoefs(x = rxis, y = y, parTable = parTable.std)

    corrtermsIncorrect <- sum(2 * lcoefsIncorrect * rcoefsIncorrect * corrs)
    
    b3Incorrect <- getCoefs(y = y, x = xz, parTable = parTable.std)
    projVarY_XZ <- b3Incorrect ^ 2 + corrtermsIncorrect # this should be the same 
                                                        # for both the correctly, and
                                                        # incorrectly standardized terms
                                                        # and is the identity which makes it
                                                        # possible to standardize the terms correctly

    # Correctly standardized terms
    b3Correct <- B3 * abs(prod(sds) / sdy) # in case some variances are negative
                                           # we want to make sure we don't flip
                                           # the sign...
    
    lcoefsCorrect <- lcoefsIncorrect
    rcoefsCorrect <- rcoefsIncorrect
    lcoefsCorrect[lxis == xz] <- b3Correct
    rcoefsCorrect[rxis == xz] <- b3Correct
    corrtermsCorrect <- sum(2 * lcoefsCorrect * rcoefsCorrect * corrs)

    # Calculate the correct standard deviation of the interaction term
    # Using the identity:
    #   projVarY_XZ = b3Correct ^ 2 * sd(xz) ^ 2 + sd(xz) * corrterms
    # Solve for sd(xz) using the quadratic formula:
    #   sd(xz) = (- corrterms +/- sqrt(corrterms^2 + 4 * b3Correct ^ 2 * projVarY_XZ)) / 2 * b3Correct ^ 2
    numerator <- -corrtermsCorrect + sign(projVarY_XZ) * 
      sqrt(corrtermsCorrect^2 + 4 * (b3Correct ^ 2) * projVarY_XZ)
    denominator <- 2 * (b3Correct ^ 2)
    sdXZ <- numerator / denominator # correctly standardized sd(xz)

    # Correct parTable
    # We apply scaling factors to the parameters and the standard errors.
    # In theory we should use the delta method to calculate the new standard errors
    # but in practice it shouldn't really matter...
    lequalY  <- parTable.std$lhs == y
    requalXZ <- parTable.std$rhs == xz
    lequalXZ <- parTable.std$lhs == xz
    isCov    <- parTable.std$op == "~~"
    isCoef   <- parTable.std$op == "~"

    isCovTerm  <- xor(lequalXZ, requalXZ) & isCov
    isCoefTerm <- lequalY & requalXZ & isCoef
    isVarTerm  <- lequalXZ & requalXZ & isCov

    # Scaling factors
    scalefVar  <- (sdXZ^2) / 1
    scalefCoef <- b3Correct / b3Incorrect
    scalefCov  <- sdXZ

    # Apply
    parTable.std[isCovTerm, cols]  <- parTable.std[isCovTerm, cols]  * scalefCov
    parTable.std[isCoefTerm, cols] <- parTable.std[isCoefTerm, cols] * scalefCoef
    parTable.std[isVarTerm, cols]  <- parTable.std[isVarTerm, cols]  * scalefVar
  }

  # Recalculate custom parameters
  constrExprs <- sortConstrExprsFinalPt(parTable.std)
  parTable.std <- parTable.std[parTable.std$op != ":=", ]

  errorWarning <- function(e) {
    warning2("Calculation of custom parameter failed: ", e)
    NA
  }

  for (i in seq_len(NROW(constrExprs))) {
    row <- constrExprs[i, , drop=FALSE]

    expr      <- parse(text=constrExprs[i, "rhs"])
    labelList <- parTableLabelsToList(parTable.std) # must be updated for each iteration
  
    oldVal <- row$est
    newVal <- tryCatch(eval(expr, envir = labelList), error = errorWarning)

    ratio  <- newVal / oldVal
    values <- row[ , cols]
   
    row[, cols] <- values * ratio

    parTable.std <- rbind(parTable.std, row)
  }

  etas <- unique(parTable.std[parTable.std$op == "~", "lhs"])
  varEtas <- calcVarParTable(etas, parTable.std)

  warnif(
    any(abs(varEtas - 1) > 1e-10), 
    "Some variances are not equal to 1! ",
    "This indicates that the solution was not standardized correctly!",
    immediate. = FALSE
  )

  attr(parTable.std, "var.etas") <- varEtas
  parTable.std
}
