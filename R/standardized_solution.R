transformedSolutionCOEFS <- function(object,
                                     monte.carlo = FALSE,
                                     mc.reps = 10000,
                                     tolerance.zero = 1e-10,
                                     delta.epsilon = 1e-8,
                                     grouping = NULL,
                                     center = TRUE,
                                     standardize = TRUE,
                                     ...) {
  stopif(!inherits(object, c("modsem_da", "modsem_pi", "lavaan", "modsem_mplus")),
         "The model must be of class `modsem_da`, `modsem_mplus`, `modsem_pi` or `lavaan`!")

  isLav   <- inherits(object, "lavaan")
  isDA    <- inherits(object, "modsem_da")
  isMplus <- inherits(object, "modsem_mplus")

  if (isLav) {
    vcov <- lavaan::vcov # load vcov and coef from lavaan if dealing with a lavaan object
    coef <- lavaan::coef
  }

  parTable <- parameter_estimates(object, colon.pi = TRUE, high.order.as.measr = FALSE)
  parTable <- subsetByGrouping(parTable, grouping = grouping) # if NULL no subsetting

  if (!NROW(parTable)) return(NULL)

  if (isDA || isMplus) {
    parTable <- parTable[c("lhs", "op", "rhs", "label", "est", "std.error")]

  } else { # modsem_pi or lavaan
    if (!"label" %in% names(parTable)) parTable$label <- ""
    if (!"se"    %in% names(parTable)) parTable$se    <- NA

    parTable <- parTable[c("lhs", "op", "rhs", "label", "est", "se")]
    parTable <- rename(parTable, se = "std.error")
  }

  if (center && (isLav || isDA || isMplus)) { # not relevant for modsem_pi
    warnif(isLav, "Replacing interaction (co-)",
           "variances when centering the model!\n", immediate. = FALSE)

    if (isDA || isMplus)
      parTable <- meanInteractions(parTable) # get means for interaction terms

    parTable <- var_interactions(parTable, ignore.means = TRUE, mc.reps = mc.reps)
  }

  lVs      <- getLVs(parTable)
  intTerms <- getIntTerms(parTable)
  etas     <- getSortedEtas(parTable, isLV = FALSE)
  xis      <- getXis(parTable, etas = etas, isLV = FALSE)
  indsLVs  <- getIndsLVs(parTable, lVs)
  allInds  <- unique(unlist(indsLVs))

  originalLabels <- parTable$label
  labels         <- getParTableLabels(parTable, labelCol="label")
  labels.clean   <- getParTableLabels(parTable, labelCol="label", replace.dup = TRUE)
  parTable$label <- labels
  parTable$std.error <- NA

  # Get vcov and coefs
  V     <- tryCatch(vcov(object), error = \(e) NULL)
  coefs <- structure(parTable$est, names = labels)

  if (is.null(V)) { # calc.se == FALSE
    k    <- length(coefs)
    pars <- names(coefs)
    V    <- matrix(0, nrow = k, ncol = k,
                   dimnames = list(pars, pars))
  }

  # Subset and expand based on unique labels
  labels.u <- unique(labels)
  V        <- expandVCOV(V, labels=labels.u)
  coefs    <- coefs[labels.u]

  # Replace labels with cleaned versions
  # E.g., if two parameters share the same label, we want to split them int two different labels
  parTable$label <- labels.clean
  legalNames     <- stringr::str_replace_all(labels.clean, OP_REPLACEMENTS)
  paramMapping   <- structure(labels, names = legalNames)

  V       <- V[labels, labels]
  coefs   <- coefs[labels]

  dimnames(V)  <- list(labels.clean, labels.clean)
  names(coefs) <- labels.clean

  if (monte.carlo) {
    COEFS <- as.data.frame(mvtnorm::rmvnorm(mc.reps, mean = coefs, sigma = V))
  } else { # delta method
    k <- length(coefs)
    COEFS <- matrix(coefs, nrow=k, ncol=k, byrow=TRUE)
    DeltaEpsilon <- diag(delta.epsilon, k)
    isFixed <- abs(diag(V)) < tolerance.zero
    DeltaEpsilon[isFixed, isFixed] <- 0
    COEFS_p <- as.data.frame(COEFS + DeltaEpsilon)
    COEFS_m <- as.data.frame(COEFS - DeltaEpsilon)
    COEFS <- rbind(COEFS_p, COEFS_m)
  }

  # Get legal parameter names
  names(COEFS) <- legalNames
  names(coefs) <- legalNames
  COEFS <- rbind(as.data.frame(as.list(coefs)), COEFS) # first row is the original values

  colnames(V) <- legalNames
  rownames(V) <- legalNames

  parTable$label <- stringr::str_replace_all(parTable$label, OP_REPLACEMENTS)
  parTable       <- parTable[c("lhs", "op", "rhs", "label")]

  # Center interactions
  if (center) {
    COEFS <- centerInteractionsCOEFS(parTable, COEFS = COEFS) # re-estimate path-coefficients
    parTable <- parTable[!parTable$op %in% c("~1", "|"), ]
  }

  # Unstandardized copies
  parTable.ustd  <- parTable
  COEFS.ustd     <- COEFS

  if (standardize) {
    # get variances
    vars <- unique(c(allInds, lVs, intTerms, xis, etas))
    varianceEquations <- structure(getCovEqExprs(
       x = vars,
       y = vars,
       parTable = parTable,
       measurement.model = TRUE
       ), names = vars)
    variances <- lapply(varianceEquations, FUN = \(eq) eval(eq, envir = COEFS))

    # Center interaction terms

    # Factor Loadings
    lambda     <- NULL
    selectRows <- NULL

    for (lV in lVs) {
      for (ind in indsLVs[[lV]]) {
        selectRows  <- parTable$lhs == lV & parTable$op == "=~" & parTable$rhs == ind
        label <- parTable[selectRows, "label"]

        # est in parTable
        scalingCoef <- sqrt(variances[[lV]]) / sqrt(variances[[ind]])
        lambda      <- COEFS[[label]] * scalingCoef

        COEFS[[label]] <- lambda
      }
    }

    # Structural Coefficients
    gamma               <- NULL
    selectStrucExprsEta <- NULL
    structExprsEta      <- NULL
    selectStrucExprs    <- parTable$op == "~" & parTable$lhs %in% etas

    for (eta in etas) {
      selectStrucExprsEta <- selectStrucExprs & parTable$lhs == eta
      structExprsEta      <- parTable[selectStrucExprsEta, ]

      for (xi in structExprsEta$rhs) {
        selectRows  <- selectStrucExprsEta & parTable$rhs == xi
        scalingCoef <- sqrt(variances[[xi]]) / sqrt(variances[[eta]])
        label       <- parTable[selectRows, "label"]
        gamma       <- COEFS[[label]] * scalingCoef

        COEFS[[label]] <- gamma
      }
    }

    # (Co-) Variances of xis
    selectCovXis <- parTable$op == "~~" &
      (parTable$lhs %in% c(xis, intTerms) | parTable$lhs %in% c(xis, intTerms))

    covRowsXis <- parTable[selectCovXis, , drop = FALSE]
    for (i in seq_len(nrow(covRowsXis))) {
      lhs         <- covRowsXis$lhs[[i]]
      rhs         <- covRowsXis$rhs[[i]]
      xis         <- c(lhs, rhs)
      selectRows  <- selectCovXis & parTable$lhs %in% xis & parTable$rhs %in% xis
      scalingCoef <- sqrt(variances[[lhs]]) * sqrt(variances[[rhs]])

      if (lhs != rhs) selectRows <- selectRows & parTable$lhs != parTable$rhs

      label <- parTable[selectRows, "label"]
      covs <- COEFS[[label]] / scalingCoef

      COEFS[[label]] <- covs
    }

    # Residual Variances etas
    selectRows <- NULL
    residual   <- NULL
    for (eta in etas) {
      selectRows <- parTable$lhs == eta & parTable$op == "~~" & parTable$rhs == eta
      label <- parTable[selectRows, "label"]
      residual <- COEFS[[label]] / variances[[eta]]

      COEFS[[label]] <- residual
    }

    # residual variances inds
    for (ind in allInds) {
      selectRows <- parTable$lhs == ind & parTable$op == "~~" & parTable$rhs == ind
      label <- parTable[selectRows, "label"]
      residual <- COEFS[[label]] / variances[[ind]]

      COEFS[[label]] <- residual
    }

    # Correct Scale of interaction terms
    COEFS <- correctStdSolutionCOEFS(parTable = parTable, # for generating equations
                                     COEFS.std = COEFS,
                                     COEFS.ustd = COEFS.ustd,
                                     variances = variances,
                                     intTerms = intTerms)

  }
  # recalculate custom parameters
  constrExprs <- sortConstrExprsFinalPt(parTable)
  for (i in seq_len(NROW(constrExprs))) {
    row <- constrExprs[i, , drop=FALSE]
    label <- row$label
    expr <- parse(text=constrExprs[i, "rhs"])
    newVals <- eval(expr, envir = COEFS)

    COEFS[[label]] <- newVals
  }

  nzeros <- round(log10(1 / tolerance.zero), 0)
  coefs <- unlist(COEFS[1, ])
  COEFS <- round(COEFS[2:(nrow(COEFS)), ],  nzeros) # skip first row

  if (monte.carlo) {
    vcov <- stats::cov(COEFS)
  } else {
    # delta method
    p <- 1:k
    m <- (k+1):(k + k)

    COEFS_p <- t(COEFS[p, ])
    COEFS_m <- t(COEFS[m, ])

    J <- (COEFS_p - COEFS_m) / (2 * delta.epsilon)
    vcov <- J %*% V %*% t(J)
  }

  # fill parTable
  std.errors <- suppressWarnings(sqrt(diag(vcov)))
  warnFunc <- function(type, row) {
    warning2("Unable to calculate standardized ", type, " for: ",
             paste0(row$lhs, row$op, row$rhs))
  }

  verboseLabels <- stringr::str_replace_all(parTable$label, OP_REPLACEMENTS)
  for (i in seq_len(nrow(parTable))) {
    row <- parTable[i, , drop = FALSE]
    label <- verboseLabels[[i]]

    if (!label %in% names(coefs)) {
      warnFunc("coefficient", row)
      est <- NA
    } else est <- coefs[[label]]

    if (!label %in% names(std.errors)) {
      warnFunc("std.error", row)
      se <- NA
    } else se <- std.errors[[label]]

    parTable[i, "est"] <- est
    parTable[i, "std.error"] <- se
  }

  # Fill in NA on zero-standard errors, and create vcov and coefs
  parTable[!is.na(parTable$std.error) &
           abs(parTable$std.error) < tolerance.zero, "std.error"] <- NA

  isFree <- !is.na(parTable$std.error)
  coefs <- structure(parTable$est[isFree], names = parTable$label[isFree])
  params <- intersect(names(coefs), rownames(vcov))

  coefs <- coefs[params]
  vcov <- vcov[params, params]

  cleanedParamLabels <- paramMapping[params]
  names(coefs) <- cleanedParamLabels
  colnames(vcov) <- rownames(vcov) <- cleanedParamLabels

  # Finalize parTable
  parTable$label    <- paramMapping[parTable$label]
  parTable$z.value  <- parTable$est / parTable$std.error
  parTable$p.value  <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - CI_WIDTH * parTable$std.error
  parTable$ci.upper <- parTable$est + CI_WIDTH * parTable$std.error

  # Remove added labels
  labelInOrig       <- parTable$label %in% originalLabels
  parTable[!labelInOrig, "label"] <- ""

  if (isDA) {
    indsHigherOrderLVs <- object$model$info$indsHigherOrderLVs
    parTable <- higherOrderStruct2Measr(parTable = parTable,
                                        indsHigherOrderLVs = indsHigherOrderLVs)
    parTable <- sortParTableDA(parTable, model = object$model)
  }

  # Reset index
  rownames(parTable) <- NULL

  list(parTable  = modsemParTable(parTable),
       coefs     = coefs,
       vcov      = vcov,
       COEFS     = COEFS)
}


correctStdSolutionCOEFS <- function(parTable,
                                    COEFS.std,
                                    COEFS.ustd,
                                    variances,
                                    intTerms) {
  for (XZ in intTerms) {
    elems <- stringr::str_split(XZ, pattern = ":")[[1L]]

    vars  <- as.data.frame(variances[elems])
    sds   <- sqrt(vars)

    rowsXZ <- parTable[parTable$rhs == XZ & parTable$op == "~", , drop = FALSE]
    Y <- rowsXZ$lhs[[1]]

    if (!length(Y)) {
      warning2("No endogenous variable found for interaction term '", XZ, "'.",
               immediate. = FALSE)
      next
    }

    # Find the relevant exogenous variables
    struct <- parTable[parTable$lhs == Y & parTable$op == "~", , drop = FALSE]
    xis    <- struct$rhs

    # Unstandardized terms
    varY  <- variances[[Y]]
    sdY   <- sqrt(varY)
    B3    <- getCOEFS(y = Y, x = XZ, COEFS = COEFS.ustd, parTable = parTable)

    # Correlations
    combosXis <- getUniqueCombos(xis)
    combosXis <- combosXis[combosXis[[1]] == XZ |
                           combosXis[[2]] == XZ, , drop = FALSE]

    lxis <- combosXis[[1]]
    rxis <- combosXis[[2]]

    eqCorrs <- getCovEqExprs(x = lxis, y = rxis, parTable = parTable)
    corrs <- matrix(NA, nrow = NROW(COEFS.std), ncol = length(lxis))

    for (i in seq_len(ncol(corrs))) {
      eqCorr <- eqCorrs[[i]]
      corr <- eval(eqCorr, envir = COEFS.std)
      corr[is.na(corr)] <- 0 # in case there is no connection, it should be zero

      corrs[, i] <- corr
    }

    # Incorrectly standardized terms
    lcoefsIncorrect <- getCOEFS(x = lxis, y = Y, COEFS = COEFS.std, parTable = parTable)
    rcoefsIncorrect <- getCOEFS(x = rxis, y = Y, COEFS = COEFS.std, parTable = parTable)

    corrtermsIncorrect <- rowSums(2 * lcoefsIncorrect * rcoefsIncorrect * corrs)

    b3Incorrect <- getCOEFS(y = Y, x = XZ, COEFS = COEFS.std, parTable = parTable)
    projVarY_XZ <- b3Incorrect ^ 2 + corrtermsIncorrect # this should be the same
                                                        # for both the correctly, and
                                                        # incorrectly standardized terms
                                                        # and is the identity which makes it
                                                        # possible to standardize the terms correctly

    # Correctly standardized terms
    rowProds <- apply(sds, MARGIN = 1, FUN = prod)
    b3Correct <- B3 * abs(rowProds / sdY) # in case some variances are negative
                                           # we want to make sure we don't flip
                                           # the sign...

    lcoefsCorrect <- lcoefsIncorrect
    rcoefsCorrect <- rcoefsIncorrect
    if (any(lxis == XZ)) lcoefsCorrect[lxis == XZ] <- b3Correct
    if (any(rxis == XZ)) rcoefsCorrect[rxis == XZ] <- b3Correct
    corrtermsCorrect <- rowSums(2 * lcoefsCorrect * rcoefsCorrect * corrs)

    # Calculate the correct standard deviation of the interaction term
    # Using the identity:
    #   projVarY_XZ = b3Correct ^ 2 * sd(xz) ^ 2 + sd(xz) * corrterms
    # Solve for sd(xz) using the quadratic formula:
    #   sd(xz) = (- corrterms +/- sqrt(corrterms^2 + 4 * b3Correct ^ 2 * projVarY_XZ)) / 2 * b3Correct ^ 2
    numerator <- -corrtermsCorrect + sign(projVarY_XZ) *
      sqrt(corrtermsCorrect^2 + 4 * (b3Correct ^ 2) * projVarY_XZ)
    denominator <- 2 * (b3Correct ^ 2)
    sdXZ <- numerator / denominator # correctly standardized sd(xz)

    # Scaling factors
    scalefVar  <- (sdXZ^2) / 1
    scalefCoef <- b3Correct / b3Incorrect
    scalefCov  <- sdXZ

    lequalY  <- parTable$lhs == Y
    requalXZ <- parTable$rhs == XZ
    lequalXZ <- parTable$lhs == XZ
    isCov    <- parTable$op == "~~"
    isCoef   <- parTable$op == "~"

    covTerms <- parTable$label[xor(lequalXZ, requalXZ) & isCov]
    coefTerm <- parTable$label[lequalY & requalXZ & isCoef]
    varTerm  <- parTable$label[lequalXZ & requalXZ & isCov]

    COEFS.std[covTerms] <- COEFS.std[covTerms] * scalefCov
    COEFS.std[coefTerm] <- COEFS.std[coefTerm] * scalefCoef
    COEFS.std[varTerm]  <- COEFS.std[varTerm]  * scalefVar
  }

  COEFS.std
}


standardizedSolutionCOEFS <- function(object,
                                      monte.carlo = FALSE,
                                      mc.reps = 10000,
                                      tolerance.zero = 1e-10,
                                      delta.epsilon = 1e-8,
                                      grouping = NULL,
                                      center = TRUE,
                                      ...) {
  transformedSolutionCOEFS(
    object = object,
    monte.carlo = monte.carlo,
    mc.reps = mc.reps,
    tolerance.zero = tolerance.zero,
    delta.epsilon = delta.epsilon,
    grouping = grouping,
    standardize = TRUE,
    ...
  )
}


centeredSolutionCOEFS <- function(object,
                                  monte.carlo = FALSE,
                                  mc.reps = 10000,
                                  tolerance.zero = 1e-10,
                                  delta.epsilon = 1e-8,
                                  grouping = NULL,
                                  ...) {
  transformedSolutionCOEFS(
    object         = object,
    monte.carlo    = monte.carlo,
    mc.reps        = mc.reps,
    tolerance.zero = tolerance.zero,
    delta.epsilon  = delta.epsilon,
    standardize    = FALSE,
    center         = TRUE,
    ...
  )
}


getMeanFormula <- function(x, parTable, label.col = "label") {
  stopif(length(x) > 1, "x must be a single string")

  meanY <- getIntercept(x, parTable = parTable, col = label.col)
  gamma <- parTable[parTable$lhs == x & parTable$op == "~", , drop = FALSE]

  if (NROW(gamma) == 0) return(meanY)
  for (i in seq_len(NROW(gamma))) {
    meanX <- getMeanFormula(gamma[i, "rhs"], parTable = parTable)
    meanY <- paste0("(", meanY, "+", gamma[i, label.col], "*", meanX, ")")
  }

  if (!length(meanY)) "0" else paste0("(", meanY, ")")
}


centerInteractionsCOEFS <- function(parTable, COEFS, center.means = TRUE,
                                    label.col = "label") {
  rows <- getIntTermRows(parTable)

  for (i in seq_len(NROW(rows))) {
    Y <- rows[i, "lhs"]
    XZ <- unlist(stringr::str_split(rows[i, "rhs"], ":"))
    X <- XZ[[1]]
    Z <- XZ[[2]]
    gamma <- parTable[parTable$lhs == Y & parTable$op == "~", , drop = FALSE]

    formulaMeanX <- parse(text = getMeanFormula(X, parTable = parTable))
    formulaMeanZ <- parse(text = getMeanFormula(Z, parTable = parTable))

    labelGammaXZ <- rows[i, label.col]
    labelGammaX  <- gamma[gamma$rhs == X, label.col] # length should always be 1, but just in case...
    labelGammaZ  <- gamma[gamma$rhs == Z, label.col]

    meanX   <- eval(formulaMeanX, envir = COEFS)
    meanZ   <- eval(formulaMeanZ, envir = COEFS)
    gammaXZ <- COEFS[[labelGammaXZ]]

    if (length(labelGammaX) == 1) {
      gammaX  <- COEFS[[labelGammaX]]
      COEFS[[labelGammaX]] <- gammaX + gammaXZ * meanZ
    }

    if (length(labelGammaZ) == 1) {
      gammaZ  <- COEFS[[labelGammaZ]]
      COEFS[[labelGammaZ]] <- gammaZ + gammaXZ * meanX
    }
  }

  if (center.means) {
    innerVars <- unique(unlist(parTable[parTable$op == "~", c("rhs", "lhs")]))
    interceptLabels <- parTable[parTable$lhs %in% innerVars &
                                parTable$op == "~1", label.col]

    for (label in interceptLabels)
      COEFS[[label]] <- 0
  }

  COEFS
}


addTransformedEstimatesPT <- function(parTable,
                                      FUN,
                                      pass.parTable = TRUE,
                                      values.to = "transformed",
                                      values.from = "est",
                                      merge.by = c("lhs", "op", "rhs"),
                                      ...) {
    if (pass.parTable) parTable.transform <- FUN(parTable = parTable, ...)
    else               parTable.transform <- FUN(...)

    parTable.transform[[values.to]] <- parTable.transform[[values.from]]
    parTable.transform <- parTable.transform[c(merge.by, values.to)]

    leftJoin(left   = parTable,
             right  = parTable.transform,
             by     = merge.by)
}


applyTransformationByGrouping <- function(parTable,
                                          FUN,
                                          groupingcols = c("block", "group"),
                                          ...) {
  if (any(groupingcols %in% colnames(parTable))) {
    groupingcols <- intersect(groupingcols, colnames(parTable))
    categories   <- unique(parTable[groupingcols])

    parTable.out <- NULL

    for (i in seq_len(NROW(categories))) {
      grouping       <- structure(unlist(categories[i, ]),
                                  names = groupingcols)
      parTable.out_i <- FUN(..., grouping = grouping)

      if (is.null(parTable.out_i))
        next

      for (group in names(grouping))
        parTable.out_i[[group]] <- grouping[[group]]

      parTable.out <- rbind(parTable.out, parTable.out_i)
    }

    colblock1 <- c("lhs", "op", "rhs")
    colblock2 <- groupingcols
    colblock3 <- setdiff(colnames(parTable.out), c(colblock1, colblock2))

    parTable.out <- parTable.out[c(colblock1, colblock2, colblock3)]

  } else parTable.out <- FUN(...)

  parTable.out
}
