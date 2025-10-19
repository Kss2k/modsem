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

  if (!"group" %in% names(parTable)) {
    parTable$group <- 1L
  }

  if (!NROW(parTable)) return(NULL)

  if (isDA || isMplus) {
    cols.keep <- c("lhs", "op", "rhs", "label", "est", "std.error")
    if ("group" %in% names(parTable)) cols.keep <- c(cols.keep, "group")
    cols.keep <- intersect(cols.keep, names(parTable))
    parTable <- parTable[cols.keep]

  } else { # modsem_pi or lavaan
    if (!"label" %in% names(parTable)) parTable$label <- ""
    if (!"se"    %in% names(parTable)) parTable$se    <- NA

    cols.keep <- c("lhs", "op", "rhs", "label", "est", "se")
    if ("group" %in% names(parTable)) cols.keep <- c(cols.keep, "group")
    cols.keep <- intersect(cols.keep, names(parTable))
    parTable <- parTable[cols.keep]
    parTable <- rename(parTable, se = "std.error")
  }

  if (center && (isLav || isDA || isMplus)) { # not relevant for modsem_pi
    if (isDA || isMplus)
      parTable <- meanInteractions(parTable) # get means for interaction terms

    # if we had to center the solution, we have to replace the existing variances
    # if we're using LMS/QML/Mplus there aren't any variances
    isNonCentered <- isNonCenteredParTable(parTable)
    missingVars   <- !hasIntTermVariances(parTable)
    addVariances  <- (isNonCentered && isLav) || isMplus || isDA || missingVars

    if (addVariances) {
      warnif(isLav, "Replacing interaction (co-)", "variances when centering the model!\n", immediate. = FALSE)
      parTable <- var_interactions(parTable, ignore.means = TRUE, mc.reps = mc.reps)
    }
  }

  originalLabels <- parTable$label
  labels         <- getParTableLabels(parTable, labelCol="label")
  labels.clean   <- getParTableLabels(parTable, labelCol="label", replace.dup = TRUE)
  parTable$label <- labels
  parTable$std.error <- NA

  groups <- getGroupsParTable(parTable)

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
  if ("group" %in% names(parTable)) {
    parTable <- parTable[c("lhs", "op", "rhs", "label", "group")]
  } else {
    parTable <- parTable[c("lhs", "op", "rhs", "label")]
  }

  # Center interactions
  if (center) {
    for (g in groups) {
      mask_g <- parTable$group == g
      if (!any(mask_g, na.rm = TRUE)) next
      parTable_g <- parTable[mask_g, , drop = FALSE]
      COEFS <- centerInteractionsCOEFS(parTable_g, COEFS = COEFS) # re-estimate path-coefficients
    }
    parTable <- parTable[!(parTable$op %in% c("~1", "|") &
                           parTable$group %in% groups), , drop = FALSE]
  }

  # Unstandardized copies
  parTable.ustd  <- parTable
  COEFS.ustd     <- COEFS

  if (standardize) for (g in groups) {
    parTable_g <- parTable[parTable$group == g, , drop = FALSE]
    if (!NROW(parTable_g)) next

    lVs_g      <- getLVs(parTable_g)
    intTerms_g <- getIntTerms(parTable_g)
    etas_g     <- getSortedEtas(parTable_g, isLV = FALSE)
    xis_g      <- getXis(parTable_g, etas = etas_g, isLV = FALSE)
    indsLVs_g  <- getIndsLVs(parTable_g, lVs_g)
    allInds_g  <- unique(unlist(indsLVs_g))

    vars_g <- unique(c(allInds_g, lVs_g, intTerms_g, xis_g, etas_g))
    varianceEquations_g <- structure(
      getCovEqExprs(
        x = vars_g,
        y = vars_g,
        parTable = parTable_g,
        measurement.model = TRUE
      ),
      names = vars_g
    )
    variances_g <- lapply(varianceEquations_g, FUN = \(eq) eval(eq, envir = COEFS))

    # Factor Loadings
    for (lV in lVs_g) {
      inds_lV <- indsLVs_g[[lV]]
      if (!length(inds_lV)) next

      for (ind in inds_lV) {
        selectRows  <- parTable_g$lhs == lV & parTable_g$op == "=~" & parTable_g$rhs == ind
        if (!any(selectRows)) next
        label <- parTable_g[selectRows, "label"]

        var_lV <- variances_g[[lV]]
        var_ind <- variances_g[[ind]]
        if (is.null(var_lV) || is.null(var_ind)) next

        scalingCoef <- sqrt(var_lV) / sqrt(var_ind)
        lambda      <- COEFS[[label]] * scalingCoef

        COEFS[[label]] <- lambda
      }
    }

    # Structural Coefficients
    selectStrucExprs <- parTable_g$op == "~" & parTable_g$lhs %in% etas_g

    for (eta in etas_g) {
      selectStrucExprsEta <- selectStrucExprs & parTable_g$lhs == eta
      structExprsEta      <- parTable_g[selectStrucExprsEta, ]

      for (xi in structExprsEta$rhs) {
        selectRows  <- selectStrucExprsEta & parTable_g$rhs == xi
        if (!any(selectRows)) next
        var_xi <- variances_g[[xi]]
        var_eta <- variances_g[[eta]]
        if (is.null(var_xi) || is.null(var_eta)) next

        scalingCoef <- sqrt(var_xi) / sqrt(var_eta)
        label       <- parTable_g[selectRows, "label"]
        gamma       <- COEFS[[label]] * scalingCoef

        COEFS[[label]] <- gamma
      }
    }

    # (Co-) Variances of xis
    selectCovXis <- parTable_g$op == "~~" &
      (parTable_g$lhs %in% c(xis_g, intTerms_g) |
         parTable_g$rhs %in% c(xis_g, intTerms_g))

    covRowsXis <- parTable_g[selectCovXis, , drop = FALSE]
    for (i in seq_len(nrow(covRowsXis))) {
      lhs <- covRowsXis$lhs[[i]]
      rhs <- covRowsXis$rhs[[i]]
      xis_pair <- c(lhs, rhs)
      selectRows  <- selectCovXis &
        parTable_g$lhs %in% xis_pair &
        parTable_g$rhs %in% xis_pair

      var_lhs <- variances_g[[lhs]]
      var_rhs <- variances_g[[rhs]]
      if (is.null(var_lhs) || is.null(var_rhs)) next

      scalingCoef <- sqrt(var_lhs) * sqrt(var_rhs)

      if (lhs != rhs) {
        selectRows <- selectRows & parTable_g$lhs != parTable_g$rhs
      }

      label <- parTable_g[selectRows, "label"]
      covs <- COEFS[[label]] / scalingCoef

      COEFS[[label]] <- covs
    }

    # Residual Variances etas
    for (eta in etas_g) {
      selectRows <- parTable_g$lhs == eta & parTable_g$op == "~~" & parTable_g$rhs == eta
      if (!any(selectRows)) next
      var_eta <- variances_g[[eta]]
      if (is.null(var_eta)) next
      label <- parTable_g[selectRows, "label"]
      residual <- COEFS[[label]] / var_eta

      COEFS[[label]] <- residual
    }

    # residual variances inds
    for (ind in allInds_g) {
      selectRows <- parTable_g$lhs == ind & parTable_g$op == "~~" & parTable_g$rhs == ind
      if (!any(selectRows)) next
      var_ind <- variances_g[[ind]]
      if (is.null(var_ind)) next

      label <- parTable_g[selectRows, "label"]
      residual <- COEFS[[label]] / var_ind

      COEFS[[label]] <- residual
    }

    # Correct Scale of interaction terms
    COEFS <- correctStdSolutionCOEFS(
      parTable = parTable_g, # for generating equations
      COEFS.std = COEFS,
      COEFS.ustd = COEFS.ustd,
      variances = variances_g,
      intTerms = intTerms_g
    )
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
                                      merge.by = c("lhs", "op", "rhs", "group"),
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
