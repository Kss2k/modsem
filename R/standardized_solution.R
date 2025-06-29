standardizedSolutionCOEFS <- function(object, 
                                      monte.carlo = FALSE, 
                                      mc.reps = 10000, 
                                      tolerance.zero = 1e-10, 
                                      seed = 123, 
                                      delta.epsilon = 1e-8, 
                                      grouping = NULL, 
                                      ...) {
  set.seed(seed) # only relevant if monte.carlo is TRUE

  stopif(!inherits(object, c("modsem_da", "modsem_pi")), 
         "The model must be of class 'modsem_da' or 'modsem_pi'!")

  parTable <- parameter_estimates(object, colon.pi = TRUE)
  parTable <- subsetByGrouping(parTable, grouping = grouping) # if NULL no subsetting

  if (!NROW(parTable)) return(NULL)

  if (inherits(object, "modsem_da")) {
    parTable <- parTable[c("lhs", "op", "rhs", "label", "est", "std.error")]
    parTable <- centerInteraction(parTable) # re-estimate path-coefficients 
    parTable <- parTable[parTable$op != "~1", ] # when intercepts are zero
    parTable <- var_interactions(removeInteractionVariances(parTable))

  } else { # modsem_pi
    if (!"label" %in% names(parTable)) parTable$label <- ""
    parTable <- parTable[c("lhs", "op", "rhs", "label", "est", "se")]
    parTable <- rename(parTable, se = "std.error")
    parTable <- parTable[parTable$op != "~1", ] # when intercepts are zero
  }

  lVs      <- getLVs(parTable)
  intTerms <- getIntTerms(parTable)
  etas     <- getSortedEtas(parTable, isLV = FALSE)
  xis      <- getXis(parTable, etas = etas, isLV = FALSE)
  indsLVs  <- getIndsLVs(parTable, lVs)
  allInds  <- unique(unlist(indsLVs))

  originalLabels <- parTable$label
  labels         <- getParTableLabels(parTable, labelCol="label")
  parTable$label <- labels
  parTable$std.error <- NA

  # Get vcov and coefs
  V     <- vcov(object)
  coefs <- structure(parTable$est, names = labels)
 
  # Get unique labels, and legal names
  labels     <- unique(labels)
  legalNames <- stringr::str_replace_all(labels, OP_REPLACEMENTS)

  # Subset and expand based on unique labels
  V     <- expandVCOV(V, labels=labels)
  coefs <- coefs[labels]

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

  # Unstandardized copies
  parTable.ustd  <- parTable
  COEFS.ustd     <- COEFS

  # get variances
  vars <- unique(c(allInds, lVs, intTerms, xis, etas))
  varianceEquations <- structure(getCovEqExprs(
    x = vars, 
    y = vars, 
    parTable = parTable,
    measurement.model = TRUE
  ), names = vars)
  variances <- lapply(varianceEquations, FUN = \(eq) eval(eq, envir = COEFS))

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
  COEFS <- correctStdSolutionCOEFS(
    parTable = parTable, # for generating equations
    COEFS.std = COEFS,
    COEFS.ustd = COEFS.ustd,
    variances = variances,
    intTerms = intTerms
  )

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

  cleanedParamLabels <- stringr::str_replace_all(params, OP_REPLACEMENTS_INV)
  names(coefs) <- cleanedParamLabels 
  colnames(vcov) <- rownames(vcov) <- cleanedParamLabels

  # Finalize parTable
  parTable[!parTable$label %in% originalLabels, "label"] <- "" 
  parTable$z.value  <- parTable$est / parTable$std.error
  parTable$p.value  <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - CI_WIDTH * parTable$std.error
  parTable$ci.upper <- parTable$est + CI_WIDTH * parTable$std.error

  list(parTable  = parTable,
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
    elems <- stringr::str_split_fixed(XZ, ":", 2)
    X     <- elems[[1]]
    Z     <- elems[[2]]

    vars  <- as.data.frame(variances[c(X, Z)])
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
