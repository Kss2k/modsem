optimizeStartingParamsDA <- function(model,
                                     args = list(orthogonal.x = FALSE,
                                                 orthogonal.y = FALSE,
                                                 auto.fix.first = TRUE,
                                                 auto.fix.sinlge = TRUE,
                                                 robust.se = FALSE),
                                     group = NULL,
                                     engine = c("pi", "sam")) {
  engine <- tolower(engine)
  engine <- match.arg(engine)

  etas     <- model$info$etas
  indsEtas <- model$info$allIndsEtas
  xis      <- model$info$xis
  numXis   <- model$info$numXis
  indsXis  <- model$info$allIndsXis
  data     <- model$data.raw
  missing  <- tolower(args$missing)

  robust.se       <- args$robust.se
  has.interaction <- model$info$has.interaction

  syntax <- paste(model$syntax, model$models[[1L]]$covModel$syntax,
                  model$info$lavOptimizerSyntaxAdditions, sep = "\n")

  acceptable.missing <- c("listwise", "ml", "direct", "fiml")
  if (has.interaction && missing %in% c("ml", "fiml", "direct")) {
    # It is not worth the extra compute to use fiml (if fiml is set)
    # for getting good starting estiamtes for the model, for non-linear models.
    # FIML can be very slow if there are a lot of product indicators.
    # However, if the model is linear, the estimates should be the
    # exact same, and product indicators aren't an issue. Thus
    # we will save time.
    missing <- "listwise"
  } else if (!missing %in% acceptable.missing) missing <- "listwise"

  if (!has.interaction && robust.se) {
    estimator <- "MLR"

    # Not worth the extra compute for non-linear models
  } else estimator <- "ML"

  if (engine == "pi") {
    estPI <- modsem_pi(
      model.syntax    = syntax,
      data            = data,
      method          = "dblcent",
      estimator       = estimator,
      meanstructure   = TRUE,
      orthogonal.x    = args$orthogonal.x,
      orthogonal.y    = args$orthogonal.y,
      auto.fix.first  = args$auto.fix.first,
      auto.fix.single = args$auto.fix.single,
      res.cov.method  = "simple.no.warn",
      res.cov.across  = TRUE,
      match           = TRUE,
      match.recycle   = TRUE,
      missing         = missing,
      suppress.warnings.match = TRUE,
      suppress.warnings.lavaan = TRUE
    )
    parTable   <- parameter_estimates(estPI, colon.pi = TRUE)
    lavaan.fit <- extract_lavaan(estPI)

  } else if (engine == "sam") {
    fitSAM   <- parameterEstimatesLavSAM(
      syntax          = syntax,
      data            = data,
      estimator       = estimator,
      missing         = missing,
      meanstructure   = TRUE,
      orthogonal.x    = args$orthogonal.x,
      orthogonal.y    = args$orthogonal.y,
      auto.fix.first  = args$auto.fix.first,
      auto.fix.single = args$auto.fix.single,
      group           = group,
      suppress.warnings.lavaan = TRUE,
    )

    parTable   <- fitSAM$parTable
    lavaan.fit <- fitSAM$fit
  }

  stopif(is.null(parTable), "lavaan failed!")

  if (isHigherOrderParTable(parTable))
    parTable <- higherOrderMeasr2Struct(parTable)

  params <- model$params
  SELECT_THETA_LAB  <- params$SELECT_THETA_LAB
  SELECT_THETA_COV  <- params$SELECT_THETA_COV
  SELECT_THETA_MAIN <- params$SELECT_THETA_MAIN
  THETA <- params$theta
    
  thetaLabel <- getLabeledParamsLavaan(parTable, params$constrExprs$fixedParams)
  THETA[SELECT_THETA_LAB[[1L]]][names(thetaLabel)] <- thetaLabel

  for (g in seq_len(model$info$n.groups)) {
    submodel <- model$models[[g]]
    
    if ("group" %in% colnames(parTable))
      parTable_g <- parTable[parTable$group == g, , drop = FALSE]
    else
      parTable_g <- parTable

    fillLabelsMatrix <- function(matNumeric, matLabel, symmetric = FALSE) {
      if (all(matLabel == ""))
        return(matNumeric)

      labels <- c(matLabel[matLabel != ""])
      labels <- labels[labels %in% names(thetaLabel)]

      for (label in labels)
        matNumeric[matLabel == label] <- thetaLabel[[label]]

      if (symmetric)
        matNumeric[upper.tri(matNumeric)] <- t(matNumeric)[upper.tri(matNumeric)]

      matNumeric
    }

    # Main Model
    matricesMain      <- submodel$matrices
    labelMatricesMain <- submodel$labelMatrices

    LambdaX <- findEstimatesParTable(matricesMain$lambdaX, parTable_g, op = "=~",
                                     rows_lhs = FALSE, fill = 0.7)
    LambdaY <- findEstimatesParTable(matricesMain$lambdaY, parTable_g, op = "=~",
                                     rows_lhs = FALSE, fill = 0.7)

    ThetaEpsilon <- findEstimatesParTable(matricesMain$thetaEpsilon, parTable_g,
                                          op = "~~", fill = 0.2)
    ThetaDelta   <- findEstimatesParTable(matricesMain$thetaDelta, parTable_g,
                                          op = "~~", fill = 0.2)

    Psi <- findEstimatesParTable(matricesMain$psi, parTable_g, op = "~~", fill = 0)
    Phi <- findEstimatesParTable(matricesMain$phi, parTable_g, op = "~~", fill = 0)
    A   <- findEstimatesParTable(matricesMain$A, parTable_g, op = "~~", fill = 0)

    # Matrices which can be corrected to ensure viable starting parameters need to
    # get filled in using labels as well, just for the checks them selves
    Psi <- fillLabelsMatrix(Psi, labelMatricesMain$psi, symmetric = TRUE)
    Phi <- fillLabelsMatrix(Phi, labelMatricesMain$phi, symmetric = TRUE)
    A   <- fillLabelsMatrix(A, labelMatricesMain$A, symmetric = FALSE)

    ThetaEpsilon <- fillLabelsMatrix(ThetaEpsilon,
                                     labelMatricesMain$thetaEpsilon,
                                     symmetric = TRUE)

    ThetaDelta <- fillLabelsMatrix(ThetaDelta,
                                   labelMatricesMain$thetaDelta,
                                   symmetric = TRUE)

    correctDiag <- function(M, fill = 1, tol = 0) {
      M[M < tol & is.diag(M)] <- fill
      M
    }

    # Check for negative diagonals
    ThetaEpsilon <- correctDiag(ThetaEpsilon, tol = 0) # no negative values
    ThetaDelta   <- correctDiag(ThetaDelta, tol = 0) # no negative values
    Psi          <- correctDiag(Psi, tol = 0) # no negative values
    Phi          <- correctDiag(Phi, tol = 0) # no negative values
    A            <- correctDiag(A, tol = 0)

    as.I <- function(M) { # If Phi/A is non-invertible we want I instead
      I <- diag(NROW(M))
      dimnames(I) <- dimnames(M)
      I
    }

    if (!is.invertible(Phi)) Phi <- as.I(Phi)
    # Residuals don't need to be invertible...
    # if (!is.invertible(Psi))          Psi <- as.I(Psi)
    # if (!is.invertible(ThetaEpsilon)) ThetaEpsilon <- as.I(ThetaEpsilon)
    # if (!is.invertible(ThetaDelta))   ThetaDelta   <- as.I(ThetaDelta)

    A[upper.tri(A)] <- t(A)[upper.tri(A)]
    A <- t(tryCatch(chol(A), error = function(x) as.I(A)))

    beta0 <- findInterceptsParTable(matricesMain$beta0, parTable_g, fill = 0)
    alpha <- findInterceptsParTable(matricesMain$alpha, parTable_g, fill = 0)

    GammaEta <- findEstimatesParTable(matricesMain$gammaEta, parTable_g, op = "~", fill = 0)
    GammaXi  <- findEstimatesParTable(matricesMain$gammaXi, parTable_g, op = "~", fill = 0)

    OmegaEtaXi <- findInteractionEstimatesParTable(matricesMain$omegaEtaXi,
                                                   parTable = parTable_g, fill = 0)
    OmegaXiXi <- findInteractionEstimatesParTable(matricesMain$omegaXiXi,
                                                  parTable = parTable_g, fill = 0)
    tauX <- findInterceptsParTable(matricesMain$tauX, parTable_g, fill = 0)
    tauY <- findInterceptsParTable(matricesMain$tauY, parTable_g, fill = 0)

    thetaMain <- unlist(list(LambdaX[is.na(matricesMain$lambdaX)],
                             LambdaY[is.na(matricesMain$lambdaY)],
                             tauX[is.na(matricesMain$tauX)],
                             tauY[is.na(matricesMain$tauY)],
                             ThetaDelta[is.na(matricesMain$thetaDelta)],
                             ThetaEpsilon[is.na(matricesMain$thetaEpsilon)],
                             Phi[is.na(matricesMain$phi)],
                             A[is.na(matricesMain$A)],
                             Psi[is.na(matricesMain$psi)],
                             alpha[is.na(matricesMain$alpha)],
                             beta0[is.na(matricesMain$beta0)],
                             GammaXi[is.na(matricesMain$gammaXi)],
                             GammaEta[is.na(matricesMain$gammaEta)],
                             OmegaXiXi[is.na(matricesMain$omegaXiXi)],
                             OmegaEtaXi[is.na(matricesMain$omegaEtaXi)]))

    # Cov Model
    matricesCov      <- model$covModel$matrices
    labelMatricesCov <- model$covModel$labelMatrices

    if (!is.null(matricesCov)) {
      PsiCovModel <- findEstimatesParTable(matricesCov$psi, parTable_g, op = "~~", fill = 0)
      PhiCovModel <- findEstimatesParTable(matricesCov$phi, parTable_g, op = "~~", fill = 0)

      GammaEtaCovModel <- findEstimatesParTable(matricesCov$gammaEta, parTable_g, op = "~", fill = 0)
      GammaXiCovModel <- findEstimatesParTable(matricesCov$gammaXi, parTable_g, op = "~", fill = 0)

      PhiCovModel <- correctDiag(PhiCovModel, tol = 0)
      PsiCovModel <- correctDiag(PsiCovModel, tol = 0)

      PhiCovModel <- fillLabelsMatrix(PhiCovModel, labelMatricesCov$phi, symmetric = TRUE)

      if (!is.invertible(PhiCovModel)) PhiCovModel <- as.I(PhiCovModel)
      # Residuals don't need to be invertible...
      # if (!is.invertible(PsiCovModel)) PsiCovModel <- as.I(PsiCovModel)

      thetaCov <- unlist(list(PhiCovModel[is.na(matricesCov$phi)],
                              PsiCovModel[is.na(matricesCov$psi)],
                              GammaXiCovModel[is.na(matricesCov$gammaXi)],
                              GammaEtaCovModel[is.na(matricesCov$gammaEta)]))
    } else thetaCov <- NULL

    selectThetaMain <- SELECT_THETA_MAIN[[g]]
    selectThetaCov  <- SELECT_THETA_COV[[g]]

    if (length(selectThetaMain) == length(thetaMain))
      THETA[selectThetaMain] <- thetaMain

    if (length(selectThetaCov) == length(thetaCov) && length(thetaCov) > 0L)
      THETA[selectThetaCov] <- thetaCov
  }

  if (length(THETA) == length(model$theta)) {
    names(THETA) <- names(model$theta)
    model$theta <- THETA
  }

  model$lavaan.fit <- lavaan.fit
  model
}


findEstimatesParTable <- function(mat, parTable, op = NULL, rows_lhs = TRUE,
                                  fill = NULL) {
  if (is.null(op)) stop("Missing operator")
  for (row in rownames(mat)) {
    for (col in colnames(mat)) {
      if (is.na(mat[row, col]))
        mat[row, col] <- extractFromParTable(row = row, op = op, col = col,
                                             parTable = parTable,
                                             rows_lhs = rows_lhs, fill = fill)
    }
  }
  mat
}


findInterceptsParTable <- function(mat, parTable, fill = NULL) {
  for (row in rownames(mat)) {
    if (is.na(mat[row, ]))
      mat[row, ] <- extractFromParTable(row = row, op = "~1", col = "",
                                        parTable = parTable, rows_lhs = TRUE,
                                        fill = fill)
  }
  mat
}


findInteractionEstimatesParTable <- function(omega, parTable, fill = NULL) {
  rows <- rownames(omega)
  cols <- colnames(omega)

  for (row in rows) for (col in cols) {
    if (!is.na(omega[row, col])) next
    eta <- getEtaRowLabelOmega(row)
    x   <- getXiRowLabelOmega(row)
    xz  <- createDoubleIntTerms(x = x, z = col, sep = ":")
    omega[row, col] <- extractFromParTable(eta, "~", xz, parTable = parTable,
                                           rows_lhs = TRUE, fill = fill)
  }
  omega
}


extractFromParTable <- function(row, op, col, parTable, rows_lhs = TRUE, fill = NULL) {
  if (rows_lhs) {
    out <- parTable[parTable$lhs == row &
                    parTable$op == op &
                    parTable$rhs %in% col, "est"]
  } else {
    out <- parTable[parTable$lhs == col &
                    parTable$op == op &
                    parTable$rhs %in% row, "est"]
  }

  if (length(out) == 0 && op == "~~") {
    out <- parTable[parTable$lhs == col & parTable$op == op &
                    parTable$rhs %in% row, "est"]
  }

  if (length(out) == 0) {
    stopif(is.null(fill), "No match found")
    out <- fill
  }

  stopif(length(out) > 1, "Incorrect length of matches")

  out
}


sortParTable <- function(parTable, lhs, op, rhs) {
  out <- NULL
  for (l in lhs) {
    for (r in rhs) {
      row <- parTable[parTable$lhs == l & parTable$op == op & parTable$rhs == r, ]
      if (NROW(row) == 0) next
      out <- rbind(out, row)
    }
  }
  out$est
}


parameterEstimatesLavSAM <- function(syntax,
                                     data,
                                     estimator       = "ml",
                                     missing         = "listwise",
                                     meanstructure   = TRUE,
                                     orthogonal.x    = FALSE,
                                     orthogonal.y    = FALSE,
                                     auto.fix.first  = TRUE,
                                     auto.fix.single = TRUE,
                                     suppress.warnings.lavaan = TRUE,
                                     group = NULL,
                                     ...) {
  parTable <- modsemify(syntax)
  higherOrderLVs <- getHigherOrderLVs(parTable)
  isHigherOrder  <- length(higherOrderLVs) > 0L
  isNonCentered  <- isNonCenteredParTable(parTable)
  lowerOrderInds <- unlist(getIndsLVs(parTable, lVs = higherOrderLVs,
                                      isOV = FALSE))

  if (suppress.warnings.lavaan) wrapper <- suppressWarnings
  else                          wrapper <- \(x) x # do nothing

  if (!any(grepl(":", parTable$rhs) | grepl(":", parTable$lhs))) {
    fitSEM <- wrapper(lavaan::sem(
      model           = syntax,
      data            = data,
      meanstructure   = meanstructure,
      estimator       = estimator,
      missing         = missing,
      orthogonal.x    = orthogonal.x,
      orthogonal.y    = orthogonal.y,
      auto.fix.first  = auto.fix.first,
      auto.fix.single = auto.fix.single,
      group           = group,
      ...
    ))

    return(list(fit = fitSEM,
                parTable = lavaan::parameterEstimates(fitSEM)))
  }

  # Get SAM structural model with measurement model from a CFA
  lVs <- getLVs(parTable)

  getCFARows <- function(pt) {
    rhs <- pt$rhs
    lhs <- pt$lhs
    op  <- pt$op

    cond1 <- op == "=~"
    cond2 <- op == "~1"
    cond3 <- op == "~~" & !lhs %in% lVs & !rhs %in% lVs

    # residual variances for higher order lvs are not returned from SAM
    # estimates
    cond4 <- op == "~~" & lhs %in% lowerOrderInds & rhs %in% lowerOrderInds

    pt[cond1 | cond2 | cond3 | cond4, , drop = FALSE]
  }

  parTableOuter <- getCFARows(parTable)

  syntaxCFA <- parTableToSyntax(parTableOuter)

  fitCFA <- wrapper(lavaan::cfa(
    model           = syntaxCFA,
    data            = data,
    meanstructure   = meanstructure,
    estimator       = estimator,
    missing         = missing,
    orthogonal.x    = orthogonal.x,
    orthogonal.y    = orthogonal.y,
    auto.fix.first  = auto.fix.first,
    auto.fix.single = auto.fix.single,
    group           = group,
    ...
  ))

  if (isHigherOrder || isNonCentered) {
    # use factor scores instead
    # using `sam.method="fsr"` doesn't work for this purpose (yet)
    # so we do it manually instead
    dataSAM <- tryCatch(lavaan::lavPredict(fitCFA, transform = TRUE),
                        error = \(e) lavaan::lavPredict(fitCFA))

    structvars <- unique(c(
      colnames(dataSAM),
      parTable[grepl(":", parTable$lhs), "lhs"],
      parTable[grepl(":", parTable$rhs), "rhs"]
    ))

    parTableInner <- parTable[parTable$lhs %in% structvars &
                              parTable$rhs %in% structvars &
                              parTable$op != "=~", , drop = FALSE]
    syntaxSAM <- parTableToSyntax(parTableInner)
    SAMFUN    <- lavaan::sem

  } else {
    syntaxSAM <- parTableToSyntax(parTable)
    dataSAM   <- data
    SAMFUN    <- lavaan::sam
  }

  fitSAM <- wrapper(SAMFUN(
    model           = syntaxSAM,
    data            = dataSAM,
    se              = "none",
    estimator       = estimator,
    missing         = missing,
    orthogonal.x    = orthogonal.x,
    orthogonal.y    = orthogonal.y,
    auto.fix.first  = auto.fix.first,
    auto.fix.single = auto.fix.single,
    group           = group,
    ...
  ))

  measr  <- getCFARows(lavaan::parameterEstimates(fitCFA))
  struct <- lavaan::parameterEstimates(fitSAM)

  addcol <- \(pt, col, val) if (!col %in% colnames(pt)) {pt[[col]] <- val; pt} else pt
  cols.x <- c("lhs", "op", "rhs")
  cols.y <- c("label", "group", "est")
  cols   <- c(cols.x, cols.y)

  measr  <- addcol(measr, col = "label", val = "")
  measr  <- addcol(measr, col = "group", val = 1L)
  struct <- addcol(struct, col = "label", val = "")
  struct <- addcol(struct, col = "group", val = 1L)

  parTableFull <- rbind(measr[cols], struct[cols])
  parTableFull <- parTableFull[!duplicated(parTableFull[cols.x]), , drop = FALSE]

  # if (!isNonCentered && !isHigherOrder) # if latent mean structure is not included
  parTableFull <- recalcInterceptsY(parTableFull)

  list(fit = fitCFA, parTable = parTableFull)
}
