optimizeStartingParamsDA <- function(model,
                                     args = list(orthogonal.x = FALSE,
                                                 orthogonal.y = FALSE,
                                                 auto.fix.first = TRUE,
                                                 auto.fix.sinlge = TRUE)) {
  etas     <- model$info$etas
  indsEtas <- model$info$allIndsEtas
  xis      <- model$info$xis
  numXis   <- model$info$numXis
  indsXis  <- model$info$allIndsXis
  data     <- model$data.raw
  missing  <- tolower(args$missing)

  syntax <- paste(model$syntax, model$covModel$syntax,
                  model$info$lavOptimizerSyntaxAdditions,
                  sep = "\n")

  if (grepl(":", syntax) && missing %in% c("ml", "fiml", "direct")) {
    # It is not worth the extra compute to use fiml (if fiml is set)
    # for getting good starting estiamtes for the model, for non-linear models.
    # FIML can be very slow if there are a lot of product indicators.
    # However, if the model is linear, the estimates should be the
    # exact same, and product indicators aren't an issue. Thus
    # we will save time.
    missing <- "listwise"
  }

  estPI <- modsem_pi(
    model.syntax    = syntax,
    data            = data,
    method          = "dblcent",
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

  parTable <- parameter_estimates(estPI, colon.pi = TRUE)

  stopif(is.null(parTable), "lavaan failed!")

  # labelled parameters
  thetaLabel <- getLabeledParamsLavaan(parTable, model$constrExprs$fixedParams)

  fillLabelsMatrix <- function(matNumeric, matLabel, symmetric = FALSE) {
    if (all(matLabel == ""))
      return(matNumeric)

    labels <- c(matLabel[matLabel != ""])
    for (label in labels)
      matNumeric[matLabel == label] <- thetaLabel[[label]]

    if (symmetric)
      matNumeric[upper.tri(matNumeric)] <- t(matNumeric)[upper.tri(matNumeric)]

    matNumeric
  }

  # Main Model
  matricesMain      <- model$matrices
  labelMatricesMain <- model$labelMatrices

  LambdaX <- findEstimatesParTable(matricesMain$lambdaX, parTable, op = "=~",
                                   rows_lhs = FALSE, fill = 0.7)
  LambdaY <- findEstimatesParTable(matricesMain$lambdaY, parTable, op = "=~",
                                   rows_lhs = FALSE, fill = 0.7)

  ThetaEpsilon <- findEstimatesParTable(matricesMain$thetaEpsilon, parTable,
                                        op = "~~", fill = 0.2)
  ThetaDelta   <- findEstimatesParTable(matricesMain$thetaDelta, parTable,
                                        op = "~~", fill = 0.2)

  Psi <- findEstimatesParTable(matricesMain$psi, parTable, op = "~~", fill = 0)
  Phi <- findEstimatesParTable(matricesMain$phi, parTable, op = "~~", fill = 0)
  A   <- findEstimatesParTable(matricesMain$A, parTable, op = "~~", fill = 0)

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

  beta0 <- findInterceptsParTable(matricesMain$beta0, parTable, fill = 0)
  alpha <- findInterceptsParTable(matricesMain$alpha, parTable, fill = 0)

  GammaEta <- findEstimatesParTable(matricesMain$gammaEta, parTable, op = "~", fill = 0)
  GammaXi  <- findEstimatesParTable(matricesMain$gammaXi, parTable, op = "~", fill = 0)

  OmegaEtaXi <- findInteractionEstimatesParTable(matricesMain$omegaEtaXi,
                                                 parTable = parTable, fill = 0)
  OmegaXiXi <- findInteractionEstimatesParTable(matricesMain$omegaXiXi,
                                                parTable = parTable, fill = 0)
  tauX <- findInterceptsParTable(matricesMain$tauX, parTable, fill = 0)
  tauY <- findInterceptsParTable(matricesMain$tauY, parTable, fill = 0)

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
    PsiCovModel <- findEstimatesParTable(matricesCov$psi, parTable, op = "~~", fill = 0)
    PhiCovModel <- findEstimatesParTable(matricesCov$phi, parTable, op = "~~", fill = 0)

    GammaEtaCovModel <- findEstimatesParTable(matricesCov$gammaEta, parTable, op = "~", fill = 0)
    GammaXiCovModel <- findEstimatesParTable(matricesCov$gammaXi, parTable, op = "~", fill = 0)

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

  # Combinging the two
  theta <- c(thetaLabel, thetaCov, thetaMain)
  if (length(theta) == length(model$theta)) {
    names(theta) <- names(model$theta)
    model$theta <- theta
  }

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


parameterEstimatesLavSAM <- function(syntax, data, ...) {
  parTable <- modsemify(syntax)

  if (!any(grepl(":", parTable$rhs) | grepl(":", parTable$lhs))) {
    fitSEM <- lavaan::sem(syntax, data = data, meanstructure = TRUE, ...)
    return(lavaan::parameterEstimates(fitSEM))
  }

  # Get SAM structural model with measurement model from a CFA
  lVs <- getLVs(parTable)

  getCFARows <- function(pt) {
    pt[pt$op == "=~" |
       pt$op == "~1" |
      (pt$op == "~~" & !pt$lhs %in% lVs & !pt$rhs %in% lVs), ]
  }

  parTableOuter <- getCFARows(parTable)

  syntaxCFA <- parTableToSyntax(parTableOuter)
  syntaxSAM <- parTableToSyntax(parTable)

  fitCFA <- suppressWarnings(lavaan::cfa(syntaxCFA, data = data, meanstructure = TRUE))
  fitSAM <- suppressWarnings(lavaan::sam(syntaxSAM, data = data, se = "none"))

  measr  <- getCFARows(lavaan::parameterEstimates(fitCFA))
  struct <- suppressWarnings(centered_estimates(fitSAM))

  addlab <- \(pt) if (!"label" %in% colnames(pt)) {pt$label <- ""; pt} else pt
  cols <- c("lhs", "op", "rhs", "label", "est")
  rbind(addlab(measr)[cols], addlab(struct)[cols])
}
