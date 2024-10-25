optimizeStartingParamsDA <- function(model) {
  etas     <- model$info$etas
  indsEtas <- model$info$allIndsEtas
  xis      <- model$info$xis
  numXis   <- model$info$numXis
  indsXis  <- model$info$allIndsXis
  data     <- model$data

  syntax <- paste(model$syntax, model$covModel$syntax,
                  model$info$lavOptimizerSyntaxAdditions,
                  sep = "\n")
  parTable <- modsem_pi(syntax, data, method = "dblcent",
                        meanstructure = TRUE,
                        suppress.warnings.lavaan = TRUE)$coefParTable
  if (is.null(parTable)) {
    modsem_pi(syntax, data, method = "dblcent", meanstructure = TRUE)
    stop2("lavaan failed")
  }
  # Main Model
  matricesMain <- model$matrices
  LambdaX <- findEstimatesParTable(matricesMain$lambdaX, parTable, op = "=~",
                                   rows_lhs = FALSE)
  LambdaY <- findEstimatesParTable(matricesMain$lambdaY, parTable, op = "=~",
                                   rows_lhs = FALSE)

  ThetaEpsilon <- findEstimatesParTable(matricesMain$thetaEpsilon, parTable, op = "~~")
  ThetaDelta   <- findEstimatesParTable(matricesMain$thetaDelta, parTable, op = "~~")

  Psi <- findEstimatesParTable(matricesMain$psi, parTable, op = "~~")
  Phi <- findEstimatesParTable(matricesMain$phi, parTable, op = "~~")

  A <- findEstimatesParTable(matricesMain$A, parTable, op = "~~")
  A[upper.tri(A)] <- t(A[lower.tri(A)])
  A <- t(tryCatch(chol(A), error = function(x) diag(ncol(A))))

  beta0 <- findInterceptsParTable(matricesMain$beta0, parTable, fill = 0)
  alpha <- findInterceptsParTable(matricesMain$alpha, parTable, fill = 0)

  GammaEta <- findEstimatesParTable(matricesMain$gammaEta, parTable, op = "~")
  GammaXi  <- findEstimatesParTable(matricesMain$gammaXi, parTable, op = "~")

  OmegaEtaXi <- findInteractionEstimatesParTable(matricesMain$omegaEtaXi,
                                                 parTable = parTable)
  OmegaXiXi <- findInteractionEstimatesParTable(matricesMain$omegaXiXi,
                                                parTable = parTable)
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
  matricesCov <- model$covModel$matrices
  if (!is.null(matricesCov)) {
    PsiCovModel <- findEstimatesParTable(matricesCov$psi, parTable, op = "~~")
    PhiCovModel <- findEstimatesParTable(matricesCov$phi, parTable, op = "~~")
    ACovModel <- findEstimatesParTable(matricesCov$A, parTable, op = "~~")
    ACovModel[upper.tri(ACovModel)] <- t(ACovModel[lower.tri(ACovModel)])
    ACovModel <- t(tryCatch(chol(ACovModel), error = function(x)
                            diag(ncol(ACovModel))))
    GammaEtaCovModel <- findEstimatesParTable(matricesCov$gammaEta, parTable, op = "~")
    GammaXiCovModel <- findEstimatesParTable(matricesCov$gammaXi, parTable, op = "~")

    thetaCov <- unlist(list(PhiCovModel[is.na(matricesCov$phi)],
                            ACovModel[is.na(matricesCov$A)],
                            PsiCovModel[is.na(matricesCov$psi)],
                            GammaXiCovModel[is.na(matricesCov$gammaXi)],
                            GammaEtaCovModel[is.na(matricesCov$gammaEta)]))
  } else thetaCov <- NULL

  # labelTheta
  thetaLabel <- getLabeledParamsLavaan(parTable, model$constrExprs$fixedParams)

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


findInteractionEstimatesParTable <- function(omega, parTable) {
  rows <- rownames(omega)
  cols <- colnames(omega)

  for (row in rows) for (col in cols) {
    if (!is.na(omega[row, col])) next
    eta <- getEtaRowLabelOmega(row)
    x   <- getXiRowLabelOmega(row)
    xz  <- createDoubleIntTerms(x = x, z = col, sep = "")
    omega[row, col] <- extractFromParTable(eta, "~", xz, parTable,
                                           rows_lhs = TRUE)
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
    if (is.null(fill)) stop("No match found")
    out <- fill
  }

  if (length(out) > 1) stop("Incorrect length of matches")
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
