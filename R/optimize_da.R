optimizeStartingParamsDA <- function(model) {
  etas <- model$info$etas
  indsEtas <- model$info$allIndsEtas
  xis <- model$info$xis
  numXis <- model$info$numXis
  indsXis <- model$info$allIndsXis
  data <- model$data

  syntax <- paste(model$syntax, model$covModel$syntax, sep = "\n")
  parTable <- modsem(syntax, data, method = "dblcent")$coefParTable

  # Main Model
  matricesMain <- model$matrices

  LambdaX <- findEstimatesParTable(matricesMain$lambdaX, parTable, op = "=~", 
                                   rows_lhs = FALSE)
  LambdaY <- findEstimatesParTable(matricesMain$lambdaY, parTable, op = "=~",
                                   rows_lhs = FALSE)

  ThetaEpsilon <- findEstimatesParTable(matricesMain$thetaEpsilon, parTable, op = "~~")
  ThetaDelta <- findEstimatesParTable(matricesMain$thetaDelta, parTable, op = "~~")

  Psi <- findEstimatesParTable(matricesMain$psi, parTable, op = "~~")
  Phi <- findEstimatesParTable(matricesMain$phi, parTable, op = "~~")

  A <- findEstimatesParTable(matricesMain$A, parTable, op = "~~")
  A[upper.tri(A)] <- t(A[lower.tri(A)])
  A <- t(tryCatch(chol(A), error = function(x) diag(ncol(A))))

  alpha <- matricesMain$alpha
  alpha[is.na(alpha)] <- 0

  GammaEta <- findEstimatesParTable(matricesMain$gammaEta, parTable, op = "~")
  GammaXi <- findEstimatesParTable(matricesMain$gammaXi, parTable, op = "~")  

  OmegaEtaXi <- findInteractionEstimatesParTable(matricesMain$omegaEtaXi, 
                                                 etas = etas, 
                                                 parTable = parTable)
  OmegaXiXi <- findInteractionEstimatesParTable(matricesMain$omegaXiXi, 
                                                etas = etas, 
                                                parTable = parTable)

  tauX <- apply(data[, indsXis, drop = FALSE], MARGIN = 2, FUN = mean)
  tauY <- apply(data[, indsEtas, drop = FALSE], MARGIN = 2, FUN = mean)

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


findEstimatesParTable <- function(mat, parTable, op = NULL, rows_lhs = TRUE) {
  if (is.null(op)) stop2("Missing operator")
  for (row in rownames(mat)) {
    for (col in colnames(mat)) {
      if (is.na(mat[row, col])) 
        mat[row, col] <- extractFromParTable(row, op, col, parTable, rows_lhs)    
    }
  }
  mat 
}


findInteractionEstimatesParTable <- function(omega, etas, parTable) {
  nrowSubOmega <- nrow(omega) %/% length(etas)
  firstRow <- 1
  lastRow <- nrowSubOmega
  for (eta in etas) {
    subOmega <- omega[firstRow:lastRow, , drop = FALSE] 
    rows <- rownames(subOmega)
    cols <- colnames(subOmega)

    for (row in seq_len(NROW(subOmega))) {
      for (col in seq_len(NCOL(subOmega))) {
        if (!is.na(subOmega[row, col])) next
        xz <- createDoubleIntTerms(x = rows[[row]], z = cols[[col]], sep = "")
        subOmega[row, col] <- extractFromParTable(eta, "~", xz, parTable, 
                                                  rows_lhs = TRUE)
      }
    }
    omega[firstRow:lastRow, ] <- subOmega
    firstRow <- lastRow + 1
    lastRow <- lastRow + nrowSubOmega
  }
  omega
}


extractFromParTable <- function(row, op, col, parTable, rows_lhs = TRUE) {
  if (rows_lhs) {
    out <- (parTable[parTable$lhs == row & 
                     parTable$op == op & 
                     parTable$rhs %in% col, "est"])
  } else {
    out <- parTable[parTable$lhs == col & 
                    parTable$op == op & 
                    parTable$rhs %in% row, "est"]
  }

  if (length(out) == 0 && op == "~~") {
    out <- parTable[parTable$lhs == col & parTable$op == op & parTable$rhs %in% row, "est"]
  } else if (length(out) != 1) stop2("Incorrect length of matches")

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
