optimizeStartingParamsLms <- function(model) {
  etas <- model$info$etas
  indsEtas <- model$info$allIndsEtas
  xis <- model$info$xis
  numXis <- model$info$numXis
  indsXis <- model$info$allIndsXis
  data <- model$data
  syntax <- paste(model$syntax, model$covModel$syntax, sep = "\n")
  pt <- modsem(syntax, data, method = "dblcent")$coefParTable

  # Main Model
  matricesMain <- model$matrices
  LambdaX <- findEstimatesParTable(matricesMain$lambdaX, pt, op = "=~", 
                                   rows_lhs = FALSE)
  LambdaY <- findEstimatesParTable(matricesMain$lambdaY, pt, op = "=~",
                                   rows_lhs = FALSE)
  ThetaEpsilon <- findEstimatesParTable(matricesMain$thetaEpsilon, pt, op = "~~")
  ThetaDelta <- findEstimatesParTable(matricesMain$thetaDelta, pt, op = "~~")
  Psi <- findEstimatesParTable(matricesMain$psi, pt, op = "~~")
  Phi <- findEstimatesParTable(matricesMain$phi, pt, op = "~~")
  A <- findEstimatesParTable(matricesMain$A, pt, op = "~~")
  A[upper.tri(A)] <- t(A[lower.tri(A)])
  A <- t(tryCatch(chol(A), error = function(x) diag(ncol(A))))
  alpha <- matricesMain$alpha
  alpha[is.na(alpha)] <- 0
  GammaEta <- findEstimatesParTable(matricesMain$gammaEta, pt, op = "~")
  GammaXi <- findEstimatesParTable(matricesMain$gammaXi, pt, op = "~")  
  OmegaEtaXi <- findInteractionEstimatesParTable(matricesMain$omegaEtaXi, lhs = etas,
                                                 rhs1 = etas, rhs2 = xis, pt = pt)
  OmegaXiXi <- findInteractionEstimatesParTable(matricesMain$omegaXiXi, lhs = etas,
                                                rhs1 = xis, rhs2 = xis, pt = pt)
  tauX <- apply(data[, indsXis], 2, mean)
  tauY <- apply(data[, indsEtas], 2, mean)
  thetaMainModel <- unlist(list(LambdaX[is.na(matricesMain$lambdaX)], 
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
    PsiCovModel <- findEstimatesParTable(matricesCov$psi, pt, op = "~~")
    PhiCovModel <- findEstimatesParTable(matricesCov$phi, pt, op = "~~")
    ACovModel <- findEstimatesParTable(matricesCov$A, pt, op = "~~")
    ACovModel[upper.tri(ACovModel)] <- t(ACovModel[lower.tri(ACovModel)])
    ACovModel <- t(tryCatch(chol(ACovModel), error = function(x) 
                            diag(ncol(ACovModel))))
    GammaEtaCovModel <- findEstimatesParTable(matricesCov$gammaEta, pt, op = "~")
    GammaXiCovModel <- findEstimatesParTable(matricesCov$gammaXi, pt, op = "~")  

    thetaCovModel <- unlist(list(PhiCovModel[is.na(matricesCov$phi)], 
                                 ACovModel[is.na(matricesCov$A)], 
                                 PsiCovModel[is.na(matricesCov$psi)], 
                                 GammaXiCovModel[is.na(matricesCov$gammaXi)], 
                                 GammaEtaCovModel[is.na(matricesCov$gammaEta)]))
  } else thetaCovModel <- NULL

  # Combinging the two
  theta <- c(thetaCovModel, thetaMainModel)
  if (length(theta) == length(model$theta)) {
    names(theta) <- names(model$theta)
    model$theta <- theta
  }

  model
}


findEstimatesParTable <- function(mat, pt, op = NULL, rows_lhs = TRUE) {
  if (is.null(op)) stop2("Missing operator")
  for (row in rownames(mat)) {
    for (col in colnames(mat)) {
      if (is.na(mat[row, col])) 
        mat[row, col] <- extractFromParTable(row, op, col, pt, rows_lhs)    
    }
  }
  mat 
}


findInteractionEstimatesParTable <- function(omega, lhs, rhs1, rhs2, pt) {
  if (length(rhs1) != length(rhs2) || !all(rhs1 != rhs2)) {
    combos <- rbind(expand.grid(rhs1, rhs2), expand.grid(rhs2, rhs1))
  } else {
    combos <- expand.grid(rhs1, rhs2)
  }
  rhs <- combos |> apply(1, stringr::str_c, collapse = "") |>
    unlist() |> unique()
  omega[is.na(omega)] <- sortParTable(pt, lhs, "~", rhs)  
  omega
}


extractFromParTable <- function(row, op, col, pt, rows_lhs = TRUE) {
  if (rows_lhs) out <- (pt[pt$lhs == row & pt$op == op & pt$rhs == col, "est"])
  else out <- pt[pt$lhs == col & pt$op == op & pt$rhs == row, "est"]
  if (length(out) == 0 && op == "~~") {
    out <- pt[pt$lhs == col & pt$op == op & pt$rhs == row, "est"]
  } else if (length(out) != 1) stop2("Incorrect length of matches")
  out
}


sortParTable <- function(pt, lhs, op, rhs) {
  out <- NULL
  for (l in lhs) {
    for (r in rhs) {
      row <- pt[pt$lhs == l & pt$op == op & pt$rhs == r, ]
      if (NROW(row) == 0) next
      out <- rbind(out, row)
    }
  }
  out$est
}
