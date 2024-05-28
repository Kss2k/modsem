optimizeStartingParamsLms <- function(model) {
  etas <- model$info$etas
  indsEtas <- model$info$allIndsEtas
  xis <- model$info$xis
  numXis <- model$info$numXis
  indsXis <- model$info$allIndsXis
  data <- model$data
  pt <- modsem(model$syntax, data, method = "dblcent")$coefParTable
  matrices <- model$matrices 
  LambdaX <- findEstimatesParTable(matrices$lambdaX, pt, op = "=~", 
                                   rows_lhs = FALSE)
  LambdaY <- findEstimatesParTable(matrices$lambdaY, pt, op = "=~",
                                   rows_lhs = FALSE)
  ThetaEpsilon <- findEstimatesParTable(matrices$thetaEpsilon, pt, op = "~~")
  ThetaDelta <- findEstimatesParTable(matrices$thetaDelta, pt, op = "~~")
  Psi <- findEstimatesParTable(matrices$psi, pt, op = "~~")
  Phi <- findEstimatesParTable(matrices$phi, pt, op = "~~")
  subPhi <- Phi[seq_len(numXis), seq_len(numXis)]
  A <- findEstimatesParTable(matrices$A, pt, op = "~~")
  subA <- A[seq_len(numXis), seq_len(numXis)]
  subA[upper.tri(subA)] <- t(subA[lower.tri(subA)])
  subA <- t(tryCatch(chol(subA), error = function(x) diag(ncol(subA))))
  A[seq_len(numXis), seq_len(numXis)] <- subA
  alpha <- matrices$alpha
  alpha[is.na(alpha)] <- 0
  GammaEta <- findEstimatesParTable(matrices$gammaEta, pt, op = "~")
  GammaXi <- findEstimatesParTable(matrices$gammaXi, pt, op = "~")  
  OmegaEtaXi <- findInteractionEstimatesParTable(matrices$omegaEtaXi, lhs = etas,
                                                 rhs1 = etas, rhs2 = xis, pt = pt)
  OmegaXiXi <- findInteractionEstimatesParTable(matrices$omegaXiXi, lhs = etas,
                                                rhs1 = xis, rhs2 = xis, pt = pt)
  tauX <- apply(data[, indsXis], 2, mean)
  tauY <- apply(data[, indsEtas], 2, mean)
  theta <- unlist(list(LambdaX[is.na(matrices$lambdaX)], 
                       LambdaY[is.na(matrices$lambdaY)], 
                       tauX[is.na(matrices$tauX)], 
                       tauY[is.na(matrices$tauY)],
                       ThetaDelta[is.na(matrices$thetaDelta)],
                       ThetaEpsilon[is.na(matrices$thetaEpsilon)],
                       Phi[is.na(matrices$phi)], 
                       A[is.na(matrices$A)], 
                       Psi[is.na(matrices$psi)], 
                       alpha[is.na(matrices$alpha)],
                       GammaXi[is.na(matrices$gammaXi)], 
                       GammaEta[is.na(matrices$gammaEta)], 
                       OmegaXiXi[is.na(matrices$omegaXiXi)], 
                       OmegaEtaXi[is.na(matrices$omegaEtaXi)]))
  if (length(theta) == length(model$theta)) {
    names(theta) <- names(model$theta)
    model$theta <- theta
  }
  model
}


findEstimatesParTable <- function(mat, pt, op = NULL, rows_lhs = TRUE) {
  if (is.null(op)) stop("Missing operator")
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
  } else if (length(out) != 1) stop("Incorrect length of matches")
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
