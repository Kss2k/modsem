optimizeStartingParamsLms <- function(model) {
  etas <- model$info$etas
  indsEtas <- model$info$allIndsEtas
  xis <- model$info$xis
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
  A <- findEstimatesParTable(matrices$A, pt, op = "~~")
  A[upper.tri(A)] <- t(A[lower.tri(A)])
  A <- t(tryCatch(chol(A), error = function(x) diag(ncol(A))))
  alpha <- matrices$alpha
  alpha[is.na(alpha)] <- 0
  GammaEta <- findEstimatesParTable(matrices$gammaEta, pt, op = "~")
  GammaXi <- findEstimatesParTable(matrices$gammaXi, pt, op = "~")  
  OmegaEtaXi <- findInteractionEstimatesParTable(matrices$omegaEtaXi, etas,
                                                 c(etas, xis), pt)
  OmegaXiXi <- findInteractionEstimatesParTable(matrices$omegaXiXi, 
                                                 etas, xis, pt)
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


createParamVector <- function(model, start = NULL) {
  set.seed(123)
  matrices <- model$matrices
  lambdaX <- as.vector(matrices$lambdaX)
  lambdaY <- as.vector(matrices$lambdaY)
  thetaDelta <- as.vector(matrices$thetaDelta)
  thetaEpsilon <- as.vector(matrices$thetaEpsilon)
  phi <- as.vector(matrices$phi)
  A <- as.vector(matrices$A)
  psi <- as.vector(matrices$psi)
  tauX <- as.vector(matrices$tauX)
  tauY <- as.vector(matrices$tauY)
  alpha <- as.vector(matrices$alpha)
  gammaXi <- as.vector(matrices$gammaXi)
  gammaEta <- as.vector(matrices$gammaEta)
  omgeaXiXi <- as.vector(matrices$omegaXiXi)
  omegaEtaXi <- as.vector(matrices$omegaEtaXi)
  theta <- c("lambdaX" = lambdaX,
             "lambdaY" = lambdaY,
             "tauX" = tauX,
             "tauY" = tauY,
             "thetaDelta" = thetaDelta,
             "thetaEpsilon" = thetaEpsilon,
             "phi" = phi,
             "A" = A,
             "psi" = psi,
             "alpha" = alpha,
             "gammaXi" = gammaXi,
             "gammaEta" = gammaEta,
             "omegaXiXi" = omgeaXiXi,
             "omegaEtaXi" = omegaEtaXi)
  theta <- theta[is.na(theta)]
  if (is.null(start)) {
   theta <- vapply(theta,
      FUN.VALUE = vector("numeric", 1L),
      FUN = function(x) runif(1)
    )
  }
  theta
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


findInteractionEstimatesParTable <- function(omega, lhs, rhs, pt) {
  rhs <- expand.grid(rhs, rhs) |> apply(1, stringr::str_c, collapse = "") |>
    unlist()
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
