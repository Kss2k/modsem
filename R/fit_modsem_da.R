# not finished yet
calcChiSqr_da <- function(model) {
  parTable <- model$parTable
  if (any(grepl(":", parTable$rhs))) {
    warning2("Fit measures for LMS and QML should be calculated for baseline model ", 
             "i.e., the model without the interaction effect")
    parTable <- centerInteraction(parTable)
  }
  
  O <- stats::cov(model$data)
  N <- NROW(model$data) 
  p <- NCOL(model$data)

  matrices <- model$model$matrices
  gammaXi <- matrices$gammaXi
  gammaEta <- matrices$gammaEta
  phi <- matrices$phi
  psi <- matrices$psi
  lambdaX <- matrices$lambdaX
  lambdaY <- matrices$lambdaY
  thetaY <- matrices$thetaEpsilon
  thetaX <- matrices$thetaDelta
  Ieta <- matrices$Ieta
  Binv <- solve(Ieta - gammaEta)
  
  covX <- lambdaX %*% phi %*% t(lambdaX) + thetaX
  covXY <- lambdaY %*% (Binv %*% gammaXi %*% phi) %*% t(lambdaX) 
  covY <- lambdaY %*% 
    (Binv %*% (gammaXi %*% phi %*% t(gammaXi) + psi) %*% t(Binv)) %*% 
    t(lambdaY) + thetaY

  E <- rbind(cbind(covX, t(covXY)),
              cbind(covXY, covY))

  chi.sq <- N * (log(det(E)) + tr(E %*% solve(O)) - log(det(O)) - p)
  chi.sq
}
