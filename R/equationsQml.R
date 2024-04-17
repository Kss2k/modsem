logLikQml <- function(theta, model) {
  modelFilled <- fillModel(model, theta)
  numEta <- model$info$numEta
  m <- modelFilled$matrices

  m$x <- model$data[, model$info$allIndsXis]
  for (i in seq_len(ncol(m$x))) {
    m$x[, i] <- m$x[, i] - m$tauX[i]
  }

  m$y <- model$data[, model$info$allIndsEtas] 
  for (i in seq_len(ncol(m$y))) m$y[, i] <- m$y[, i] - m$tauY[i]

  t <- NROW(m$x)
  if (ncol(m$y) > 1) {
    m$R <- m$emptyR
    m$R[is.na(m$R)] <- -m$lambdaY[!m$selectScalingY] # fill R with -Beta
    m$u <- m$y %*% t(m$R)
    m$Beta <- m$lambdaY[m$selectBetaRows, ]
  } else m$u <- 0
  
  m$subThetaEpsilon <- m$subThetaEpsilon
  m$subThetaEpsilon[is.na(m$subThetaEpsilon)] <- 
    m$thetaEpsilon[m$selectThetaEpsilon]

  m$subPhi <- m$phi[seq_len(model$info$numXis), seq_len(model$info$numXis)]
  m$RER <- m$R %*% m$thetaEpsilon %*% t(m$R)
  invRER <- solve(m$RER)
  m$LXPLX <- m$lambdaX %*% m$subPhi %*% t(m$lambdaX) + m$thetaDelta
  invLXPLX <- solve(m$LXPLX)

  m$L1 <- m$subPhi %*% t(m$lambdaX) %*% invLXPLX
  #m$L1 <- diagPartitioned(m$subL1, numEta) 

  m$L2 <- - m$subThetaEpsilon %*% t(m$Beta) %*% invRER
  #m$L2 <- diagPartitioned(m$subL2, numEta)

  m$Sigma1 <- m$subPhi - m$subPhi %*% t(m$lambdaX) %*% 
    invLXPLX %*% m$lambdaX %*% m$subPhi
  #m$Sigma1 <- diagPartitioned(m$subSigma1, numEta)

  m$Sigma2 <- m$psi + m$subThetaEpsilon -
    m$subThetaEpsilon ^ 2 %*% t(m$Beta) %*%  
    invRER %*% m$Beta
  #m$Sigma2 <- diagPartitioned(m$subSigma2, numEta)
  Ey <- muQmlCpp(m, t)
  sigmaEpsilon <- sigmaQmlCpp(m, t)

  sigmaXU <- rbind(cbind(m$LXPLX, matrix(0, ncol = ncol(m$RER),
                                      nrow = nrow(m$LXPLX))), 
                   cbind(matrix(0, ncol = ncol(m$LXPLX),
                                nrow = nrow(m$RER)), m$RER))
  mean <- rep(0, ncol(sigmaXU))
  f2 <- dMvn(cbind(m$x, m$u), mean = mean, sigma = sigmaXU, log = TRUE)
  # original implementation: produces NaN when sds are negative
  f3 <- dnormCpp(m$y[,1], mu = Ey, sigma = sqrt(sigmaEpsilon))
  -sum(f2 + f3)
}


mstepQml <- function(model, theta, negHessian = TRUE,
                     maxIter = 150, verbose = FALSE,
                     convergence = 1e-2,
                     control = list(), ...) {
  control$iter.max <- maxIter
  control$eval.max <- maxIter * 2
  control$rel.tol <- convergence
  if (verbose) cat("Starting M-step\n")
  est <- stats::nlminb(start = theta, objective = logLikQml,
                model = model,
                upper = model$info$bounds$upper,
                lower = model$info$bounds$lower, control = control, ...)

  if (negHessian){
    if (verbose) cat("Calculating Hessian\n") 
    est$hessian <- nlme::fdHess(pars=est$par, fun = logLikQml,
                                model = model, 
                                .relStep = .Machine$double.eps^(1/5))$Hessian
  }
  est
}
