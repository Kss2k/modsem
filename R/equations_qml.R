logLikQml <- function(theta, model) {
  modelFilled <- fillModel(model, theta, method = "qml")
  numEta <- model$info$numEta
  numXi <- model$info$numXi
  kOmegaEta <- model$info$kOmegaEta

  m <- modelFilled$matrices
  m$numEta <- numEta
  m$numXi <- numXi
  m$kOmegaEta <- kOmegaEta

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

  m$RER <- m$R %*% m$thetaEpsilon %*% t(m$R)
  invRER <- solve(m$RER)
  m$LXPLX <- m$lambdaX %*% m$phi %*% t(m$lambdaX) + m$thetaDelta
  invLXPLX <- solve(m$LXPLX)

  m$L1 <- m$phi %*% t(m$lambdaX) %*% invLXPLX
  m$L2 <- - m$subThetaEpsilon %*% t(m$Beta) %*% invRER

  m$kronXi <- calcKronXi(m, t)
  m$Binv <- calcBinvCpp(m, t)

  m$Sigma1 <- m$phi - m$phi %*% t(m$lambdaX) %*% 
    invLXPLX %*% m$lambdaX %*% m$phi

  m$Sigma2ThetaEpsilon <-  # m$Binv %*% m$psi %*% t(m$Binv) + 
    m$subThetaEpsilon -
    m$subThetaEpsilon ^ 2 %*% t(m$Beta) %*%  
    invRER %*% m$Beta

  Ey <- muQmlCpp(m, t)
  sigmaEpsilon <- sigmaQmlCpp(m, t)

  sigmaXU <- rbind(cbind(m$LXPLX, matrix(0, ncol = ncol(m$RER),
                                      nrow = nrow(m$LXPLX))), 
                   cbind(matrix(0, ncol = ncol(m$LXPLX),
                                nrow = nrow(m$RER)), m$RER))
  mean <- rep(0, ncol(sigmaXU))

  f2 <- dMvn(cbind(m$x, m$u), mean = mean, sigma = sigmaXU, log = TRUE)
  if (numEta <= 1) {
    f3 <- dnormCpp(m$y[,1], mu = Ey, sigma = sqrt(sigmaEpsilon))
  } else {
    f3 <- rep_dmvnorm(m$y[, colnames(m$subThetaEpsilon)], expected = Ey, 
                      sigma = sigmaEpsilon, t = t, cores = 2)
  }
  -sum(f2 + f3)
}


mstepQml <- function(model, theta, hessian = TRUE,
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

  if (hessian){
    if (verbose) cat("Calculating Hessian\n") 
    est$hessian <- nlme::fdHess(pars=est$par, fun = logLikQml,
                                model = model, 
                                .relStep = .Machine$double.eps^(1/5))$Hessian
  }
  est
}
