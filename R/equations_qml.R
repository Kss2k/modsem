logLikQml <- function(theta, model, sum = TRUE, sign = -1) {
  modelFilled <- fillModel(model, theta, method = "qml")
  numEta <- model$info$numEta
  numXi <- model$info$numXi
  kOmegaEta <- model$info$kOmegaEta

  m <- modelFilled$matrices
  m$numEta <- numEta
  m$numXi <- numXi
  m$kOmegaEta <- kOmegaEta

  m$x <- model$data[, model$info$allIndsXis, drop = FALSE]
  for (i in seq_len(ncol(m$x))) {
    m$x[, i] <- m$x[, i] - m$tauX[i]
  }

  m$y <- model$data[, model$info$allIndsEtas, drop = FALSE]
  for (i in seq_len(ncol(m$y))) m$y[, i] <- m$y[, i] - m$tauY[i]

  t <- NROW(m$x)
  if (ncol(m$y) > 1) {
    m$R <- m$emptyR
    m$R[is.na(m$R)] <- -m$lambdaY[!m$selectScalingY] # fill R with -Beta
    m$u <- m$y %*% t(m$R)
    m$Beta <- m$lambdaY[m$selectBetaRows, ]
  } else {
    m$u <- 0
  }


  m$subThetaEpsilon <- m$subThetaEpsilon
  m$subThetaEpsilon[is.na(m$subThetaEpsilon)] <-
    m$thetaEpsilon[m$selectThetaEpsilon]

  m$RER <- m$R %*% m$thetaEpsilon %*% t(m$R)
  invRER <- solve(m$RER)
  m$LXPLX <- m$lambdaX %*% m$phi %*% t(m$lambdaX) + m$thetaDelta
  invLXPLX <- solve(m$LXPLX)

  m$L1 <- m$phi %*% t(m$lambdaX) %*% invLXPLX
  m$L2 <- -m$subThetaEpsilon %*% t(m$Beta) %*% invRER

  m$kronXi <- calcKronXi(m, t)
  m$Binv <- calcBinvCpp(m, t)

  m$Sigma1 <- m$phi - m$phi %*% t(m$lambdaX) %*%
    invLXPLX %*% m$lambdaX %*% m$phi

  m$Sigma2ThetaEpsilon <- # m$Binv %*% m$psi %*% t(m$Binv) +
    m$subThetaEpsilon -
    m$subThetaEpsilon^2 %*% t(m$Beta) %*%
    invRER %*% m$Beta

  Ey <- muQmlCpp(m, t)
  sigmaEpsilon <- sigmaQmlCpp(m, t)

  sigmaXU <- rbind(
    cbind(m$LXPLX, matrix(0,
      ncol = ncol(m$RER),
      nrow = nrow(m$LXPLX)
    )),
    cbind(matrix(0,
      ncol = ncol(m$LXPLX),
      nrow = nrow(m$RER)
    ), m$RER)
  )
  mean <- rep(0, ncol(sigmaXU))

  f2 <- dmvn(cbind(m$x, m$u), mean = mean, sigma = sigmaXU, log = TRUE)
  if (numEta <= 1) {
    f3 <- dnormCpp(m$y[, 1], mu = Ey, sigma = sqrt(sigmaEpsilon))
  } else {
    f3 <- rep_dmvnorm(m$y[, colnames(m$subThetaEpsilon)],
      expected = Ey,
      sigma = sigmaEpsilon, t = t, cores = 2
    )
  }
  if (sum) {
    return(sign * sum(f2 + f3))
  }
  sign * (f2 + f3)
}


gradientLogLikQml <- function(theta, model, dt = 1e-8) {
  baseLL <- logLikQml(theta, model)
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + dt
    (logLikQml(theta, model) - baseLL) / dt
  })
}


# log likelihood for each observation -- not all
# wrapper fro logLikQml(sum = FALSE)
logLikQml_i <- function(theta, model, sign = -1) {
  logLikQml(theta, model, sum = FALSE, sign = sign)
}


# gradient function of logLikQml_i
gradientLogLikQml_i <- function(theta, model, sign = -1, dt = 1e-8) {
  baseLL <- logLikQml_i(theta, model, sign = sign)
  lapplyMatrix(seq_along(theta), FUN = function(i) {
    theta[[i]] <- theta[[i]] + dt
    (logLikQml_i(theta, model, sign = sign) - baseLL) / dt
  }, FUN.VALUE = numeric(nrow(model$data)))
}


mstepQml <- function(model, theta,
                     max.iter = 150, verbose = FALSE,
                     convergence = 1e-2,
                     control = list(),
                     optimizer = "nlminb", ...) {
  if (verbose) cat("Starting M-step\n")
  if (optimizer == "nlminb") {
    control$iter.max <- max.iter
    control$eval.max <- max.iter * 2
    control$rel.tol <- convergence
    est <- stats::nlminb(
      start = theta, objective = logLikQml,
      model = model,
      gradient = gradientLogLikQml,
      upper = model$info$bounds$upper,
      lower = model$info$bounds$lower,
      control = control, ...
    )

  } else if (optimizer == "optim") {
    control$maxit <- max.iter
    control$reltol <- convergence
    est <- stats::optim(
      par = theta, fn = logLikQml,
      model = model,
      gr = gradientLogLikQml,
      method = "BFGS",
      control = control, ...
    )
  } else stop2("Unknown optimizer")
  
  est
}
