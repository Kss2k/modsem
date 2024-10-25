logLikQml <- function(theta, model, sum = TRUE, sign = -1) {
  modelFilled <- fillModel(model, theta, method = "qml")
  numXi       <- model$info$numXis
  numEta      <- model$info$numEtas
  kOmegaEta   <- model$info$kOmegaEta
  latentEtas  <- model$info$latentEtas

  m <- modelFilled$matrices
  m$numEta    <- numEta
  m$numXi     <- numXi
  m$kOmegaEta <- kOmegaEta

  m$tauX      <- m$tauX + m$lambdaX %*% m$beta0
  m$x <- model$data[, model$info$allIndsXis, drop = FALSE]
  m$y <- model$data[, model$info$allIndsEtas, drop = FALSE]
  m$x <- centerIndicators(m$x, tau = m$tauX)
  m$y <- centerIndicators(m$y, tau = m$tauY)

  t <- NROW(m$x)
  if (!is.null(m$emptyR)) {
    m$R <- m$emptyR
    m$R[is.na(m$R)] <- -m$lambdaY[!m$selectScalingY] # fill R with -Beta
    m$fullR[m$rowsR, m$colsR] <- m$R
    m$u <- m$y %*% t(m$fullR)
    m$fullU[, m$colsU] <- m$u
    m$Beta <- m$lambdaY[m$selectBetaRows, latentEtas]
    m$subThetaEpsilon <- m$subThetaEpsilon
    m$subThetaEpsilon[is.na(m$subThetaEpsilon)] <-
      m$thetaEpsilon[m$selectThetaEpsilon]

    m$RER <- m$R %*% m$thetaEpsilon[m$colsR, m$colsR] %*% t(m$R)
    invRER <- solve(m$RER)
    m$L2 <- -m$subThetaEpsilon %*% t(m$Beta) %*% invRER
    m$fullL2[m$selectSubL2] <- m$L2

    m$Sigma2ThetaEpsilon <- m$subThetaEpsilon - m$subThetaEpsilon ^ 2 %*%
      t(m$Beta) %*% invRER %*% m$Beta

    m$fullSigma2ThetaEpsilon[m$selectSubSigma2ThetaEpsilon] <-
      m$Sigma2ThetaEpsilon
  }

  m$Sigma2ThetaEpsilon <- m$fullSigma2ThetaEpsilon
  m$L2     <- m$fullL2
  m$u      <- m$fullU
  m$LXPLX  <- m$lambdaX %*% m$phi %*% t(m$lambdaX) + m$thetaDelta
  invLXPLX <- solve(m$LXPLX)
  m$L1     <- m$phi %*% t(m$lambdaX) %*% invLXPLX
  m$Sigma1 <- m$phi - m$phi %*% t(m$lambdaX) %*% invLXPLX %*% m$lambdaX %*% m$phi
  m$kronXi <- calcKronXi(m, t)
  m$Binv   <- calcBinvCpp(m, t)

  Ey            <- muQmlCpp(m, t)
  sigmaEpsilon  <- sigmaQmlCpp(m, t)
  sigmaXU       <- calcSigmaXU(m)

  normalInds    <- colnames(sigmaXU)
  indsY         <- colnames(m$y)
  nonNormalInds <- indsY[!indsY %in% normalInds]

  f2 <- probf2(matrices = m, normalInds = normalInds, sigma = sigmaXU)
  f3 <- probf3(matrices = m, nonNormalInds = nonNormalInds, expected = Ey,
               sigma = sigmaEpsilon, t = t, numEta = numEta)

  if (sum) return(sign * sum(f2 + f3))
  sign * (f2 + f3)
}


calcSigmaXU <- function(matrices) {
  if (is.null(matrices$emptyR)) return(matrices$LXPLX)
  diagBindSquareMatrices(matrices$LXPLX, matrices$RER)
}


probf2 <- function(matrices, normalInds, sigma) {
  mu <- rep(0, ncol(sigma))
  X  <- cbind(matrices$x, matrices$u)[ , normalInds]
  dmvn(X, mean = mu, sigma = sigma, log = TRUE)
}


probf3 <- function(matrices, nonNormalInds, expected, sigma, t, numEta) {
  if (numEta == 1) {
    return(dnormCpp(matrices$y[, 1], mu = expected, sigma = sqrt(sigma)))
  }
  rep_dmvnorm(matrices$y[, nonNormalInds], expected = expected,
              sigma = sigma, t = t)
}


centerIndicators <- function(X, tau) {
  for (i in seq_len(ncol(X))) X[, i] <- X[, i] - tau[[i]]
  X
}


gradientLogLikQml <- function(theta, model, epsilon = 1e-8, sign = -1) {
  baseLL <- logLikQml(theta, model, sign = sign)
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (logLikQml(theta, model, sign = sign) - baseLL) / epsilon
  })
}


# log likelihood for each observation -- not all
# wrapper fro logLikQml(sum = FALSE)
logLikQml_i <- function(theta, model, sign = -1) {
  logLikQml(theta, model, sum = FALSE, sign = sign)
}


# gradient function of logLikQml_i
gradientLogLikQml_i <- function(theta, model, sign = -1, epsilon = 1e-8) {
  baseLL <- logLikQml_i(theta, model, sign = sign)
  lapplyMatrix(seq_along(theta), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (logLikQml_i(theta, model, sign = sign) - baseLL) / epsilon
  }, FUN.VALUE = numeric(nrow(model$data)))
}


mstepQml <- function(model,
                     theta,
                     max.iter = 500,
                     verbose = FALSE,
                     convergence = 1e-6,
                     control = list(),
                     optimizer = "nlminb",
                     epsilon = 1e-8,
                     ...) {
  gradient <- function(theta, model, sign)
    gradientLogLikQml(theta = theta, model = model, epsilon = epsilon, sign = sign)

  if (verbose) cat("Starting M-step\n")

  if (optimizer == "nlminb") {
    control$iter.max <- max.iter
    control$eval.max <- max.iter * 2
    control$rel.tol  <- convergence

    est <- stats::nlminb(start = theta, objective = logLikQml, model = model,
                         gradient = gradient, sign = -1,
                         upper = model$info$bounds$upper,
                         lower = model$info$bounds$lower, control = control, ...)

  } else if (optimizer == "L-BFGS-B") {
    control$factr <- convergence
    control$maxit <- max.iter

    est <- stats::optim(par = theta, fn = logLikQml, model = model,
                        gr = gradient, method = optimizer, sign = -1,
                        control = control, ...)

    est$objective  <- est$value
    est$iterations <- est$counts[["function"]]
  } else stop2("Unrecognized optimizer, must be either 'nlminb' or 'L-BFGS-B'")

  est
}
