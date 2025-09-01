OptimizerInfoQML <- rlang::env(eval = 0, logLiks = 0)


resetOptimizerInfoQML <- function() {
  OptimizerInfoQML$eval    <- 0
  OptimizerInfoQML$logLiks <- -Inf
}


incrementIterations <- function(logLik) {
  eval       <- OptimizerInfoQML$eval + 1
  logLiks    <- c(OptimizerInfoQML$logLiks, logLik)
  diffLL     <- getDiffTwoMax(logLiks)
  deltaLL    <- diffLL$abs
  relDeltaLL <- diffLL$rel

  clearConsoleLine() # clear before printing
  printf("\rEval=%d LogLik=%.2f \u0394LL=%.2g rel\u0394LL=%.2g",
         eval, logLik, abs(deltaLL), abs(relDeltaLL))

  OptimizerInfoQML$eval    <- eval
  OptimizerInfoQML$logLiks <- logLiks
}


logLikQml <- function(theta, model, sum = TRUE, sign = -1, verbose = FALSE) {
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
  m$x <- model$data$data.full[, model$info$allIndsXis, drop = FALSE]
  m$y <- model$data$data.full[, model$info$allIndsEtas, drop = FALSE]
  m$x <- centerIndicators(m$x, tau = m$tauX)
  m$y <- centerIndicators(m$y, tau = m$tauY)

  # Fill in residual variances for latent etas with only a single indicator
  # only needed if the resiudal variances are non-zero
  if (any(m$selectThetaEpsilon2)) {
    m$subThetaEpsilon2[is.na(m$subThetaEpsilon2)] <-
      m$thetaEpsilon[m$selectThetaEpsilon2]

    selectRows <- apply(m$selectSubSigma2ThetaEpsilon, MARGIN = 1, FUN = \(x) all(!x))
    selectCols <- apply(m$selectSubSigma2ThetaEpsilon, MARGIN = 2, FUN = \(x) all(!x))
    m$fullSigma2ThetaEpsilon[selectRows, selectCols] <- m$subThetaEpsilon2
  }

  t <- NROW(m$x)
  if (!is.null(m$emptyR)) {
    m$R <- m$emptyR
    m$R[is.na(m$R)] <- -m$lambdaY[!m$selectScalingY] # fill R with -Beta
    m$fullR[m$rowsR, m$colsR] <- m$R
    m$u <- m$y %*% t(m$fullR)
    m$fullU[, m$colsU] <- m$u
    m$Beta <- m$lambdaY[m$selectBetaRows, latentEtas]

    m$subThetaEpsilon1[is.na(m$subThetaEpsilon1)] <-
      m$thetaEpsilon[m$selectThetaEpsilon1]

    m$RER <- m$R %*% m$thetaEpsilon[m$colsR, m$colsR] %*% t(m$R)
    invRER <- solve(m$RER)
    m$L2 <- -m$subThetaEpsilon1 %*% t(m$Beta) %*% invRER
    m$fullL2[m$selectSubL2] <- m$L2

    m$Sigma2ThetaEpsilon <- m$subThetaEpsilon1 - m$subThetaEpsilon1 ^ 2 %*%
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
  m$kronXi <- calcKronXi(m, t, ncores = ThreadEnv$n.threads)
  m$Binv   <- calcBinvCpp(m, t, ncores = ThreadEnv$n.threads)

  Ey            <- muQmlCpp(m, t, ncores = ThreadEnv$n.threads)
  sigmaEpsilon  <- sigmaQmlCpp(m, t, ncores = ThreadEnv$n.threads)
  sigmaXU       <- calcSigmaXU(m)

  normalInds    <- colnames(sigmaXU)
  indsY         <- colnames(m$y)
  nonNormalInds <- indsY[!indsY %in% normalInds]

  f2 <- probf2(matrices = m, normalInds = normalInds, sigma = sigmaXU, sum = sum)
  f3 <- probf3(matrices = m, nonNormalInds = nonNormalInds, expected = Ey,
               sigma = sigmaEpsilon, t = t, numEta = numEta, sum = sum)

  logLik <- (f2 + f3)

  if (verbose) incrementIterations(logLik)

  sign * logLik
}


calcSigmaXU <- function(matrices) {
  if (is.null(matrices$emptyR)) return(matrices$LXPLX)
  diagBindSquareMatrices(matrices$LXPLX, matrices$RER)
}


probf2 <- function(matrices, normalInds, sigma, sum=FALSE) {
  browser()
  mu <- rep(0, ncol(sigma))
  X  <- cbind(matrices$x, matrices$u)[ , normalInds]

  p <- dmvn(X, mean = mu, sigma = sigma, log = TRUE)

  if (sum) sum(p) else p
}


probf3 <- function(matrices, nonNormalInds, expected, sigma, t, numEta, sum = FALSE) {
  if (numEta == 1)
    p <- dnormCpp(matrices$y[, 1], mu = expected, sigma = sqrt(sigma), ncores = ThreadEnv$n.threads)
  else
    p <- rep_dmvnorm(matrices$y[, nonNormalInds], expected = expected,
                     sigma = sigma, t = t, ncores = ThreadEnv$n.threads)
  if (sum) sum(p) else p
}


centerIndicators <- function(X, tau) {
  for (i in seq_len(ncol(X))) X[, i] <- X[, i] - tau[[i]]
  X
}


gradientLogLikQml <- function(theta, model, epsilon = 1e-8, sign = -1, data=NULL) {
  if (!is.null(data)) model$data <- data
  baseLL <- logLikQml(theta, model, sign = sign, verbose = FALSE)
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (logLikQml(theta, model, sign = sign, verbose = FALSE) - baseLL) / epsilon
  })
}


# log likelihood for each observation -- not all
# wrapper fro logLikQml(sum = FALSE)
logLikQml_i <- function(theta, model, sign = -1) {
  logLikQml(theta, model, sum = FALSE, sign = sign, verbose = FALSE)
}


# gradient function of logLikQml_i
gradientLogLikQml_i <- function(theta, model, sign = -1, epsilon = 1e-8) {
  baseLL <- logLikQml_i(theta, model, sign = sign)
  lapplyMatrix(seq_along(theta), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (logLikQml_i(theta, model, sign = sign) - baseLL) / epsilon
  }, FUN.VALUE = numeric(model$data$n))
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
  resetOptimizerInfoQML()

  gradient <- function(theta, model, sign, ...) # use ... to capture unused args
    gradientLogLikQml(theta = theta, model = model, epsilon = epsilon, sign = sign)

  if (optimizer == "nlminb") {
    control$iter.max <- max.iter
    control$eval.max <- max.iter * 2
    control$rel.tol  <- convergence

    est <- stats::nlminb(start = theta, objective = logLikQml, model = model,
                         gradient = gradient, sign = -1, verbose = verbose,
                         upper = model$info$bounds$upper,
                         lower = model$info$bounds$lower,
                         control = control, ...) |> suppressWarnings()

  } else if (optimizer == "L-BFGS-B") {
    control$factr <- convergence
    control$maxit <- max.iter

    est <- stats::optim(par = theta, fn = logLikQml, model = model,
                        gr = gradient, method = optimizer, sign = -1,
                        control = control, ...)

    est$objective  <- est$value
    est$iterations <- est$counts[["function"]]
  } else stop2("Unrecognized optimizer, must be either 'nlminb' or 'L-BFGS-B'")

  if (verbose) cat("\n")

  est
}
