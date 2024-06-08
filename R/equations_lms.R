muLms <- function(model, z1) {
  matrices <- model$matrices
  A <- matrices$A
  Oxx <- matrices$omegaXiXi
  Oex <- matrices$omegaEtaXi
  Ie <- matrices$Ieta
  lY <- matrices$lambdaY
  lX <- matrices$lambdaX
  tY <- matrices$tauY 
  tX <- matrices$tauX
  Gx <- matrices$gammaXi
  Ge <- matrices$gammaEta
  a <- matrices$alpha
  psi <- matrices$psi

  k <- model$quad$k
  zVec <- c(z1[0:k], rep(0, model$info$numXis - k))
  kronZ <- kronecker(Ie, A %*% zVec)
  if (ncol(Ie) == 1) Binv <- Ie else Binv <- solve(Ie - Ge - t(kronZ) %*% Oex)

  muX <- tX + lX %*% A %*% zVec
  muY <- tY + 
    lY %*% (Binv %*% (a + 
            Gx %*% A %*% zVec + 
            t(kronZ) %*% Oxx %*% A %*% zVec))
  rbind(muX, muY)
}


sigmaLms <- function(model, z1) {
  matrices <- model$matrices
  Oxx <- matrices$omegaXiXi
  Oex <- matrices$omegaEtaXi
  Ie <- matrices$Ieta
  A <- matrices$A
  lY <- matrices$lambdaY
  lX <- matrices$lambdaX
  Gx <- matrices$gammaXi
  Ge <- matrices$gammaEta
  dX <- matrices$thetaDelta
  dY <- matrices$thetaEpsilon
  psi <- matrices$psi
  k <- model$quad$k
  zVec <- c(z1[0:k], rep(0, model$info$numXis - k))
  kronZ <- kronecker(Ie, A %*% zVec)
  if (ncol(Ie) == 1) Binv <- Ie else Binv <- solve(Ie - Ge - t(kronZ) %*% Oex)
  
  OI <- diag(1, model$info$numXis)
  diag(OI) <- c(rep(0, k), rep(1, model$info$numXis - k))

  Sxx <- lX %*% A %*% OI %*%  
    t(A) %*% t(lX) + dX
  Sxy <- lX %*% A %*% OI %*% 
                 t(Binv %*% (Gx %*% A + t(kronZ) %*% Oxx %*% A)) %*%  t(lY)
  Syy <- lY %*% 
    (Binv %*% (Gx %*% A + t(kronZ) %*% Oxx %*% A)) %*%
                 OI %*%
    t(Binv %*% (Gx %*% A + t(kronZ) %*% Oxx %*% A)) %*% t(lY) +
    lY %*% (Binv %*% psi %*% t(Binv)) %*% t(lY) + dY
  rbind(cbind(Sxx,  Sxy),
      cbind(t(Sxy), Syy))
}


estepLms <- function(model, theta, data, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  V <- modFilled$quad$n       # matrix of node vectors m x k
  w <- modFilled$quad$w       # weights
  # the probability of each observation is derived as the sum of the probabilites 
  # of observing the data given the parameters at each point in the k dimensional
  # space of the distribution of the non-linear latent variables. To approximate
  # this we use individual nodes sampled from the k-dimensional space with a 
  # corresponding probability (weight w) for each node appearing. I.e., whats the 
  # probability of observing the data given the nodes, and what is the probability
  # of observing the given nodes (i.e., w). Sum the row probabilities = 1
  P <- matrix(0, nrow = nrow(data), ncol = length(w))
  sapply(seq_along(w), FUN = function(i) {
      P[,i] <<- w[[i]] * dMvn(data, mean = muLmsCpp(model = modFilled, z = V[i,]),
                              sigma = sigmaLmsCpp(model = modFilled, z = V[i,]))
  })
  P / rowSums(P)
}



derivativeMuSigmaTheta <- function(model, theta) {

}
stochasticGradient <- function(theta, model, data, P, 
                               sampleGrad = NULL, ...) {
  baseline <- logLikLms(theta, model, data, P)
  grad <- rep(0, length(theta))
  if (!is.null(sampleGrad)) params <- sample(seq_along(theta), sampleGrad)
  else params <- seq_along(theta)
  for (i in params) {
    theta[i] <- theta[i] + 1e-12
    newLik <- logLikLms(theta, model, data, P) 
    grad[i] <- (newLik - baseline) / 1e-12
  }
  grad
} 


logLikLms <- function(theta, model, data, P, sampleGrad = NULL, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  k <- model$quad$k 
  V <- modFilled$quad$n
  # summed log probability of observing the data given the parameters
  # weighted my the posterior probability calculated in the E-step
  r <- vapply(seq_len(nrow(V)), FUN.VALUE = numeric(1L), FUN = function(i){
    lls <- sum(dMvn(data, mean = muLmsCpp(model = modFilled, z = V[i,]),
                    sigma = sigmaLmsCpp(model = modFilled, z = V[i,]),
                    log = TRUE) * P[,i])
    lls
  }) |> sum()
  -r
}


# Maximization step of EM-algorithm (see Klein & Moosbrugger, 2000)
mstepLms <- function(theta, model, data, P, hessian = FALSE,
                     maxstep,
                     verbose = FALSE,
                     control = list(), sampleGrad,...) {
  if (is.null(sampleGrad)) stochasticGradient <- NULL
  if (is.null(control$iter.max)) control$iter.max <- maxstep
  est <- stats::nlminb(start = theta, objective = logLikLms, data = data,
                       gradient = stochasticGradient, sampleGrad = sampleGrad,
                       model = model, P = P, 
                       upper = model$info$bounds$upper,
                       lower = model$info$bounds$lower, control = control,
                       ...) |> suppressWarnings()
  if (hessian) {
    if (verbose) cat("Calculating Hessian\n")
    est$hessian <- nlme::fdHess(pars = est$par, fun = logLikLms, 
                                model = model, data = data, P = P,
                                .relStep = .Machine$double.eps^(1/5))$Hessian
  }
  est
}
