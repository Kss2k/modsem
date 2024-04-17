muLms <- function(model, z1) {
  matrices <- model$matrices
  A <- matrices$A
  subA <- A[seq_len(model$info$numXis), seq_len(model$info$numXis)]
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
  collapseEta <- rep(1, model$info$numEtas)
  cOex <- matrices$selectionMatrixOmegaEtaXi

  k <- model$quad$k
  zVec <- c(z1, rep(0, model$info$numXis - k))
  zMat <- zToMatrix(zVec, model$info$numEtas)
  if (ncol(Ie) == 1) Binv <- Ie else Binv <- matlib::inv(Ie - Ge - (t(zMat) %*% t(A) %*% Oex) %*% cOex)

  muX <- tX + lX %*% subA %*% zVec
  muY <- tY + 
    lY %*% (Binv %*% (a + 
            Gx %*% subA %*% zVec + 
            (t(zMat) %*% t(A) %*% Oxx %*% A %*% zMat) %*%
            collapseEta))
  rbind(muX, muY)
}


sigmaLms <- function(model, z1) {
  matrices <- model$matrices
  Oxx <- matrices$omegaXiXi
  Oex <- matrices$omegaEtaXi
  Ie <- matrices$Ieta
  A <- matrices$A
  subA <- A[seq_len(model$info$numXis), seq_len(model$info$numXis)]
  O <- matrices$omega
  lY <- matrices$lambdaY
  lX <- matrices$lambdaX
  Gx <- matrices$gammaXi
  Ge <- matrices$gammaEta
  dX <- matrices$thetaDelta
  dY <- matrices$thetaEpsilon
  psi <- matrices$psi
  k <- model$quad$k
  zVec <- c(z1, rep(0, model$info$numXis - k))
  zMat <- zToMatrix(zVec, model$info$numEtas)
  cOex <- matrices$selectionMatrixOmegaEtaXi
  if (ncol(Ie) == 1) Binv <- Ie else Binv <- matlib::inv(Ie - Ge - (t(zMat) %*% t(A) %*% Oex) %*% cOex)
  
  OI <- diag(1, model$info$numXis)
  diag(OI) <- c(rep(0, k), rep(1, model$info$numXis - k))

  select <- matrices$selectionMatrixOmega

  phiEta <- 
    ((Binv %*% (Gx %*% subA + (t(zMat) %*% t(A) %*% Oxx %*% A) %*% select)) %*%
                 OI %*%
    t(Binv %*% (Gx %*% subA + (t(zMat) %*% t(A) %*% Oxx %*% A) %*% select))) +
    (Binv %*% psi %*% t(Binv))

  Sxx <- lX %*% subA %*% OI %*%  
    t(subA) %*% t(lX) + dX
  Sxy <- lX %*% (subA %*% OI %*% 
                 t(Binv %*% (Gx %*% subA + (t(zMat) %*% t(A) %*% Oxx %*% A) %*% 
                             select))) %*% t(lY)
  Syy <- lY %*% 
    ((Binv %*% (Gx %*% subA + (t(zMat) %*% t(A) %*% Oxx %*% A) %*% select)) %*%
                 OI %*%
    t(Binv %*% (Gx %*% subA + (t(zMat) %*% t(A) %*% Oxx %*% A) %*% select))) %*% t(lY) +
    lY %*% (Binv %*% psi %*% t(Binv)) %*% t(lY) + dY
  rbind(cbind(Sxx,  Sxy),
      cbind(t(Sxy), Syy))
}


estepLms <- function(model, theta, dat, ...) {
  if (countFreeParams(model) != length(theta))
    stop("length paramaters does not match free parameters in model")
  modFilled <- fillModel(model = model, theta = theta)
  V <- modFilled$quad$n       # matrix of node vectors m x k
  w <- modFilled$quad$w       # weights
  # the probability of each observation is derived as the sum of the probabilites 
  # of observing the data given the parameters at each point in the k dimensional
  # space of the distribution of the non-linear latent variables. To approximate
  # this we use individual nodes sampled from the k-dimensional space with a 
  # corresponding probability (weight w) for each node appearing. I.e., whats the 
  # probability of observing the data given the nodes, and what is the probability
  # of observing the given nodes (i.e., w). Sum the row probabilities = 1
  P <- matrix(0, nrow = nrow(dat), ncol = length(w))
  sapply(seq_along(w), FUN = function(i) {
      P[,i] <<- w[[i]] * dMvn(dat, mean = muLmsCpp(model = modFilled, z = V[i,]),
                              sigma = sigmaLmsCpp(model = modFilled, z = V[i,]))
  })
  P / rowSums(P)
}


stochasticGradient <- function(theta, model, dat, P, 
                               sampleGrad = NULL, ...) {
  baseline <- logLikLms(theta, model, dat, P)
  grad <- rep(0, length(theta))
  if (!is.null(sampleGrad)) params <- sample(seq_along(theta), sampleGrad)
  else params <- seq_along(theta)
  for (i in params) {
    theta[i] <- theta[i] + 1e-12
    newLik <- logLikLms(theta, model, dat, P) 
    grad[i] <- (newLik - baseline) / 1e-12
  }
  grad
} 


logLikLms <- function(theta, model, dat, P, sampleGrad = NULL, ...) {
  modFilled <- fillModel(model = model, theta = theta)
  k <- model$quad$k 
  V <- modFilled$quad$n
  # summed log probability of observing the data given the parameters
  # weighted my the posterior probability calculated in the E-step
  r <- vapply(seq_len(nrow(V)), FUN.VALUE = numeric(1L), FUN = function(i){
    lls <- sum(dMvn(dat, mean = muLmsCpp(model = modFilled, z = V[i,]),
                    sigma = sigmaLmsCpp(model = modFilled, z = V[i,]),
                    log = TRUE) * P[,i])
    lls
  }) |> sum()
  -r
}


# Maximization step of EM-algorithm (see Klein & Moosbrugger, 2000)
mstepLms <- function(theta, model, dat, P, negHessian = FALSE,
                      maxstep,
                      verbose = FALSE,
                      control=list(), sampleGrad,...) {
  if (is.null(sampleGrad)) stochasticGradient <- NULL
  if (is.null(control$iter.max)) control$iter.max <- maxstep
  est <- stats::nlminb(start = theta, objective = logLikLms, dat = dat,
                gradient = stochasticGradient, sampleGrad = sampleGrad,
                model = model, P = P, 
                upper = model$info$bounds$upper,
                lower = model$info$bounds$lower, control = control,
               ...) |> suppressWarnings()
  if (negHessian) {
    if (verbose) cat("Calculating Hessian\n")
    est$hessian <- nlme::fdHess(pars = est$par, fun=logLikLms, 
                                model = model, dat=dat, P=P,
                                .relStep = .Machine$double.eps^(1/5))$Hessian
  }
  est
}


# This is probably inneficient, but i'll fix it later
zToMatrix <- function(zVec, nEta) {
  mat <- matrix(0, nrow = nEta * length(zVec), ncol = nEta)
  for (i in seq_len(nEta)) 
    mat[seq_len(length(zVec)) + (i - 1) * length(zVec), i] <- zVec
  mat
}


collapsePartitionedMatrixRow <- function(x, nEtas) {
  out <- matrix(0, nrow = nrow(x) / nEtas, ncol = ncol(x))
  rowOffset <- 0
  for (i in seq_len(nEtas)) {
    out <- out + x[seq_len(nrow(out)) + rowOffset, ]
    rowOffset <- rowOffset + nrow(out) * nEtas
  }
  out
}
