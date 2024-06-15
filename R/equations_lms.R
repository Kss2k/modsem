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


function(model, N) {
  V <- model$quad$n       # matrix of node vectors m x k
  w <- model$quad$w       # weights
  P <- matrix(rep(w, N), nrow = N, ncol = length(w),
              byrow = TRUE)
  P
}


logLikLms <- function(theta, model, data, P, sign = -1, ...) {
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
  sign * r
}


gradientLogLikLms <- function(theta, model, data, P, sign = -1, dt = 1e-6) {
  baseLL <- logLikLms(theta, model = model, data = data, P = P, sign = sign) 
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
     theta[[i]] <- theta[[i]] + dt 
     (logLikLms(theta, model = model, data = data, P = P, sign = sign) - baseLL) / dt
  })
}


# log likelihood for each observation -- not all
logLikLms_i <- function(theta, model, data, P, sign = -1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  k <- model$quad$k 
  V <- modFilled$quad$n
  # summed log probability of observing the data given the parameters
  # weighted my the posterior probability calculated in the E-step
  r <- lapplyMatrix(seq_len(nrow(V)), FUN = function(i){
    lls <- dMvn(data, mean = muLmsCpp(model = modFilled, z = V[i,]),
                sigma = sigmaLmsCpp(model = modFilled, z = V[i,]),
                log = TRUE) * P[,i]
    lls
  }, FUN.VALUE = numeric(nrow(data)))

  sign * apply(r, MARGIN = 1, FUN = sum)
}


# gradient function of logLikLms_i
gradientLogLikLms_i <- function(theta, model, data, P, sign = -1, dt = 1e-6) {
  baseLL <- logLikLms_i(theta, model, data = data, P = P, sign = sign) 
  lapplyMatrix(seq_along(theta), FUN = function(i) {
     theta[[i]] <- theta[[i]] + dt 
     (logLikLms_i(theta, model, data = data, P = P, sign = sign) - baseLL) / dt
  }, FUN.VALUE = numeric(nrow(data)))
}


# Maximization step of EM-algorithm (see Klein & Moosbrugger, 2000)
mstepLms <- function(theta, model, data, P, 
                     max.step,
                     verbose = FALSE,
                     control = list(), ...) {
  if (is.null(control$iter.max)) control$iter.max <- max.step
  est <- stats::nlminb(start = theta, objective = logLikLms, data = data,
                       model = model, P = P, 
                       upper = model$info$bounds$upper,
                       lower = model$info$bounds$lower, control = control,
                       ...) |> suppressWarnings()
  est
}
