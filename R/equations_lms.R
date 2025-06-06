muLms <- function(model, z1) { # for testing purposes
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


sigmaLms <- function(model, z1) { # for testing purposes
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
    t(Binv %*% (Gx %*% A + t(kronZ) %*% Oxx %*% A)) %*% t(lY)
  Syy <- lY %*%
    (Binv %*% (Gx %*% A + t(kronZ) %*% Oxx %*% A)) %*%
    OI %*%
    t(Binv %*% (Gx %*% A + t(kronZ) %*% Oxx %*% A)) %*% t(lY) +
    lY %*% (Binv %*% psi %*% t(Binv)) %*% t(lY) + dY
  rbind(
    cbind(Sxx, Sxy),
    cbind(t(Sxy), Syy)
  )
}


llSingleLms <- function(z, modFilled, data) {
  mu <- muLmsCpp(model = modFilled, z = z)
  sigma <- sigmaLmsCpp(model = modFilled, z = z)
  sum(dmvn(data, mean = mu, sigma = sigma, log = TRUE))
}


flatLogLikLms <- function(z, modFilled, data) {
  vapply(seq_len(nrow(z)), FUN.VALUE = numeric(1L), FUN = function(i) {
    llSingleLms(z = z[i, , drop=FALSE], modFilled = modFilled, data = data)
  })
}


estepLms <- function(model, theta, data, last.integral = NULL, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  if (model$quad$adaptive) {
    quad <- adaptiveGaussQuadrature(
      flatLogLikLms, 
      a = model$quad$a, 
      b = model$quad$b, 
      modFilled = modFilled, 
      data = data,
      m.start = model$quad$m.start, 
      m.max = model$quad$m, 
      k = model$quad$k,
      tol = model$quad$adaptive.quad.tol,
      total.integral = last.integral
    )
  } else quad <- model$quad

  V <- as.matrix(quad$n)
  w <- quad$w

  P <- matrix(0, nrow = nrow(data), ncol = length(w))

  for (i in seq_along(w)) {
    P[, i] <- dmvn(data, mean = muLmsCpp(model = modFilled, z = V[i, ]),
                   sigma = sigmaLmsCpp(model = modFilled, z = V[i, ]),
                   log = FALSE) * w[[i]]
  }

  density        <- rowSums(P)
  observedLogLik <- sum(log(density))
  P              <- P / density
  
  wMeans <- vector("list", length = length(w))
  wCovs  <- vector("list", length = length(w))
  tGamma <- vector("list", length = length(w))

  for (i in seq_along(w)) {
    p    <- P[, i]
    wm   <- colSums(data * p) / sum(p)
    X    <- data - matrix(wm, nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
    wcov <- t(X) %*% (X * p)

    P[, i]      <- p
    wMeans[[i]] <- wm
    wCovs[[i]]  <- wcov
    tGamma[[i]] <- sum(p)
  }

  list(P = P, mean = wMeans, cov = wCovs, tgamma = tGamma, V = V, w = w, 
       obsLL = observedLogLik, integral.quad = quad$integral)
}


logLikLms <- function(theta, model, data, P, sign = -1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  k <- model$quad$k
  V <- P$V
  n <- nrow(data)
  d <- ncol(data)
  # summed log probability of observing the data given the parameters
  # weighted my the posterior probability calculated in the E-step
  r <- vapply(seq_len(nrow(V)), FUN.VALUE = numeric(1L), FUN = function(i) {
    if (P$tgamma[[i]] < .Machine$double.xmin) return(0)

    totalDmvnWeightedCpp(mu=muLmsCpp(model=modFilled, z=V[i, ]),
                         sigma=sigmaLmsCpp(model=modFilled, z=V[i, ]), 
                         nu=P$mean[[i]], S=P$cov[[i]], tgamma=P$tgamma[[i]], 
                         n = n, d = d)

  })

  
  sign * sum(r)
}


gradientLogLikLms <- function(theta, model, data, P, sign = -1, epsilon = 1e-4) {
  baseLL <- logLikLms(theta, model = model, data = data, P = P, sign = sign)
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (logLikLms(theta, model = model, data = data, P = P, sign = sign) - baseLL) / epsilon
  })
}


# log likelihood for each observation -- not all
logLikLms_i <- function(theta, model, data, P, sign = -1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  k <- model$quad$k
  V <- P$V

  # summed log probability of observing the data given the parameters
  # weighted my the posterior probability calculated in the E-step
  r <- lapplyMatrix(seq_len(nrow(V)), FUN = function(i) {
    dmvn(data,
      mean = muLmsCpp(model = modFilled, z = V[i, ]),
      sigma = sigmaLmsCpp(model = modFilled, z = V[i, ]),
      log = TRUE
    ) * P$P[, i]
  }, FUN.VALUE = numeric(nrow(data)))

  sign * apply(r, MARGIN = 1, FUN = sum)
}


# Maximization step of EM-algorithm (see Klein & Moosbrugger, 2000)
mstepLms <- function(theta, model, data, P,
                     max.step,
                     verbose = FALSE,
                     control = list(),
                     optimizer = "nlminb",
                     optim.method = "L-BFGS-B",
                     epsilon = 1e-6,
                     ...) {
  gradient <- function(theta, model, data, P, sign) {
    gradientLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      data = data, epsilon = epsilon)
  }

  if (optimizer == "nlminb") {
    if (is.null(control$iter.max)) control$iter.max <- max.step
    preObjective <- logLikLms(theta=theta, model=model, P=P, sign=1, data=data)
    est <- stats::nlminb(start = theta, objective = logLikLms, data = data,
                         model = model, P = P, gradient = gradient,
                         sign = -1,
                         upper = model$info$bounds$upper,
                         lower = model$info$bounds$lower, control = control,
                         ...) |> suppressWarnings()
    postObjective <- logLikLms(theta=est$par, model=model, P=P, sign=1, data=data)

  } else if (optimizer == "L-BFGS-B") {
    if (is.null(control$maxit)) control$maxit <- max.step
    est <- stats::optim(par = theta, fn = logLikLms, data = data,
                        model = model, P = P, gr = gradient,
                        method = optimizer, control = control,
                        sign = -1, lower = model$info$bounds$lower,
                        upper = model$info$bounds$upper, ...)

    est$objective  <- est$value
    est$iterations <- est$counts[["function"]]
  } else {
    stop2("Unrecognized optimizer, must be either 'nlminb' or 'L-BFGS-B'")
  }

  est
}
