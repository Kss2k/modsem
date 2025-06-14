estepLms <- function(model, theta, data, lastQuad = NULL, recalcQuad = FALSE, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  if (model$quad$adaptive && (recalcQuad || is.null(lastQuad))) {
    m <- model$quad$m
    a <- model$quad$a
    b <- model$quad$b
    m <- model$quad$m
    k <- model$quad$k

    if (!is.null(lastQuad)) m.ceil <- lastQuad$m.ceil 
    else if (k > 1) m.ceil <- m
    else m.ceil <- round(estMForNodesInRange(m, a = -5, b = 5))

    quad <- adaptiveGaussQuadrature(
      fun = densityLms, collapse = \(x) sum(log(rowSums(x))),
      modFilled = modFilled, data = data, a = a, b = b, m = m, 
      k = k, m.ceil = m.ceil
    )

  } else if (model$quad$adaptive) {
    quad <- lastQuad

  } else quad <- model$quad

  V <- quad$n
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
       obsLL = observedLogLik, quad = quad)
}


# Maximization step of EM-algorithm (see Klein & Moosbrugger, 2000)
mstepLms <- function(theta, model, P,
                     max.step,
                     verbose = FALSE,
                     control = list(),
                     optimizer = "nlminb",
                     optim.method = "L-BFGS-B",
                     epsilon = 1e-6,
                     ...) {
  gradient <- function(theta) {
    gradientLogLikLms(theta = theta, model = model, P = P, sign = -1,
                      epsilon = epsilon)
  }

  objective <- function(theta) {
    logLikLms(theta = theta, model = model, P = P, sign = -1, epsilon = epsilon)
  }

  if (optimizer == "nlminb") {
    if (is.null(control$iter.max)) control$iter.max <- max.step
    est <- stats::nlminb(start = theta, objective = objective, 
                         gradient = gradient,
                         upper = model$info$bounds$upper,
                         lower = model$info$bounds$lower, control = control,
                         ...) |> suppressWarnings()

  } else if (optimizer == "L-BFGS-B") {
    if (is.null(control$maxit)) control$maxit <- max.step
    est <- stats::optim(par = theta, fn = objective, gr = gradient,
                        method = optim.method, control = control,
                        lower = model$info$bounds$lower,
                        upper = model$info$bounds$upper, ...)

    est$objective  <- est$value
    est$iterations <- est$counts[["function"]]
  } else {
    stop2("Unrecognized optimizer, must be either 'nlminb' or 'L-BFGS-B'")
  }

  est
}


logLikLms <- function(theta, model, P, sign = -1, ...) {
  tryCatch({
    modFilled <- fillModel(model = model, theta = theta, method = "lms")
    sign * completeLogLikLmsCpp(modelR=modFilled, P=P, quad=P$quad)
  }, error = \(e) NA)
}


gradientLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6) {
  hasCovModel <- model$gradientStruct$hasCovModel

  if (hasCovModel) gradient <- complicatedGradientLogLikLms
  else             gradient <- simpleGradientLogLikLms

  gradient(theta = theta, model = model, P = P, sign = sign, epsilon = epsilon)
}


simpleGradientLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6) {
  # simple gradient which should work if constraints are well-behaved functions 
  # which can be derivated by Deriv::Deriv, and there is no covModel
  modelR     <- fillModel(model=model, theta=theta, method="lms")
  locations  <- model$gradientStruct$locations
  Jacobian   <- model$gradientStruct$Jacobian
  nlinDerivs <- model$gradientStruct$nlinDerivs

  block     <- locations$block
  row       <- locations$row
  col       <- locations$col
  param     <- locations$param

  grad <- gradLogLikLmsCpp(modelR, P = P, block = block, 
                           row = row, col = col, eps=epsilon)
  # grad <- structure(c(grad), names = param)

  if (length(nlinDerivs)) {
    evalTheta  <- model$gradientStruct$evalTheta
    param.full <- colnames(Jacobian)
    param.part <- rownames(Jacobian)
    THETA      <- list2env(as.list(evalTheta(theta)))

    for (dep in names(nlinDerivs)) {
      derivs <- eval(expr = nlinDerivs[[dep]], envir = THETA)

      for (indep in names(derivs)) {
        match.full <- param.full == dep
        match.part <- param.part == indep
        Jacobian[match.part, match.full] <- derivs[[indep]]
      }
    }
  }

  c(sign * Jacobian %*% grad)
}


complicatedGradientLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-4) {
  # when we cannot simplify gradient using simpleGradientLogLikLms, we use 
  # good old forward difference...
  baseLL <- logLikLms(theta, model = model, P = P, sign = sign)
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (logLikLms(theta, model = model, P = P, sign = sign) - baseLL) / epsilon
  })
}


obsLogLikLms <- function(theta, model, data, P, sign = 1, ...) {
  sum(obsLogLikLms_i(theta, model = model, data = data, P = P, sign = sign))
}


gradientObsLogLikLms <- function(theta, model, data, P, sign = -1, epsilon = 1e-4) {
  baseLL <- logLikLms(theta, model = model, data = data, P = P, sign = sign)

  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (obsLogLikLms(theta, model = model, data = data, P = P, sign = sign) - baseLL) / epsilon
  })
}


obsLogLikLms_i <- function(theta, model, data, P, sign = 1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  V <- P$V
  w <- P$w
  N <- nrow(data)
  m <- nrow(V)
  px <- numeric(N)

  for (i in seq_len(m)) {
    z_i     <- V[i, ]
    mu_i    <- muLmsCpp(  model = modFilled, z = z_i)
    sigma_i <- sigmaLmsCpp(model = modFilled, z = z_i)
    dens_i  <- dmvn(data, mean = mu_i, sigma = sigma_i, log = FALSE)
    px <- px + w[i] * dens_i
  }
 
  sign * log(px)
}


# gradient function of logLikLms_i
gradientObsLogLikLms_i <- function(theta, model, data, P, sign = 1, epsilon = 1e-4) {
  baseLL <- obsLogLikLms_i(theta, model, data = data, P = P, sign = sign)

  lapplyMatrix(seq_along(theta), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (obsLogLikLms_i(theta, model, data = data, P = P, sign = sign) - baseLL) / epsilon
  }, FUN.VALUE = numeric(nrow(data)))
}


densitySingleLms <- function(z, modFilled, data) {
  mu <- muLmsCpp(model = modFilled, z = z)
  sigma <- sigmaLmsCpp(model = modFilled, z = z)
  dmvn(data, mean = mu, sigma = sigma)
}


densityLms <- function(z, modFilled, data) {
  if (is.null(dim(z))) z <- matrix(z, ncol = modFilled$quad$k)

  lapplyMatrix(seq_len(nrow(z)), FUN.VALUE = numeric(NROW(data)), FUN = function(i) {
    densitySingleLms(z = z[i, , drop=FALSE], modFilled = modFilled, data = data)
  })
}


