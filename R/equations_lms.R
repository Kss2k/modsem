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
  tGamma <- vector("numeric", length = length(w))

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
    gradientCompLogLikLms(theta = theta, model = model, P = P, sign = -1,
                      epsilon = epsilon)
  }

  objective <- function(theta) {
    compLogLikLms(theta = theta, model = model, P = P, sign = -1, epsilon = epsilon)
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


compLogLikLms <- function(theta, model, P, sign = -1, ...) {
  tryCatch({
    modFilled <- fillModel(model = model, theta = theta, method = "lms")
    sign * completeLogLikLmsCpp(modelR=modFilled, P=P, quad=P$quad)
  }, error = \(e) NA)
}


gradientCompLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6) {
  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign, 
                       epsilon = epsilon, FGRAD = gradLogLikLmsCpp, FOBJECTIVE = compLogLikLms)
}


gradientAllLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6,
                                 FGRAD, FOBJECTIVE) {
  hasCovModel <- model$gradientStruct$hasCovModel

  if (hasCovModel) gradient <- \(...) complicatedGradientAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE)
  else             gradient <- \(...) simpleGradientAllLogLikLms(..., FGRAD = FGRAD)

  c(gradient(theta = theta, model = model, P = P, sign = sign, epsilon = epsilon))
}


complicatedGradientAllLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-4, FOBJECTIVE) {
  # when we cannot simplify gradient using simpleGradientLogLikLms, we use 
  # good old forward difference...
  baseLL <- FOBJECTIVE(theta, model = model, P = P, sign = sign)
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (FOBJECTIVE(theta, model = model, P = P, sign = sign) - baseLL) / epsilon
  })
}


simpleGradientAllLogLikLms <- function(theta, model, P, data, sign = -1, epsilon = 1e-6, FGRAD) {
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

  grad <- FGRAD(modelR = modelR, 
                P = P, 
                block = block, 
                row = row, 
                col = col, 
                eps=epsilon)

  if (length(nlinDerivs)) {
    evalTheta  <- model$gradientStruct$evalTheta
    param.full <- colnames(Jacobian)
    param.part <- rownames(Jacobian)
    THETA      <- list2env(as.list(evalTheta(theta)))

    for (dep in names(nlinDerivs)) {
      derivs  <- nlinDerivs[[dep]]

      for (indep in names(derivs)) {
        deriv <- eval(expr = derivs[[indep]], envir = THETA)

        match.full <- param.full == dep
        match.part <- param.part == indep
        Jacobian[match.part, match.full] <- deriv
      }
    }
  }

  sign * Jacobian %*% grad
}


obsLogLikLms <- function(theta, model, data, P, sign = 1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  ll <- observedLogLikLmsCpp(modFilled, data = data, P = P, ncores = ThreadEnv$n.threads)
  ll * sign
}


gradientObsLogLikLms <- function(theta, model, data, P, sign = 1, epsilon = 1e-6) {
  FGRAD <- function(modelR, P, block, row, col, eps) {
    gradObsLogLikLmsCpp(modelR = modelR, data = data, P = P,
                        block = block, row = row, col = col,
                        eps = eps, ncores = ThreadEnv$n.threads)
  }

  FOBJECTIVE <- function(theta, model, P, sign) {
    obsLogLikLms(theta = theta, model = model, data = data, P = P, sign = sign)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign, 
                       epsilon = epsilon, FGRAD = FGRAD, FOBJECTIVE = FOBJECTIVE)
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


# gradient function of obsLogLikLms_i
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


hessianAllLogLikLms <- function(theta, model, P, sign = -1,
                                 FHESS, FOBJECTIVE, .relStep = Machine$double.eps ^ (1/5)) {
  hasCovModel <- model$gradientStruct$hasCovModel

  if (hasCovModel) hessian <- \(...) complicatedHessianAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
  else             hessian <- \(...) simpleHessianAllLogLikLms(..., FHESS = FHESS, .relStep = .relStep)

  hessian(theta = theta, model = model, P = P, sign = sign)
}


compHessianAllLogLikLms <- function(theta, model, P, data, sign = -1, 
                                    .relStep = .Machine$double.eps ^ (1/5), 
                                    FOBJECTIVE) {
  nlme::fdHess(pars = theta, FOBJECTIVE, model = model, P = P, 
               data = data, sign = sign, .relStep = .relStep)
}


simpleHessianAllLogLikLms <- function(theta, model, P, data, sign = -1, 
                                      .relStep = .Machine$double.eps ^ (1/5), 
                                      FHESS) {
  # simple gradient which should work if constraints are well-behaved functions 
  # which can be derivated by Deriv::Deriv, and there is no covModel
  modelR      <- fillModel(model=model, theta=theta, method="lms")
  locations   <- model$gradientStruct$locations
  Jacobian    <- model$gradientStruct$Jacobian
  Jacobian2   <- model$gradientStruct$Jacobian2
  nlinDerivs  <- model$gradientStruct$nlinDerivs
  nlinDerivs2 <- model$gradientStruct$nlinDerivs2

  block     <- locations$block
  row       <- locations$row
  col       <- locations$col
  param     <- locations$param

  HESS <- FHESS(modelR = modelR, 
                P = P, 
                block = block, 
                row = row, 
                col = col, 
                .relStep = .relStep)

  H    <- HESS$Hessian
  grad <- HESS$gradient

  if (length(nlinDerivs)) {
    evalTheta  <- model$gradientStruct$evalTheta
    param.full <- colnames(Jacobian)
    param.part <- rownames(Jacobian)
    THETA      <- list2env(as.list(evalTheta(theta)))

    for (dep in names(nlinDerivs)) {
      derivs1 <- nlinDerivs[[dep]]
      derivs2 <- nlinDerivs2[[dep]]

      for (indep in names(derivs1)) {
        deriv1 <- eval(expr = derivs1[[indep]], envir = THETA)
        deriv2 <- eval(expr = derivs2[[indep]], envir = THETA)

        match.full <- param.full == dep
        match.part <- param.part == indep

        Jacobian[match.part, match.full] <- deriv1
        Jacobian2[match.part, match.full] <- deriv2
      }
    }
  }

  term1 <- Jacobian %*% H %*% t(Jacobian)          
  term2 <- diag(drop(Jacobian2 %*% grad), nrow = nrow(Jacobian))

  sign * (term1 + term2)
}


complicatedHessianAllLogLikLms <- function(theta, model, P, sign = -1, FOBJECTIVE, 
                                           .relStep = .Machine$double.eps ^ (1/5)) {
  nlme::fdHess(pars = theta, fun = FOBJECTIVE, model = model, P = P, sign = sign,
               .relStep = .relStep)$Hessian
}


hessianObsLogLikLms <- function(theta, model, data, P, sign = -1, 
                                .relStep = .Machine$double.eps ^ (1/5)) {

  FHESS <- function(modelR, P, block, row, col, eps, .relStep) {
    hessObsLogLikLmsCpp(modelR = modelR, data = data, P = P,
                        block = block, row = row, col = col,
                        relStep = .relStep, minAbs = 0.0,
                        ncores = ThreadEnv$n.threads)
  }

  FOBJECTIVE <- function(theta, model, P, sign) {
    obsLogLikLms(theta = theta, model = model, data = data, P = P, sign = sign)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign, 
                      FHESS = FHESS, FOBJECTIVE = FOBJECTIVE,
                      .relStep = .relStep)
}


hessianCompLogLikLms <- function(theta, model, P, sign = -1, 
                                 .relStep = .Machine$double.eps ^ (1/5)) {

  FHESS <- function(modelR, P, block, row, col, eps, .relStep) {
    hessCompLogLikLmsCpp(modelR = modelR, P = P,
                         block = block, row = row, col = col,
                         relStep = .relStep, minAbs = 0.0,
                         ncores = ThreadEnv$n.threads)
  }

  FOBJECTIVE <- function(theta, model, P, sign) {
    compLogLikLms(theta = theta, model = model, P = P, sign = sign)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign, 
                      FHESS = FHESS, FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
}
