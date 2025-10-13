estepLms <- function(model, theta, lastQuad = NULL, recalcQuad = FALSE,
                     adaptive.quad.tol = 1e-12, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  P_GROUPS <- vector("list", length = model$info$n.groups)
  for (g in seq_len(model$info$n.groups)) {
    P_GROUPS[[g]] <- estepLmsGroup(
      submodel = modFilled$models[[g]], lastQuad = lastQuad,
      recalcQuad = recalcQuad, adaptive.quad.tol = adaptive.quad.tol, ...
    )
  }

  P_GROUPS
}


estepLmsGroup <- function(submodel, lastQuad = NULL, recalcQuad = FALSE,
                          adaptive.quad.tol = 1e-12, ...) {
  data <- submodel$data

  if (submodel$quad$adaptive && (recalcQuad || is.null(lastQuad))) {
    m <- submodel$quad$m
    a <- submodel$quad$a
    b <- submodel$quad$b
    m <- submodel$quad$m
    k <- submodel$quad$k

    if (!is.null(lastQuad)) m.ceil <- lastQuad$m.ceil
    else if (k > 1) m.ceil <- m
    else m.ceil <- round(estMForNodesInRange(m, a = -5, b = 5))

    quad <- tryCatch({
        adaptiveGaussQuadrature(
          fun = densityLms, collapse = \(x) sum(log(rowSums(x))),
          modFilled = submodel, data = data, a = a, b = b, m = m,
          k = k, m.ceil = m.ceil, tol = adaptive.quad.tol,
        )
      }, error = function(e) {
        warning2("Calculation of adaptive quadrature failed!\n", e,
                 immediate. = FALSE)
        NULL
      }
    )

    if (is.null(quad)) {
      estep.fixed <- estepLms(
        submodel = submodel,
        theta = theta,
        data  = data,
        lastQuad = if (!is.null(lastQuad)) lastQuad else submodel$quad,
        recalcQuad = FALSE,
        ...
      )

      return(estep.fixed)
    }

    P <- quad$W * quad$F # P is already calculated
    V <- quad$n
    w <- quad$w

  } else {
    quad <- if (submodel$quad$adaptive) lastQuad else submodel$quad
    V    <- quad$n
    w    <- quad$w
    W    <- matrix(w, nrow = data$n, ncol = length(w), byrow = TRUE)
    P    <- W * densityLms(V, modFilled = subModel, data = data)
  }

  density        <- rowSums(P)
  observedLogLik <- sum(log(density))
  P              <- P / density

  wMeans <- vector("list", length = length(w))
  wCovs  <- vector("list", length = length(w))
  tGamma <- vector("list", length = length(w))

  for (i in seq_along(w)) {
    p <- P[, i]
    tGamma[[i]] <- sum(p)

    offset <- 1L

    wMeans[[i]] <- vector("list", length = length(data$ids))
    wCovs[[i]]  <- vector("list", length = length(data$ids))
    tGamma[[i]] <- numeric(length = length(data$ids))

    # wmean <- colSums(data$data.full * p, na.rm = TRUE) / sum(p)
    for (j in data$ids) {
      n.pattern <- data$n.pattern[[j]]
      end       <- offset + n.pattern - 1L

      data.id <- data$data.split[[j]]
      colidx  <- data$colidx[[j]]

      pj   <- p[offset:end]
      # wm   <- wmean[colidx]
      wm   <- colSums(data.id * pj) / sum(pj)
      X    <- data.id - matrix(wm, nrow=nrow(data.id), ncol=ncol(data.id), byrow=TRUE)
      wcov <- t(X) %*% (X * pj)

      wMeans[[i]][[j]] <- wm
      wCovs[[i]][[j]]  <- wcov
      tGamma[[i]][[j]] <- sum(pj)

      offset <- end + 1L
    }
  }

  list(P = P, mean = wMeans, cov = wCovs, tgamma = tGamma, V = V, w = w,
       obsLL = observedLogLik, quad = quad)
}


rawGradientGroupLms <- function(submodel, theta, P, data, epsilon) {
  GS <- submodel$gradientStruct
  if (is.null(GS$locations)) {
    stop("Gradient structure for subgroup is not available; cannot compute raw gradient.")
  }

  modelR     <- fillModel(model = submodel, theta = theta, method = "lms")
  locations  <- GS$locations

  gradLogLikLmsCpp(modelR    = modelR,
                   P         = P,
                   block     = locations$block,
                   row       = locations$row,
                   col       = locations$col,
                   symmetric = locations$symmetric,
                   colidxR   = data$colidx0,
                   npatterns = data$p,
                   n         = data$n.pattern,
                   d         = data$d.pattern,
                   eps       = epsilon,
                   ncores    = ThreadEnv$n.threads)
}


rawHessianGroupLms <- function(submodel, theta, P, data, .relStep) {
  GS <- submodel$gradientStruct
  if (is.null(GS$locations)) {
    stop("Hessian structure for subgroup is not available; cannot compute raw Hessian.")
  }

  modelR     <- fillModel(model = submodel, theta = theta, method = "lms")
  locations  <- GS$locations

  hessCompLogLikLmsCpp(modelR = modelR,
                       P = P,
                       block = locations$block,
                       row = locations$row,
                       col = locations$col,
                       symmetric = locations$symmetric,
                       colidxR = data$colidx0,
                       n = data$n.pattern,
                       d = data$d.pattern,
                       npatterns = data$p,
                       relStep = .relStep,
                       minAbs = 0.0,
                       ncores = ThreadEnv$n.threads)
}


rawHessianObsGroupLms <- function(submodel, theta, P, data, .relStep) {
  GS <- submodel$gradientStruct
  if (is.null(GS$locations)) {
    stop("Observed Hessian structure for subgroup is not available; cannot compute raw Hessian.")
  }

  modelR    <- fillModel(model = submodel, theta = theta, method = "lms")
  locations <- GS$locations

  hessObsLogLikLmsCpp(modelR = modelR,
                      dataR = data$data.split,
                      P = P,
                      block = locations$block,
                      row = locations$row,
                      col = locations$col,
                      symmetric = locations$symmetric,
                      colidxR = data$colidx0,
                      n = data$n.pattern,
                      npatterns = data$p,
                      relStep = .relStep,
                      minAbs = 0.0,
                      ncores = ThreadEnv$n.threads)
}


# Maximization step of EM-algorithm (see Klein & Moosbrugger, 2000)
mstepLms <- function(theta, model, P, data,
                     max.step,
                     verbose = FALSE,
                     control = list(),
                     optimizer = "nlminb",
                     optim.method = "L-BFGS-B",
                     epsilon = 1e-6,
                     ...) {
  gradient <- function(theta) {
    gradientCompLogLikLms(theta = theta, model = model, P = P, sign = -1,
                          data = data, epsilon = epsilon)
  }

  objective <- function(theta) {
    compLogLikLms(theta = theta, model = model, P = P, sign = -1, data = data,
                  epsilon = epsilon)
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


compLogLikLms <- function(theta, model, P, data, sign = -1, ...) {
  if (isMultiGroupModelDA(model)) {
    submodels <- model$groupModels
    P_groups <- P$groups
    stopifnot(length(P_groups) == length(submodels))
    total <- 0

    for (g in seq_along(submodels)) {
      submodel <- submodels[[g]]
      theta_g <- getThetaGroupDA(model, theta, g)
      P_g <- P_groups[[g]]
      data_g <- submodel$data

      total <- total + compLogLikLms(theta = theta_g, model = submodel,
                                     P = P_g, data = data_g, sign = 1, ...)
    }

    return(sign * total)
  }

  tryCatch({
    modFilled <- fillModel(model = model, theta = theta, method = "lms")
    sign * completeLogLikLmsCpp(modelR=modFilled, P=P, quad=P$quad,
                                colidxR = data$colidx0, n = data$n.pattern,
                                d = data$d.pattern, npatterns = data$p)
  }, error = \(e) NA)
}


gradientCompLogLikLms <- function(theta, model, P, data, sign = -1, epsilon = 1e-6) {
  if (isMultiGroupModelDA(model)) {
    GS <- model$gradientStruct
    submodels <- model$groupModels
    P_groups <- P$groups
    stopifnot(length(P_groups) == length(submodels))

    if (isTRUE(GS$hasCovModel)) {
      return(complicatedGradientAllLogLikLms(theta = theta, model = model, P = P,
                                             data = data, sign = sign, epsilon = epsilon,
                                             FOBJECTIVE = compLogLikLms))
    }

    raw_grad <- numeric(length(GS$param.full))

    for (g in seq_along(submodels)) {
      submodel <- submodels[[g]]
      theta_g  <- getThetaGroupDA(model, theta, g)
      P_g      <- P_groups[[g]]
      data_g   <- submodel$data

      grad_g <- rawGradientGroupLms(submodel, theta_g, P_g, data_g, epsilon)
      idx <- GS$groupColumns[[g]]
      raw_grad[idx] <- raw_grad[idx] + grad_g
    }

    Jacobian <- GS$Jacobian
    if (length(GS$nlinDerivs)) {
      evalTheta  <- GS$evalTheta
      param.full <- colnames(Jacobian)
      param.part <- rownames(Jacobian)
      THETA      <- list2env(as.list(evalTheta(theta)))

      for (dep in names(GS$nlinDerivs)) {
        derivs <- GS$nlinDerivs[[dep]]
        for (indep in names(derivs)) {
          deriv <- eval(expr = derivs[[indep]], envir = THETA)

          match.full <- param.full == dep
          match.part <- param.part == indep
          if (any(match.full) && any(match.part)) {
            Jacobian[match.part, match.full] <- deriv
          }
        }
      }
    }

    grad_theta <- as.vector(Jacobian %*% raw_grad)
    names(grad_theta) <- rownames(Jacobian)
    return(sign * grad_theta)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign, data = data,
                       epsilon = epsilon, FGRAD = gradLogLikLmsCpp, FOBJECTIVE = compLogLikLms)
}


gradientAllLogLikLms <- function(theta, model, P, data, sign = -1, epsilon = 1e-6,
                                 FGRAD, FOBJECTIVE) {
  hasCovModel <- model$gradientStruct$hasCovModel

  if (hasCovModel) gradient <- \(...) complicatedGradientAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE)
  else             gradient <- \(...) simpleGradientAllLogLikLms(..., FGRAD = FGRAD)

  c(gradient(theta = theta, model = model, P = P, data = data, sign = sign, epsilon = epsilon))
}


complicatedGradientAllLogLikLms <- function(theta, model, P, data, sign = -1, epsilon = 1e-4, FOBJECTIVE) {
  # when we cannot simplify gradient using simpleGradientLogLikLms, we use
  # good old forward difference...
  baseLL <- FOBJECTIVE(theta, model = model, P = P, data = data, sign = sign)
  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (FOBJECTIVE(theta, model = model, P = P, data = data, sign = sign) - baseLL) / epsilon
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
  symmetric <- locations$symmetric

  grad <- FGRAD(modelR    = modelR,
                P         = P,
                block     = block,
                row       = row,
                col       = col,
                symmetric = symmetric,
                colidxR   = data$colidx0,
                npatterns = data$p,
                n         = data$n.pattern,
                d         = data$d.pattern,
                eps       = epsilon,
                ncores    = ThreadEnv$n.threads)

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
  ll <- observedLogLikLmsCpp(modFilled, dataR = data$data.split, P = P,
                             colidxR = data$colidx0, n = data$n.pattern,
                             npatterns = data$p, ncores = ThreadEnv$n.threads)
  ll * sign
}


gradientObsLogLikLms <- function(theta, model, data, P, sign = 1, epsilon = 1e-6) {
  FGRAD <- function(modelR, P, block, row, col, symmetric, colidxR, npatterns,
                    eps, ncores, n, ...) {
    gradObsLogLikLmsCpp(modelR = modelR, dataR = data$data.split, P = P,
                        block = block, row = row, col = col,
                        symmetric = symmetric, colidxR = colidxR,
                        n = n, npatterns = npatterns, eps = eps,
                        ncores = ncores)
  }

  FOBJECTIVE <- function(theta, model, P, colidxR, npatterns, sign, ...) {
    obsLogLikLms(theta = theta, model = model, data = data$data.split,
                 P = P, colidxR = colidxR, npatterns = npatterns,
                 sign = sign)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, data = data,
                       FGRAD = FGRAD, FOBJECTIVE = FOBJECTIVE)
}


obsLogLikLms_i <- function(theta, model, data, P, sign = 1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  V <- P$V
  w <- P$w
  m <- nrow(V)
  px <- numeric(data$n)

  for (i in seq_len(m)) {
    z_i     <- V[i, ]
    mu_i    <- muLmsCpp(  model = modFilled, z = z_i)
    sigma_i <- sigmaLmsCpp(model = modFilled, z = z_i)

    dens_i <- numeric(data$n)

    offset <- 1L
    for (id in data$ids) { # go along patterns
      n.pattern <- data$n.pattern[[id]]

      colidx <- data$colidx[[id]]
      dataid <- data$data.split[[id]]

      end <- offset + n.pattern - 1L
      dens_i[offset:end] <- dmvn(data$data.split[[id]],
                                 mean = mu_i[colidx],
                                 sigma = sigma_i[colidx, colidx],
                                 log = FALSE)
      offset <- end + 1L
    }

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
  }, FUN.VALUE = numeric(data$n))
}


densitySingleLms <- function(z, modFilled, data) {
  mu <- muLmsCpp(model = modFilled, z = z)
  sigma <- sigmaLmsCpp(model = modFilled, z = z)

  density <- numeric(data$n)

  offset <- 1L
  for (id in data$ids) { # go along patterns
    n.pattern <- data$n.pattern[[id]]

    colidx <- data$colidx[[id]]
    dataid <- data$data.split[[id]]

    end <- offset + n.pattern - 1L
    density[offset:end] <- dmvn(data$data.split[[id]],
                                mean = mu[colidx],
                                sigma = sigma[colidx, colidx])
    offset <- end + 1L
  }

  density
}


densityLms <- function(z, modFilled, data) {
  if (is.null(dim(z))) z <- matrix(z, ncol = modFilled$quad$k)

  lapplyMatrix(seq_len(nrow(z)), FUN.VALUE = numeric(data$n), FUN = function(i) {
    densitySingleLms(z = z[i, , drop=FALSE], modFilled = modFilled, data = data)
  })
}


hessianAllLogLikLms <- function(theta, model, P, data, sign = -1,
                                FHESS, FOBJECTIVE, .relStep = .Machine$double.eps ^ (1/5)) {
  hasCovModel <- model$gradientStruct$hasCovModel

  if (hasCovModel) hessian <- \(...) complicatedHessianAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
  else             hessian <- \(...) simpleHessianAllLogLikLms(..., FHESS = FHESS, .relStep = .relStep)

  hessian(theta = theta, model = model, P = P, data = data, sign = sign)
}


compHessianAllLogLikLms <- function(theta, model, P, data, sign = -1,
                                    .relStep = .Machine$double.eps ^ (1/6),
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
  symmetric <- locations$symmetric

  HESS <- FHESS(modelR    = modelR,
                P         = P,
                block     = block,
                row       = row,
                col       = col,
                colidxR   = data$colidx0,
                n         = data$n.pattern,
                d         = data$d.pattern,
                npatterns = data$p,
                symmetric = symmetric,
                .relStep  = .relStep,
                ncores    = ThreadEnv$n.threads)

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


complicatedHessianAllLogLikLms <- function(theta, model, P, data, sign = -1, FOBJECTIVE,
                                           .relStep = .Machine$double.eps ^ (1/5)) {
  nlme::fdHess(pars = theta, fun = FOBJECTIVE, model = model, P = P, sign = sign,
               .relStep = .relStep, data = data)$Hessian
}


hessianObsLogLikLms <- function(theta, model, data, P, sign = -1,
                                .relStep = .Machine$double.eps ^ (1/5)) {
  if (isMultiGroupModelDA(model)) {
    GS <- model$gradientStruct
    if (isTRUE(GS$hasCovModel)) {
      return(complicatedHessianAllLogLikLms(theta = theta, model = model, P = P,
                                           data = data, sign = sign,
                                           FOBJECTIVE = obsLogLikLms,
                                           .relStep = .relStep))
    }

    submodels <- model$groupModels
    P_groups  <- P$groups
    stopifnot(length(P_groups) == length(submodels))

    k <- length(GS$param.full)
    raw_H <- matrix(0, k, k, dimnames = list(GS$param.full, GS$param.full))
    raw_grad <- numeric(k)

    for (g in seq_along(submodels)) {
      submodel <- submodels[[g]]
      theta_g  <- getThetaGroupDA(model, theta, g)
      P_g      <- P_groups[[g]]
      data_g   <- submodel$data

      res <- rawHessianObsGroupLms(submodel, theta_g, P_g, data_g, .relStep)
      idx <- GS$groupColumns[[g]]

      raw_H[idx, idx] <- raw_H[idx, idx] + res$Hessian
      raw_grad[idx]   <- raw_grad[idx] + res$gradient
    }

    Jacobian  <- GS$Jacobian
    Jacobian2 <- GS$Jacobian2

    if (length(GS$nlinDerivs)) {
      evalTheta  <- GS$evalTheta
      param.full <- colnames(Jacobian)
      param.part <- rownames(Jacobian)
      THETA      <- list2env(as.list(evalTheta(theta)))

      for (dep in names(GS$nlinDerivs)) {
        derivs1 <- GS$nlinDerivs[[dep]]
        derivs2 <- GS$nlinDerivs2[[dep]]

        for (indep in names(derivs1)) {
          deriv1 <- eval(expr = derivs1[[indep]], envir = THETA)
          deriv2 <- eval(expr = derivs2[[indep]], envir = THETA)

          match.full <- param.full == dep
          match.part <- param.part == indep
          if (any(match.full) && any(match.part)) {
            Jacobian[match.part, match.full]  <- deriv1
            Jacobian2[match.part, match.full] <- deriv2
          }
        }
      }
    }

    term1 <- Jacobian %*% raw_H %*% t(Jacobian)
    term2 <- diag(drop(Jacobian2 %*% raw_grad), nrow = nrow(Jacobian))

    return(sign * (term1 + term2))
  }

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR, n,
                    npatterns, ncores, ...) {
    hessObsLogLikLmsCpp(modelR = modelR, dataR = data$data.split, P = P,
                        block = block, row = row, col = col,
                        symmetric = symmetric, npatterns = npatterns,
                        colidxR = colidxR, n = n, relStep = .relStep,
                        minAbs = 0.0, ncores = ncores)
  }

  FOBJECTIVE <- function(theta, model, P, sign, data, ...) {
    obsLogLikLms(theta = theta, model = model, P = P, data = data, sign = sign)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      data = data, FHESS = FHESS, FOBJECTIVE = FOBJECTIVE,
                      .relStep = .relStep)
}


hessianCompLogLikLms <- function(theta, model, P, data, sign = -1,
                                 .relStep = .Machine$double.eps ^ (1/5)) {
  if (isMultiGroupModelDA(model)) {
    GS <- model$gradientStruct
    if (isTRUE(GS$hasCovModel)) {
      return(complicatedHessianAllLogLikLms(theta = theta, model = model, P = P,
                                           data = data, sign = sign,
                                           FOBJECTIVE = compLogLikLms,
                                           .relStep = .relStep))
    }

    submodels <- model$groupModels
    P_groups  <- P$groups
    stopifnot(length(P_groups) == length(submodels))

    k <- length(GS$param.full)
    raw_H <- matrix(0, k, k, dimnames = list(GS$param.full, GS$param.full))
    raw_grad <- numeric(k)

    for (g in seq_along(submodels)) {
      submodel <- submodels[[g]]
      theta_g  <- getThetaGroupDA(model, theta, g)
      P_g      <- P_groups[[g]]
      data_g   <- submodel$data

      res <- rawHessianGroupLms(submodel, theta_g, P_g, data_g, .relStep)
      idx <- GS$groupColumns[[g]]

      raw_H[idx, idx] <- raw_H[idx, idx] + res$Hessian
      raw_grad[idx]   <- raw_grad[idx] + res$gradient
    }

    Jacobian  <- GS$Jacobian
    Jacobian2 <- GS$Jacobian2

    if (length(GS$nlinDerivs)) {
      evalTheta  <- GS$evalTheta
      param.full <- colnames(Jacobian)
      param.part <- rownames(Jacobian)
      THETA      <- list2env(as.list(evalTheta(theta)))

      for (dep in names(GS$nlinDerivs)) {
        derivs1 <- GS$nlinDerivs[[dep]]
        derivs2 <- GS$nlinDerivs2[[dep]]

        for (indep in names(derivs1)) {
          deriv1 <- eval(expr = derivs1[[indep]], envir = THETA)
          deriv2 <- eval(expr = derivs2[[indep]], envir = THETA)

          match.full <- param.full == dep
          match.part <- param.part == indep
          if (any(match.full) && any(match.part)) {
            Jacobian[match.part, match.full]  <- deriv1
            Jacobian2[match.part, match.full] <- deriv2
          }
        }
      }
    }

    term1 <- Jacobian %*% raw_H %*% t(Jacobian)
    term2 <- diag(drop(Jacobian2 %*% raw_grad), nrow = nrow(Jacobian))

    return(sign * (term1 + term2))
  }

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR,
                    n, d, npatterns, ncores) {
    hessCompLogLikLmsCpp(modelR = modelR, P = P, block = block,
                         row = row, col = col, symmetric = symmetric,
                         colidxR = colidxR, n = n, d = d, relStep = .relStep,
                         npatterns = npatterns, minAbs = 0.0, ncores = ncores)
  }

  FOBJECTIVE <- function(theta, model, P, data, sign) {
    compLogLikLms(theta = theta, model = model, P = P, data = data, sign = sign)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, data = data, sign = sign,
                      FHESS = FHESS, FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
}


.logdensAllObsNode <- function(theta, model, data, z) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  # densitySingleLms returns densities for all rows at given z (not log)
  dens <- densitySingleLms(z = z, modFilled = modFilled, data = data)
  # guard against underflow/zeros
  log(pmax(dens, .Machine$double.xmin))
}


# per-node, per-observation complete-data score via finite difference
# Returns an n x p matrix S_j with row i = s_{ij}^T = grad_theta log p(y_i, z_j | theta)
.completeScoresNodeFD <- function(theta, model, data, z,
                                  epsilon = 1e-6, scheme = c("forward","central")) {
  scheme <- match.arg(scheme)
  p <- length(theta)
  n <- data$n
  S <- matrix(0.0, nrow = n, ncol = p)
  if (scheme == "forward") {
    f0 <- .logdensAllObsNode(theta, model, data, z)
    for (k in seq_len(p)) {
      th1 <- theta; th1[k] <- th1[k] + epsilon
      f1 <- .logdensAllObsNode(th1, model, data, z)
      S[, k] <- (f1 - f0) / epsilon
    }
  } else { # central
    for (k in seq_len(p)) {
      thp <- theta; thp[k] <- thp[k] + epsilon
      thm <- theta; thm[k] <- thm[k] - epsilon
      fp <- .logdensAllObsNode(thp, model, data, z)
      fm <- .logdensAllObsNode(thm, model, data, z)
      S[, k] <- (fp - fm) / (2*epsilon)
    }
  }
  dimnames(S) <- list(NULL, names(theta))
  S
}

# I_obs = I_com - I_mis using Louis' identity
observedInfoFromLouisLms <- function(model,
                                     theta,
                                     data,
                                     P = NULL,
                                     recompute.P = is.null(P),
                                     adaptive.quad.tol = 1e-12,
                                     fd.epsilon = 1e-6,
                                     fd.scheme = c("forward","central"),
                                     symmetrize = TRUE,
                                     jitter = 0.0,
                                     ...) {
  fd.scheme <- match.arg(fd.scheme)

  if (isMultiGroupModelDA(model)) {
    if (recompute.P || is.null(P)) {
      P <- estepLms(model = model, theta = theta, data = data,
                    lastQuad = NULL, recalcQuad = FALSE,
                    adaptive.quad.tol = adaptive.quad.tol, ...)
    }

    submodels <- model$groupModels
    P_groups <- P$groups
    stopifnot(length(P_groups) == length(submodels))

    p <- length(theta)
    lbl <- names(theta)
    Iobs_total <- matrix(0, p, p, dimnames = list(lbl, lbl))
    Icom_total <- matrix(0, p, p, dimnames = list(lbl, lbl))
    Imis_total <- matrix(0, p, p, dimnames = list(lbl, lbl))

    for (g in seq_along(submodels)) {
      submodel <- submodels[[g]]
      theta_g <- getThetaGroupDA(model, theta, g)
      data_g <- submodel$data
      P_g <- P_groups[[g]]

      res_g <- observedInfoFromLouisLms(
        model = submodel,
        theta = theta_g,
        data = data_g,
        P = P_g,
        recompute.P = FALSE,
        adaptive.quad.tol = adaptive.quad.tol,
        fd.epsilon = fd.epsilon,
        fd.scheme = fd.scheme,
        symmetrize = FALSE,
        jitter = 0.0,
        ...
      )

      Iobs_total <- embedGroupMatrixDA(Iobs_total, model, g, res_g$I.obs)
      Icom_total <- embedGroupMatrixDA(Icom_total, model, g, res_g$I.com)
      Imis_total <- embedGroupMatrixDA(Imis_total, model, g, res_g$I.mis)
    }

    if (symmetrize) {
      sym <- function(A) 0.5 * (A + t(A))
      Iobs_total <- sym(Iobs_total)
      Icom_total <- sym(Icom_total)
      Imis_total <- sym(Imis_total)
    }

    if (jitter > 0) {
      diag_add <- diag(jitter, nrow = p)
      Iobs_total <- Iobs_total + diag_add
      Icom_total <- Icom_total + diag_add
      Imis_total <- Imis_total + diag_add
    }

    return(list(I.obs = Iobs_total, I.com = Icom_total, I.mis = Imis_total, P = P))
  }

  # E-step (if needed)
  if (recompute.P) {
    P <- estepLms(model = model, theta = theta, data = data,
                  lastQuad = NULL, recalcQuad = FALSE,
                  adaptive.quad.tol = adaptive.quad.tol, ...)
  }

  # labels and sizes
  p <- length(theta)
  n <- data$n
  J <- length(P$w)             # number of quadrature nodes
  lbl <- names(theta)

  # Complete Information
  Icom <- hessianCompLogLikLms(theta = theta, model = model, P = P, data = data,
                               sign = -1)

  # Imis = sum_i ( E[s_ij s_ij^T] - E[s_ij]E[s_ij]^T ), expectations over j with weights r_ij = P$P[i,j]
  # Accumulate efficiently:
  #   total_M = sum_j S_j^T diag(r_.j) S_j, where S_j is n x p score matrix at node j
  #   Sbar (n x p) = sum_j r_.j * S_j  -> then sum_i sbar_i sbar_i^T = t(Sbar) %*% Sbar
  total_M <- matrix(0.0, p, p, dimnames = list(lbl, lbl))
  Sbar    <- matrix(0.0, n, p)   # will hold per-observation E[s_ij]; no need to keep per-i outer products

  for (j in seq_len(J)) {
    z_j <- P$V[j, , drop = FALSE]       # node
    S_j <- .completeScoresNodeFD(theta, model, data, z_j,
                                       epsilon = fd.epsilon, scheme = fd.scheme)  # n x p
    r_j <- P$P[, j]                     # length-n weights (posterior r_ij)

    # total_M += t( S_j * sqrt(r_j) ) %*% ( S_j * sqrt(r_j) )
    # (avoid forming diag(r_j) explicitly)
    Rhalf <- sqrt(pmax(r_j, 0))
    X <- S_j * Rhalf       # row-wise scaling
    total_M <- total_M + crossprod(X)   # t(X) %*% X  -> p x p

    # Sbar += r_j * S_j   (row-wise scaling)
    Sbar <- Sbar + (S_j * r_j)
  }

  sbar_outer <- crossprod(Sbar) # sum_i sbar_i sbar_i^T  -> p x p
  Imis <- total_M - sbar_outer

  # Louis Identity
  Iobs <- Icom - Imis

  # hygiene
  if (symmetrize) {
    sym <- function(A) 0.5*(A + t(A))
    Icom <- sym(Icom); Imis <- sym(Imis); Iobs <- sym(Iobs)
  }

  if (jitter > 0) {
    Icom <- Icom + diag(jitter, p)
    Imis <- Imis + diag(jitter, p)
    Iobs <- Iobs + diag(jitter, p)
  }

  list(I.obs = Iobs, I.com = Icom, I.mis = Imis, P = P)
}
