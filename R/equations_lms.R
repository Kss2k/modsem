estepLms <- function(model, theta, lastQuad = NULL, recalcQuad = FALSE,
                     adaptive.quad.tol = 1e-12, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  P <- list(P_GROUPS = vector("list", length = model$info$n.groups),
            quad = NULL, obsLL = NULL)

  for (g in seq_len(model$info$n.groups)) {
    lastQuad.g <- if (!is.null(lastQuad)) lastQuad[[g]] else NULL
    submodel   <- modFilled$models[[g]]

    P$P_GROUPS[[g]] <- estepLmsGroup(
      submodel          = submodel,
      lastQuad          = lastQuad.g,
      recalcQuad        = recalcQuad,
      adaptive.quad.tol = adaptive.quad.tol,
      ...
    )
  }

  P$quad  <- lapply(P$P_GROUPS, FUN = \(P) P$quad)
  P$obsLL <- sum(vapply(P$P_GROUPS, FUN.VALUE = numeric(1L), FUN = \(P) P$obsLL))

  P
}


estepLmsGroup <- function(submodel, lastQuad = NULL, recalcQuad = FALSE,
                          adaptive.quad.tol = 1e-12, ...) {
  data             <- submodel$data
  sampling.weights <- data$weights

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
      estep.fixed <- estepLmsGroup(
        submodel = submodel,
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
    P    <- W * densityLms(V, modFilled = submodel, data = data)
  }

  density        <- rowSums(P)
  observedLogLik <- sum(log(density))
  P              <- P / density

  # The sampling weights are already incorporated in `densityLms()`, so
  # `observedLogLik` is correct. But the P/density correction is not.
  # Here we correct P/density (if needed).
  if (!is.null(sampling.weights))
    P <- sampling.weights * P

  wMeans <- vector("list", length = length(w))
  wCovs  <- vector("list", length = length(w))
  tGamma <- vector("list", length = length(w))

  for (i in seq_along(w)) {
    p <- P[, i]
    offset <- 1L

    wMeans[[i]] <- vector("list", length = length(data$ids))
    wCovs[[i]]  <- vector("list", length = length(data$ids))
    tGamma[[i]] <- numeric(length = length(data$ids))

    for (j in data$ids) {
      n.pattern <- data$n.pattern[[j]]
      end       <- offset + n.pattern - 1L

      data.id <- data$data.split[[j]]
      colidx  <- data$colidx[[j]]

      pj   <- p[offset:end]
      wm   <- colSums(data.id * pj) / sum(pj)
      X    <- data.id - matrix(wm, nrow=nrow(data.id), ncol=ncol(data.id), byrow=TRUE)
      wcov <- t(X) %*% (X * pj)

      wMeans[[i]][[j]] <- wm
      wCovs[[i]][[j]]  <- wcov
      tGamma[[i]][[j]] <- sum(pj)

      offset <- end + 1L
    }
  }

  # Create a vector for sampling weights, needed in some C++ code
  if (!is.null(sampling.weights)) sampling.weights.vec <- sampling.weights
  else                            sampling.weights.vec <- rep(1, NROW(P))

  list(P = P, mean = wMeans, cov = wCovs, tgamma = tGamma, V = V, w = w,
       obsLL = observedLogLik, quad = quad, sampling.weights = sampling.weights.vec)
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
    compLogLikLms(theta = theta, model = model, P = P, sign = -1,
                  epsilon = epsilon)
  }

  if (optimizer == "nlminb") {
    if (is.null(control$iter.max)) control$iter.max <- max.step
    est <- stats::nlminb(start = theta, objective = objective,
                         gradient = gradient,
                         upper = model$params$bounds$upper,
                         lower = model$params$bounds$lower, control = control,
                         ...) |> suppressWarnings()

  } else if (optimizer == "L-BFGS-B") {
    if (is.null(control$maxit)) control$maxit <- max.step
    est <- stats::optim(par = theta, fn = objective, gr = gradient,
                        method = optim.method, control = control,
                        lower = model$params$bounds$lower,
                        upper = model$params$bounds$upper, ...)

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

    ll <- 0
    for (g in seq_len(model$info$n.groups)) {
      submodel <- modFilled$models[[g]]
      data <- submodel$data

      ll <- ll + completeLogLikLmsCpp(
        modelR=submodel, P=P$P_GROUPS[[g]], quad=P$quad[[g]],
        colidxR = data$colidx0, n = data$n.pattern,
        d = data$d.pattern, npatterns = data$p
      )
    }

    sign * ll

  }, error = \(e) NA)
}


compLogLikLmsGroup <- function(submodel, P, sign = -1, ...) {
  tryCatch({
    data <- submodel$data

    ll <- completeLogLikLmsCpp(
      modelR = submodel, P = P, quad = P$quad,
      colidxR = data$colidx0, n = data$n.pattern,
      d = data$d.pattern, npatterns = data$p
    )

    sign * ll

  }, error = \(e) NA)
}



gradientCompLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6) {
  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, FGRAD = gradLogLikLmsCpp, FOBJECTIVE = compLogLikLmsGroup)
}


gradientAllLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6,
                                 FGRAD, FOBJECTIVE) {
  hasCovModel <- model$params$gradientStruct$hasCovModel

  if (hasCovModel) gradient <- \(...) complicatedGradientAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE)
  else             gradient <- \(...) simpleGradientAllLogLikLms(..., FGRAD = FGRAD)

  c(gradient(theta = theta, model = model, P = P, sign = sign, epsilon = epsilon))
}


complicatedGradientAllLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-4, FOBJECTIVE, sum = TRUE, ...) {
  params <- model$params

  SELECT_THETA_LAB  <- params$SELECT_THETA_LAB
  SELECT_THETA_COV  <- params$SELECT_THETA_COV
  SELECT_THETA_MAIN <- params$SELECT_THETA_MAIN

  N <- vapply(model$models, FUN.VALUE = numeric(1L), FUN = \(sub) sub$data$n)
  N.start <- c(1, cumsum(N)[-length(N)] + 1L)
  N.end   <- cumsum(N)

  if (sum) n <- 1L
  else     n <- sum(N)

  k  <- length(theta)

  grad <- matrix(0, nrow = n, ncol = k, dimnames = list(NULL, names(theta)))

  FOBJECTIVE_GROUP <- function(theta, g) {
    modFilled <- fillModel(theta = theta, model = model, method = "lms")
    FOBJECTIVE(submodel = modFilled$models[[g]], P = P$P_GROUPS[[g]], sign = sign, ...)
  }

  for (g in seq_len(model$info$n.groups)) {
    f0 <- FOBJECTIVE_GROUP(theta = theta, g = g)

    if (sum) J <- 1L
    else     J <- seq(N.start[[g]], N.end[[g]], by = 1L)

    indices <- c(
      SELECT_THETA_LAB[[g]],
      SELECT_THETA_COV[[g]],
      SELECT_THETA_MAIN[[g]]
    )

    for (i in indices) {
      theta_i <- theta
      theta_i[i] <- theta_i[i] + epsilon

      fi <- FOBJECTIVE_GROUP(theta_i, g = g)
      grad[J, i] <- grad[J, i] + (fi - f0) / epsilon
    }
  }

  if (n == 1L) as.vector(grad) else grad
}


simpleGradientAllLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6, FGRAD, N.FGRAD = 1L) {
  # simple gradient which should work if constraints are well-behaved functions
  # which can be derivated by Deriv::Deriv, and there is no covModel
  modelR     <- fillModel(model=model, theta=theta, method="lms")
  locations  <- model$params$gradientStruct$locations
  Jacobian   <- model$params$gradientStruct$Jacobian
  nlinDerivs <- model$params$gradientStruct$nlinDerivs

  # grad <- stats::setNames(numeric(NROW(locations)), nm = locations$param)
  grad <- matrix(0, nrow = NROW(locations), ncol = N.FGRAD,
                 dimnames = list(locations$param, NULL))

  for (g in seq_len(modelR$info$n.groups)) {
    submodelR   <- modelR$models[[g]]
    locations.g <- locations[locations$group == g, , drop = FALSE]

    data      <- submodelR$data
    block     <- locations.g$block
    row       <- locations.g$row
    col       <- locations.g$col
    param     <- locations.g$param
    symmetric <- locations.g$symmetric

    grad.g <- FGRAD(modelR    = submodelR,
                    P         = P$P_GROUPS[[g]],
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

    grad[locations.g$param, ] <- grad.g
  }

  if (length(nlinDerivs)) {
    evalTheta  <- model$params$gradientStruct$evalTheta
    param.full <- stringr::str_split_i(colnames(Jacobian), pattern = "#", i = 1L) # non-unique
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


obsLogLikLms <- function(theta, model, P, sign = 1, ...) {
  modFilled  <- fillModel(model = model, theta = theta, method = "lms")

  ll <- 0
  for (g in seq_len(modFilled$info$n.groups)) {
    submodel <- modFilled$models[[g]]
    data.g   <- submodel$data
    P.g      <- P$P_GROUPS[[g]]

    ll.g <- observedLogLikLmsCpp(submodel,
                                 dataR = data.g$data.split,
                                 colidxR = data.g$colidx0,
                                 P = P.g,
                                 n = data.g$n.pattern,
                                 npatterns = data.g$p,
                                 ncores = ThreadEnv$n.threads)
    ll <- ll + ll.g
  }

  sign * ll
}


obsLogLikLmsGroup <- function(submodel, P, sign = -1, ...) {
  tryCatch({
    data.g <- submodel$data

    ll <- observedLogLikLmsCpp(submodel,
                               dataR = data.g$data.split,
                               colidxR = data.g$colidx0,
                               P = P,
                               n = data.g$n.pattern,
                               npatterns = data.g$p,
                               ncores = ThreadEnv$n.threads)
    sign * ll

  }, error = \(e) NA)
}


gradientObsLogLikLms <- function(theta, model, P, sign = 1, epsilon = 1e-6) {
  FGRAD <- function(modelR, P, block, row, col, symmetric, colidxR, npatterns,
                    eps, ncores, n, ...) {
    dataR <- modelR$data
    gradObsLogLikLmsCpp(modelR = modelR, dataR = dataR$data.split, P = P,
                        block = block, row = row, col = col,
                        symmetric = symmetric, colidxR = colidxR,
                        n = n, npatterns = npatterns, eps = eps,
                        ncores = ncores)
  }

  FOBJECTIVE <- function(theta, model, P, colidxR, npatterns, sign, ...) {
    obsLogLikLmsGroup(theta = theta, model = model, P = P, sign = sign)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, FGRAD = FGRAD,
                       FOBJECTIVE = FOBJECTIVE)
}


obsLogLikLmsGroup_i <- function(submodel, P, sign = 1) {
  data <- submodel$data
  V  <- P$V
  w  <- P$w
  m  <- nrow(V)
  px <- numeric(data$n)

  for (i in seq_len(m)) {
    z_i     <- V[i, ]
    mu_i    <- muLmsCpp(model = submodel, z = z_i)
    sigma_i <- sigmaLmsCpp(model = submodel, z = z_i)

    dens_i <- numeric(data$n)
    offset <- 1L

    for (id in data$ids) {
      n.pattern <- data$n.pattern[[id]]
      colidx    <- data$colidx[[id]]

      end <- offset + n.pattern - 1L
      dens_i[offset:end] <- dmvn(data$data.split[[id]],
                                 mean = mu_i[colidx],
                                 sigma = sigma_i[colidx, colidx],
                                 log = FALSE)
      offset <- end + 1L
    }

    px <- px + w[i] * dens_i
  }

  logdens <- log(px)

  sampling.weights <- data$weights
  if (!is.null(sampling.weights))
    logdens <- sampling.weights * logdens

  sign * logdens
}


gradientObsLogLikLms_i <- function(theta, model, P, sign = 1, epsilon = 1e-4) {
  complicatedGradientAllLogLikLms(theta = theta, model = model, P = P,
                                  sign = sign, epsilon = epsilon, sum = FALSE,
                                  FOBJECTIVE = obsLogLikLmsGroup_i)
}


densitySingleLms <- function(z, modFilled, data) {
  mu <- muLmsCpp(model = modFilled, z = z)
  sigma <- sigmaLmsCpp(model = modFilled, z = z)

  cov.lv <- covLvLmsCpp(model = modFilled, z = z)
  mu.lv  <- muLvLmsCpp(model = modFilled, z = z)

  k <- length(z)
  diag(cov.lv)[seq_len(k)] <- 1

  A.lv <- t(chol(cov.lv))

  inds         <- modFilled$info$allInds
  quad.ordered <- modFilled$quad.ordered
  V.o          <- quad.ordered$n
  w.o          <- quad.ordered$weights

  Theta        <- modFilled$matrices$thetaDelta
  Thresholds   <- modFilled$matrices$thresholds
  Lambda       <- modFilled$matrices$lambdaX

  density <- numeric(data$n)
  offset <- 1L

  for (id in data$ids) { # go along patterns
    n.pattern <- data$n.pattern[[id]]

    colidx      <- data$colidx[[id]]
    is.ordered  <- modFilled$info$is.ordered[colidx]
    colidx.cont <- colidx[!is.ordered]
    colidx.ord  <- colidx[is.ordered]
    dataid      <- data$data.split[[id]]

    end <- offset + n.pattern - 1L
    density[offset:end] <- dmvn(data$data.split[[id]][, colidx.cont],
                                mean = mu[colidx.cont],
                                sigma = sigma[colidx.cont, colidx.cont])

    for (q in seq_len(NROW(V.o))) {
      v <- mu.lv + t(A.lv) %*% V.o[q, ]
      response <- Lambda %*% v

      for (j in colidx.ord) {
        t <- data$data.split[[id]][, j]
        lower <- Thresholds[j, t]
        upper <- Thresholds[j, t+1L]
        sd.j  <- sqrt(Theta[j, j])
        v.j   <- response[j, ]

        density[offset:end] <- density[offset:end] * (
          pnorm(upper, mean = v.j, sd = sd.j) - pnorm(lower, mean = v.j, sd = sd.j)
        )
      }
    }

    density[offset:end] 
    offset <- end + 1L
  }

  sampling.weights <- data$weights
  if (!is.null(sampling.weights))
    density <- exp(log(density) * sampling.weights) # can this be simplified?

  density
}


densityLms <- function(z, modFilled, data) {
  if (is.null(dim(z))) z <- matrix(z, ncol = modFilled$quad$k)

  lapplyMatrix(seq_len(nrow(z)), FUN.VALUE = numeric(data$n), FUN = function(i) {
    densitySingleLms(z = z[i, , drop=FALSE], modFilled = modFilled, data = data)
  })
}


hessianAllLogLikLms <- function(theta, model, P, sign = -1,
                                FHESS, FOBJECTIVE, .relStep = .Machine$double.eps ^ (1/5)) {
  hasCovModel <- model$params$gradientStruct$hasCovModel

  if (hasCovModel) hessian <- \(...) complicatedHessianAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
  else             hessian <- \(...) simpleHessianAllLogLikLms(..., FHESS = FHESS, .relStep = .relStep)

  hessian(theta = theta, model = model, P = P, sign = sign)
}


complicatedHessianAllLogLikLms <- function(theta, model, P, sign = -1,
                                          .relStep = .Machine$double.eps ^ (1/6),
                                          FOBJECTIVE) {
  params <- model$params

  SELECT_THETA_LAB  <- params$SELECT_THETA_LAB
  SELECT_THETA_COV  <- params$SELECT_THETA_COV
  SELECT_THETA_MAIN <- params$SELECT_THETA_MAIN

  k <- length(theta)
  H <- matrix(0, nrow = k, ncol = k, dimnames = list(names(theta), names(theta)))

  for (g in seq_len(model$info$n.groups)) {
    indices <- c(
      SELECT_THETA_LAB[[g]],
      SELECT_THETA_COV[[g]],
      SELECT_THETA_MAIN[[g]]
    )

    FOBJECTIVE_GROUP <- function(theta.g) {
      theta[indices] <- theta.g # local copy of theta
      modFilled <- fillModel(theta = theta, model = model, method = "lms")
      FOBJECTIVE(submodel = modFilled$models[[g]], sign = sign, P = P$P_GROUPS[[g]])
    }

    theta.g <- theta[indices]
    Hg <- fdHESS(pars = theta.g, fun = FOBJECTIVE_GROUP, .relStep = .Machine$double.eps^(1/5))

    H[indices, indices] <- H[indices, indices] + Hg
  }

  H
}


simpleHessianAllLogLikLms <- function(theta, model, P, sign = -1,
                                      .relStep = .Machine$double.eps ^ (1/5),
                                      FHESS) {
  # simple Hessian which should work if constraints are well-behaved functions
  # which can be derivated by Deriv::Deriv, and there is no covModel
  modelR      <- fillModel(model = model, theta = theta, method = "lms")
  locations   <- model$params$gradientStruct$locations
  Jacobian    <- model$params$gradientStruct$Jacobian
  Jacobian2   <- model$params$gradientStruct$Jacobian2
  nlinDerivs  <- model$params$gradientStruct$nlinDerivs
  nlinDerivs2 <- model$params$gradientStruct$nlinDerivs2

  n.loc <- NROW(locations)
  H <- matrix(0.0, nrow = n.loc, ncol = n.loc,
              dimnames = list(locations$param, locations$param))
  grad <- numeric(n.loc)
  names(grad) <- locations$param

  for (g in seq_len(modelR$info$n.groups)) {
    locations.g <- locations[locations$group == g, , drop = FALSE]
    if (!NROW(locations.g)) next

    submodelR <- modelR$models[[g]]
    data.g    <- submodelR$data

    HESS.g <- FHESS(modelR    = submodelR,
                    P         = P$P_GROUPS[[g]],
                    block     = locations.g$block,
                    row       = locations.g$row,
                    col       = locations.g$col,
                    colidxR   = data.g$colidx0,
                    n         = data.g$n.pattern,
                    d         = data.g$d.pattern,
                    npatterns = data.g$p,
                    symmetric = locations.g$symmetric,
                    .relStep  = .relStep,
                    ncores    = ThreadEnv$n.threads)

    H.g    <- HESS.g$Hessian
    grad.g <- HESS.g$gradient

    dimnames(H.g) <- list(locations.g$param, locations.g$param)
    names(grad.g) <- locations.g$param

    H[locations.g$param, locations.g$param] <- H[locations.g$param, locations.g$param] + H.g
    grad[locations.g$param] <- grad[locations.g$param] + grad.g
  }

  if (length(nlinDerivs)) {
    evalTheta  <- model$params$gradientStruct$evalTheta
    param.full <- stringr::str_split_i(colnames(Jacobian), pattern = "#", i = 1L) # non-unique
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


hessianObsLogLikLms <- function(theta, model, P, sign = -1,
                                data = NULL,
                                .relStep = .Machine$double.eps ^ (1/5)) {

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR, n,
                    npatterns, ncores, ...) {
    dataR <- modelR$data
    hessObsLogLikLmsCpp(modelR = modelR, dataR = dataR$data.split, P = P,
                        block = block, row = row, col = col,
                        symmetric = symmetric, npatterns = npatterns,
                        colidxR = colidxR, n = n, relStep = .relStep,
                        minAbs = 0.0, ncores = ncores)
  }

  FOBJECTIVE <- function(theta, submodel, P, sign, ...) {
    obsLogLikLmsGroup(theta = theta, submodel = model, P = P, sign = sign)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      FHESS = FHESS, FOBJECTIVE = FOBJECTIVE,
                      .relStep = .relStep)
}


hessianCompLogLikLms <- function(theta, model, P, sign = -1,
                                 .relStep = .Machine$double.eps ^ (1/5)) {

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR,
                    n, d, npatterns, ncores) {
    hessCompLogLikLmsCpp(modelR = modelR, P = P, block = block,
                         row = row, col = col, symmetric = symmetric,
                         colidxR = colidxR, n = n, d = d, relStep = .relStep,
                         npatterns = npatterns, minAbs = 0.0, ncores = ncores)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      FHESS = FHESS, FOBJECTIVE = compLogLikLmsGroup,
                      .relStep = .relStep)
}


.logdensAllObsNode <- function(theta, model, data, z, group = NULL) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  submodel <- if (!is.null(modFilled$models)) {
    if (is.null(group)) modFilled$models[[1L]] else modFilled$models[[group]]
  } else modFilled
  # densitySingleLms returns densities for all rows at given z (not log)
  dens <- densitySingleLms(z = z, modFilled = submodel, data = data)
  # guard against underflow/zeros
  log(pmax(dens, .Machine$double.xmin))
}


.activeThetaIndicesLms <- function(model, group, p) {
  select_lab  <- model$params$SELECT_THETA_LAB[[group]]
  select_cov  <- model$params$SELECT_THETA_COV[[group]]
  select_main <- model$params$SELECT_THETA_MAIN[[group]]

  if (is.null(select_lab) && is.null(select_cov) && is.null(select_main)) {
    return(seq_len(p))
  }

  idx <- c(select_lab,
           select_cov,
           select_main)

  idx <- sort(unique(as.integer(idx)))
  idx <- idx[idx >= 1L & idx <= p]

  if (length(idx)) idx else seq_len(p)
}


# per-node, per-observation complete-data score via finite difference
# Returns an n x p matrix S_j with row i = s_{ij}^T = grad_theta log p(y_i, z_j | theta)
.completeScoresNodeFD <- function(theta, model, data, z,
                                  epsilon = 1e-6, scheme = c("forward","central"),
                                  group = NULL, active = NULL) {
  scheme <- match.arg(scheme)
  p <- length(theta)

  if (is.null(active)) {
    active <- seq_len(p)
  } else {
    active <- unique(as.integer(active))
    active <- active[active >= 1L & active <= p]
  }

  n <- data$n
  n_active <- length(active)

  if (!n_active) {
    return(matrix(0.0, nrow = n, ncol = 0L))
  }

  col_names <- names(theta)
  if (!length(col_names)) col_names <- rep("", p)

  S <- matrix(0.0, nrow = n, ncol = n_active,
              dimnames = list(NULL, col_names[active]))

  if (scheme == "forward") {
    f0 <- .logdensAllObsNode(theta, model, data, z, group = group)
    for (pos in seq_len(n_active)) {
      k <- active[pos]
      th1 <- theta; th1[k] <- th1[k] + epsilon
      f1 <- .logdensAllObsNode(th1, model, data, z, group = group)
      S[, pos] <- (f1 - f0) / epsilon
    }
  } else { # central
    for (pos in seq_len(n_active)) {
      k <- active[pos]
      thp <- theta; thp[k] <- thp[k] + epsilon
      thm <- theta; thm[k] <- thm[k] - epsilon
      fp <- .logdensAllObsNode(thp, model, data, z, group = group)
      fm <- .logdensAllObsNode(thm, model, data, z, group = group)
      S[, pos] <- (fp - fm) / (2 * epsilon)
    }
  }
  S
}

# I_obs = I_com - I_mis using Louis' identity
observedInfoFromLouisLms <- function(model,
                                     theta,
                                     P = NULL,
                                     recompute.P = is.null(P),
                                     adaptive.quad.tol = 1e-12,
                                     fd.epsilon = 1e-6,
                                     fd.scheme = c("forward","central"),
                                     symmetrize = TRUE,
                                     jitter = 0.0,
                                     ...) {
  fd.scheme <- match.arg(fd.scheme)

  # E-step (if needed)
  if (recompute.P) {
    P <- estepLms(model = model, theta = theta,
                  lastQuad = NULL, recalcQuad = FALSE,
                  adaptive.quad.tol = adaptive.quad.tol, ...)
  }

  p   <- length(theta)
  lbl <- names(theta)

  n_total   <- sum(vapply(model$models, function(sub) sub$data$n, numeric(1L)))

  Icom <- hessianCompLogLikLms(theta = theta, model = model, P = P, sign = -1)

  total_M <- matrix(0.0, p, p, dimnames = list(lbl, lbl))
  Sbar    <- matrix(0.0, n_total, p)

  row_offset <- 0L
  for (g in seq_len(model$info$n.groups)) {
    submodel <- model$models[[g]]
    data.g   <- submodel$data
    P.g      <- P$P_GROUPS[[g]]
    rows     <- seq_len(data.g$n) + row_offset

    active_idx <- .activeThetaIndicesLms(model, g, p)
    Jg <- length(P.g$w)
    for (j in seq_len(Jg)) {
      z_j <- P.g$V[j, , drop = FALSE]
      S_j <- .completeScoresNodeFD(theta, model, data.g, z_j,
                                   epsilon = fd.epsilon, scheme = fd.scheme,
                                   group = if (model$info$n.groups > 1L) g else NULL,
                                   active = active_idx)
      r_j <- P.g$P[, j]

      Rhalf <- sqrt(pmax(r_j, 0))
      if (NROW(S_j) && NCOL(S_j)) {
        X <- S_j * Rhalf
        block <- total_M[active_idx, active_idx, drop = FALSE]
        block <- block + crossprod(X)
        total_M[active_idx, active_idx] <- block

        Sbar_block <- Sbar[rows, active_idx, drop = FALSE]
        Sbar_block <- Sbar_block + (S_j * r_j)
        Sbar[rows, active_idx] <- Sbar_block
      }
    }

    row_offset <- row_offset + data.g$n
  }

  sbar_outer <- crossprod(Sbar)
  Imis <- total_M - sbar_outer
  Iobs <- Icom - Imis

  if (symmetrize) {
    sym <- function(A) 0.5 * (A + t(A))
    Icom <- sym(Icom); Imis <- sym(Imis); Iobs <- sym(Iobs)
  }

  if (jitter > 0) {
    Icom <- Icom + diag(jitter, p)
    Imis <- Imis + diag(jitter, p)
    Iobs <- Iobs + diag(jitter, p)
  }

  list(I.obs = Iobs, I.com = Icom, I.mis = Imis, P = P)
}
