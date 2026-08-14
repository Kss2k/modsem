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

  P$quad.err <- sum(vapply(
    X = P$P_GROUPS,
    FUN.VALUE = numeric(1L),
    FUN = function(P) {
      err <- P$quad$error.abs
      if (!length(err) || !is.finite(err[[1L]])) return(0.0)
      abs(err[[1L]])
    }
  ))

  P$quad  <- lapply(
    X = P$P_GROUPS,
    FUN = \(P) P$quad
  )

  P$obsLL <- sum(vapply(
    X = P$P_GROUPS,
    FUN.VALUE = numeric(1L),
    FUN = \(P) P$obsLL
  ))

  P
}


estepLmsGroup <- function(submodel, lastQuad = NULL, recalcQuad = FALSE,
                          adaptive.quad.tol = 1e-12, ...) {
  data             <- submodel$data
  sampling.weights <- data$weights
  adaptiveMode     <- submodel$quad$adaptive

  if (adaptiveMode %in% c("aghq", "fixed-n"))
    return(estepLmsGroupAghq(submodel = submodel, lastQuad = lastQuad,
                             adapt = adaptiveMode == "aghq"))

  if (adaptiveMode == "quasi" && (recalcQuad || is.null(lastQuad))) {
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
          fun = densityLms, modFilled = submodel,
          data = data, a = a, b = b, m = m,
          k = k, m.ceil = m.ceil, tol = adaptive.quad.tol,
        )
      }, error = function(e) {
        mod_msg_warn_immediate(
          paste0("Calculation of adaptive quadrature failed!\n", e),
          .newline = TRUE
        )
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

    V <- quad$n
    w <- quad$w
    P <- sweep(quad$F, MARGIN = 2, STATS = w, FUN = "*")

  } else {
    quad <- if (adaptiveMode == "quasi") lastQuad else submodel$quad
    V    <- quad$n
    w    <- quad$w
    densityVals <- densityLms(V, modFilled = submodel, data = data)
    P    <- sweep(densityVals, MARGIN = 2, STATS = w, FUN = "*")
  }

  density <- rowSums(P)
  P <- P / density
  if (is.null(sampling.weights)) observedLogLik <- sum(log(density))
  else observedLogLik <- sum(sampling.weights * log(density))


  # Sampling weights multiply each observation's log likelihood and posterior
  # contribution; they must not alter the posterior node probabilities.
  if (!is.null(sampling.weights))
    P <- sampling.weights * P

  stats <- estepSuffStatLmsCpp(P = P, dataR = data$data.split,
                               n = data$n.pattern, npatterns = data$p,
                               ncores = ThreadEnv$n.threads)

  # Create a vector for sampling weights, needed in some C++ code
  if (!is.null(sampling.weights)) sampling.weights.vec <- sampling.weights
  else                            sampling.weights.vec <- rep(1, NROW(P))

  list(P = P, mean = stats$mean, cov = stats$cov, tgamma = stats$tgamma,
       V = V, w = w, obsLL = observedLogLik, quad = quad,
       sampling.weights = sampling.weights.vec)
}


# E-step for `adaptive.quad %in% c("aghq", "fixed-n")`: every row gets its own
# adapted (or, for "fixed-n", deliberately non-adapted) node/weight grid.
# Unlike estepLmsGroup(), there is no shared-node sufficient-statistic
# aggregation (estepSuffStatLmsCpp) here -- see the note atop the analogous
# C++ section in equations_lms.cpp for why that's unnecessary for this path.
# The per-row modes are warm-started from `lastQuad$Zstar` (previous
# iteration), and always recomputed (regardless of `recalcQuad`/
# `adaptive.frequency`, which only govern the "quasi" path): AGHQ's benefit
# comes from staying tightly adapted to the current theta.
estepLmsGroupAghq <- function(submodel, lastQuad = NULL, adapt = TRUE, ...) {
  data             <- submodel$data
  sampling.weights <- data$weights
  k <- submodel$quad$k
  m <- submodel$quad$m

  base <- quadrature(m = m, k = k, adaptive = "fixed")

  Z0 <- if (!is.null(lastQuad$Zstar) &&
            NROW(lastQuad$Zstar) == data$n && NCOL(lastQuad$Zstar) == k) {
    lastQuad$Zstar
  } else {
    matrix(0, nrow = data$n, ncol = k)
  }

  grids <- aghqRowGridsLmsCpp(
    modelR = submodel, dataR = data$data.split, colidxR = data$colidx0,
    n = data$n.pattern, baseN = base$n, baseW = base$w, Z0 = Z0,
    adapt = adapt, npatterns = data$p, ncores = ThreadEnv$n.threads
  )

  est <- aghqEstepLmsCpp(
    modelR = submodel, Z = grids$Z, W = grids$W, dataR = data$data.split,
    colidxR = data$colidx0,
    samplingWeights = if (!is.null(sampling.weights)) sampling.weights else numeric(0L),
    n = data$n.pattern, npatterns = data$p, ncores = ThreadEnv$n.threads
  )

  mod_warnif_immediate(
    mean(grids$ok) < 0.9 && adapt,
    sprintf(paste0("AGHQ mode-finding fell back to the non-adapted node for ",
                   "%.1f%% of observations in this group (usually harmless, ",
                   "but a high fraction can indicate a poorly identified model)."),
            100 * (1 - mean(grids$ok))),
    .newline = TRUE
  )

  if (!is.null(sampling.weights)) sampling.weights.vec <- sampling.weights
  else                            sampling.weights.vec <- rep(1, data$n)

  quad <- list(adaptive = submodel$quad$adaptive, k = k, m = m,
              Zstar = grids$Zstar, ok = grids$ok)

  # Z/W are kept (not just the posterior Gamma) so that obsLogLikLmsGroup()/
  # gradientObsLogLikLms() can recompute the posterior at a *different* trial
  # theta using these SAME (E-step-time) adapted nodes -- the AGHQ analogue of
  # how observedGradientReverseFromModel() reuses the shared-node V/w for the
  # observed-data gradient via Fisher's identity.
  list(P = est$Gamma, Z = grids$Z, W = grids$W,
       mean = NULL, cov = NULL, tgamma = NULL, V = NULL, w = NULL,
       obsLL = est$obsLL, quad = quad,
       sampling.weights = sampling.weights.vec)
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
  thetaInput <- theta

  gradient <- function(theta) {
    grad <- gradientCompLogLikLms(
      theta = theta,
      model = model,
      P = P,
      sign = -1,
      epsilon = epsilon
    )

    non.finite <- !is.finite(grad)
    if (any(non.finite)) {
      # This might happen during a bad line search in optim()
      # as long as obsLogLikLms() is non finite as well, it should be ok,
      # as it will be rejected...
      ll <- compLogLikLms(
        theta = theta,
        model = model,
        P = P, sign = -1,
        epsilon = epsilon
      )

      mod_warnif(is.finite(ll),
        "Some gradient values are NaN, but objective is finite!"
      )

      dir <- sign(theta - thetaInput)
      dir[!is.finite(dir)] <- 0

      grad[non.finite] <- dir[non.finite] * .Machine$double.xmax^(1/10)
    }

    grad
  }

  objective <- function(theta) {
    ll <- compLogLikLms(
      theta = theta,
      model = model,
      P = P, sign = -1,
      epsilon = epsilon
    )

    if (!is.finite(ll)) .Machine$double.xmax^(1/10) else ll # optim doesn't handle Inf well
  }

  if (optimizer == "nlminb") {
    if (is.null(control$iter.max)) {
      # max.step is usually 1
      control$iter.max <- max.step
    }

    if (is.null(control$eval.max)) {
      # defaults to 200, so we keep the default unless iter.max > 100
      control$eval.max <- max(200L, 2 * control$iter.max)
    }

    est <- suppressWarnings(stats::nlminb(
      start     = theta,
      objective = objective,
      gradient  = gradient,
      upper     = model$params$bounds$upper,
      lower     = model$params$bounds$lower,
      control   = control,
      ...
    ))

  } else if (optimizer == "L-BFGS-B") {
    if (is.null(control$maxit)) {
      # max.step is usually 1
      control$maxit <- max.step
    }

    est <- stats::optim(
      par    = theta,
      fn     = objective,
      gr     = gradient,
      method = optim.method, control = control,
      lower  = model$params$bounds$lower,
      upper  = model$params$bounds$upper,
      ...
    )

    est$objective  <- est$value
    est$iterations <- est$counts[["function"]]

  } else {
    mod_msg_stop("Unrecognized optimizer, must be either 'nlminb' or 'L-BFGS-B'")
  }

  est
}


completeLogLikLmsGroupCpp <- function(submodel, P) {
  data <- submodel$data

  if (submodel$quad$adaptive %in% c("aghq", "fixed-n")) {
    completeLogLikLmsAghqCpp(
      modelR = submodel, Z = P$Z, Gamma = P$P,
      dataR = data$data.split, colidxR = data$colidx0,
      n = data$n.pattern, npatterns = data$p, ncores = ThreadEnv$n.threads
    )
  } else {
    completeLogLikLmsCpp(
      modelR = submodel, P = P, quad = P$quad,
      colidxR = data$colidx0, n = data$n.pattern,
      d = data$d.pattern, npatterns = data$p
    )
  }
}


compLogLikLms <- function(theta, model, P, sign = -1, ...) {
  tryCatch({
    modFilled <- fillModel(model = model, theta = theta, method = "lms")

    ll <- 0
    for (g in seq_len(model$info$n.groups)) {
      submodel <- modFilled$models[[g]]

      ll <- ll + completeLogLikLmsGroupCpp(submodel = submodel, P = P$P_GROUPS[[g]])
    }

    sign * ll

  }, error = \(e) NA)
}


compLogLikLmsGroup <- function(submodel, P, sign = -1, ...) {
  tryCatch({
    ll <- completeLogLikLmsGroupCpp(submodel = submodel, P = P)

    sign * ll

  }, error = \(e) NA)
}


obsLogLikLms <- function(theta, model, P, sign = -1, ...) {
  tryCatch({
    modFilled <- fillModel(model = model, theta = theta, method = "lms")
    ll <- 0

    for (g in seq_len(model$info$n.groups)) {
      ll <- ll + obsLogLikLmsGroup(
        submodel = modFilled$models[[g]], P = P$P_GROUPS[[g]], sign = 1
      )
    }

    sign * ll
  }, error = \(e) NA_real_)
}


# Recomputes the AGHQ posterior (and marginal density/obsLL) at `submodel`'s
# CURRENT theta, reusing the FIXED (E-step-time) adapted nodes/weights in
# `P$Z`/`P$W` -- i.e. Fisher's identity applied to the per-row grid: the
# observed-data quantities at a *trial* theta are obtained by re-weighting
# the same nodes, not by re-adapting them. Mirrors how observedGradientReverseFromModel()
# reuses the shared-node V/w for the non-aghq path.
aghqPosteriorAtTheta <- function(submodel, P) {
  data <- submodel$data
  aghqEstepLmsCpp(
    modelR = submodel, Z = P$Z, W = P$W, dataR = data$data.split,
    colidxR = data$colidx0, samplingWeights = P$sampling.weights,
    n = data$n.pattern, npatterns = data$p, ncores = ThreadEnv$n.threads
  )
}


obsLogLikLmsGroup <- function(submodel, P, sign = -1, ...) {
  tryCatch({
    data <- submodel$data

    ll <- if (submodel$quad$adaptive %in% c("aghq", "fixed-n")) {
      aghqPosteriorAtTheta(submodel = submodel, P = P)$obsLL
    } else {
      observedLogLikLmsCpp(
        modelR = submodel, dataR = data$data.split, colidxR = data$colidx0,
        P = P, n = data$n.pattern, npatterns = data$p,
        ncores = ThreadEnv$n.threads
      )
    }

    sign * ll
  }, error = \(e) NA_real_)
}



gradientCompLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6) {
  FGRAD <- function(modelR, P, block, row, col, symmetric, colidxR,
                    n, d, npatterns, eps, ncores) {
    if (modelR$quad$adaptive %in% c("aghq", "fixed-n")) {
      gradLogLikLmsAghqCpp(
        modelR = modelR, Z = P$Z, Gamma = P$P,
        dataR = modelR$data$data.split, colidxR = colidxR,
        n = n, block = block, row = row, col = col, symmetric = symmetric,
        npatterns = npatterns, ncores = ncores
      )
    } else {
      gradLogLikLmsCpp(
        modelR = modelR, P = P, block = block, row = row, col = col,
        symmetric = symmetric, colidxR = colidxR, n = n, d = d,
        npatterns = npatterns, eps = eps, ncores = ncores
      )
    }
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, FGRAD = FGRAD, FOBJECTIVE = compLogLikLmsGroup)
}


gradientObsLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6) {
  FGRAD <- function(modelR, P, block, row, col, symmetric, colidxR,
                    n, d, npatterns, eps, ncores) {
    if (modelR$quad$adaptive %in% c("aghq", "fixed-n")) {
      Gamma <- aghqPosteriorAtTheta(submodel = modelR, P = P)$Gamma
      gradLogLikLmsAghqCpp(
        modelR = modelR, Z = P$Z, Gamma = Gamma,
        dataR = modelR$data$data.split, colidxR = colidxR, n = n,
        block = block, row = row, col = col, symmetric = symmetric,
        npatterns = npatterns, ncores = ncores
      )
    } else {
      gradObsLogLikLmsCpp(
        modelR = modelR, dataR = modelR$data$data.split, colidxR = colidxR,
        P = P, block = block, row = row, col = col, symmetric = symmetric,
        n = n, eps = eps, npatterns = npatterns, ncores = ncores
      )
    }
  }

  gradientAllLogLikLms(
    theta = theta, model = model, P = P, sign = sign, epsilon = epsilon,
    FGRAD = FGRAD, FOBJECTIVE = obsLogLikLmsGroup
  )
}


gradientAllLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6,
                                 FGRAD, FOBJECTIVE) {
  useFDGradient <- model$params$gradientStruct$useFDGradient

  if (useFDGradient) gradient <- \(...) complicatedGradientAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE)
  else               gradient <- \(...) simpleGradientAllLogLikLms(..., FGRAD = FGRAD)

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
  Jacobian   <- refreshCovModelJacobian(theta, model, Jacobian)$J
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


obsLogLikLmsGroup_i <- function(submodel, P, sign = 1) {
  data <- submodel$data
  V  <- P$V
  w  <- P$w
  m  <- nrow(V)
  px <- numeric(data$n)

  for (i in seq_len(m)) {
    z_i        <- V[i, ]
    ms_i       <- muSigmaLmsCpp(model = submodel, z = z_i)
    mu_i    <- ms_i$mu
    sigma_i <- ms_i$sigma

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
  ms    <- muSigmaLmsCpp(model = modFilled, z = z)
  mu    <- ms$mu
  sigma <- ms$sigma

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
  densityMatrixLmsCpp(modelR = modFilled, V = z,
                      dataR = data$data.split, colidxR = data$colidx0,
                      n = data$n.pattern, samplingWeights = numeric(0L),
                      npatterns = data$p, ncores = ThreadEnv$n.threads)
}


hessianAllLogLikLms <- function(theta, model, P, sign = -1,
                                FHESS, FOBJECTIVE, .relStep = .Machine$double.eps ^ (1/5)) {
  useFDGradient <- model$params$gradientStruct$useFDGradient

  if (useFDGradient) hessian <- \(...) complicatedHessianAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
  else               hessian <- \(...) simpleHessianAllLogLikLms(..., FHESS = FHESS, .relStep = .relStep)

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
  jac         <- refreshCovModelJacobian(theta, model, Jacobian, Jacobian2)
  Jacobian    <- jac$J
  Jacobian2   <- jac$J2
  nlinDerivs  <- model$params$gradientStruct$nlinDerivs
  nlinDerivs2 <- model$params$gradientStruct$nlinDerivs2

  n.loc <- NROW(locations)
  nm <- locations$param
  H <- matrix(0.0, nrow = n.loc, ncol = n.loc, dimnames = list(nm, nm))
  grad <- stats::setNames(numeric(n.loc), nm = nm)

  for (g in seq_len(modelR$info$n.groups)) {
    locations.g <- locations[
      locations$group == g, , drop = FALSE
    ]

    if (!NROW(locations.g)) next

    submodelR <- modelR$models[[g]]
    data.g    <- submodelR$data

    HESS.g <- FHESS(
      modelR    = submodelR,
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
      ncores    = ThreadEnv$n.threads
    )

    H.g    <- HESS.g$Hessian
    grad.g <- HESS.g$gradient
    nm.g   <- locations.g$param

    dimnames(H.g) <- list(nm.g, nm.g)
    names(grad.g) <- nm.g

    H[nm.g, nm.g] <- H[nm.g, nm.g] + H.g
    grad[nm.g] <- grad[nm.g] + grad.g
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


hessianCompLogLikLms <- function(theta, model, P, sign = -1,
                                 .relStep = .Machine$double.eps ^ (1/5)) {

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR,
                    n, d, npatterns, ncores) {
    if (modelR$quad$adaptive %in% c("aghq", "fixed-n")) {
      hessCompLogLikLmsAghqCpp(
        modelR = modelR, Z = P$Z, Gamma = P$P, dataR = modelR$data$data.split,
        block = block, row = row, col = col, symmetric = symmetric,
        colidxR = colidxR, n = n, relStep = .relStep, npatterns = npatterns,
        minAbs = 0.0, ncores = ncores
      )
    } else {
      hessCompLogLikLmsCpp(modelR = modelR, P = P, block = block,
                           row = row, col = col, symmetric = symmetric,
                           colidxR = colidxR, n = n, d = d, relStep = .relStep,
                           npatterns = npatterns, minAbs = 0.0, ncores = ncores)
    }
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      FHESS = FHESS, FOBJECTIVE = compLogLikLmsGroup,
                      .relStep = .relStep)
}


hessianObsLogLikLms <- function(theta, model, P, sign = -1,
                                .relStep = .Machine$double.eps ^ (1/5)) {
  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep,
                    colidxR, n, d, npatterns, ncores) {
    hessObsLogLikLmsCpp(
      modelR = modelR, dataR = modelR$data$data.split, P = P,
      block = block, row = row, col = col, symmetric = symmetric,
      colidxR = colidxR, n = n, npatterns = npatterns,
      relStep = .relStep, minAbs = 0.0, ncores = ncores
    )
  }

  hessianAllLogLikLms(
    theta = theta, model = model, P = P, sign = sign,
    FHESS = FHESS, FOBJECTIVE = obsLogLikLmsGroup, .relStep = .relStep
  )
}


logdensAllObsNode <- function(theta, model, data, z, group = NULL) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  submodel <- if (!is.null(modFilled$models)) {
    if (is.null(group)) modFilled$models[[1L]] else modFilled$models[[group]]
  } else modFilled
  # densitySingleLms returns densities for all rows at given z (not log)
  dens <- densitySingleLms(z = z, modFilled = submodel, data = data)
  # guard against underflow/zeros
  log(pmax(dens, .Machine$double.xmin))
}


activeThetaIndicesLms <- function(model, group, p) {
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
completeScoresNodeFD <- function(theta, model, data, z,
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
    f0 <- logdensAllObsNode(theta, model, data, z, group = group)
    for (pos in seq_len(n_active)) {
      k <- active[pos]
      th1 <- theta; th1[k] <- th1[k] + epsilon
      f1 <- logdensAllObsNode(th1, model, data, z, group = group)
      S[, pos] <- (f1 - f0) / epsilon
    }
  } else { # central
    for (pos in seq_len(n_active)) {
      k <- active[pos]
      thp <- theta; thp[k] <- thp[k] + epsilon
      thm <- theta; thm[k] <- thm[k] - epsilon
      fp <- logdensAllObsNode(thp, model, data, z, group = group)
      fm <- logdensAllObsNode(thm, model, data, z, group = group)
      S[, pos] <- (fp - fm) / (2 * epsilon)
    }
  }
  S
}


lmsFirstDerivativeJacobian <- function(theta, model) {
  gradientStruct <- model$params$gradientStruct
  Jacobian <- refreshCovModelJacobian(
    theta, model, gradientStruct$Jacobian
  )$J
  nlinDerivs <- gradientStruct$nlinDerivs

  if (length(nlinDerivs)) {
    evalTheta  <- gradientStruct$evalTheta
    param.full <- stringr::str_split_i(
      colnames(Jacobian), pattern = "#", i = 1L
    )
    param.part <- rownames(Jacobian)
    THETA <- list2env(as.list(evalTheta(theta)))

    for (dep in names(nlinDerivs)) {
      for (indep in names(nlinDerivs[[dep]])) {
        deriv <- eval(expr = nlinDerivs[[dep]][[indep]], envir = THETA)
        Jacobian[param.part == indep, param.full == dep] <- deriv
      }
    }
  }

  Jacobian
}


completeScoresNodeAnalytical <- function(submodel, data, z, locations,
                                         Jacobian, active) {
  rawScores <- completeScoresNodeAnalyticalLmsCpp(
    modelR = submodel,
    dataR = data$data.split,
    z = c(z),
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric,
    colidxR = data$colidx0,
    n = data$n.pattern,
    npatterns = data$p,
    ncores = ThreadEnv$n.threads
  )

  mapping <- Jacobian[active, locations$param, drop = FALSE]
  rawScores %*% t(mapping)
}

# I_obs = I_com - I_mis using Louis' identity
observedInfoFromLouisLms <- function(model,
                                     theta,
                                     P = NULL,
                                     recompute.P = is.null(P),
                                     adaptive.quad.tol = 1e-12,
                                     fd.epsilon = 1e-6,
                                     fd.scheme = c("forward","central"),
                                     score.method = c("auto", "analytical", "finite.difference"),
                                     symmetrize = TRUE,
                                     jitter = 0.0,
                                     ...) {
  fd.scheme <- match.arg(fd.scheme)
  score.method <- match.arg(score.method)
  if (score.method == "auto") {
    score.method <- if (isTRUE(model$params$gradientStruct$useFDGradient)) {
      "finite.difference"
    } else {
      "analytical"
    }
  }

  aghqGroups <- vapply(model$models, FUN.VALUE = logical(1L),
                       FUN = \(sub) sub$quad$adaptive %in% c("aghq", "fixed-n"))
  mod_stopif(any(aghqGroups) && score.method != "analytical",
    "Louis' identity for `adaptive.quad %in% c(\"aghq\", \"fixed-n\")` currently ",
    "requires `score.method = \"analytical\"` (the default, unless the model's ",
    "gradient structure forces `useFDGradient = TRUE`).")

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

  if (score.method == "analytical") {
    modelFilled <- fillModel(model = model, theta = theta, method = "lms")
    locations <- model$params$gradientStruct$locations
    Jacobian <- lmsFirstDerivativeJacobian(theta, model)
  }

  total_M <- matrix(0.0, p, p, dimnames = list(lbl, lbl))
  Sbar    <- matrix(0.0, n_total, p)

  row_offset <- 0L
  for (g in seq_len(model$info$n.groups)) {
    submodel <- model$models[[g]]
    data.g   <- submodel$data
    P.g      <- P$P_GROUPS[[g]]
    rows     <- seq_len(data.g$n) + row_offset

    active_idx <- activeThetaIndicesLms(model, g, p)
    if (score.method == "analytical") {
      submodelFilled <- modelFilled$models[[g]]
      locations.g <- locations[locations$group == g, , drop = FALSE]
    }
    sampling_weights <- P.g$sampling.weights
    mod_stopif(length(sampling_weights) != data.g$n,
               "Invalid sampling-weight vector in LMS posterior object.")
    sqrt_weights <- sqrt(pmax(sampling_weights, 0))

    if (aghqGroups[[g]]) {
      # Every row has its own adapted nodes, so there's no shared z_j to loop
      # over -- louisRawScoresAghqCpp() accumulates total_M/Sbar directly in
      # one pass over (row, node) pairs, in raw "locations" space. Map through
      # the Jacobian once here (equivalent to completeScoresNodeAnalytical()
      # doing it per-node, since the mapping is linear and commutes with the
      # row/node sums).
      raw <- louisRawScoresAghqCpp(
        modelR = submodelFilled, Z = P.g$Z, Gamma = P.g$P,
        sqrtSamplingWeights = sqrt_weights,
        dataR = data.g$data.split, colidxR = data.g$colidx0,
        n = data.g$n.pattern, block = locations.g$block,
        row = locations.g$row, col = locations.g$col,
        symmetric = locations.g$symmetric, npatterns = data.g$p,
        ncores = ThreadEnv$n.threads
      )

      mapping <- Jacobian[active_idx, locations.g$param, drop = FALSE]

      total_M[active_idx, active_idx] <- total_M[active_idx, active_idx] +
        mapping %*% raw$total_M %*% t(mapping)
      Sbar[rows, active_idx] <- Sbar[rows, active_idx] +
        raw$Sbar %*% t(mapping)

      row_offset <- row_offset + data.g$n
      next
    }

    Jg <- length(P.g$w)
    for (j in seq_len(Jg)) {
      z_j <- P.g$V[j, , drop = FALSE]
      S_j <- if (score.method == "analytical") {
        completeScoresNodeAnalytical(
          submodel = submodelFilled, data = data.g, z = z_j,
          locations = locations.g, Jacobian = Jacobian, active = active_idx
        )
      } else {
        completeScoresNodeFD(
          theta, model, data.g, z_j,
          epsilon = fd.epsilon, scheme = fd.scheme,
          group = if (model$info$n.groups > 1L) g else NULL,
          active = active_idx
        )
      }
      # P.g$P contains w_i * gamma_ij. This is directly suitable for
      # E_w[s s']; for E_w[s]E_w[s]' we need sqrt(w_i) * E[s], not
      # w_i * E[s], otherwise crossprod() squares the sampling weights.
      r_j <- P.g$P[, j]
      mean_score_weights <- numeric(length(r_j))
      positive_weights <- sqrt_weights > 0
      mean_score_weights[positive_weights] <-
        r_j[positive_weights] / sqrt_weights[positive_weights]

      Rhalf <- sqrt(pmax(r_j, 0))
      if (NROW(S_j) && NCOL(S_j)) {
        X <- S_j * Rhalf
        block <- total_M[active_idx, active_idx, drop = FALSE]
        block <- block + crossprod(X)
        total_M[active_idx, active_idx] <- block

        Sbar_block <- Sbar[rows, active_idx, drop = FALSE]
        Sbar_block <- Sbar_block + (S_j * mean_score_weights)
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
    Icom <- sym(Icom)
    Imis <- sym(Imis)
    Iobs <- sym(Iobs)
  }

  if (jitter > 0) {
    Icom <- Icom + diag(jitter, p)
    Imis <- Imis + diag(jitter, p)
    Iobs <- Iobs + diag(jitter, p)
  }

  list(I.obs = Iobs, I.com = Icom, I.mis = Imis, P = P)
}
