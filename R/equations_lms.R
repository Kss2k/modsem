estepLms <- function(model, theta, lastQuad = NULL, recalcQuad = FALSE,
                     adaptive.quad.tol = 1e-12, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  P <- list(P_GROUPS = vector("list", length = model$info$n.groups),
            quad = NULL, obsLL = NULL)

  for (g in seq_len(model$info$n.groups)) {
    P$P_GROUPS[[g]] <- estepLmsGroup(
      submodel = modFilled$models[[g]], lastQuad = lastQuad[[g]],
      recalcQuad = recalcQuad, adaptive.quad.tol = adaptive.quad.tol, ...
    )
  }

  P$quad  <- lapply(P$P_GROUPS, FUN = \(P) P$quad)
  P$obsLL <- sum(vapply(P$P_GROUPS, FUN.VALUE = numeric(1L), FUN = \(P) P$obsLL))

  P
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
    P    <- W * densityLms(V, modFilled = submodel, data = data)
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


gradientCompLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6) {
  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, FGRAD = gradLogLikLmsCpp, FOBJECTIVE = compLogLikLms)
}


gradientAllLogLikLms <- function(theta, model, P, sign = -1, epsilon = 1e-6,
                                 FGRAD, FOBJECTIVE) {
  hasCovModel <- model$params$gradientStruct$hasCovModel

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
    locations_g <- locations[locations$group == g, , drop = FALSE]

    data      <- submodelR$data
    block     <- locations_g$block
    row       <- locations_g$row
    col       <- locations_g$col
    param     <- locations_g$param
    symmetric <- locations_g$symmetric

    grad_g <- FGRAD(modelR    = submodelR,
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

    grad[locations_g$param, ] <- grad_g
  }

  if (length(nlinDerivs)) {
    evalTheta  <- model$params$gradientStruct$evalTheta
    param.full <- stringr::str_split_i(colnames(Jacobian), i = 1L) # non-unique
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
    data_g   <- submodel$data
    P_g      <- P$P_GROUPS[[g]]

    ll_g <- observedLogLikLmsCpp(submodel,
                                 dataR = data_g$data.split,
                                 colidxR = data_g$colidx0,
                                 P = P_g,
                                 n = data_g$n.pattern,
                                 npatterns = data_g$p,
                                 ncores = ThreadEnv$n.threads)
    ll <- ll + ll_g
  }

  sign * ll
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
    obsLogLikLms(theta = theta, model = model, P = P, sign = sign)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, FGRAD = FGRAD,
                       FOBJECTIVE = FOBJECTIVE)
}


obsLogLikLms_i_single <- function(submodel, data, P, sign = 1) {
  V <- P$V
  w <- P$w
  m <- nrow(V)
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

  sign * log(px)
}


obsLogLikLms_i <- function(theta, model, P, sign = 1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  contribs <- vector("list", length = modFilled$info$n.groups)

  for (g in seq_len(modFilled$info$n.groups)) {
    submodel <- modFilled$models[[g]]
    data_g   <- submodel$data
    P_g      <- P$P_GROUPS[[g]]

    contribs[[g]] <- obsLogLikLms_i_single(submodel, data_g, P_g, sign = sign)
  }

  unlist(contribs, use.names = FALSE)
}


# gradient function of obsLogLikLms_i
gradientObsLogLikLms_i <- function(theta, model, P, sign = 1, epsilon = 1e-4) {
  group_sizes <- vapply(model$models, function(sub) sub$data$n, numeric(1L))
  n_total     <- sum(group_sizes)
  p           <- length(theta)

  baseLL <- obsLogLikLms_i(theta, model, P = P, sign = sign)

  offsets <- c(0L, cumsum(group_sizes))
  rows_by_group <- lapply(seq_along(group_sizes), function(g) {
    idx_start <- offsets[g] + 1L
    idx_end   <- offsets[g + 1L]
    seq.int(idx_start, idx_end)
  })

  base_by_group <- Map(function(rows) baseLL[rows], rows_by_group)

  groups_by_param <- vector("list", length = p)
  for (g in seq_along(group_sizes)) {
    active_idx <- .activeThetaIndicesLms(model, g, p)
    for (idx in active_idx) {
      groups_by_param[[idx]] <- unique(c(groups_by_param[[idx]], g))
    }
  }

  grad_mat <- matrix(0.0, nrow = n_total, ncol = p,
                     dimnames = list(NULL, names(theta)))

  for (idx in seq_len(p)) {
    param_groups <- groups_by_param[[idx]]
    if (!length(param_groups)) next

    theta_pert <- theta
    theta_pert[[idx]] <- theta_pert[[idx]] + epsilon

    modFilled <- fillModel(model = model, theta = theta_pert, method = "lms")

    for (g in param_groups) {
      rows    <- rows_by_group[[g]]
      submodel <- modFilled$models[[g]]
      data_g   <- submodel$data
      P_g      <- P$P_GROUPS[[g]]

      ll_new <- obsLogLikLms_i_single(submodel, data_g, P_g, sign = sign)
      grad_mat[rows, idx] <- (ll_new - base_by_group[[g]]) / epsilon
    }
  }

  grad_mat
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


hessianAllLogLikLms <- function(theta, model, P, sign = -1,
                                FHESS, FOBJECTIVE, .relStep = .Machine$double.eps ^ (1/5)) {
  hasCovModel <- model$params$gradientStruct$hasCovModel

  if (hasCovModel) hessian <- \(...) complicatedHessianAllLogLikLms(..., FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
  else             hessian <- \(...) simpleHessianAllLogLikLms(..., FHESS = FHESS, .relStep = .relStep)

  hessian(theta = theta, model = model, P = P, sign = sign)
}


compHessianAllLogLikLms <- function(theta, model, P, sign = -1,
                                    .relStep = .Machine$double.eps ^ (1/6),
                                    FOBJECTIVE) {
  nlme::fdHess(pars = theta, FOBJECTIVE, model = model, P = P,
               sign = sign, .relStep = .relStep)
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
    locations_g <- locations[locations$group == g, , drop = FALSE]
    if (!NROW(locations_g)) next

    submodelR <- modelR$models[[g]]
    data_g    <- submodelR$data

    HESS_g <- FHESS(modelR    = submodelR,
                    P         = P$P_GROUPS[[g]],
                    block     = locations_g$block,
                    row       = locations_g$row,
                    col       = locations_g$col,
                    colidxR   = data_g$colidx0,
                    n         = data_g$n.pattern,
                    d         = data_g$d.pattern,
                    npatterns = data_g$p,
                    symmetric = locations_g$symmetric,
                    .relStep  = .relStep,
                    ncores    = ThreadEnv$n.threads)

    H_g    <- HESS_g$Hessian
    grad_g <- HESS_g$gradient

    dimnames(H_g) <- list(locations_g$param, locations_g$param)
    names(grad_g) <- locations_g$param

    H[locations_g$param, locations_g$param] <- H[locations_g$param, locations_g$param] + H_g
    grad[locations_g$param] <- grad[locations_g$param] + grad_g
  }

  if (length(nlinDerivs)) {
    evalTheta  <- model$params$gradientStruct$evalTheta
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
               .relStep = .relStep)$Hessian
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

  FOBJECTIVE <- function(theta, model, P, sign, data, ...) {
    obsLogLikLms(theta = theta, model = model, P = P, sign = sign)
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

  FOBJECTIVE <- function(theta, model, P, data, sign) {
    compLogLikLms(theta = theta, model = model, P = P, sign = sign)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      FHESS = FHESS, FOBJECTIVE = FOBJECTIVE, .relStep = .relStep)
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
    data_g   <- submodel$data
    P_g      <- P$P_GROUPS[[g]]
    rows     <- seq_len(data_g$n) + row_offset

    active_idx <- .activeThetaIndicesLms(model, g, p)
    Jg <- length(P_g$w)
    for (j in seq_len(Jg)) {
      z_j <- P_g$V[j, , drop = FALSE]
      S_j <- .completeScoresNodeFD(theta, model, data_g, z_j,
                                   epsilon = fd.epsilon, scheme = fd.scheme,
                                   group = if (model$info$n.groups > 1L) g else NULL,
                                   active = active_idx)
      r_j <- P_g$P[, j]

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

    row_offset <- row_offset + data_g$n
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
