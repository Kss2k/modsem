estepLms <- function(model, theta, data, lastQuad = NULL, recalcQuad = FALSE,
                     adaptive.quad.tol = 1e-12, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  ## --- Quadrature setup (unchanged) ---
  if (model$quad$adaptive && (recalcQuad || is.null(lastQuad))) {
    m <- model$quad$m; a <- model$quad$a; b <- model$quad$b; k <- model$quad$k
    if (!is.null(lastQuad)) m.ceil <- lastQuad$m.ceil
    else if (k > 1) m.ceil <- m
    else m.ceil <- round(estMForNodesInRange(m, a = -5, b = 5))

    quad <- tryCatch({
      adaptiveGaussQuadrature(
        fun = densityLms,                                     # OK for ML; unused for PML below
        collapse = \(x) sum(log(rowSums(x))),                 # idem
        modFilled = modFilled, data = data, a = a, b = b,
        m = m, k = k, m.ceil = m.ceil, tol = adaptive.quad.tol
      )
    }, error = function(e) {
      warning2("Calculation of adaptive quadrature failed!\n", e, immediate. = FALSE)
      NULL
    })

    if (is.null(quad)) {
      return(estepLms(model, theta, data, lastQuad = if (!is.null(lastQuad)) lastQuad else model$quad,
                      recalcQuad = FALSE, ...))
    }
    V <- quad$n; w <- quad$w
  } else {
    quad <- if (model$quad$adaptive) lastQuad else model$quad
    V    <- quad$n
    w    <- quad$w
  }

  estimator <- tolower(modFilled$info$estimator)

  ## ================= PML BRANCH (pair-level responsibilities) =================
  if (estimator == "pml") {
    Q <- nrow(V)

    # For each node q: list over patterns of (n_j x npairs_j) log-densities
    log_list_per_q <- vector("list", Q)
    for (q in seq_len(Q)) {
      zq <- V[q, , drop = FALSE]
      log_list_per_q[[q]] <- .pairwise_logdens_node(z = zq, modFilled = modFilled, data = data)
    }

    # Stack by *columns* (pair-major) within each pattern, then concatenate patterns.
    # This matches C++ 'vectorise(Lq, dim=0)' used in completeLogLikFromModelPML.
    stack_by_pattern <- function(L_list) {
      unlist(lapply(L_list, function(M) as.vector(M)), use.names = FALSE)
    }

    # Build a big R_total x Q matrix of log f^{(q)}_{pair}(y)
    L_big <- do.call(
      cbind,
      lapply(log_list_per_q, stack_by_pattern)
    )  # rows: (pattern blocks), each block = c(vec(L_{pair1}), vec(L_{pair2}), ...)

    # Add log quadrature weights inside the pair mixture
    LW <- sweep(L_big, 2L, log(w), FUN = "+")

    # Row-wise log-sum-exp across q gives the observed composite loglik contributions
    log_row_sum <- .rowLogSumExps(LW)
    obsLL <- sum(log_row_sum)

    # Pair-level responsibilities (R_total x Q)
    R_big <- exp(LW - log_row_sum)

    # Return a P list that C++ knows how to consume:
    # - P$P is the pair-level responsibilities (R_total x Q)
    # - P$V, P$w, P$quad as usual
    # - keep placeholders for ML fields to avoid breaking existing wrappers

    getEmptyList <- \(x) 
      lapply(seq_along(data$ids), FUN = \(j) x)
    getEmptyListList <- \(x) 
      lapply(seq_along(data$ids), FUN = \(j) lapply(seq_len(Q), FUN = \(i) x))

    Pout <- list(
      P    = R_big,
      V    = V,
      w    = w,
      quad = quad,
      mean = getEmptyListList(0),
      cov = getEmptyListList(matrix(0)),
      tgamma = getEmptyList(0)
    )

    return(list(
      P     = R_big,
      mean  = Pout$mean,     # not used by PML C++ path
      cov   = Pout$cov,      # idem
      tgamma= Pout$tgamma,   # idem
      V     = V,
      w     = w,
      obsLL = obsLL,
      quad  = quad
  ))
}

  ## ================= ML BRANCH (unchanged) ===================================
  # Your original ML E-step (observation-level responsibilities) stays as-is:
  W <- matrix(w, nrow = data$n, ncol = length(w), byrow = TRUE)
  P  <- W * densityLms(V, modFilled = modFilled, data = data)  # n x Q
  density        <- rowSums(P)
  observedLogLik <- sum(log(density))
  P              <- P / density

  # (Your weighted means/covs code unchanged)
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

  list(P = P, V = V, w = w, quad = quad,
       mean = wMeans, cov = wCovs, tgamma = tGamma, obsLL = observedLogLik)
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

  gradient2 <- function(theta) {
    theta0 <- theta
    f0 <- objective(theta0)
    eps <- 1e-4
    grad <- numeric(length(theta0))
    for (i in seq_along(theta0)) {
      theta1 <- theta0
      theta1[i] <- theta0[i] + eps
      
      f1 <- objective(theta1)
      grad[i] <- (f1 - f0) / eps
    }
    grad
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
  # tryCatch({
    modFilled <- fillModel(model = model, theta = theta, method = "lms")
    sign * completeLogLikLmsCpp(modelR=modFilled, dataR = data$data.split, P=P, quad=P$quad,
                                colidxR = data$colidx0, n = data$n.pattern,
                                d = data$d.pattern, npatterns = data$p)
  # }, error = \(e) NA)
}


gradientCompLogLikLms <- function(theta, model, P, data, sign = -1, epsilon = 1e-6) {
  FGRAD <- function(modelR, P, block, row, col, symmetric, colidxR, npatterns,
                    eps, ncores, n, ...) {
    gradLogLikLmsCpp(modelR = modelR, dataR = data$data.split, P = P, 
                     block = block, row = row, col = col,
                     symmetric = symmetric, colidxR = colidxR,
                     n = n,
                     d = data$d.pattern,
                     npatterns = npatterns,
                     eps = eps,
                     ncores = ncores)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign, data = data,
                       epsilon = epsilon, FGRAD = FGRAD, FOBJECTIVE = compLogLikLms)
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


# Returns an n x V matrix of log-likelihoods:
# - ML:   full MVN log-densities per node (dmvn(..., log=TRUE))
# - PML:  pairwise composite *log*-likelihoods per node (probPML_Fast returns logs)
densityLms <- function(z, modFilled, data) {
  if (is.null(dim(z))) z <- matrix(z, ncol = modFilled$quad$k)

  # lapplyMatrix should combine the returned numeric vectors column-wise
  lapplyMatrix(seq_len(nrow(z)), FUN.VALUE = numeric(data$n), FUN = function(i) {
    densitySingleLms(z = z[i, , drop = FALSE], modFilled = modFilled, data = data)
  })
}


densitySingleLms <- function(z, modFilled, data) {
  estimator <- tolower(modFilled$info$estimator)

  mu    <- muLmsCpp(model = modFilled, z = z)      # length p
  Sigma <- sigmaLmsCpp(model = modFilled, z = z)   # p x p

  # This function returns a length-n *log*-likelihood vector
  out_log <- numeric(data$n)

  offset <- 1L
  for (id in data$ids) {  # iterate missingness patterns
    n.pattern <- data$n.pattern[[id]]
    end       <- offset + n.pattern - 1L

    colidx <- data$colidx[[id]]
    Xid    <- data$data.split[[id]]               # n_id x p_id

    if (estimator == "ml") {
      # Full multivariate normal, LOG scale
      out_log[offset:end] <-
        dmvn(
          Xid,
          mean  = mu[colidx],
          sigma = Sigma[colidx, colidx, drop = FALSE],
          log   = FALSE
        )

    } else if (estimator == "pml") {
      # Pairwise composite, LOG scale via your Rcpp probPML_Fast
      # probPML_Fast signature (Rcpp): (data, mu, Sigma, isOrderedEnum, thresholds)
      # - isOrderedEnum: pass as integer (0=continuous, >0 = row index into thresholds)
      # - thresholds:    matrix; rows aligned with those indices
      isOrd <- as.integer(modFilled$info$isOrderedEnum[colidx])

      out_log[offset:end] <-
        exp(probPML_Fast(
          data          = as.matrix(Xid),
          mu            = mu[colidx],
          Sigma         = Sigma[colidx, colidx, drop = FALSE],
          isOrderedEnum = isOrd,
          thresholds    = modFilled$matrices$thresholds
        ))

    } else {
      stop("Unknown estimator '", estimator, "'. Expected 'ml' or 'pml'.")
    }

    offset <- end + 1L
  }

  out_log
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

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR,
                    n, d, npatterns, ncores) {
    hessCompLogLikLmsCpp(modelR = modelR, dataR = data$data.split, P = P, block = block,
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


probPMLr <- function(z, modFilled, data) {
  M <- modFilled$matrices
  T <- M$thresholds

 
  X <- data$data.split[[1L]]
  
  n <-NROW(X)
  p <- NCOL(X)
  ldensity <- rep(0, n)
  getlower <- \(t, i) sapply(t, \(ti) T[i, ti])
  getupper <- \(t, i) sapply(t, \(ti) T[i, ti+1L])
  PMVNORM <- function(lower.r, upper.r, lower.v, upper.v, mean, Sigma) {
    lower <- cbind(lower.r, lower.v)
    upper <- cbind(upper.r, upper.v)

    log(vapply(
      seq_along(lower.r), FUN.VALUE = numeric(1L), FUN = \(i)
      mvtnorm::pmvnorm(lower=lower[i, ], upper=upper[i, ], mean = mean, sigma = Sigma)
    ))
  }
  
  mu    <- as.vector(muLmsCpp(model = modFilled, z = z))      # length p
  Sigma <- sigmaLmsCpp(model = modFilled, z = z)   # p x p
  for (i in seq_len(p)) for (j in seq_len(i - 1)) {
    r <- as.integer(X[, i])
    v <- as.integer(X[, j])
    lower.r <- getlower(r, i)
    lower.v <- getlower(v, j)
    upper.r <- getupper(r, i)
    upper.v <- getupper(v, j)

    if (any(lower.r > upper.r)) browser()
    ldensity_ij <- PMVNORM(lower.r, upper.r, lower.v, upper.v, mu[c(i, j)], Sigma[c(i, j), c(i, j)])
    ldensity <- ldensity + ldensity_ij
  }

  ldensity_ij2 <- probPML_Fast(cbind(r, v), mu = mu[c(i,j)],
                          Sigma = Sigma[c(i,j), c(i,j)], isOrderedEnum = c(i, j), thresholds = T)
  ldensity2 <- probPML_Fast(as.matrix(X), mu = mu, Sigma = Sigma, isOrderedEnum = modFilled$info$isOrderedEnum,
                       thresholds = T)

  ldensity
}


# Enumerate column pairs r<s for a p-column block; returns a 2xK matrix (combn order)
.pair_index_mat <- function(p) utils::combn(p, 2, simplify = TRUE)

# row-wise log-sum-exp for a numeric matrix
.rowLogSumExps <- function(M) {
  # M: R x Q
  m <- matrixStats::rowMaxs(M)                 # length R
  m + log(rowSums(exp(M - m)))
}

# Compute per-observation log bivariate density for a TWO-COLUMN matrix under mean & 2x2 Sigma.
# Handles continuous or ordinal (via your probPML_Fast for 2 columns).
# - X2: n x 2 (numeric or integer categories if ordinal)
# - mu2: length-2 numeric
# - Sig2: 2x2 numeric
# - isOrd2: length-2 integer flags (0=cont, >0=ordinal id into thresholds)
# - thresholds: your thresholds matrix
.log_bvn_pair <- function(X2, mu2, Sig2, isOrd2, thresholds) {
  ord_any <- any(isOrd2 != 0L)
  if (!ord_any) {
    # continuous-continuous -> log N2(x; mu2, Sig2) for each row
    mvtnorm::dmvnorm(X2, mean = mu2, sigma = Sig2, log = TRUE)
  } else {
    # ordinal/ordinal or mixed -> rectangle prob (log) using your fast routine
    # probPML_Fast returns *log* prob for a matrix when given isOrderedEnum != 0
    as.numeric(probPML_Fast(
      data          = as.matrix(X2),
      mu            = mu2,
      Sigma         = Sig2,
      isOrderedEnum = as.integer(isOrd2),
      thresholds    = thresholds
    ))
  }
}


# Returns a list over patterns j; each entry is an n_j x npairs_j matrix of *log* bivariate densities
# under the given node z.
.pairwise_logdens_node <- function(z, modFilled, data) {
  mu    <- muLmsCpp(   model = modFilled, z = z)
  Sigma <- sigmaLmsCpp(model = modFilled, z = z)
  Tmat  <- modFilled$matrices$thresholds

  out <- vector("list", length(data$ids))

  for (jj in seq_along(data$ids)) {
    id     <- data$ids[[jj]]
    Xj     <- data$data.split[[id]]      # n_j x p_j
    idx    <- data$colidx[[id]]          # indices into (X,Y)
    isOrd  <- as.integer(modFilled$info$isOrderedEnum[idx])
    n_j    <- NROW(Xj)
    p_j    <- NCOL(Xj)
    if (p_j < 2L || n_j == 0L) { out[[jj]] <- matrix(0, n_j, 0); next }

    pairs  <- .pair_index_mat(p_j)       # 2 x npairs
    npairs <- NCOL(pairs)
    Lj     <- matrix(NA_real_, n_j, npairs)

    # subset mean/cov for this pattern once
    mu_sub  <- mu[idx]
    Sig_sub <- Sigma[idx, idx, drop = FALSE]

    for (k in seq_len(npairs)) {
      r <- pairs[1L, k]; s <- pairs[2L, k]         # local (1..p_j)
      glob <- idx[c(r, s)]
      mu2  <- mu_sub[c(r, s)]
      Sig2 <- Sig_sub[c(r, s), c(r, s), drop = FALSE]
      X2   <- Xj[, c(r, s), drop = FALSE]
      is2  <- isOrd[c(r, s)]
      Lj[, k] <- .log_bvn_pair(X2, mu2, Sig2, is2, Tmat)  # length n_j
    }
    out[[jj]] <- Lj  # n_j x npairs
  }
  out
}
