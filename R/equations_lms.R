estepLms <- function(model, theta, data, lastQuad = NULL, recalcQuad = FALSE,
                     adaptive.quad.tol = 1e-12, ...) {

  # Fill model once (we'll reuse for mu at each node)
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  # ---------------------------
  # Quadrature weights/nodes
  # ---------------------------
  if (model$quad$adaptive && (recalcQuad || is.null(lastQuad))) {
    m <- model$quad$m
    a <- model$quad$a
    b <- model$quad$b
    k <- model$quad$k

    if (!is.null(lastQuad)) m.ceil <- lastQuad$m.ceil
    else if (k > 1) m.ceil <- m
    else m.ceil <- round(estMForNodesInRange(m, a = -5, b = 5))

    quad <- tryCatch({
      adaptiveGaussQuadrature(
        fun = densityLms, collapse = \(x) sum(log(rowSums(x))),
        modFilled = modFilled, data = data, a = a, b = b, m = m,
        k = k, m.ceil = m.ceil, tol = adaptive.quad.tol
      )
    }, error = function(e) {
      warning2("Calculation of adaptive quadrature failed!\n", e, immediate. = FALSE)
      NULL
    })

    if (is.null(quad)) {
      return(estepLms(model = model, theta = theta, data = data,
                      lastQuad = if (!is.null(lastQuad)) lastQuad else model$quad,
                      recalcQuad = FALSE, ...))
    }

    # Posterior numerators already computed (quad$W * quad$F)
    P <- quad$W * quad$F          # n x J
    V <- quad$n                   # J x k (nodes)
    w <- quad$w                   # length J
  } else {
    quad <- if (model$quad$adaptive) lastQuad else model$quad
    V    <- quad$n
    w    <- quad$w
    W    <- matrix(w, nrow = data$n, ncol = length(w), byrow = TRUE) # n x J
    P    <- W * densityLms(V, modFilled = modFilled, data = data)    # n x J
  }

  # ---------------------------
  # Normalize posteriors
  # ---------------------------
  density        <- pmax(rowSums(P), .Machine$double.xmin)
  observedLogLik <- sum(log(density))
  P              <- P / density  # each row sums to 1

  # ---------------------------
  # Precompute names & masks
  # ---------------------------
  allX <- model$info$allIndsXis
  allY <- model$info$allIndsEtas
  pX   <- length(allX)

  ordX_idx <- model$info$ordinalX_idx %||% setNames(rep(FALSE, length(allX)), allX)
  ordY_idx <- model$info$ordinalY_idx %||% setNames(rep(FALSE, length(allY)), allY)

  thrX <- model$matrices$thresholdsX %||% list()
  thrY <- model$matrices$thresholdsY %||% list()

  # ---------------------------
  # Allocate outputs
  # ---------------------------
  J <- length(w)
  wMeans <- vector("list", J)
  wCovs  <- vector("list", J)
  tGamma <- vector("list", J)

  # ---------------------------
  # Loop over nodes
  # ---------------------------
  for (i in seq_len(J)) {
    p_ij <- P[, i]                 # length n posterior for node i
    tGamma[[i]] <- sum(p_ij)       # overall weight at this node (across all obs)

    # node-specific implied means (for all observed: X followed by Y)
    z_i    <- V[i, , drop = FALSE]
    mu_all <- as.numeric(muLmsCpp(modFilled, as.numeric(z_i)))  # length = pX + pY
    muX    <- setNames(mu_all[seq_len(pX)], allX)
    muY    <- setNames(mu_all[seq_along(allY) + pX], allY)

    # per-pattern containers
    wMeans[[i]] <- vector("list", length = length(data$ids))
    wCovs [[i]] <- vector("list", length = length(data$ids))
    tGamma[[i]] <- numeric(length = length(data$ids))

    offset <- 1L
    for (jj in seq_along(data$ids)) {
      id        <- data$ids[[jj]]
      n.pattern <- data$n.pattern[[id]]
      end       <- offset + n.pattern - 1L

      data.id <- data$data.split[[id]]         # n_i x q_i matrix (observed, by pattern)
      cols    <- colnames(data.id)
      pj      <- p_ij[offset:end]              # length n_i

      # build effective means and per-row var adds
      eff    <- data.id
      addVar <- matrix(0, nrow = nrow(eff), ncol = ncol(eff))
      colnames(addVar) <- cols

      # masks for this pattern's columns
      ordX_cols <- names(ordX_idx)[ordX_idx]
      ordX_mask <- cols %in% ordX_cols
      ordY_cols <- names(ordY_idx)[ordY_idx]
      ordY_mask <- cols %in% ordY_cols

      # ----- ORDINAL X -----
      if (any(ordX_mask, na.rm = TRUE)) {
        for (nm in cols[ordX_mask]) {
          xcat <- as.integer(eff[, nm])
          if (min(xcat, na.rm = TRUE) == 0L) xcat <- xcat + 1L
          C <- max(xcat, na.rm = TRUE)
          tau <- thrX[[nm]]
          if (is.null(tau) || length(tau) != (C - 1L))
            stop(sprintf("Thresholds for ordinal X '%s' missing/wrong length.", nm))
          mu_k <- muX[[nm]]
          tL <- c(-Inf, tau); tU <- c(tau, Inf)
          a <- tL[xcat] - mu_k
          b <- tU[xcat] - mu_k
          tm <- .trunc_moments_std(a, b)
          eff[, nm]    <- mu_k + tm$mean
          addVar[, nm] <- tm$var
        }
      }

      # ----- ORDINAL Y -----
      if (any(ordY_mask, na.rm = TRUE)) {
        for (nm in cols[ordY_mask]) {
          ycat <- as.integer(eff[, nm])
          if (min(ycat, na.rm = TRUE) == 0L) ycat <- ycat + 1L
          C <- max(ycat, na.rm = TRUE)
          tau <- thrY[[nm]]
          if (is.null(tau) || length(tau) != (C - 1L))
            stop(sprintf("Thresholds for ordinal Y '%s' missing/wrong length.", nm))
          mu_k <- muY[[nm]]
          tL <- c(-Inf, tau); tU <- c(tau, Inf)
          a <- tL[ycat] - mu_k
          b <- tU[ycat] - mu_k
          tm <- .trunc_moments_std(a, b)
          eff[, nm]    <- mu_k + tm$mean
          addVar[, nm] <- tm$var
        }
      }

      # weighted mean and covariance for this pattern at node i
      wsum <- sum(pj)
      wm   <- colSums(eff * pj) / wsum
      Xc   <- sweep(eff, 2, wm, "-")
      wcov <- t(Xc) %*% (Xc * pj)

      # add truncated-normal within-row variances on the diagonal
      if (any(ordX_mask, na.rm = TRUE) || any(ordY_mask, na.rm = TRUE)) {
        diag(wcov) <- diag(wcov) + colSums(addVar * pj)
      }

      wMeans[[i]][[jj]] <- wm
      wCovs [[i]][[jj]] <- wcov
      tGamma[[i]][[jj]] <- wsum

      offset <- end + 1L
    }
  }

  list(P = P, mean = wMeans, cov = wCovs, tgamma = tGamma, V = V, w = w,
       obsLL = observedLogLik, quad = quad)
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
  tryCatch({
    modFilled <- fillModel(model = model, theta = theta, method = "lms")
    sign * completeLogLikLmsCpp(modelR=modFilled, P=P, quad=P$quad,
                                colidxR = data$colidx0, n = data$n.pattern,
                                d = data$d.pattern, npatterns = data$p)
  }, error = \(e) NA)
}


gradientCompLogLikLms <- function(theta, model, P, data, sign = -1, epsilon = 1e-6) {
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


# Robust truncated standard-normal moments for interval (a, b).
# Vectorised over a,b. Handles +/-Inf cleanly and is stable in the tails.
.trunc_moments_std <- function(a, b) {
  stopifnot(length(a) == length(b))
  n <- length(a)

  # Guard against accidental a >= b due to threshold ties; jitter minimally
  bad <- !(a < b)
  if (any(bad)) {
    bump <- 1e-10
    b[bad] <- a[bad] + bump
  }

  # helpers
  logPhi <- function(x, lower = TRUE) pnorm(x, log.p = TRUE, lower.tail = lower)
  phi    <- function(x) dnorm(x)
  # log(exp(u) - exp(v)) for u >= v
  logspace_sub <- function(u, v) {
    out <- u + log1p(-exp(v - u))
    # when u==v, result is log(0) = -Inf; that's fine (will clamp later)
    out
  }

  # Cases: both finite, a=-Inf, b=+Inf, both infinite
  is_minf <- is.infinite(a) & a < 0
  is_pinf <- is.infinite(b) & b > 0
  both_inf <- is_minf & is_pinf
  mid <- !(is_minf | is_pinf)         # both finite
  left <- is_minf & !is_pinf          # (-Inf, b]
  right <- !is_minf & is_pinf         # [a, +Inf)

  # Allocate outputs
  m <- numeric(n)
  v <- numeric(n)

  # -------- both infinite: (-Inf, +Inf) -> standard normal --------
  if (any(both_inf)) {
    m[both_inf] <- 0
    v[both_inf] <- 1
  }

  # -------- finite a,b: use log-space for Z --------
  if (any(mid)) {
    ai <- a[mid]; bi <- b[mid]
    # Z = Phi(b) - Phi(a) in log-space
    lPhi_b <- logPhi(bi, lower = TRUE)
    lPhi_a <- logPhi(ai, lower = TRUE)
    # For valid intervals, Phi(b) >= Phi(a). Still be defensive:
    swap <- lPhi_b < lPhi_a
    tmp  <- lPhi_b; lPhi_b[swap] <- lPhi_a[swap]; lPhi_a[swap] <- tmp

    logZ <- logspace_sub(lPhi_b, lPhi_a)
    Z    <- pmax(exp(logZ), .Machine$double.xmin)

    phia <- phi(ai); phib <- phi(bi)
    mi   <- (phia - phib) / Z
    vi   <- 1 + (ai * phia - bi * phib) / Z - mi^2

    m[mid] <- mi
    v[mid] <- pmax(vi, 1e-12)
  }

  # -------- left-open: (-Inf, b] --------
  if (any(left)) {
    bi <- b[left]
    # Z = Phi(b); use logPhi to remain stable in far left tail
    lZ  <- logPhi(bi, lower = TRUE)
    Z   <- pmax(exp(lZ), .Machine$double.xmin)
    phib <- phi(bi)
    mi   <- - phib / Z                    # -(phi(b)/Phi(b))
    vi   <- 1 - (bi * phib) / Z - mi^2    # 1 + (0 - b*phi(b))/Z - m^2

    m[left] <- mi
    v[left] <- pmax(vi, 1e-12)
  }

  # -------- right-open: [a, +Inf) --------
  if (any(right)) {
    ai <- a[right]
    # Z = 1 - Phi(a); compute with upper-tail log CDF for stability
    lZ  <- logPhi(ai, lower = FALSE)
    Z   <- pmax(exp(lZ), .Machine$double.xmin)
    phia <- phi(ai)
    mi   <-  phia / Z                     #  phi(a)/(1-Phi(a))
    vi   <- 1 + (ai * phia) / Z - mi^2    # 1 + (a*phi(a) - 0)/Z - m^2

    m[right] <- mi
    v[right] <- pmax(vi, 1e-12)
  }

  list(mean = m, var = v)
}
