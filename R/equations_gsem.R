emGsem <- function(model,
                   algorithm = c("EMA", "EM"),
                   em.control = list(),
                   verbose = FALSE,
                   convergence.abs = 1e-4,
                   convergence.rel = 1e-10,
                   max.iter = 500,
                   max.step = 1,
                   control = list(),
                   calc.se = TRUE,
                   FIM = "observed",
                   OFIM.hessian = FALSE,
                   EFIM.S = 3e4,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   epsilon = 1e-6,
                   optimizer = "nlminb",
                   R.max = 1e6,
                   adaptive.quad = FALSE,
                   quad.range = -Inf,
                   adaptive.quad.tol = 1e-12,
                   nodes = 24,
                   cr1s = TRUE,
                   # new knobs for FS/Louis
                   fs.matrix = c("Iobs","Icom"),        # use Iobs (Louis) or fall back to Icom
                   fs.fd.scheme = c("forward","central"),
                   fs.fd.epsilon = 1e-6,
                   fs.jitter.mult = sqrt(.Machine$double.eps),
                   ...) {

  theta.lower  <- model$params$bounds$lower
  theta.upper  <- model$params$bounds$upper
  bounds.all   <- c(theta.lower, theta.upper)

  tryCatch({
    logLikNew <- -Inf
    logLikOld <- -Inf
    thetaNew  <- model$theta

    iterations <- 0L
    run <- TRUE

    testSimpleGradient <- !model$params$gradientStruct$hasCovModel

    while (run) {
      iterations <- iterations + 1L
      logLikOld  <- logLikNew
      thetaOld   <- thetaNew

      # E-step at thetaOld
      P <- estepGsem(model = model, theta = thetaOld)

      logLikNew  <- P$obsLL
      deltaLL    <- logLikNew - logLikOld
      relDeltaLL <- if (is.finite(logLikOld)) deltaLL / abs(logLikOld) else Inf

      updateStatusLog(iterations, "EM", logLikNew, deltaLL, relDeltaLL, verbose)

      converged <- (abs(deltaLL) < convergence.abs) ||
                   (abs(relDeltaLL) < convergence.rel)

      if (iterations >= max.iter || converged) break

      if (deltaLL < -1e-8) {
        if (verbose) cat("\n")
        warning2(sprintf("Loglikelihood decreased by %.2g", deltaLL))
      }

      mstep <- mstepGsem(model = model, P = P, theta = thetaOld,
                         max.step = max.step, control = control, ...)
      if (!any(is.na(mstep$par))) thetaNew <- mstep$par
    }

    if (verbose) cat("\n")
    warnif(iterations >= max.iter, "Maximum iterations reached!\n",
           "Consider tweaking these parameters:\n",
           formatParameters(convergence.abs, convergence.rel, algorithm,
                            max.step, max.iter, nodes, adaptive.quad,
                            adaptive.quad.tol, quad.range))

  })


  browser()
  return(mstep)
}


estepGsem <- function(model, theta, lastQuad = NULL, recalcQuad = FALSE, adaptive.quad.tol = 1e-10) {
  modFilled <- fillModelGsem(model = model, theta = theta)

  # First step: Compute interaction term quadrature points
  # This should probably be moved inside a new funciton
  # So that we can apply it post AGHQ
  for (g in seq_along(modFilled$models)) {
    submodel <- modFilled$models[[g]]

    quad.g   <- submodel$quad
    etas.g   <- submodel$info$etas 
    xis.g    <- submodel$info$xis
    OMEGA.g  <- submodel$matrices$OMEGA
    gamma.g  <- submodel$matrices$gamma
    psi.g    <- submodel$matrices$psi
    psi.g    <- psi.g[c(xis.g, etas.g), c(xis.g, etas.g), drop = FALSE]
    alpha.g  <- submodel$matrices$alpha
    A.g      <- chol(psi.g)

    colnames_gamma <- colnames(gamma.g)
    latent_names   <- colnames(psi.g)

    # Pre-compute index maps
    idx_latent_in_gamma <- match(latent_names, colnames_gamma)
    idx_xis   <- match(xis.g, colnames_gamma)
    idx_etas  <- match(etas.g, colnames_gamma)
    xwith_names <- rownames(OMEGA.g)
    idx_xwith   <- match(xwith_names, colnames_gamma)

    alpha_full <- numeric(length(colnames_gamma))
    if (!is.null(alpha.g) && nrow(alpha.g)) {
      alpha_full[match(rownames(alpha.g), colnames_gamma)] <- alpha.g[, 1L]
    }

    gamma_mat <- as.matrix(gamma.g)

    eta_parent_idx <- vector("list", length(idx_etas))
    if (length(idx_etas)) {
      for (k in seq_along(idx_etas)) {
        row_idx <- idx_etas[k]
        parents <- which(abs(gamma_mat[row_idx, ]) > 0)
        eta_parent_idx[[k]] <- parents
      }
    }

    xwith_parent_idx <- vector("list", length(idx_xwith))
    if (length(idx_xwith)) {
      parent_names <- colnames_gamma
      for (k in seq_along(idx_xwith)) {
        row <- OMEGA.g[k, , drop = TRUE]
        xwith_parent_idx[[k]] <- which(row)
      }
    }

    row_prod <- function(mat) {
      if (!ncol(mat)) return(rep(1, nrow(mat)))
      out <- mat[, 1L]
      if (ncol(mat) > 1L) {
        for (j in 2L:ncol(mat))
          out <- out * mat[, j]
      }
      out
    }

    Z_list <- quad.g$n
    G      <- length(Z_list)
    Zx.g   <- vector("list", length = G)

    for (i in seq_len(G)) {
      Z_gi <- Z_list[[i]]
      n_i  <- nrow(Z_gi)

      ZETA <- Z_gi %*% A.g
      colnames(ZETA) <- latent_names

      Zx_gi <- matrix(0, nrow = n_i, ncol = length(colnames_gamma))
      colnames(Zx_gi) <- colnames_gamma

      if (length(idx_latent_in_gamma))
        Zx_gi[, idx_latent_in_gamma] <- Z_gi

      Yvals <- matrix(0, nrow = n_i, ncol = length(colnames_gamma))
      if (length(idx_latent_in_gamma)) {
        base <- ZETA
        base <- sweep(base, 2L, alpha_full[idx_latent_in_gamma], FUN = "+")
        Yvals[, idx_latent_in_gamma] <- base
      }

      defined <- rep(FALSE, length(colnames_gamma))
      if (length(idx_xis))
        defined[idx_xis] <- TRUE

      pending_etas <- idx_etas
      pending_xwith <- seq_along(idx_xwith)

      while (length(pending_etas) || length(pending_xwith)) {
        progress <- FALSE

        if (length(pending_xwith)) {
          to_remove <- integer(0L)
          for (k in pending_xwith) {
            target_idx <- idx_xwith[k]
            parents <- xwith_parent_idx[[k]]
            if (!length(parents)) {
              Yvals[, target_idx] <- 1
              Zx_gi[, target_idx] <- 1
              defined[target_idx] <- TRUE
              to_remove <- c(to_remove, k)
              progress <- TRUE
            } else if (all(defined[parents])) {
              prod_vals <- row_prod(Yvals[, parents, drop = FALSE])
              Yvals[, target_idx] <- prod_vals
              Zx_gi[, target_idx] <- prod_vals
              defined[target_idx] <- TRUE
              to_remove <- c(to_remove, k)
              progress <- TRUE
            }
          }
          if (length(to_remove))
            pending_xwith <- setdiff(pending_xwith, to_remove)
        }

        if (length(pending_etas)) {
          eta_idx <- pending_etas[1L]
          eta_pos <- match(eta_idx, idx_etas)
          parents <- eta_parent_idx[[eta_pos]]

          vals_eta <- Yvals[, eta_idx, drop = FALSE]
          vals_eta <- vals_eta[, 1L]
          if (length(parents))
            vals_eta <- vals_eta + as.vector(Yvals[, parents, drop = FALSE] %*% gamma_mat[eta_idx, parents])

          Yvals[, eta_idx] <- vals_eta
          defined[eta_idx] <- TRUE
          pending_etas <- pending_etas[-1L]
          progress <- TRUE
        }

        if (!progress) break
      }

      Zx.g[[i]] <- Zx_gi
    }

    modFilled$models[[g]]$quad$n <- Zx.g
  }

  QUAD <- lapply(modFilled$models, FUN = \(submod) submod$quad)
  P    <- P_Step_GSEM(modFilled, normalized = FALSE)

  density <- rowSums(P)
  P       <- P / density
  obsLL   <- sum(log(density))

  list(P = P, QUAD = QUAD, obsLL = obsLL)
}


Q_Gsem <- function(theta, model, P_Step) {
  modFilled <- fillModelGsem(model = model, theta = theta)
  # This should probably be refactored to C++ code
  for (g in seq_along(modFilled$models)) { # This should probably be refactored to C++ code
    quad.g <- P_Step$QUAD[[g]]
    modFilled$models[[g]]$quad <- quad.g
  }

  Q_GSEM(modFilled, P = P_Step$P)
}


# fdgrad <- function(x, .f, ...,  eps = 1e-5) {
#   g <- numeric(length(x))
# 
#   f0 <- .f(x, ...)
#   for (i in seq_along(x)) {
#     xi <- x
#     xi[i] <- x[i] + eps
# 
#     fi <- .f(xi, ...)
#     g[i] <- (fi - f0) / eps
#   }
# 
#   g
# }
# fastOneStepNLMINB <- function(start, objective, gradient = fdgrad, lower = NULL, upper = NULL, ...) {
#   direction <- -fdgrad(x = start, .f = objective, ...)
# 
#   f0 <- objective(start, ...)
# 
#   success <- FALSE
# 
#   if (!is.null(lower) || !is.null(upper)) {
#     if (is.null(lower)) lower <- -Inf
#     if (is.null(upper)) upper <- Inf
# 
#     truncatePars <- function(x) {
#       x[x < lower] <- lower[x < lower]
#       x[x > upper] <- upper[x > upper]
#       x
#     }
# 
#   } else truncatePars <- \(x) x
# 
#   alpha <- 1
#   while (alpha > 1e-5) {
#     test <- truncatePars(start + alpha * direction)
#     f1 <- objective(test, ...)
# 
#     if (!is.nan(f1) && !is.na(f1) && f1 < f0) {
#       success <- TRUE
#       break
#     }
# 
#     alpha <- alpha / 2L
#   }
# 
#   list(
#     objective = if (success) f1   else f0,
#     par       = if (success) test else start,
#     success   = success
#   )
# }

gradQ_Gsem <- function(theta, model, P_Step) {
  modFilled <- fillModelGsem(model = model, theta = theta)
  for (g in seq_along(modFilled$models)) { # This should probably be refactored to C++ code
    quad.g <- P_Step$QUAD[[g]]
    modFilled$models[[g]]$quad <- quad.g
  }

  Grad_Q_GSEM(modelR = modFilled, P = P_Step$P)[names(theta)]
}

mstepGsem <- function(theta, model, P_Step, max.step = 1L, control = list(),
                      lower = NULL, upper = NULL) {
  control$iter.max <- max.step

  gradient  <- \(theta) -gradQ_Gsem(theta = theta, model = model, P_Step = P_Step)
  objective <- \(theta) -Q_Gsem(theta = theta, model = model, P_Step = P_Step)

  # timeExpr(
  # fit2 <- fastOneStepNLMINB(start = theta, objective = objective,
  #                           lower = lower, upper = upper)
  # )

  fit <- nlminb(start = theta, objective = objective, control = control,
                lower = lower, upper = upper, gradient = gradient)

  fit
}
