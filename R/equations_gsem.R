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

      updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)

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

    ## Pack quad.g$n (list of matrices) into a 3D array
    Z_list   <- quad.g$n
    G        <- length(Z_list)
    rcz      <- dim(Z_list[[1L]])  # c(r, c)
    r        <- rcz[1L]
    cZ       <- rcz[2L]

    # Combine without names overhead
    Z_arr <- array(
                   data = unlist(Z_list, use.names = FALSE),
                   dim  = c(r, cZ, G)
    )

    ## Batch compute ZETA for all slices: r*G by cZ times A.g
    # ZETA_arr[,,i] = Z_arr[,,i] %*% A.g
    ZETA_mat <- matrix(Z_arr, nrow = r * G, ncol = cZ) %*% A.g   # (rG × cZ) %*% (cZ × p)
    p        <- ncol(psi.g)
    ZETA_arr <- array(ZETA_mat, dim = c(r, p, G))
    colnames_ZETA <- colnames(psi.g)  # reuse once

    ## Prepare constants reused in the inner loop
    colnames_gamma <- colnames(gamma.g)
    undef_all_XWITH <- rownames(OMEGA.g)

    Zx.g <- vector("list", length = G)

    for (i in seq_len(G)) {
      # pull slice i once
      Z_gi      <- Z_arr[,,i]
      ZETA      <- ZETA_arr[,,i, drop = FALSE][,,1L]
      colnames(ZETA) <- colnames_ZETA

      # Zx.gi scaffold: zeros for XWITH columns, then original Z_gi
      Zx.gi <- cbind(matrix(0, nrow = r, ncol = nrow(OMEGA.g)), Z_gi)
      colnames(Zx.gi) <- colnames_gamma

      definedLVs     <- xis.g
      undefinedLVs   <- etas.g
      undefinedXWITH <- undef_all_XWITH

      # Initialize Y.gi with xis: alpha + ZETA[, xi]
      # (build as data.frame with named columns in xis order)
      Y.gi <- as.data.frame(stats::setNames(
         lapply(xis.g, \(xi.g) alpha.g[xi.g, 1L] + ZETA[, xi.g]),
         xis.g
      ))

      while (length(undefinedLVs)) {
        lv.i <- undefinedLVs[[1L]]

        undefInXWITH <- apply(OMEGA.g[, undefinedLVs, drop = FALSE], 1L, any)
        OMEGA_DEF.g  <- OMEGA.g[!undefInXWITH, , drop = FALSE]

        # Resolve XWITH rows that now have only defined LVs
        if (nrow(OMEGA_DEF.g)) {
          for (xwith in intersect(undefinedXWITH, rownames(OMEGA_DEF.g))) {
            row  <- OMEGA_DEF.g[xwith, , drop = TRUE]
            vars <- names(row)[row]

            # product over the specific vars (more precise & avoids unnecessary cols)
            # Base R row-wise product across selected columns:
            prod <- Reduce(`*`, Y.gi[vars])

            Y.gi[[xwith]]  <- prod
            Zx.gi[, xwith] <- prod

            # mark as defined
            undefinedXWITH <- setdiff(undefinedXWITH, xwith)
          }
        }

        # latent lv.i = alpha + residual(ZETA) + gamma * parents
        vals.eta  <- alpha.g[lv.i, 1L] + ZETA[, lv.i]
        gamma.eta <- gamma.g[lv.i, definedLVs, drop = TRUE]
        if (length(gamma.eta)) {
          # accumulate effects from already-defined parents
          for (j in seq_along(gamma.eta)) {
            indep.eta.j <- names(gamma.eta)[[j]]
            vals.eta    <- vals.eta + gamma.eta[[j]] * Y.gi[[indep.eta.j]]
          }
        }

        Y.gi[[lv.i]] <- vals.eta

        definedLVs   <- c(definedLVs, lv.i)
        undefinedLVs <- undefinedLVs[-1L]
      }

      Zx.g[[i]] <- Zx.gi
    }

    modFilled$models[[g]]$quad$n <- Zx.g
  }

  QUAD <- lapply(modFilled$models, FUN = \(submod) submod$quad)
  browser()
  P    <- P_Step_GSEM(modFilled, normalized = FALSE)

  density <- rowSums(P)
  P       <- P / density
  obsLL   <- sum(log(density))

  browser()
  list(P = P, QUAD = QUAD, obsLL = sum(log(P)))
}


Q_Gsem <- function(theta, model, P_Step) {
  modFilled <- fillModelGsem(model = model, theta = theta)
  # This should probably be refactored to C++ code
  for (g in seq_along(modFilled)) { # This should probably be refactored to C++ code
    quad.g <- P_Step$QUAD[[g]]
    modFilled$models[[g]] <- quad.g
  }
}


mstepGsem <- function(theta, model, P_Step) {


}
