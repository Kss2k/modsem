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

  algorithm    <- toupper(match.arg(algorithm))
  fs.matrix    <- match.arg(fs.matrix)
  fs.fd.scheme <- match.arg(fs.fd.scheme)

  theta.lower  <- model$params$bounds$lower
  theta.upper  <- model$params$bounds$upper
  bounds.all   <- c(theta.lower, theta.upper)

  if (all(is.infinite(bounds.all))) {
    boundedTheta <- \(theta) theta # don't do anything

  } else {
    boundedTheta <- function(theta) {
      underflow <- theta < theta.lower
      overflow  <- theta > theta.upper

      theta[underflow] <- theta.lower[underflow]
      theta[overflow]  <- theta.upper[overflow]

      theta
    }
  }

  tryCatch({
    tau  <- convergence.rel
    tau1 <- if (is.null(em.control$tau1)) tau*1e6 else em.control$tau1
    tau2 <- if (is.null(em.control$tau2)) tau*2   else em.control$tau2
    tau3 <- if (is.null(em.control$tau3)) tau     else em.control$tau3

    logLikNew <- -Inf
    logLikOld <- -Inf
    thetaNew  <- model$theta
    direction <- NULL

    lastQuad     <- NULL
    adaptiveQuad <- model$models[[1L]]$quad$adaptive # fixed across groups
    adaptiveFreq <- model$models[[1L]]$quad$adaptive.frequency

    qn_env <- new.env(parent = emptyenv())
    qn_env$LBFGS_M <- 5
    qn_env$s_list <- list()
    qn_env$y_list <- list()

    mode <- "EM"
    iterations <- 0L
    run <- TRUE

    testSimpleGradient <- !model$params$gradientStruct$hasCovModel

    while (run) {
      iterations <- iterations + 1L
      logLikOld  <- logLikNew
      thetaOld   <- thetaNew
      recalcQuad <- adaptiveQuad && iterations %% adaptiveFreq == 0L

      # E-step at thetaOld
      P <- estepGsem(model = model, theta = thetaOld,
                    lastQuad = lastQuad, recalcQuad = recalcQuad,
                    adaptive.quad.tol = adaptive.quad.tol, ...)

      if (testSimpleGradient) {
        tryCatch({
          gradientCompLogLikGsem(theta = thetaNew, model = model, P = P)
        }, error = \(e) {
          warning2("Optimized computation of gradient failed! Switching gradient type.\n",
                   "Message: ", conditionMessage(e))
          model$gradientStruct$hasCovModel <<- TRUE
          model$gradientStruct$isNonLinear <<- TRUE
        })
        testSimpleGradient <- FALSE
      }

      lastQuad <- P$quad

      logLikNew  <- P$obsLL
      deltaLL    <- logLikNew - logLikOld
      relDeltaLL <- if (is.finite(logLikOld)) deltaLL / abs(logLikOld) else Inf

      updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)

      converged <- (abs(deltaLL) < convergence.abs) ||
                   (abs(relDeltaLL) < convergence.rel)
      converged.em <- converged && mode == "EM"

      if (iterations >= max.iter || converged.em) break

      if (deltaLL < -1e-8) {
        if (verbose) cat("\n")
        warning2(sprintf("Loglikelihood decreased by %.2g", deltaLL))
      }

      # EMA controller
      if (algorithm == "EMA") {
        previousMode <- mode
        dl <- abs(relDeltaLL)
        mode <- switch(mode,
          EM = if (dl < tau1) "QN" else "EM",
          QN = if (dl < tau2) "FS" else if (dl >= tau1) "EM" else "QN",
          FS = if (dl < tau3) "STOP" else if (dl >= tau2) "QN" else "FS",
          "EM"
        )

        if (converged) # converged but not in EM mode, switch to EM
          mode <- "EM"

        if (mode == "STOP") break
        if (mode != previousMode) {
          updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)
        }
      }

      # accelerated step
      if (algorithm != "EM" && mode != "EM") {
        grad <- computeGradient(theta = thetaOld, model = model, P = P, epsilon = epsilon)

        if (mode == "QN") {
          direction <- if (length(qn_env$s_list)) {
            lbfgs_two_loop(-grad, qn_env$s_list, qn_env$y_list)
          } else -grad

        } else if (mode == "FS") {
          # FS direction using Iobs = Icom - Imis (Louis), with fallback to Icom
          I_fs <- NULL
          if (fs.matrix == "Iobs") {
            I_fs <- tryCatch({
              L <- observedInfoFromLouisGsem(model = model, theta = thetaOld,
                                            P = P, recompute.P = FALSE,
                                            fd.scheme = fs.fd.scheme,
                                            fd.epsilon = fs.fd.epsilon,
                                            symmetrize = TRUE)
              L$I.obs
            }, error = function(e) NULL)
          }
          if (is.null(I_fs)) {
            I_fs <- computeFullIcom(theta = thetaOld, model = model, P = P)
          } else {
            I_fs <- 0.5 * (I_fs + t(I_fs))
          }

          # conditioning & scaled jitter
          rc <- suppressWarnings(rcond(I_fs))
          if (!is.finite(rc) || rc < 1e-10) {
            jj <- fs.jitter.mult * max(1, max(diag(I_fs), na.rm = TRUE))
            diag(I_fs) <- diag(I_fs) + jj
          }
          I_fs <- 0.5 * (I_fs + t(I_fs))
          direction <- -tryCatch(solve(I_fs, grad), error = function(e) NULL)
        }

        # line search on observed LL (and weakly on Q)
        if (!is.null(direction)) {
          alpha     <- 1
          success   <- FALSE
          refQ      <- compLogLikGsem(theta = thetaOld, model = model, P = P, sign = 1)

          while (alpha > 1e-5) {
            thetaTrial  <- boundedTheta(thetaOld + alpha * direction)
            llQTrial    <- suppressWarnings(compLogLikGsem(theta = thetaTrial,  model = model, P = P, sign = 1))
            ok <- !is.na(llQTrial) && (llQTrial >= refQ)
            if (ok) { success <- TRUE; break }
            alpha <- alpha / 2
          }

          if (success) {
            thetaNew <- thetaTrial
            if (mode == "QN") {
              # refresh P and grad at thetaNew before adding curvature pair
              P_new <- estepGsem(model = model, theta = thetaNew,
                                lastQuad = lastQuad, recalcQuad = FALSE,
                                adaptive.quad.tol = adaptive.quad.tol, ...)
              gradNew <- computeGradient(theta = thetaNew, model = model,
                                         P = P_new, epsilon = epsilon)
              s_vec <- thetaNew - thetaOld
              y_vec <- gradNew - grad
              if (sum(s_vec * y_vec) > 1e-8) {
                qn_env$s_list <- c(qn_env$s_list, list(s_vec))
                qn_env$y_list <- c(qn_env$y_list, list(y_vec))
                if (length(qn_env$s_list) > qn_env$LBFGS_M) {
                  qn_env$s_list <- qn_env$s_list[-1]
                  qn_env$y_list <- qn_env$y_list[-1]
                }
              }
              P <- P_new; lastQuad <- P_new$quad
            }
          } else {
            mode <- "EM"
          }
        } else {
          mode <- "EM"
        }
      }

      # EM M-step (plain or fallback)
      if (algorithm == "EM" || mode == "EM") {
        mstep <- mstepGsem(model = model, P = P, theta = thetaOld,
                          max.step = max.step, epsilon = epsilon,
                          optimizer = optimizer, control = control, ...)
        if (!any(is.na(mstep$par))) thetaNew <- mstep$par
      }
    } # while

    if (verbose) cat("\n")
    warnif(iterations >= max.iter, "Maximum iterations reached!\n",
           "Consider tweaking these parameters:\n",
           formatParameters(convergence.abs, convergence.rel, algorithm,
                            max.step, max.iter, nodes, adaptive.quad,
                            adaptive.quad.tol, quad.range))

    # final E-step
    P <- estepGsem(model = model, theta = thetaNew,
                  lastQuad = lastQuad, recalcQuad = FALSE,
                  adaptive.quad.tol = adaptive.quad.tol, ...)
  })
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
    Z.g      <- quad.g$n
    Zx.g     <- vector("list", length = length(Z.g))

    for (i in seq_along(Z.g)) {
      Z.gi  <- Z.g[[i]]
      Zx.gi <- cbind(matrix(0, nrow = nrow(Z.gi), ncol = NROW(OMEGA.g)), Z.gi)

      colnames(Zx.gi) <- colnames(gamma.g)

      ZETA <- Z.gi %*% A.g
      colnames(ZETA) <- colnames(psi.g)


      definedLVs     <- xis.g
      undefinedLVs   <- etas.g
      undefinedXWITH <- rownames(OMEGA.g)

      # Initialize Y.g
      Y.gi <- as.data.frame(
        lapplyNamed(xis.g, FUN = \(xi.g) alpha.g[xi.g, 1L] + ZETA[, xi.g],
                    names = xis.g)
      )

      while (length(undefinedLVs)) {
        lv.i <- undefinedLVs[[1L]]

        undefInXWITH <- apply(OMEGA.g[, undefinedLVs, drop = FALSE], MARGIN = 1L, FUN = any)
        OMEGA_DEF.g <- OMEGA.g[!undefInXWITH, , drop = FALSE]

        for (xwith in intersect(undefinedXWITH, rownames(OMEGA_DEF.g))) {
          row  <- OMEGA_DEF.g[xwith, , drop = TRUE]
          vars <- colnames(OMEGA_DEF.g)[row]
          prod <- apply(Y.gi, MARGIN = 1, FUN = prod)

          Y.gi[[xwith]]  <- prod
          Zx.gi[, xwith] <- prod

          undefinedXWITH <- setdiff(xwith, undefinedXWITH)
        }

        vals.eta <- alpha.g[lv.i, 1L] + ZETA[, lv.i] # Start with residual + intercept
        gamma.eta <- gamma.g[lv.i, definedLVs, drop = TRUE]
        for (j in seq_along(gamma.eta)) {
          indep.eta.j  <- names(gamma.eta)[[j]]
          gamma.eta.j  <- gamma.eta[[j]]

          vals.eta <- vals.eta + gamma.eta.j * Y.gi[[indep.eta.j]]
        }
        
        Y.gi[[lv.i]] <- vals.eta

        definedLVs   <- c(definedLVs, undefinedLVs[1L])
        undefinedLVs <- undefinedLVs[-1L]
      }

      Zx.g[[i]] <- Zx.gi
    }

    modFilled$models[[g]]$quad$n <- Zx.g
  }

 
  browser()
  P <- P_Step_GSEM(modFilled)
  Q <- Q_GSEM(modFilled, P = P)
  Qi <- Qi_GSEM(modFilled, P = P)

  browser()

  P
}
