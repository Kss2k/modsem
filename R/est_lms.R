computeGradient <- function(theta, model, P, epsilon) {
  gradientCompLogLikLms(theta = theta, model = model, P = P, sign = -1,
                        epsilon = epsilon)
}


lbfgs_two_loop <- function(grad, s_list, y_list) {
  q <- grad
  m <- length(s_list)
  if (m == 0L) return(q)

  alpha <- numeric(m)
  rho   <- numeric(m)

  for (i in rev(seq_len(m))) {
    s_i <- s_list[[i]]; y_i <- y_list[[i]]
    sy  <- sum(y_i * s_i)
    if (!is.finite(sy) || sy <= 1e-12) { alpha[i] <- 0; next }
    rho[i]   <- 1 / sy
    alpha[i] <- rho[i] * sum(s_i * q)
    q <- q - alpha[i] * y_i
  }

  s_last <- s_list[[m]]; y_last <- y_list[[m]]
  denom  <- max(sum(s_last * y_last), 1e-12)
  gamma0 <- sum(s_last * s_last) / denom
  r <- gamma0 * q

  for (i in seq_len(m)) {
    s_i <- s_list[[i]]; y_i <- y_list[[i]]
    sy  <- sum(y_i * s_i); if (!is.finite(sy) || sy <= 1e-12) next
    beta_i <- (sum(y_i * r)) / sy
    r <- r + s_i * (alpha[i] - beta_i)
  }
  r
}


computeFullIcom <- function(theta, model, P) {
  Ic <- hessianCompLogLikLms(theta = theta, model = model, P = P, sign = -1)
  0.5 * (Ic + t(Ic))
}


emLms <- function(model,
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
      P <- estepLms(model = model, theta = thetaOld,
                    lastQuad = lastQuad, recalcQuad = recalcQuad,
                    adaptive.quad.tol = adaptive.quad.tol, ...)

      if (FALSE && testSimpleGradient) {
        tryCatch({
          gradientCompLogLikLms(theta = thetaNew, model = model, P = P)
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
              L <- observedInfoFromLouisLms(model = model, theta = thetaOld,
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
          refQ      <- compLogLikLms(theta = thetaOld, model = model, P = P, sign = 1)

          while (alpha > 1e-5) {
            thetaTrial  <- boundedTheta(thetaOld + alpha * direction)
            llQTrial    <- suppressWarnings(compLogLikLms(theta = thetaTrial,  model = model, P = P, sign = 1))
            ok <- !is.na(llQTrial) && (llQTrial >= refQ)
            if (ok) { success <- TRUE; break }
            alpha <- alpha / 2
          }

          if (success) {
            thetaNew <- thetaTrial
            if (mode == "QN") {
              # refresh P and grad at thetaNew before adding curvature pair
              P_new <- estepLms(model = model, theta = thetaNew,
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
        mstep <- mstepLms(model = model, P = P, theta = thetaOld,
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
    P <- estepLms(model = model, theta = thetaNew,
                  lastQuad = lastQuad, recalcQuad = FALSE,
                  adaptive.quad.tol = adaptive.quad.tol, ...)

    finalizeModelEstimatesDA(
      model             = model,
      theta             = thetaNew,
      method            = "lms",
      data              = lapply(model$models, FUN = \(submodel) submodel$data),
      logLik            = P$obsLL,
      iterations        = iterations,
      converged         = iterations < max.iter,
      optimizer         = paste(algorithm, optimizer, sep = "-"),
      calc.se           = calc.se,
      FIM               = FIM,
      OFIM.hessian      = OFIM.hessian,
      EFIM.S            = EFIM.S,
      EFIM.parametric   = EFIM.parametric,
      robust.se         = robust.se,
      epsilon           = epsilon,
      cr1s              = cr1s,
      R.max             = R.max,
      verbose           = verbose,
      P                 = P,
      includeStartModel = TRUE,
      startModel        = model
    )

  }, error = function(e) {
    if (verbose) cat("\n")
    warning2(paste0(
      "Model estimation failed, returning starting values!\n",
      "Message: ", conditionMessage(e)
    ))
    P0 <- tryCatch(
      estepLms(model = model, theta = model$theta,
               lastQuad = NULL, recalcQuad = FALSE,
               adaptive.quad.tol = adaptive.quad.tol, ...),
      error = function(e2) NULL
    )
    ll0 <- if (!is.null(P0)) P0$obsLL else NA_real_

    finalizeModelEstimatesDA(
      model             = model,
      theta             = model$theta,
      method            = "lms",
      data              = lapply(model$models, FUN = \(submodel) submodel$data),
      logLik            = ll0,
      iterations        = 0L,
      converged         = FALSE,
      optimizer         = paste(algorithm, optimizer, sep = "-"),
      calc.se           = FALSE,
      FIM               = FIM,
      OFIM.hessian      = OFIM.hessian,
      EFIM.S            = EFIM.S,
      EFIM.parametric   = EFIM.parametric,
      robust.se         = robust.se,
      epsilon           = epsilon,
      cr1s              = cr1s,
      R.max             = R.max,
      verbose           = verbose,
      P                 = P0,
      includeStartModel = TRUE,
      startModel        = model
    )
  })
}
