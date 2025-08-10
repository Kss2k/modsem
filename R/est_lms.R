## emLms with optional EMA (Quasi-Newton + Fisher Scoring) acceleration Simplified and cleaned logic: hybrid switching based on gradient norm and relative ΔLL.


# Compute complete-data gradient
computeGradient <- function(theta, model, P, data, epsilon) {
  gradientCompLogLikLms(theta = theta, model = model, P = P, sign = -1,
                        epsilon = epsilon, data = data)
}


# L-BFGS two-loop recursion to approximate H * gradient
lbfgs_two_loop <- function(grad, s_list, g_list) {
  q <- grad; m <- length(s_list)
  alpha <- numeric(m);
  rho <- rep(NA, m)
  for (i in m:1) {
    s_i <- s_list[[i]]; g_i <- g_list[[i]]
    rho[i] <- 1 / sum(g_i * s_i)
    alpha[i] <- rho[i] * sum(s_i * q)
    q <- q - alpha[i] * g_i
  }
  gamma0 <- if (m > 0) {
    last_s <- s_list[[m]]; last_g <- g_list[[m]]
    sum(last_s * last_g) / sum(last_g * last_g)
  } else 1
  r <- gamma0 * q
  for (i in seq_len(m)) {
    s_i <- s_list[[i]]; g_i <- g_list[[i]]
    beta_i <- rho[i] * sum(g_i * r)
    r <- r + s_i * (alpha[i] - beta_i)
  }
  r
}


# Compute complete-data Fisher information via inverse Hessian
computeFullIcom <- function(theta, model, data, P) {
  # fLogLik <- function(par) compLogLikLms(theta = par, model = model, P = P)
  # Ic <- fdHESS(pars = theta, fun = fLogLik)

  Ic <- hessianCompLogLikLms(theta = theta, model = model, P = P,
                             sign = -1, data = data)
  diag(Ic) <- diag(Ic) + 1e-8  # ensure invertibility
  Ic
}


updateStatusLog <- function(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose = FALSE) {
  if (verbose) {
    clearConsoleLine()
    printf("\rIter=%d Mode=%s LogLik=%.2f \u0394LL=%.2g rel\u0394LL=%.2g",
           iterations, mode, logLikNew, deltaLL, relDeltaLL)
  }
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
                  ...) {

  algorithm <- toupper(algorithm)
  algorithm <- match.arg(algorithm)
  data <- model$data

  tryCatch({

    # thresholds (overridable)
    tau  <- convergence.rel
    tau1 <- if (is.null(em.control$tau1)) tau*1e6 else em.control$tau1
    tau2 <- if (is.null(em.control$tau2)) tau*1e1 else em.control$tau2
    tau3 <- if (is.null(em.control$tau3)) tau*2   else em.control$tau3

    logLikNew <- -Inf
    logLikOld <- -Inf
    thetaNew  <- model$theta
    logLiks   <- NULL
    direction <- NULL

    lastQuad     <- NULL
    adaptiveQuad <- model$quad$adaptive
    adaptiveFreq <- model$quad$adaptive.frequency

    qn_env <- new.env(parent = emptyenv())
    qn_env$LBFGS_M <- 5
    qn_env$s_list <- list()
    qn_env$g_list <- list()

    mode <- "EM"
    iterations <- 0
    run <- TRUE

    testSimpleGradient <- !model$gradientStruct$hasCovModel

    while (run) {
      iterations <- iterations + 1
      logLikOld  <- logLikNew
      thetaOld   <- thetaNew
      recalcQuad <- adaptiveQuad && iterations %% adaptiveFreq == 0

      # E-step
      P <- estepLms(model = model, theta = thetaOld, data = data,
                    lastQuad = lastQuad, recalcQuad = recalcQuad,
                    adaptive.quad.tol = adaptive.quad.tol, ...)

      if (testSimpleGradient) {
        tryCatch({
          gradientCompLogLikLms(theta = thetaNew, model = model, P = P, data = data)
        }, error = \(e) {
          warning2("Optimized computation of gradient failed! Switching gradient type.")
          model$gradientStruct$hasCovModel <<- TRUE
          model$gradientStruct$isNonLinear <<- TRUE
        })
        testSimpleGradient <- FALSE
      }

      lastQuad <- P$quad

      logLikNew  <- P$obsLL
      logLiks    <- c(logLiks, logLikNew)
      deltaLL    <- logLikNew - logLikOld
      relDeltaLL <- ifelse(is.finite(deltaLL), deltaLL / abs(logLikOld), Inf)

      updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)

      converged <- abs(deltaLL) < convergence.abs ||
                   abs(relDeltaLL) < convergence.rel
      converged.em <- converged && mode == "EM" # optimallly we want to finish in EM mode

      if (iterations >= max.iter || converged.em) break

      if (deltaLL < 0) {
        if (verbose) cat("\n")
        warning2("Loglikelihood is decreasing!")
      }

      if (algorithm == "EMA") {
        previousMode <- mode
        dl <- abs(relDeltaLL)

        mode <- switch(mode,
          EM = if (dl < tau1) "QN",
          QN = if (dl < tau2) "FS"   else if (dl >= tau1) "EM",
          FS = if (dl < tau3) "STOP" else if (dl >= tau2) "QN",
          NULL
        )

        if (converged)
          mode <- "EM" # override to finish in EM mode

        mode <- ifelse(is.null(mode), previousMode, mode)

        if (mode == "STOP") { run <- FALSE; break }
        if (mode != previousMode) {
          updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)
        }
      }

      if (algorithm != "EM" && mode != "EM") {
        grad <- computeGradient(theta = thetaOld, model = model, data = data,
                                P = P, epsilon = epsilon)

        if (mode == "QN") {
          if (length(qn_env$s_list)) {
            direction <- lbfgs_two_loop(-grad, qn_env$s_list, qn_env$g_list)
          } else direction <- -grad

        } else if (mode == "FS") {
          Ic <- computeFullIcom(theta = thetaOld, model = model, P = P, data = data)
          direction <- -tryCatch(solve(Ic, grad), error = function(e) NULL)
        }

        if (!is.null(direction)) {
          alpha     <- 1
          success   <- FALSE
          refLogLik <- compLogLikLms(theta = thetaOld, model = model, P = P,
                                     data = data, sign = 1)

          while (alpha > 1e-5) {
            thetaTrial  <- thetaOld + alpha * direction
            logLikTrial <- suppressWarnings({
              compLogLikLms(theta = thetaTrial, model = model, P = P, data = data, sign = 1)
            })
            if (!is.na(logLikTrial) && logLikTrial >= refLogLik) {
              success <- TRUE; break
            }
            alpha <- alpha / 2
          }

          if (success) {
            thetaNew <- thetaTrial
            if (mode == "QN") {
              gradNew <- computeGradient(theta = thetaNew, model = model, data = data, P = P, epsilon = epsilon)
              s_vec <- thetaNew - thetaOld
              y_vec <- gradNew - grad
              if (sum(s_vec * y_vec) > 1e-8) {
                qn_env$s_list <- c(qn_env$s_list, list(s_vec))
                qn_env$g_list <- c(qn_env$g_list, list(y_vec))
                if (length(qn_env$s_list) > qn_env$LBFGS_M) {
                  qn_env$s_list <- qn_env$s_list[-1]
                  qn_env$g_list <- qn_env$g_list[-1]
                }
              }
            }
          } else mode <- "EM"
        } else mode <- "EM"
      }

      if (algorithm == "EM" || mode == "EM") {
        mstep <- mstepLms(model = model, P = P, theta = thetaOld,
                          max.step = max.step, epsilon = epsilon,
                          data = data, optimizer = optimizer, control = control, ...)
        thetaNew <- mstep$par
      }
    }

    if (verbose) cat("\n")
    warnif(iterations >= max.iter, "Maximum iterations reached!\n",
           "Consider tweaking these parameters:\n",
           formatParameters(convergence.abs, convergence.rel, algorithm,
                            max.step, max.iter, nodes, adaptive.quad,
                            adaptive.quad.tol, quad.range))

    # one last E-step to ensure P/logLik are consistent with thetaNew
    P <- estepLms(model = model, theta = thetaNew, data = data,
                  lastQuad = lastQuad, recalcQuad = FALSE,
                  adaptive.quad.tol = adaptive.quad.tol, ...)

    finalizeModelEstimatesDA(
      model             = model,
      theta             = thetaNew,
      method            = "lms",
      data              = data,
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
    # Fallback: return starting parameters without SEs, convergence = FALSE
    if (verbose) cat("\n")

    warning2(paste0(
      "Model estimation failed, returning starting values!\n",
      "Message: ", conditionMessage(e)
    ))

    P0 <- tryCatch(
      estepLms(model = model, theta = model$theta, data = data,
               lastQuad = NULL, recalcQuad = FALSE,
               adaptive.quad.tol = adaptive.quad.tol, ...),
      error = function(e2) NULL
    )
    ll0 <- if (!is.null(P0)) P0$obsLL else NA_real_

    finalizeModelEstimatesDA(
      model             = model,
      theta             = model$theta,     # start values
      method            = "lms",
      data              = data,
      logLik            = ll0,             # best effort; NA if unavailable
      iterations        = 0L,
      converged         = FALSE,           # explicit
      optimizer         = paste(algorithm, optimizer, sep = "-"),
      calc.se           = FALSE,           # <- force no standard errors
      FIM               = FIM,             # ignored when calc.se = FALSE
      OFIM.hessian      = OFIM.hessian,
      EFIM.S            = EFIM.S,
      EFIM.parametric   = EFIM.parametric,
      robust.se         = robust.se,
      epsilon           = epsilon,
      cr1s              = cr1s,
      R.max             = R.max,
      verbose           = verbose,
      P                 = P0,              # may be NULL
      includeStartModel = TRUE,
      startModel        = model
    )
  })
}
