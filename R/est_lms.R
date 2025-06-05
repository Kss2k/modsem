## emLms with optional EMA (Quasi-Newton + Fisher Scoring) acceleration
# Simplified and cleaned logic: hybrid switching based on gradient norm and relative ΔLL.

# 1. Helper: compute observed-data gradient
computeGradient <- function(theta, model, data, P, epsilon) {
  gradientLogLikLms(theta = theta, model = model, P = P, sign = -1,
                    data = data, epsilon = epsilon)
}

# 2. L-BFGS two-loop recursion to approximate H * gradient
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

# 3. Compute complete-data Fisher information via inverse Hessian
# (requires fdHESS or numDeriv)
computeFullIcom <- function(theta, model, data, P) {
  fLogLik <- function(par) logLikLms(theta = par, model = model, P = P, data = data)
  H <- fdHESS(pars = theta, fun = fLogLik)
  I_obs <- H
  diag(I_obs) <- diag(I_obs) + 1e-8  # ensure invertibility
  I_obs
}


updateStatusLog <- function(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose = FALSE) {
  if (verbose) {
    clearConsoleLine()
    cat(sprintf("\rIter=%d Mode=%s LogLik=%.2f \u0394LL=%.2g rel\u0394LL=%.2g",
                iterations, mode, logLikNew, deltaLL, relDeltaLL))
  }
}


# 5. emLms function
emLms <- function(model,
                  algorithm = c("EMA", "EM"),    # "EM" or "EMA"
                  em.control = list(),              # overrides for tau thresholds
                  verbose = FALSE,
                  convergence = 1e-2,
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
                  ...) {

  algorithm <- toupper(match.arg(algorithm))
  data <- model$data
  stopif(anyNA(data), "Remove or replace missing values from data")

  # Default thresholds (can be overridden via em.control)
  tau1 <- if (is.null(em.control$tau1)) 1e-2 else em.control$tau1 # EM→QN if gradNorm < tau1
  tau2 <- if (is.null(em.control$tau2)) 1e-6 else em.control$tau2 # EM→QN if relΔLL < tau2
  tau3 <- if (is.null(em.control$tau3)) 1e-8 else em.control$tau3 # QN→FS if gradNorm < tau3
  tau4 <- if (is.null(em.control$tau4)) 1e-9 else em.control$tau4 # FS→stop if gradNorm < tau4

  # Initialization
  logLikNew <- 0
  logLikOld <- 0
  warnedConvergence <- FALSE
  thetaNew  <- model$theta
  bestLogLik <- -Inf
  bestP <- NULL 
  bestTheta <- NULL
  logLikChanges <- NULL
  logLiks <- NULL

  qn_env <- new.env(parent = emptyenv())
  qn_env$LBFGS_M <- 5
  qn_env$s_list <- list()
  qn_env$g_list <- list()

  mode <- "EM"                # start in EM mode
  iterations <- 0
  run <- TRUE

  while (run) {
    logLikOld <- logLikNew; thetaOld <- thetaNew

    # 1. E-step
    P <- estepLms(model = model, theta = thetaOld, data = data, ...)

    # 2. Decide update based on mode and algorithm
    if (algorithm == "EM" || mode == "EM") {
      # Plain EM M-step
      mstep <- mstepLms(model = model, P = P, data = data, theta = thetaOld,
                        max.step = max.step, epsilon = epsilon,
                        optimizer = optimizer, control = control, ...)
      thetaNew <- unlist(mstep$par)
      logLikNew <- -mstep$objective
      iterations <- iterations + 1

    } else {
      # EMA: QN or FS update attempt
      # 2a. Compute gradient once
      grad <- computeGradient(thetaOld, model, data, P, epsilon)
      if (mode == "QN") {
        direction <- if (length(qn_env$s_list) > 0) 
          lbfgs_two_loop(-grad, qn_env$s_list, qn_env$g_list) else -grad
      } else if (mode == "FS") {
        I_obs <- computeFullIcom(thetaOld, model, data, P)
        direction <- tryCatch(solve(I_obs, grad), error = function(e) NULL)
      }

      # 2b. Line search if a direction is available
      if (exists("direction") && !is.null(direction)) {
        alpha <- 1; success <- FALSE
        refLogLik <- logLikOld

        while (alpha > 1e-5) {
          thetaTrial <- thetaOld + alpha * direction
          trialOK <- TRUE; logLikTrial <- NA

          tryCatch({
            P_trial <- estepLms(model = model, theta = thetaTrial, data = data, ...)
            logLikTrial <- logLikLms(theta = thetaTrial, model = model, P = P_trial, 
                                     data = data, sign = 1)

          }, error = function(e) trialOK <<- FALSE)

          if (trialOK && logLikTrial >= refLogLik) { 
            success <- TRUE
            break 
          }
          
          alpha <- alpha / 2
        }

        if (success) {
          thetaNew <- thetaTrial; logLikNew <- logLikTrial
          iterations <- iterations + 1
          # update BFGS memory if in QN mode
          if (mode == "QN") {
            gradNew <- computeGradient(thetaNew, model, data, P_trial, epsilon)
            s_vec <- thetaNew - thetaOld; y_vec <- gradNew - grad

            if (sum(s_vec * y_vec) > 1e-8) {
              qn_env$s_list <- c(qn_env$s_list, list(s_vec))
              qn_env$g_list <- c(qn_env$g_list, list(y_vec))

              if (length(qn_env$s_list) > qn_env$LBFGS_M) {
                qn_env$s_list <- qn_env$s_list[-1]
                qn_env$g_list <- qn_env$g_list[-1]
              }
            }
          }
        } else {
          # Fallback to one EM step on failure
          mode <- "EM"
          mstep <- mstepLms(model = model, P = P, data = data, theta = thetaOld,
                            max.step = max.step, epsilon = epsilon,
                            optimizer = optimizer, control = control, ...)
          thetaNew <- unlist(mstep$par)
          logLikNew <- -mstep$objective
          iterations <- iterations + 1
        }
      } else {
        # No valid direction: fallback to EM
        mode <- "EM"
        mstep <- mstepLms(model = model, P = P, data = data, theta = thetaOld,
                          max.step = max.step, epsilon = epsilon,
                          optimizer = optimizer, control = control, ...)
        thetaNew <- unlist(mstep$par)
        logLikNew <- -mstep$objective
        iterations <- iterations + 1
      }
    }


    # Monitor progress
    deltaLL <- abs(logLikNew - logLikOld)
    relDeltaLL <- if (abs(logLikNew) > 0) deltaLL / abs(logLikNew) else deltaLL

    updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)

    # 5. Mode-switch logic when EMA is active
    if (algorithm == "EMA") {
      # if (mode == "EM" && (gradNorm < tau1 || relDeltaLL < tau2)) {
      if (mode == "EM" && relDeltaLL < tau2) {
        mode <- "QN"
      } else if (mode == "QN" && relDeltaLL < tau3) {
        mode <- "FS"
      } else if (mode == "FS" && relDeltaLL < tau4) {
        run <- FALSE
      }
      # Log mode change if verbose
      # updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, gradNorm, verbose)
      updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)
    }

    # 6. Termination checks
    if (iterations >= max.iter || deltaLL < convergence || relDeltaLL < convergence.rel) break

    # 7. Warning if stuck
    if (runningAverage(logLikChanges, n = 10) < 0 && iterations > 100 && !warnedConvergence) {

      if (verbose) cat("\n")
      warning2(
        "EM appears stuck: consider tweaking these parameters:\n",
        formatParameters(convergence, algorithm, max.step, quad.range, 
                         adaptive.quad) 
      )

      if (max.step < 100) {
        warning2("Increasing max.step...")
        max.step <- 100
      }
      
      warnedConvergence <- TRUE
    }

    # Record history
    logLiks <- c(logLiks, logLikNew)
    logLikChanges <- c(logLikChanges, logLikNew - logLikOld)
    if (logLikNew > bestLogLik) {
      bestLogLik <- logLikNew; bestP <- P; bestTheta <- thetaOld
    }
  }

  if (verbose) cat("\n")
  warnif(iterations >= max.iter, "Maximum iterations reached!\n",
         "Consider a tweaking these parameters:\n", 
         formatParameters(convergence, algorithm, max.step, 
                          max.iter, quad.range, adaptive.quad))

  # Final E- and M-step for output
  P <- estepLms(model = model, theta = thetaNew, data = data, ...)
  final <- mstepLms(model = model, P = P, data = data,
                    theta = thetaNew, max.step = max.step,
                    epsilon = epsilon, optimizer = optimizer,
                    verbose = verbose, control = control, ...)

  coefficients <- final$par
  lavCoefs     <- getLavCoefs(model = model, theta = coefficients, method = "lms")
  finalModel   <- fillModel(model, coefficients, fillPhi = TRUE, method = "lms")
  info         <- model$info

  emptyModel <- getEmptyModel(parTable = model$parTable,
                              cov.syntax = model$cov.syntax,
                              parTableCovModel = model$covModel$parTable,
                              method = "lms")
  finalModel$matricesNA <- emptyModel$matrices
  finalModel$covModelNA <- emptyModel$covModel

  # Compute standard errors
  typeSE <- ifelse(!calc.se, "none", ifelse(robust.se, "robust", "standard"))
  FIM <- calcFIM_da(model = model, finalModel = finalModel, theta = coefficients,
                    data = data, method = "lms", EFIM.S = EFIM.S,
                    hessian = OFIM.hessian, calc.se = calc.se,
                    EFIM.parametric = EFIM.parametric, verbose = verbose,
                    FIM = FIM, robust.se = robust.se, epsilon = epsilon,
                    R.max = R.max, NA__ = -999)
  SE <- calcSE_da(calc.se = calc.se, FIM$vcov, rawLabels = FIM$raw.labels,
                  NA__ = -999)
  modelSE <- getSE_Model(model, se = SE, method = "lms",
                         n.additions = FIM$n.additions)
  finalModel$matricesSE <- modelSE$matrices
  finalModel$covModelSE <- modelSE$covModel

  parTable <- modelToParTable(finalModel, coefs = lavCoefs,
                              se = SE, method = "lms")
  parTable$z.value  <- parTable$est / parTable$std.error
  parTable$p.value  <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - CI_WIDTH * parTable$std.error
  parTable$ci.upper <- parTable$est + CI_WIDTH * parTable$std.error

  convergence_flag <- ifelse(iterations < max.iter & deltaLL >= convergence, TRUE, deltaLL < convergence)
  out <- list(
    model         = finalModel,
    start.model   = model,
    method        = "lms",
    optimizer     = paste(algorithm, optimizer, sep = "-"),
    data          = data,
    theta         = coefficients,
    coefs         = lavCoefs,
    parTable      = parTable,
    originalParTable = model$parTable,
    logLik        = -final$objective,
    iterations    = iterations,
    convergence   = convergence_flag,
    type.se       = typeSE,
    info.quad     = getInfoQuad(model$quad),
    type.estimates = "unstandardized",
    FIM           = FIM$FIM,
    vcov          = FIM$vcov,
    information   = FIM$type
  )
  out
}
