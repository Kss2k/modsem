## emLms with optional EMA (Quasi-Newton + Fisher Scoring) acceleration Simplified and cleaned logic: hybrid switching based on gradient norm and relative ΔLL.

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
                  ...) {
  algorithm <- toupper(match.arg(algorithm))
  data <- model$data
  stopif(anyNA(data), "Remove or replace missing values from data")

  # Default thresholds (can be overridden via em.control)
  tau1 <- if (is.null(em.control$tau1)) 1e-2 else em.control$tau1 # EM->QN if gradNorm < tau1
  tau2 <- if (is.null(em.control$tau2)) 1e-6 else em.control$tau2 # EM->QN if relΔLL < tau2
  tau3 <- if (is.null(em.control$tau3)) 1e-8 else em.control$tau3 # QN->FS if gradNorm < tau3
  tau4 <- if (is.null(em.control$tau4)) 1e-9 else em.control$tau4 # FS->stop if gradNorm < tau4

  # Initialization
  logLikNew <- -Inf
  logLikOld <- -Inf
  thetaNew  <- model$theta
  logLiks   <- NULL
  direction <- NULL

  qn_env <- new.env(parent = emptyenv())
  qn_env$LBFGS_M <- 5
  qn_env$s_list <- list()
  qn_env$g_list <- list()

  mode <- "EM" # start in EM mode
  iterations <- 0
  run <- TRUE

  while (run) {
    # New iteration
    iterations <- iterations + 1
    logLikOld <- logLikNew 
    thetaOld <- thetaNew

    # E-step
    P <- estepLms(model = model, theta = thetaOld, data = data, ...)

    # Convergence Checking
    logLikNew  <- P$obsLL
    logLiks    <- c(logLiks, logLikNew)
    deltaLL    <- logLikNew - logLikOld
    relDeltaLL <- abs(deltaLL / logLikNew)

    updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)

    if (iterations >= max.iter || 
        abs(deltaLL) < convergence.abs ||
        abs(relDeltaLL) < convergence.rel) break

    if (deltaLL < 0) {
      if (verbose) cat("\n")
      warning2("Loglikelihood is increasing!")

      if (max.step < 100) {
        message("Increasing max.step...")
        max.step <- 100
      }
    }

    # Determine Mode
    if (algorithm == "EMA") {
      previousMode <- mode

      if      (mode == "EM" && abs(relDeltaLL) < 1e-4) mode <- "QN"
      else if (mode == "QN" && abs(relDeltaLL) < tau3) mode <- "FS"
      else if (mode == "FS" && abs(relDeltaLL) < tau4) run <- FALSE

      if (mode != previousMode) { # Log mode change if verbose
        updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)
      }
    }


    if (algorithm != "EM" && mode != "EM") {
      # EMA: QN or FS update attempt
      grad <- computeGradient(theta = thetaOld, model = model, 
                              data = data, P = P, epsilon = epsilon)

      if (mode == "QN") {
        if (length(qn_env$s_list)) {
          direction <- lbfgs_two_loop(-grad, qn_env$s_list, qn_env$g_list) 
        } else direction <- -grad

      } else if (mode == "FS") {
        I_obs     <- computeFullIcom(theta = thetaOld, model = model, 
                                     data = data, P = P)
        direction <- tryCatch(solve(I_obs, grad), error = function(e) NULL)

      }

      # Line search if a direction is available
      if (!is.null(direction)) {
        alpha     <- 1 
        success   <- FALSE
        refLogLik <- logLikLms(theta = thetaOld, model = model, P = P, 
                               data = data, sign = 1)

        while (alpha > 1e-5) {
          thetaTrial  <- thetaOld + alpha * direction

          logLikTrial <- logLikLms(theta = thetaTrial, model = model, P = P,
                                   data = data, sign = 1)

          if (!is.na(logLikTrial) && logLikTrial >= refLogLik) { 
            success <- TRUE
            break 
          }
          
          alpha <- alpha / 2
        }

        if (success) {
          thetaNew <- thetaTrial
          # update BFGS memory if in QN mode

          if (mode == "QN") {
            gradNew <- computeGradient(theta = thetaNew, model = model, 
                                       data = data, P = P, epsilon = epsilon)

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
      # Plain EM M-step
      mstep <- mstepLms(model = model, P = P, data = data, theta = thetaOld,
                        max.step = max.step, epsilon = epsilon,
                        optimizer = optimizer, control = control, ...)
      thetaNew <- mstep$par
    }

  }

  if (verbose) cat("\n")
  warnif(iterations >= max.iter, "Maximum iterations reached!\n",
         "Consider a tweaking these parameters:\n", 
         formatParameters(convergence.abs, convergence.rel, algorithm, 
                          max.step, max.iter, quad.range, adaptive.quad))

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

  convergence_flag <- iterations < max.iter

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
