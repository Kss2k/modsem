## emLms with EMA (Quasi-Newton + Fisher Scoring)
# This document outlines how to modify the existing emLms() function to implement
# an EMA (accelerated EM) procedure that switches between: EM, Quasi-Newton (QN), and
# Fisher Scoring (FS) steps based on convergence thresholds, using log-likelihood changes
# rather than parameter changes as the convergence criterion.

# 1. Define thresholds for switching modes and for final convergence
epsilon_EM_SWITCH_ll <- 1e-3    # when logLik change falls below this, switch to QN
epsilon_QN_SWITCH_ll <- 1e-6    # when logLik change falls below this, switch to FS
epsilon_FS_SWITCH_ll <- 1e-9    # when logLik change falls below this, switch back or finish
epsilon_STOP_ll      <- 1e-12   # final convergence threshold on logLik change

# 2. Helper: compute observed-data gradient
#    We assume gradientLogLikLms() already returns gradient of log-likelihood.
computeGradient <- function(theta, model, data, P, epsilon) {
  # gradientLogLikLms expects arguments: theta, model, P, sign, data, epsilon
  # sign = +1 for gradient of log-likelihood
  gradientLogLikLms(theta = theta, model = model, P = P, sign = +1,
                    data = data, epsilon = epsilon)
}

# 3. L-BFGS two-loop recursion to approximate H * gradient
#    This function is unchanged and does not rely on locked globals.
lbfgs_two_loop <- function(grad, s_list, g_list) {
  q <- grad
  m <- length(s_list)
  alpha <- numeric(m)
  rho <- rep(NA, m)
  # Compute rho_i = 1 / (y_i^T s_i)
  for (i in m:1) {
    s_i <- s_list[[i]]; g_i <- g_list[[i]]
    rho[i] <- 1 / sum(g_i * s_i)
    alpha[i] <- rho[i] * sum(s_i * q)
    q <- q - alpha[i] * g_i
  }
  # Use initial H0 = (s_{m}^T y_{m})/(y_{m}^T y_{m}) * I
  if (m > 0) {
    last_s <- s_list[[m]]; last_g <- g_list[[m]]
    gamma0 <- sum(last_s * last_g) / sum(last_g * last_g)
  } else {
    gamma0 <- 1
  }
  r <- gamma0 * q
  # Loop to compute final vector via two-loop recursion
  for (i in seq_len(m)) {
    s_i <- s_list[[i]]; g_i <- g_list[[i]]
    beta_i <- rho[i] * sum(g_i * r)
    r <- r + s_i * (alpha[i] - beta_i)
  }
  return(r)
}

# 3a. Compute complete-data Fisher information via inverse Hessian
#     of the observed-data log-likelihood (logLikLms).
#     We use numDeriv::hessian to compute the Hessian and invert it.
computeFullIcom <- function(theta, model, data, P) {
  # logLikLms returns log-likelihood given theta, model, P, data
  # We need the negative Hessian of logLikLms wrt theta
  # Use numDeriv::hessian for numeric differentiation
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for computeFullIcom.")
  }
  f_loglik <- function(par) {
    logLikLms(theta = par, model = model, P = P, data = data)
  }
  H <- numDeriv::hessian(func = f_loglik, x = theta)
  I_obs <- -H
  # Add small ridge if not invertible
  diag(I_obs) <- diag(I_obs) + 1e-8
  return(I_obs)
}

# 4. Modified emLms with EMA logic using logLik-based convergence
emLms <- function(model,
                  verbose = FALSE,
                  convergence = 1e-2,
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
                  fix.estep = TRUE,
                  R.max = 1e6,
                  ...) {
  data <- model$data
  stopif(anyNA(data), "Remove or replace missing values from data")

  # Initialization
  logLikNew  <- 0; logLikOld  <- 0; iterations <- 0
  thetaNew   <- model$theta
  bestLogLik <- -Inf; bestP <- NULL; bestTheta <- NULL
  logLiks <- numeric(0); logLikChanges <- numeric(0)

  # EMA mode flag: "EM", "QN", or "FS"
  mode <- "EM"
  # Create a local environment to store L-BFGS pairs
  qn_env <- new.env(parent = emptyenv())
  qn_env$LBFGS_M <- 5
  qn_env$s_list <- list()
  qn_env$g_list <- list()

  # Flag for whether to run E-step
  doEstep <- TRUE
  nFalseConvergence <- 0; nNegCheck <- 20; pNegCheck <- 0.5
  run <- TRUE

  while (run) {
    logLikOld <- logLikNew
    thetaOld  <- thetaNew

    if (doEstep) {
      P <- estepLms(model = model, theta = thetaOld, data = data, ...)
    }

    if (mode == "EM") {
      # Standard EM M-step
      mstep <- mstepLms(model = model, P = P, data = data, theta = thetaOld,
                        max.step = max.step, epsilon = epsilon,
                        optimizer = optimizer, control = control, ...)
      logLikNew  <- -mstep$objective
      thetaNew   <- unlist(mstep$par)
      iterations <- iterations + 1

    } else if (mode == "QN") {
      # Quasi-Newton step
      grad <- computeGradient(thetaOld, model, data, P, epsilon)
      if (length(qn_env$s_list) > 0) {
        direction <- lbfgs_two_loop(-grad, qn_env$s_list, qn_env$g_list)
      } else {
        direction <- -grad
      }
      alpha <- 1; success <- FALSE
      while (alpha > 1e-5) {
        thetaTrial <- thetaOld + alpha * direction
        trialOK <- TRUE; logLikTrial <- NA
        tryCatch({
          P_trial <- estepLms(model = model, theta = thetaTrial, data = data, ...)
          mstep_res <- mstepLms(model = model, P = P_trial, data = data,
                                theta = thetaTrial, max.step = 0,
                                epsilon = epsilon, optimizer = optimizer,
                                control = control, ...)
          logLikTrial <- -mstep_res$objective
        }, error = function(e) {
          trialOK <<- FALSE
        })
        if (trialOK && logLikTrial >= logLikOld) { success <- TRUE; break }
        alpha <- alpha / 2
      }
      if (!success) {
        mode <- "EM"; thetaNew <- thetaOld; logLikNew <- logLikOld
      } else {
        thetaNew  <- thetaTrial; logLikNew <- logLikTrial
        gradNew <- computeGradient(thetaNew, model, data, P_trial, epsilon)
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
        iterations <- iterations + 1
      }

    } else if (mode == "FS") {
      # Fisher Scoring step
      I_com_mat <- computeFullIcom(thetaOld, model, data, P)
      I_obs_approx <- I_com_mat + diag(length(thetaOld)) * 1e-8
      grad <- computeGradient(thetaOld, model, data, P, epsilon)
      direction <- tryCatch(solve(I_obs_approx, grad), error = function(e) NULL)
      if (is.null(direction)) {
        mode <- "EM"; thetaNew <- thetaOld; logLikNew <- logLikOld
      } else {
        alpha <- 1; success <- FALSE
        while (alpha > 1e-5) {
          thetaTrial <- thetaOld + alpha * direction
          trialOK <- TRUE; logLikTrial <- NA
          tryCatch({
            P_trial <- estepLms(model = model, theta = thetaTrial, data = data, ...)
            mstep_res <- mstepLms(model = model, P = P_trial, data = data,
                                  theta = thetaTrial, max.step = 0,
                                  epsilon = epsilon, optimizer = optimizer,
                                  control = control, ...)
            logLikTrial <- -mstep_res$objective
          }, error = function(e) {
            trialOK <<- FALSE
          })
          if (trialOK && logLikTrial >= logLikOld) { success <- TRUE; break }
          alpha <- alpha / 2
        }
        if (!success) {
          mode <- "EM"; thetaNew <- thetaOld; logLikNew <- logLikOld
        } else {
          thetaNew  <- thetaTrial; logLikNew <- logLikTrial
          iterations <- iterations + 1
        }
      }
    }

    # Record history
    logLiks <- c(logLiks, logLikNew)
    logLikChanges <- c(logLikChanges, logLikNew - logLikOld)
    if (logLikNew > bestLogLik) {
      bestLogLik <- logLikNew; bestP <- P; bestTheta <- thetaOld
    }

    # Compute logLik change and print
    deltaLL <- abs(logLikNew - logLikOld)
    if (verbose) {
      clearConsoleLine()
      cat(sprintf("\rIter=%d Mode=%s LogLik=%.4f ΔLL=%.4g", iterations, mode, logLikNew, deltaLL))
    }

    # Mode-switch logic based on logLik change
    if (deltaLL < epsilon_STOP_ll) {
      run <- FALSE
    } else if (mode == "EM" && deltaLL < epsilon_EM_SWITCH_ll) {
      mode <- "QN"
    } else if (mode == "QN" && deltaLL < epsilon_QN_SWITCH_ll) {
      mode <- "FS"
    } else if (mode == "FS" && deltaLL < epsilon_FS_SWITCH_ll) {
      mode <- "EM"
    }

    # Check for convergence or rescue logic using logLik
    converged_ll <- abs(logLikOld - logLikNew) < convergence
    if (iterations >= max.iter || (converged_ll && (nFalseConvergence >= 3 || doEstep))) {
      run <- FALSE
    } else if (converged_ll && !doEstep) {
      nFalseConvergence <- nFalseConvergence + 1
      doEstep <- TRUE
    } else if (
      doEstep && runningAverage(logLikChanges, n = 5) < 0 && iterations > max.iter/2 &&
      nNegativeLast(logLikChanges, n = nNegCheck) >= nNegCheck * pNegCheck
    ) {
      cat("\n")
      warning("EM algorithm is not converging. Might be at a saddle point!")
      doEstep  <- FALSE
      P        <- bestP
      thetaNew <- bestTheta
    }
  }

  if (verbose) cat("\n")
  warnif(iterations >= max.iter, "Maximum iterations reached. EM might not have converged.")
  if (!doEstep) P <- estepLms(model = model, theta = thetaNew, data = data, ...)

  # Final M-step
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

  # Calculate information matrix (I) and standard errors (SE)
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

  if (iterations == max.iter) convergence <- FALSE else convergence <- TRUE

  out <- list(model     = finalModel,
              start.model = model,
              method    = "lms",
              optimizer = paste0("EM-", optimizer),
              data      = data,
              theta     = coefficients,
              coefs     = lavCoefs,
              parTable  = parTable,
              originalParTable = model$parTable,
              logLik      = -final$objective,
              iterations  = iterations,
              convergence = convergence,
              type.se     = typeSE,
              info.quad   = getInfoQuad(model$quad),
              type.estimates = "unstandardized",
              FIM         = FIM$FIM,
              vcov        = FIM$vcov,
              information = FIM$type)

  out
}

# Note: computeFullIcom() must be implemented to compute complete-data Fisher information.
