## emLms with optional EMA (Quasi-Newton + Fisher Scoring) acceleration
# Adds `algorithm` argument: "EM" or "EMA" (case-insensitive), and
# `em.control` list for thresholds (overridable).

# 2. Helper: compute observed-data gradient
computeGradient <- function(theta, model, data, P, epsilon) {
  gradientLogLikLms(theta = theta, model = model, P = P, sign = -1,
                    data = data, epsilon = epsilon)
}

# 3. L-BFGS two-loop recursion to approximate H * gradient
lbfgs_two_loop <- function(grad, s_list, g_list) {
  q <- grad; m <- length(s_list)
  alpha <- numeric(m); rho <- rep(NA, m)
  for (i in m:1) {
    s_i <- s_list[[i]]; g_i <- g_list[[i]]; rho[i] <- 1 / sum(g_i * s_i)
    alpha[i] <- rho[i] * sum(s_i * q); q <- q - alpha[i] * g_i
  }
  if (m > 0) {
    last_s <- s_list[[m]]; last_g <- g_list[[m]]
    gamma0 <- sum(last_s * last_g) / sum(last_g * last_g)
  } else {
    gamma0 <- 1
  }
  r <- gamma0 * q
  for (i in seq_len(m)) {
    s_i <- s_list[[i]]; g_i <- g_list[[i]]
    beta_i <- rho[i] * sum(g_i * r)
    r <- r + s_i * (alpha[i] - beta_i)
  }
  r
}


# 3a. Compute complete-data Fisher information via inverse Hessian (logLikLms)
computeFullIcom <- function(theta, model, data, P) {
  fLogLik <- function(par) logLikLms(theta = par, model = model, P = P, data = data, sign = -1)
  H <- fdHESS(pars=theta, fun=fLogLik)
  I_obs <- -H; diag(I_obs) <- diag(I_obs) + 1e-8
  I_obs
}


updateStatusLog <- function(iterations, mode, logLikNew, deltaLL, verbose = FALSE) {
  if (verbose) {
    clearConsoleLine()
    cat(sprintf("\rIter=%d Mode=%s LogLik=%.2f \u0394LL=%.2g", iterations, mode, logLikNew, deltaLL))
  }
}


# 4. emLms with optional EMA acceleration and overridable thresholds
emLms <- function(model,
                  algorithm = c("EMA", "EM"),  # choose algorithm
                  em.control = list(),             # list of threshold overrides
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
                  R.max = 1e6,
                  ...) {
  algorithm <- toupper(match.arg(algorithm))
  data <- model$data
  stopif(anyNA(data), "Remove or replace missing values from data")

  # Set default thresholds and override with em.control entries
  thr <- list(
    epsilon_EM_SWITCH_ll = 1e-2, # convergence * 1e3,  # 1e-3,
    epsilon_QN_SWITCH_ll = 1e-4, # convergence * 1e2, # 1e-6,
    epsilon_FS_SWITCH_ll = 1e-5, # convergence * 1e1, # 1e-9,
    epsilon_STOP_ll      = convergence        # 1e-12
  )

  for (nm in intersect(names(em.control), names(thr))) thr[[nm]] <- em.control[[nm]]

  # Initialization
  logLikNew  <- 0; logLikOld  <- 0; iterations <- 0
  thetaNew   <- model$theta
  bestLogLik <- -Inf; bestP <- NULL; bestTheta <- NULL
  logLiks <- numeric(0); logLikChanges <- numeric(0)

  # Setup for EMA
  mode <- "EM"
  qn_env <- new.env(parent = emptyenv()); qn_env$LBFGS_M <- 5
  qn_env$s_list <- list(); qn_env$g_list <- list()
  nNegCheck <- 20; pNegCheck <- 0.5
  run <- TRUE; warnedStuck <- FALSE

  while (run) {
    logLikOld <- logLikNew; thetaOld <- thetaNew

    P <- estepLms(model = model, theta = thetaOld, data = data, ...)
    mstep <- mstepLms(model = model, P = P, data = data, theta = thetaOld,
                      max.step = max.step, epsilon = epsilon,
                      optimizer = optimizer, control = control, ...)

    if (algorithm != "EM" && mode != "EM") {
      # EMA logic: compute search direction based on mode
      direction <- NULL

      if (mode == "QN") {
        grad <- computeGradient(thetaOld, model, data, P, epsilon)

        if (length(qn_env$s_list) > 0) 
          direction <- lbfgs_two_loop(-grad, qn_env$s_list, qn_env$g_list)
        else direction <- -grad

      } else if (mode == "FS") {
        grad <- computeGradient(thetaOld, model, data, P, epsilon)

        I_com_mat <- computeFullIcom(thetaOld, model, data, P)
        I_obs_approx <- I_com_mat + diag(length(thetaOld)) * 1e-8
        direction <- tryCatch(solve(I_obs_approx, grad), error = function(e) NULL)
      }
      
      if (is.null(direction)) {
        nextStep <- "fallback" 
        mode <- "EM"
      } else nextStep <- "search"

      if (!is.null(nextStep) && nextStep == "search") {
        # Line search for both QN and FS
        alpha <- 1; success <- FALSE; logLikEM <- -mstep$objective

        while (alpha > 1e-5) {
          thetaTrial <- thetaOld + alpha * direction
          trialOK <- TRUE 
          logLikTrial <- NA

          tryCatch({
            P_trial <- estepLms(model = model, theta = thetaTrial, data = data, ...)
            mstep_res <- mstepLms(model = model, P = P_trial, data = data,
                                  theta = thetaTrial, max.step = 0,
                                  epsilon = epsilon, optimizer = optimizer,
                                  control = control, ...)
            logLikTrial <- -mstep_res$objective
          }, error = function(e) trialOK <<- FALSE)
          

          if (trialOK && logLikTrial >= logLikEM) { 
            success <- TRUE
            break 
          }

          alpha <- alpha / 2
        }

        if (success) {
          thetaNew  <- thetaTrial
          logLikNew <- logLikTrial
          iterations <- iterations + 1
          if (mode == "QN") {
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
          }
        } else mode <- "EM"
      }
    }

    if (mode == "EM") {
      # Standard EM step
      thetaNew <- mstep$par
      logLikNew <- -mstep$objective
      iterations <- iterations + 1
    }
   

    # Plain EM convergence check
    deltaLL <- abs(logLikNew - logLikOld)
    updateStatusLog(iterations, mode, logLikNew, deltaLL, verbose)

    # Mode-switch logic based solely on deltaLL
    if (algorithm == "EMA") {
      if (deltaLL >= thr$epsilon_EM_SWITCH_ll) {
        mode <- "EM"
      } else if (deltaLL >= thr$epsilon_QN_SWITCH_ll) {
        mode <- "QN"
      } else if (mode == "FS" && deltaLL < thr$epsilon_FS_SWITCH_ll) {
        mode <- "EM"
      } else mode <- "FS"
   
      # update status log, in case mode has changed
      updateStatusLog(iterations, mode, logLikNew, deltaLL, verbose)
    }

    # Simplified rescue logic
    if (iterations >= max.iter || deltaLL < convergence) {
      run <- FALSE
      break

    } else if (
      runningAverage(logLikChanges, n = 5) < 0 && iterations > max.iter/2 &&
      nNegativeLast(logLikChanges, n = nNegCheck) >= nNegCheck * pNegCheck &&
      !warnedStuck
    ) {
      warning2("EM algorithm appears stuck: consider using a larger `convergence` tolerance.",
               "Consider using a larger `convergence` tolerance!")
      warnedStuck <- TRUE
    }
  }

  if (verbose) cat("\n")
  warnif(iterations >= max.iter, "Maximum iterations reached. EM might not have converged.\n",
         "Consider using a larger `convergence` tolerance!")
  P <- estepLms(model = model, theta = thetaNew, data = data, ...)

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
              optimizer = paste(algorithm, optimizer, sep="-"),
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
