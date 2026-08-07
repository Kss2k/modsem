computeGradient <- function(theta, model, P, epsilon) {
  gradientObsLogLikLms(theta = theta, model = model, P = P, sign = -1,
                       epsilon = epsilon)
}


computeFullIobs <- function(theta, model, P, louis = TRUE) {

  if (louis) {
    louisInfo <- observedInfoFromLouisLms(
      model = model,
      theta = theta,
      P = P
    )

    I <- louisInfo$I.obs

  } else {

    # negative hessian (sign = -1)
    I <- hessianObsLogLikLms(
      theta = theta,
      model = model,
      P     = P,
      sign  = -1
    )

  }

  # make sure it's symmetric
  0.5 * (I + t(I))
}


lmsIterationHistory <- function() {
  data.frame(
    iteration = integer(),
    expected.remaining.iterations = numeric(),
    expected.final.loglik = numeric()
  )
}


normalizeLmsLogLikHistory <- function(history) {
  if (is.null(history)) {
    return(data.frame(
      iteration = integer(),
      logLik = numeric(),
      deltaLL = numeric(),
      relDeltaLL = numeric()
    ))
  }

  if (is.numeric(history) && is.null(dim(history))) {
    logLik <- as.numeric(history)
    if (!length(logLik)) {
      return(data.frame(
        iteration = integer(),
        logLik = numeric(),
        deltaLL = numeric(),
        relDeltaLL = numeric()
      ))
    }
    iteration <- seq_along(logLik)
    deltaLL <- c(NA_real_, diff(logLik))
    relDeltaLL <- c(NA_real_, deltaLL[-1L] / abs(logLik[-length(logLik)]))

    return(data.frame(
      iteration = iteration,
      logLik = logLik,
      deltaLL = deltaLL,
      relDeltaLL = relDeltaLL
    ))
  }

  history <- as.data.frame(history)
  names.history <- names(history)

  findColumn <- function(candidates) {
    match <- match(tolower(candidates), tolower(names.history), nomatch = 0L)
    match <- match[match > 0L]
    if (length(match)) names.history[match[[1L]]] else NULL
  }

  logLikCol <- findColumn(c("logLik", "loglik", "ll", "obsLL"))
  mod_stopif(is.null(logLikCol), "`history` must contain a log-likelihood column.")

  iterCol <- findColumn(c("iteration", "iterations", "iter", "i"))
  deltaCol <- findColumn(c("deltaLL", "dLL", "delta", "gain"))
  relDeltaCol <- findColumn(c("relDeltaLL", "reldLL", "relDelta", "relGain"))

  logLik <- as.numeric(history[[logLikCol]])
  iteration <- if (is.null(iterCol)) seq_along(logLik) else as.integer(history[[iterCol]])
  deltaLL <- if (is.null(deltaCol)) c(NA_real_, diff(logLik)) else as.numeric(history[[deltaCol]])

  if (is.null(relDeltaCol)) {
    relDeltaLL <- rep(NA_real_, length(logLik))
    if (length(logLik) > 1L) {
      denom <- abs(logLik[-length(logLik)])
      relDeltaLL[-1L] <- deltaLL[-1L] / denom
    }
  } else {
    relDeltaLL <- as.numeric(history[[relDeltaCol]])
  }

  data.frame(
    iteration = iteration,
    logLik = logLik,
    deltaLL = deltaLL,
    relDeltaLL = relDeltaLL
  )
}


forecastLogLikLms <- function(history,
                              i = NULL,
                              convergence.abs = 1e-4,
                              convergence.rel = 1e-10,
                              min.gains = 3L,
                              window = 10L) {
  empty <- function() {
    data.frame(
      expected.remaining.iterations = NA_real_,
      expected.final.loglik = NA_real_
    )
  }

  history <- normalizeLmsLogLikHistory(history)
  if (!NROW(history)) return(empty())

  if (!is.null(i)) {
    history <- history[history$iteration <= i, , drop = FALSE]
    if (!NROW(history)) return(empty())
  }

  okLL <- is.finite(history$logLik)
  history <- history[okLL, , drop = FALSE]
  if (!NROW(history)) return(empty())

  currentLL <- history$logLik[NROW(history)]
  currentGain <- history$deltaLL[NROW(history)]
  gain <- history$deltaLL
  gainIteration <- history$iteration
  okGain <- is.finite(gain) & gain > 0
  gain <- gain[okGain]
  gainIteration <- gainIteration[okGain]

  if (length(gain) < min.gains) return(empty())

  if (is.finite(window) && length(gain) > window) {
    keep <- seq.int(length(gain) - window + 1L, length(gain))
    gain <- gain[keep]
    gainIteration <- gainIteration[keep]
  }

  fitData <- data.frame(iteration = gainIteration, logGain = log(gain))
  fit <- stats::lm(logGain ~ iteration, data = fitData)
  coefs <- stats::coef(fit)
  slope <- unname(coefs[["iteration"]])
  intercept <- unname(coefs[["(Intercept)"]])
  decay <- exp(slope)

  if (!is.finite(decay) || decay <= 0 || decay >= 1) {
    ratios <- gain[-1L] / gain[-length(gain)]
    ratios <- ratios[is.finite(ratios) & ratios > 0 & ratios < 1]
    if (!length(ratios)) return(empty())
    decay <- stats::median(ratios)
    expectedNextGain <- gain[[length(gain)]] * decay
  } else {
    expectedNextGain <- exp(intercept + slope * (history$iteration[NROW(history)] + 1L))
  }

  if (!is.finite(expectedNextGain) || expectedNextGain < 0) return(empty())

  decay <- max(0, min(decay, 1 - sqrt(.Machine$double.eps)))
  expectedRemainingGain <- if (decay == 0) expectedNextGain else expectedNextGain / (1 - decay)
  expectedFinalLogLik <- currentLL + expectedRemainingGain

  relTol <- convergence.rel * max(1, abs(currentLL), abs(expectedFinalLogLik), na.rm = TRUE)
  stopGain <- max(convergence.abs, relTol, na.rm = TRUE)
  currentAbsGain <- abs(currentGain)

  if (is.finite(currentAbsGain) && currentAbsGain < stopGain) {
    expectedRemainingIterations <- 0L
  } else if (expectedNextGain <= stopGain || decay == 0) {
    expectedRemainingIterations <- 1L
  } else {
    expectedRemainingIterations <- ceiling(log(stopGain / expectedNextGain) / log(decay) + 1)
  }

  data.frame(
    expected.remaining.iterations = expectedRemainingIterations,
    expected.final.loglik = expectedFinalLogLik
  )
}


appendLmsIterationHistory <- function(history, iteration, loglik.history,
                                      convergence.abs, convergence.rel,
                                      forecast.min.gains = 3L,
                                      forecast.window = 10L) {
  forecast <- forecastLogLikLms(
    history = loglik.history,
    i = iteration,
    convergence.abs = convergence.abs,
    convergence.rel = convergence.rel,
    min.gains = forecast.min.gains,
    window = forecast.window
  )

  rbind(history, cbind(iteration = iteration, forecast))
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
                  ema.min.iter = 10L,
                  ema.trend.min = 0,
                  ema.forecast.min.gains = 3L,
                  ema.forecast.window = 10L,
                  qn.remaining.min = 10,
                  qn.remaining.max = 200,
                  fs.remaining.min = 200,
                  fs.remaining.max = 500,
                  fallback.max.step.min = 50,
                  fallback.max.step.mult = 50,
                  qn.maxit.start = 1L,
                  qn.maxit.mult = 2,
                  qn.maxit.max = 25L,
                  qn.factr = NULL,
                  qn.pgtol = 0,
                  ema.final.qn = TRUE,
                  ema.final.qn.maxit = qn.maxit.max,
                  ema.final.qn.factr = qn.factr,
                  ema.final.qn.rounds.max = 25L,
                  fs.maxit.start = 1L,
                  ...) {

  algorithm <- toupper(match.arg(algorithm))
  ema.min.iter <- as.integer(ema.min.iter)
  ema.forecast.min.gains <- as.integer(ema.forecast.min.gains)
  qn.maxit.start <- as.integer(qn.maxit.start)
  qn.maxit.max <- as.integer(qn.maxit.max)
  ema.final.qn.maxit <- as.integer(ema.final.qn.maxit)
  ema.final.qn.rounds.max <- as.integer(ema.final.qn.rounds.max)
  fs.maxit.start <- as.integer(fs.maxit.start)

  mod_stopif(ema.min.iter < 2L, "`ema.min.iter` must be at least 2.")
  mod_stopif(ema.forecast.min.gains < 2L, "`ema.forecast.min.gains` must be at least 2.")
  mod_stopif(ema.forecast.window < ema.forecast.min.gains,
             "`ema.forecast.window` must be at least `ema.forecast.min.gains`.")
  mod_stopif(qn.remaining.min >= qn.remaining.max,
             "`qn.remaining.min` must be smaller than `qn.remaining.max`.")
  mod_stopif(fs.remaining.min >= fs.remaining.max,
             "`fs.remaining.min` must be smaller than `fs.remaining.max`.")
  mod_stopif(qn.maxit.start < 1L || qn.maxit.max < qn.maxit.start,
             "`qn.maxit.max` must be at least `qn.maxit.start`.")
  mod_stopif(qn.maxit.mult < 1, "`qn.maxit.mult` must be at least 1.")
  mod_stopif(!is.null(qn.factr) && qn.factr <= 0, "`qn.factr` must be positive.")
  mod_stopif(!is.null(ema.final.qn.factr) && ema.final.qn.factr <= 0,
             "`ema.final.qn.factr` must be positive.")
  mod_stopif(!is.null(qn.pgtol) && qn.pgtol < 0, "`qn.pgtol` must be non-negative.")
  mod_stopif(ema.final.qn.maxit < 1L, "`ema.final.qn.maxit` must be at least 1.")
  mod_stopif(ema.final.qn.rounds.max < 1L,
             "`ema.final.qn.rounds.max` must be at least 1.")
  mod_stopif(fs.maxit.start < 1L, "`fs.maxit.start` must be at least 1.")

  default.optim.factr <- if (is.finite(convergence.rel) && convergence.rel > 0) {
    max(1, convergence.rel / .Machine$double.eps)
  } else {
    NULL
  }
  qn.factr <- if (is.null(qn.factr)) default.optim.factr else qn.factr
  ema.final.qn.factr <- if (is.null(ema.final.qn.factr)) qn.factr else ema.final.qn.factr

  history <- lmsIterationHistory()
  loglik.history <- numeric()

  theta.lower  <- model$params$bounds$lower
  theta.upper  <- model$params$bounds$upper

  max.iter.fs <- fs.maxit.start
  max.iter.qn <- qn.maxit.start

  tryCatch({
    logLikNew <- -Inf
    logLikOld <- -Inf
    thetaNew  <- model$theta

    lastQuad     <- NULL
    adaptiveQuad <- model$models[[1L]]$quad$adaptive # fixed across groups
    adaptiveFreq <- model$models[[1L]]$quad$adaptive.frequency

    n.em.iter  <- 0
    mode       <- "EM"
    iterations <- 0L
    run        <- TRUE
    ema.final.phase <- "none"
    ema.final.qn.rounds <- 0L

    testSimpleGradient <- !model$params$gradientStruct$useFDGradient

    while (run) {
      iterations <- iterations + 1L
      logLikOld  <- logLikNew
      thetaOld   <- thetaNew
      recalcQuad <- adaptiveQuad && iterations %% adaptiveFreq == 0L
      tmp.max.step <- 0L

      mod_stopif(any(!is.finite(thetaOld)),
        "Missing/Non finite values in estimated parameters at iteration",
        paste0(iterations, "!")
      )

      # E-step at thetaOld
      P <- estepLms(
        model             = model,
        theta             = thetaOld,
        lastQuad          = lastQuad,
        recalcQuad        = recalcQuad,
        adaptive.quad.tol = adaptive.quad.tol,
        ...
      )

      if (testSimpleGradient) {
        tryCatch({
          gradientObsLogLikLms(theta = thetaNew, model = model, P = P)
        }, error = \(e) {
          mod_msg_warn_immediate(
            paste0("Optimized computation of gradient failed! Switching gradient type.\n",
                   "Message: ", conditionMessage(e)))
          model$params$gradientStruct$useFDGradient <<- TRUE
          model$params$gradientStruct$isNonLinear  <<- TRUE
        })
        testSimpleGradient <- FALSE
      }

      lastQuad <- P$quad

      logLikNew  <- P$obsLL
      deltaLL    <- logLikNew - logLikOld
      relDeltaLL <- if (is.finite(logLikOld)) deltaLL / abs(logLikOld) else Inf
      loglik.history <- c(loglik.history, logLikNew)

      updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)

      history <- appendLmsIterationHistory(
        history         = history,
        iteration       = iterations,
        loglik.history  = loglik.history,
        convergence.abs = convergence.abs,
        convergence.rel = convergence.rel,
        forecast.min.gains = ema.forecast.min.gains,
        forecast.window = ema.forecast.window
      )

      converged <- (abs(deltaLL) < convergence.abs) ||
                   (abs(relDeltaLL) < convergence.rel)
      converged.em <- converged && mode == "EM"
      run.final.qn <- FALSE

      if (iterations >= max.iter) break
      if (algorithm == "EMA" && isTRUE(ema.final.qn) &&
          ema.final.phase == "qn") {
        if (converged || ema.final.qn.rounds >= ema.final.qn.rounds.max) {
          ema.final.phase <- "em"
        } else {
          run.final.qn <- TRUE
        }

      } else if (converged.em) {
        if (algorithm == "EMA" && isTRUE(ema.final.qn) &&
            ema.final.phase == "none") {
          ema.final.phase <- "qn"
          run.final.qn <- TRUE
        } else {
          break
        }
      }

      # If quadrature has just been recalculated, we should be less strict
      quad.err <- as.integer(recalcQuad) * P$quad.err
      if (deltaLL < -max(1e-8, quad.err, na.rm = TRUE)) {
        if (verbose) cat("\n")
        mod_msg_warn_immediate(sprintf("Loglikelihood decreased by %.2g", deltaLL))
      }

      if (mode != "EM") { # switch back to EM
        mode <- "EM"
        n.em.iter <- 0
      } else n.em.iter <- n.em.iter + 1L

      # EMA controller
      if (run.final.qn) {
        mode <- "QN"
        updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)

      } else if (ema.final.phase == "em") {
        mode <- "EM"

      } else if (algorithm == "EMA" && n.em.iter > ema.min.iter) {
        # Check if expected number of iterations is increasing
        change <- diff(tail(history$expected.remaining.iterations, ema.min.iter))
       
        if (all(is.finite(change)) && mean(change) > ema.trend.min) {
          remaining <- tail(history$expected.remaining.iterations, 1L)[[1L]]

          # Increase the max step in case remaining > 500 or 
          # QN/FS fails, and we fall back to an EM step
          tmp.max.step <- max(fallback.max.step.min,
                            fallback.max.step.mult * max.step,
                            na.rm = TRUE)

          if (remaining > fs.remaining.min && remaining <= fs.remaining.max) {
            mode <- "FS"
          } else if (remaining > qn.remaining.min && remaining <= qn.remaining.max) {
            mode <- "QN"
          }

          updateStatusLog(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose)
        }
      }

      .obsfn <- function(x) {
        ll <- obsLogLikLms(theta = x, model = model, P = P, sign = -1)
        if (!is.finite(ll)) .Machine$double.xmax^(1/10) else ll # optim doesn't handle Inf well
      }

      .obsgr <- function(x) {
        grad <- gradientObsLogLikLms(theta = x, model = model, P = P, sign = -1, epsilon = epsilon)
   
        non.finite <- !is.finite(grad)
        if (any(non.finite)) {
          # This might happen during a bad line search in optim()
          # as long as obsLogLikLms() is non finite as well, it should be ok,
          # as it will be rejected...
          ll <- obsLogLikLms(theta = x, model = model, P = P, sign = -1)
          mod_stopif(is.finite(ll),
            "Some gradient values are NaN, but objective is finite!"
          )

          dir <- sign(x - thetaOld)
          dir[!is.finite(dir)] <- 0

          grad[non.finite] <- dir[non.finite] * .Machine$double.xmax^(1/10)
        }

        grad
      }

      .obshs <- function(x, louis = TRUE) {
        computeFullIobs(theta = x, model = model, P = P, louis = louis)
      }

      .error <- function(e) {
        mod_msg_warn(mode, "mode failed! Message:", conditionMessage(e))

        mstepLms(
          model = model,
          P = P,
          theta = thetaOld,
          max.step = max(max.step, tmp.max.step, na.rm = TRUE),
          epsilon = epsilon,
          optimizer = optimizer,
          control = control,
          ...
        )
      }

      # accelerated step
      switch(mode,
        EM = {
          mstep <- mstepLms(
            model = model,
            P = P,
            theta = thetaOld,
            max.step = max(max.step, tmp.max.step, na.rm = TRUE),
            epsilon = epsilon,
            optimizer = optimizer,
            control = control,
            ...
          )

        },

        QN = {
          # Quasi Newton
          mstep <- tryCatch({
            qn.control <- list(maxit = max.iter.qn)
            if (!is.null(qn.factr)) qn.control$factr <- qn.factr
            if (!is.null(qn.pgtol)) qn.control$pgtol <- qn.pgtol

            if (ema.final.phase == "qn") {
              qn.control$maxit <- ema.final.qn.maxit
              if (!is.null(ema.final.qn.factr)) {
                qn.control$factr <- ema.final.qn.factr
              }
            }

            optim(
              par = thetaOld,
              fn = .obsfn,
              gr = .obsgr,
              lower = theta.lower,
              upper = theta.upper,
              control = qn.control,
              method = "L-BFGS-B"
            )
          }, error = .error)

          if (ema.final.phase == "qn") {
            ema.final.qn.rounds <- ema.final.qn.rounds + 1L
          } else {
            max.iter.qn <- min(qn.maxit.max, qn.maxit.mult * max.iter.qn)
          }
        },

        FS = {
          # Fischer Scoring
          mstep <- tryCatch({
            fischerScoring(
              par = thetaOld,
              fn = .obsfn,
              gr = .obsgr,
              hs = .obshs,
              lower = theta.lower,
              upper = theta.upper,
              max.iter = max.iter.fs
            )
          }, error = .error)
        }

        )

      if (!is.null(mstep$par) &&
          length(mstep$par) == length(thetaOld) &&
          all(is.finite(mstep$par))) {
        thetaNew <- mstep$par
      } else {
        thetaNew <- thetaOld
      }
    } # while

    if (verbose) cat("\n")
    mod_warnif(iterations >= max.iter, paste0("Maximum iterations reached!\n",
           "Consider tweaking these parameters:\n",
           formatParameters(convergence.abs, convergence.rel, algorithm,
                            max.step, max.iter, nodes, adaptive.quad,
                            adaptive.quad.tol, quad.range)))

    # final E-step
    P <- estepLms(model = model, theta = thetaNew,
                  lastQuad = lastQuad, recalcQuad = FALSE,
                  adaptive.quad.tol = adaptive.quad.tol, ...)

    fit <- finalizeModelEstimatesDA(
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
    fit$iteration.history <- history
    fit

  }, error = function(e) {
    if (verbose) cat("\n")
    mod_msg_warn(paste0(
      "Model estimation failed, returning starting values!\n",
      "Message: ", conditionMessage(e)
    ), newline. = verbose)
    P0 <- tryCatch(
      estepLms(model = model, theta = model$theta,
               lastQuad = NULL, recalcQuad = FALSE,
               adaptive.quad.tol = adaptive.quad.tol, ...),
      error = function(e2) NULL
    )
    ll0 <- if (!is.null(P0)) P0$obsLL else NA_real_

    fit <- finalizeModelEstimatesDA(
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
    fit$iteration.history <- history
    fit
  })
}


fischerScoring <- function(par,
                           fn,
                           gr,
                           hs, ...,
                           lower = -Inf,
                           upper = Inf,
                           max.iter = 1L,
                           max.alpha = 1,
                           min.alpha = 1e-12,
                           epsilon = 1e-6,
                           abs.tol = 1e-4,
                           jitter.mult = sqrt(.Machine$double.eps)) {

  truncate <- function(x) {
    underflow <- x < lower
    overflow  <- x > upper 

    x[underflow] <- lower[underflow]
    x[overflow]  <- upper[overflow]

    x 
  }

  f0 <- fn(par, ...)

  for (i in seq_len(max.iter)) {
    g0 <- gr(par, ...)
    I  <- hs(par, ...)

    # conditioning & scaled jitter
    rc <- suppressWarnings(rcond(I))
    if (!is.finite(rc) || rc < 1e-10) {
      jj <- jitter.mult * max(1, max(diag(I), na.rm = TRUE))
      diag(I) <- diag(I) + jj
    }

    I <- 0.5 * (I + t(I))
    direction <- tryCatch(-solve(I, g0), error = function(e) NULL)

    if (is.null(direction)) {
      return(list(
        convergence = -1,
        par = par,
        objective = f0,
        iterations = i
      ))
    }
        
    alpha   <- max.alpha
    success <- FALSE

    while (alpha > min.alpha) {
      parTrial  <- truncate(par + alpha * direction)
      f1  <- suppressWarnings(fn(parTrial, ...))

      ok <- is.finite(f1) && f1 <= f0 
      if (ok) { success <- TRUE; break }
      alpha <- alpha / 2
    }

    # convergence
    if (success && abs(f1 - f0) <= abs.tol) {
      return(list(
        convergence = 0,
        iterations = i,
        par = parTrial,
        objective = f1
      ))

    } else if (success) {
      par <- parTrial
      f0 <- f1

    } else {
      return(list(
        convergence = 10,
        iterations = i,
        par = par,
        objective = f0
      ))
    }
  }

  list(
    convergence = 1,
    iterations  = i,
    par = par,
    objective = f0
  )
}
