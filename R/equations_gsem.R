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

    # Fill with known values for testing
    # psi.g[1, 1] <- 0.98
    # psi.g[2, 2] <- 1.04
    # psi.g[1, 2] <- psi.g[2, 1] <- 0.198
    # psi.g[3, 3] <- 0.982

    cholPsi <- chol(psi.g)
    latent_names <- colnames(psi.g)
    alpha_latent <- if (nrow(alpha.g)) alpha.g[latent_names, 1L, drop = TRUE] else rep(0, length(latent_names))

    modFilled$models[[g]]$quad$n <- expand_quadrature_cpp(
      Z_list      = quad.g$n,
      cholPsi     = cholPsi,
      gamma       = gamma.g,
      omega       = OMEGA.g,
      alphaLatent = alpha_latent,
      xiNames     = xis.g,
      etaNames    = etas.g,
      latentNames = latent_names
    )
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


fdgrad <- function(x, .f, ...,  eps = 1e-5) {
  g <- numeric(length(x))

  f0 <- .f(x, ...)
  for (i in seq_along(x)) {
    xi <- x
    xi[i] <- x[i] + eps

    fi <- .f(xi, ...)
    g[i] <- (fi - f0) / eps
  }

  g
}


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
  gradient2  <- \(theta) setNames(fdgrad(x = theta, .f = objective, eps = 1e-7), names(theta))

  g1 <- gradient(theta)
  g2 <- gradient2(theta)

  cat("\nAnalytical Gradient:\n")
  print(g1)
  cat("Numerical Gradient:\n")
  print(g2)
  cat("Difference:\n")
  print(round(g2 - g1, 4))
  stop("Stopping after grad check")
  # timeExpr(
  # fit2 <- fastOneStepNLMINB(start = theta, objective = objective,
  #                           lower = lower, upper = upper)
  # )

  fit <- nlminb(start = theta, objective = objective, control = control,
                lower = lower, upper = upper, gradient = gradient)

  fit
}
