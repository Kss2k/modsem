estQml <- function(model,
                   convergence = 1e-2,
                   verbose = FALSE,
                   max.iter = 500,
                   calc.se = TRUE,
                   FIM = "observed",
                   OFIM.hessian = FALSE,
                   EFIM.S = 3e4,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   epsilon = 1e-6,
                   optimizer = "nlminb",
                   R.max = 1e6,
                   cr1s = TRUE,
                   ...) {
  data <- model$data

  tryCatch({

    startTheta <- model$theta
    final <- mstepQml(model = model, theta = startTheta, max.iter = max.iter,
                      convergence = convergence, epsilon = epsilon,
                      verbose = verbose, optimizer = optimizer, ...)

    finalizeModelEstimatesDA(
      model             = model,
      theta             = final$par,
      method            = "qml",
      data              = data,
      logLik            = -final$objective,
      iterations        = final$iterations,
      converged         = (final$convergence == 0L) &&
                          (final$iterations <= max.iter),
      optimizer         = optimizer,
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
      includeStartModel = TRUE,
      startModel        = model
    )

  }, error = function(e) {
    warning2(paste0(
      "Model estimation failed, returning starting values!\n",
      "Message: ", conditionMessage(e)
    ))

    finalizeModelEstimatesDA(
      model             = model,
      theta             = model$theta,   # start values
      method            = "qml",
      data              = data,
      logLik            = NA_real_,      # best-effort; unknown here
      iterations        = 0L,
      converged         = FALSE,
      optimizer         = optimizer,
      calc.se           = FALSE,         # <- ensure no SEs
      FIM               = FIM,
      OFIM.hessian      = OFIM.hessian,
      EFIM.S            = EFIM.S,
      EFIM.parametric   = EFIM.parametric,
      robust.se         = robust.se,
      epsilon           = epsilon,
      cr1s              = cr1s,
      R.max             = R.max,
      verbose           = verbose,
      includeStartModel = TRUE,
      startModel        = model
    )
  })
}
