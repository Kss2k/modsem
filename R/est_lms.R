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
  logLikNew  <- 0
  logLikOld  <- 0
  iterations <- 0
  thetaNew   <- model$theta

  bestLogLik    <- -Inf
  bestP         <- NULL
  bestTheta     <- NULL
  logLiks       <- NULL
  logLikChanges <- NULL

  run     <- TRUE
  doEstep <- TRUE

  nFalseConvergence <- 0
  nNegCheck <- 20
  pNegCheck <- 0.5

  while(run) {
    logLikOld <- logLikNew
    thetaOld  <- thetaNew

    if (doEstep) P <- estepLms(model = model, theta = thetaOld, data = data, ...)

    mstep <- mstepLms(model = model, P = P, data = data, theta = thetaOld,
                      max.step = max.step, epsilon = epsilon,
                      optimizer = optimizer, control = control, ...)

    logLikNew     <- -mstep$objective
    thetaNew      <- unlist(mstep$par)
    iterations    <- iterations + 1
    logLiks       <- c(logLiks, logLikNew)
    logLikChanges <- c(logLikChanges, logLikNew - logLikOld)

    if (verbose) {
      clearConsoleLine()
      printf("\rEM: Iteration = %d, LogLik = %.2f, Change = %.3f",
             iterations, logLikNew, logLikNew - logLikOld)
    }

    if (logLikNew > bestLogLik) {
      bestLogLik <- logLikNew
      bestP      <- P
      bestTheta  <- thetaOld
    }

    converged <- abs(logLikOld - logLikNew) < convergence

    if (iterations >= max.iter || converged && (nFalseConvergence >= 3 || doEstep)) {
      run <- FALSE

    } else if (converged & !doEstep) {
      nFalseConvergence <- nFalseConvergence + 1
      doEstep <- TRUE

    } else if (
        doEstep && runningAverage(logLikChanges, n = 5) < 0 && iterations > max.iter/2 &&
        nNegativeLast(logLikChanges, n = nNegCheck) >= nNegCheck * pNegCheck
      ) {

      cat("\n")
      warning2("EM algorithm is not converging. Might be at a saddle point!")

      doEstep  <- FALSE
      P        <- bestP
      thetaNew <- bestTheta
    }
  }


  if (verbose) cat("\n")
  
  warnif(iterations >= max.iter, "Maximum number of iterations was reached. ",
         "EM algorithm might not have converged.")

  if (!doEstep) P <- estepLms(model = model, theta = thetaNew, data = data, ...)

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

  # Caclulate information matrix (I) and standard errors (SE)
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

  # convergence of em
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
