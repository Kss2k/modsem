emLms <- function(model, verbose = FALSE,
                  convergence = 1e-02, maxiter = 500,
                  maxstep = 1,
                  breakOnLogLikIncrease = FALSE,
                  sampleGrad = NULL,
                  control = list(),
                  ...) {
  data <- model$data
  model$data <- NULL # not needed in the model anymore
  if (anyNA(data)) stop("Remove or replace missing values from data")
    
  # Initialization
  logLikNew <- 0
  logLikOld <- 1.0e+10
  iterations <- 0     
  thetaNew <- model$theta
  nLogIncreased <- 0
  run <- TRUE
  while(run) { # as long as no convergence is reached
    if (logLikNew - logLikOld > 0 && iterations > 3) {
      if (breakOnLogLikIncrease) {
        warning("Loglikelihood is increasing. EM algorithm will be stopped.")
        logLikNew <- logLikOld
        thetaNew <- thetaOld
        break
      }
      nLogIncreased <- nLogIncreased + 1
    }
    # Update loglikelihood
    logLikOld <- logLikNew
    thetaOld <- thetaNew

    P <- estepLms(model = model, theta = thetaOld, data = data, ...)
    mstep <- mstepLms(model = model, P = P, data = data,
                      theta = thetaOld, maxstep = maxstep, ...,
                      sampleGrad = sampleGrad,
                      control = control)

    logLikNew <- mstep$objective
    thetaNew <- unlist(mstep$par)
    iterations <- iterations + 1

    if (verbose) {
      cat(sprintf("EM: Iteration = %5d, LogLik = %11.2f, Change = %10.3f\n",
            iterations, -logLikNew, logLikOld - logLikNew))
    }
    if(iterations >= maxiter){
      warning("Maximum number of iterations was reached. ",
              "EM algorithm might not have converged.")
      break
    }
    if (abs(logLikOld - logLikNew) < convergence) run <- FALSE
  }
  final <- mstepLms(model = model, P = P, data = data,
                    theta = thetaNew, negHessian = TRUE,
                    maxstep = maxstep, sampleGrad = NULL, 
                    verbose = verbose, control = control,
                    ...)

  coefficients <- final$par
  finalModel <- fillModel(model, coefficients, fillPhi = TRUE)
  finalModel$matricesNA <- model$matrices 
  finalModel$matricesSE <- fillModel(model, calcSE(final$hessian))$matrices 
  parTable <- finalModelToParTable(finalModel, method = "lms")
  parTable$tvalue <- parTable$est / parTable$se
  parTable$pvalue <- 2 * stats::pnorm(-abs(parTable$tvalue))
  parTable$ciLower <- parTable$est - 1.96 * parTable$se
  parTable$ciUpper <- parTable$est + 1.96 * parTable$se

  # convergence of em
  if (iterations == maxiter) convergence <- FALSE else convergence <- TRUE

  out <- list(model = finalModel, 
              data = data,
              theta = coefficients,
              parTable = parTable,
              originalParTable = model$parTable,
              logLik = -final$objective, 
              iterations = iterations,
              convergence = convergence,
              negHessian = final$hessian)

  class(out) <- "modsem_lms"
  out
}
