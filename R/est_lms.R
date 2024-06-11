emLms <- function(model, 
                  verbose = FALSE,
                  convergence = 1e-02, 
                  max.iter = 500,
                  max.step = 1,
                  sampleGrad = NULL,
                  control = list(),
                  hessian = TRUE,
                  robust.se = FALSE,
                  ...) {
  data <- model$data
  model$data <- NULL # not needed in the model anymore
  if (anyNA(data)) stop2("Remove or replace missing values from data")

  # Initialization
  logLikNew <- 0
  logLikOld <- 1.0e+10
  iterations <- 0     
  thetaNew <- model$theta
  nLogIncreased <- 0
  run <- TRUE
  while(run) { # as long as no convergence is reached
    # Update loglikelihood
    logLikOld <- logLikNew
    thetaOld <- thetaNew

    P <- estepLms(model = model, theta = thetaOld, data = data, ...)
    mstep <- mstepLms(model = model, P = P, data = data,
                      theta = thetaOld, max.step = max.step, 
                      sampleGrad = sampleGrad, hessian = FALSE,
                      control = control, ...)

    logLikNew <- mstep$objective
    thetaNew <- unlist(mstep$par)
    iterations <- iterations + 1

    if (verbose) {
      cat(sprintf("EM: Iteration = %5d, LogLik = %11.2f, Change = %10.3f\n",
            iterations, -logLikNew, logLikOld - logLikNew))
    }
    if(iterations >= max.iter){
      warning2("Maximum number of iterations was reached. ",
              "EM algorithm might not have converged.")
      break
    }
    if (abs(logLikOld - logLikNew) < convergence) run <- FALSE
  }
  final <- mstepLms(model = model, P = P, data = data,
                    theta = thetaNew, hessian = hessian,
                    max.step = max.step, sampleGrad = NULL, 
                    verbose = verbose, control = control,
                    ...)
    
  coefficients <- final$par
  finalModel <- fillModel(model, coefficients, fillPhi = TRUE,
                          method = "lms")

  emptyModel <-  getEmptyModel(parTable = model$parTable, 
                               cov.syntax = model$cov.syntax,
                               parTableCovModel = model$covModel$parTable,
                               method = "lms")

  finalModel$matricesNA <- emptyModel$matrices
  finalModel$covModelNA <- emptyModel$covModel

  if (hessian) {
    SE <- calcSE(final$hessian) 
  } else if (robust.se) {
    SE <- calcRobustSE(model, theta = coefficients, verbose = verbose, 
                       method = "lms")
  } else {
    SE <- rep(-999, length(coefficients))
  }
  modelSE <- fillModel(replaceNonNaModelMatrices(model, value = -999), 
                       SE, method = "lms")
  finalModel$matricesSE <- modelSE$matrices
  finalModel$covModelSE <- modelSE$covModel

  parTable <- rbind(finalModelToParTable(finalModel, method = "lms"),
                    covModelToParTable(finalModel, method = "lms"))

  parTable$z.value <- parTable$est / parTable$std.error
  parTable$p.value <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - 1.96 * parTable$std.error
  parTable$ci.upper <- parTable$est + 1.96 * parTable$std.error

  # convergence of em
  if (iterations == max.iter) convergence <- FALSE else convergence <- TRUE

  out <- list(model = finalModel, 
              data = data,
              theta = coefficients,
              parTable = parTable,
              originalParTable = model$parTable,
              logLik = -final$objective, 
              iterations = iterations,
              convergence = convergence,
              hessian = final$hessian)

  class(out) <- "modsem_lms"
  out
}
