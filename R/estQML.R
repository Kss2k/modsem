estQml <- function(model, 
                  convergence = 1e-2,
                  verbose = FALSE, 
                  maxIter = 1000,
                  negHessian = TRUE,
                  ...) {
  if (model$info$numEtas > 1) {
    stop("Only one eta allowed in QML estimation")
  }
  startTheta <- model$theta
  mstep <- mstepQml(model = model, theta = startTheta, maxIter = maxIter, 
                    negHessian = negHessian, verbose = verbose, ...)
  coefficients <- mstep$par
  finalModel <- fillModel(model, coefficients, fillPhi = TRUE)


  info <- model$info

  out <- list(model = finalModel, emptyModel = model,
              coefficients=coefficients, startTheta = startTheta,
              objective=-mstep$objective,
              negHessian=mstep$hessian, mstep = mstep, info = info)
  out$finalLogLik <- mstep$objective
  out$iterations <- mstep$iterations

  class(out) <- "modsemQML"
  out
}


