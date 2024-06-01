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
  final <- mstepQml(model = model, theta = startTheta, maxIter = maxIter, 
                    convergence = convergence,
                    negHessian = negHessian, verbose = verbose, ...)
  coefficients <- final$par
  finalModel <- fillModel(model, coefficients)

  info <- model$info

  coefficients <- final$par
  finalModel <- fillModel(model, coefficients)
  finalModel$matricesNA <- model$matrices 
  finalModel$matricesSE <- fillModel(model, calcSE(final$hessian))$matrices 
  finalModel$data <- NULL # not needed in the model anymore 
  parTable <- finalModelToParTable(finalModel, method = "qml")

  parTable$tvalue <- parTable$est / parTable$se
  parTable$pvalue <- 2 * stats::pnorm(-abs(parTable$tvalue))
  parTable$ciLower <- parTable$est - 1.96 * parTable$se
  parTable$ciUpper <- parTable$est + 1.96 * parTable$se

  out <- list(model = finalModel, 
              data  = model$data,
              theta = coefficients,
              parTable = parTable,
              originalParTable = model$parTable,
              logLik = -final$objective, 
              iterations = final$iterations,
              convergence = final$convergence,
              negHessian = final$hessian)

  class(out) <- "modsem_qml"
  out
}


