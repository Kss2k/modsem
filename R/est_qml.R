estQml <- function(model, 
                   convergence = 1e-2,
                   verbose = FALSE, 
                   max.iter = 1000,
                   hessian = TRUE,
                   robust.se = FALSE,
                   ...) {
  startTheta <- model$theta
  final <- mstepQml(model = model, theta = startTheta, max.iter = max.iter, 
                    convergence = convergence,
                    hessian = hessian, verbose = verbose, ...)
  coefficients <- final$par
  finalModel <- fillModel(model, coefficients)

  info <- model$info

  finalModel <- fillModel(model, coefficients, method = "qml")
  emptyModel <-  getEmptyModel(parTable = model$parTable, 
                               cov.syntax = model$cov.syntax,
                               parTableCovModel = model$covModel$parTable,
                               method = "qml")
  finalModel$matricesNA <- emptyModel$matrices
  finalModel$covModelNA <- emptyModel$covModel

  # Caclulate standard errors
  if (hessian) {
    SE <- calcSE(final$hessian) 
  } else if (robust.se) {
    SE <- calcRobustSE_qml(model, theta = coefficients, verbose = verbose)
  } else {
    SE <- rep(-999, length(coefficients))
  }
  modelSE <- fillModel(replaceNonNaModelMatrices(model, value = -999), 
                       SE, method = "qml")
  finalModel$matricesSE <- modelSE$matrices
  finalModel$covModelSE <- modelSE$covModel

  parTable <- rbind(finalModelToParTable(finalModel, method = "qml"),
                    covModelToParTable(finalModel, method = "qml"))

  parTable$z.value <- parTable$est / parTable$std.error
  parTable$p.value <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - 1.96 * parTable$std.error
  parTable$ci.upper <- parTable$est + 1.96 * parTable$std.error

  out <- list(model = finalModel, 
              data  = model$data,
              theta = coefficients,
              parTable = parTable,
              originalParTable = model$parTable,
              logLik = -final$objective, 
              iterations = final$iterations,
              convergence = final$convergence,
              hessian = final$hessian)

  class(out) <- "modsem_qml"
  out
}


