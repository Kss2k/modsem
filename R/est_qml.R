estQml <- function(model, 
                   convergence = 1e-2,
                   verbose = FALSE, 
                   maxIter = 1000,
                   hessian = TRUE,
                   ...) {
  startTheta <- model$theta
  final <- mstepQml(model = model, theta = startTheta, maxIter = maxIter, 
                    convergence = convergence,
                    hessian = hessian, verbose = verbose, ...)
  coefficients <- final$par
  finalModel <- fillModel(model, coefficients)

  info <- model$info

  finalModel <- fillModel(model, coefficients, method = "qml")

  finalModel$matricesNA <- model$matrices 
  finalModel$covModelNA <- model$covModel

  if (hessian) {
    modelSE <- fillModel(model, calcSE(final$hessian), method = "qml")
  } else {
    modelSE <- fillModel(model, rep(NA, length(coefficients)), method = "qml")
  }
  finalModel$matricesSE <- modelSE$matrices
  finalModel$covModelSE <- modelSE$covModel

  parTable <- rbind(finalModelToParTable(finalModel, method = "qml"),
                    covModelToParTable(finalModel, method = "qml"))

  parTable$t.value <- parTable$est / parTable$std.error
  parTable$p.value <- 2 * stats::pnorm(-abs(parTable$t.value))
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


