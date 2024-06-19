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
                   ...) {
  startTheta <- model$theta
  final <- mstepQml(model = model, theta = startTheta, max.iter = max.iter, 
                    convergence = convergence,
                    verbose = verbose, ...)
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

  # Caclulate information matrix (I) and standard errors (SE)
  FIM <- calcFIM_da(model = model, finalModel = finalModel, theta = coefficients, 
                    data = model$data, method = "qml", EFIM.S = EFIM.S,
                    hessian = OFIM.hessian, calc.se = calc.se, 
                    EFIM.parametric = EFIM.parametric, verbose = verbose,
                    FIM = FIM, robust.se = robust.se, NA__ = -999)
  SE <- calcSE_da(calc.se = calc.se, FIM$vcov, theta = coefficients, NA__ = -999)

  modelSE <- fillModel(replaceNonNaModelMatrices(model, value = -999), 
                       theta = SE, method = "lms")

  modelSE <- fillModel(replaceNonNaModelMatrices(model, value = -999), 
                       SE, method = "qml")
  finalModel$matricesSE <- modelSE$matrices
  finalModel$covModelSE <- modelSE$covModel

  parTable <- modelToParTable(finalModel, method = "qml")

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
              FIM = FIM$FIM,
              vcov = FIM$vcov)

  class(out) <- "modsem_qml"
  out
}


