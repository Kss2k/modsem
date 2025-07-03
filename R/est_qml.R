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
                   ...) {
  startTheta <- model$theta
  final <- mstepQml(model = model, theta = startTheta, max.iter = max.iter,
                    convergence = convergence, epsilon = epsilon,
                    verbose = verbose, optimizer = optimizer, ...)

  coefficients <- final$par
  lavCoefs     <- getLavCoefs(model = model, theta = coefficients, method = "lms")
  finalModel      <- fillModel(model, coefficients)
  info            <- model$info

  finalModel <- fillModel(model, coefficients, method = "qml")
  emptyModel <- getEmptyModel(parTable = model$parTable,
                              cov.syntax = model$cov.syntax,
                              parTableCovModel = model$covModel$parTable,
                              mean.observed = model$info$mean.observed,
                              method = "qml")
  finalModel$matricesNA <- emptyModel$matrices
  finalModel$covModelNA <- emptyModel$covModel

  # Caclulate information matrix (I) and standard errors (SE)
  typeSE <- ifelse(!calc.se, "none", ifelse(robust.se, "robust", "standard"))
  FIM <- calcFIM_da(model = model, finalModel = finalModel, theta = coefficients,
                    data = model$data, method = "qml", EFIM.S = EFIM.S,
                    hessian = OFIM.hessian, calc.se = calc.se,
                    EFIM.parametric = EFIM.parametric, verbose = verbose,
                    FIM = FIM, robust.se = robust.se, NA__ = -999,
                    epsilon = epsilon, R.max = R.max)
  SE <- calcSE_da(calc.se = calc.se, FIM$vcov.all, rawLabels = FIM$raw.labels,
                  NA__ = -999)
  modelSE <- getSE_Model(model, se = SE, method = "qml",
                         n.additions = FIM$n.additions)

  finalModel$matricesSE <- modelSE$matrices
  finalModel$covModelSE <- modelSE$covModel

  parTable <- modelToParTable(finalModel, coefs = lavCoefs$all,
                              se = SE, method = "qml")

  parTable$z.value  <- parTable$est / parTable$std.error
  parTable$p.value  <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - CI_WIDTH * parTable$std.error
  parTable$ci.upper <- parTable$est + CI_WIDTH * parTable$std.error

  warnif(final$iterations >= max.iter,
         "Maximum number of iterations was reached, ",
         "model estimation might not have converged.")

  out <- list(model      = finalModel,
              method     = "qml",
              optimizer  = optimizer,
              data       = model$data,
              theta      = coefficients,
              coefs.all  = lavCoefs$all,
              coefs.free = lavCoefs$free,
              parTable   = modsemParTable(parTable),

              originalParTable = model$parTable,

              logLik      = -final$objective,
              iterations  = final$iterations,
              convergence = final$convergence,
              type.se     = typeSE,

              type.estimates = "unstandardized",

              info.quad   = NULL,
              FIM         = FIM$FIM,
              vcov.all    = FIM$vcov.all,
              vcov.free   = FIM$vcov.free,
              information = FIM$type)

  out
}
