modsem.LMS <- function(modelSpec,
                       data,
                       qml = TRUE,
                       verbose = TRUE, 
                       optimize = TRUE,
                       m = 16, 
                       suppressWarnings = FALSE,
                       convergence = 1e-2,
                       ...) {
  model <- specifyLmsModel(modelSpec$modelSyntax, data = data, m = m)
  model <- optimizeStartingParamsLms(model)
  emLms(model, verbose = verbose, convergence = convergence, ...)
}
