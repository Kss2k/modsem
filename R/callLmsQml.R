modsem.LMS <- function(modelSpec,
                       method,
                       data,
                       qml = TRUE,
                       verbose = FALSE, 
                       optimize = TRUE,
                       m = 16, 
                       convergence = 1e-2,
                       ...) {
  model <- specifyLmsModel(modelSpec$modelSyntax, data = data, 
                           method = method, m = m)
  if (optimize) model <- optimizeStartingParamsLms(model)
  switch(method, 
         "qml" = estQml(model, verbose = verbose, 
                        convergence = convergence, ...),
         "lms" = emLms(model, verbose = verbose, 
                       convergence = convergence, ...))
}
