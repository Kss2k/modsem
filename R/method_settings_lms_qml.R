getMethodSettingsLmsQml <- function(method, args) {
    settingNames <- c("verbose", "optimize", "nodes", "convergence",
                      "center.data", "standardize.data", "standardize.out",
                      "standardize", "mean.observed", "double", "hessian")
    args <- args[settingNames]
    isMissing <- vapply(args, FUN.VALUE = logical(1L), FUN = is.null)
    missingArgs <- settingNames[isMissing]

    settings <- list(
        lms = list(verbose = FALSE, 
                   optimize = TRUE,
                   nodes = 16, 
                   convergence = 1e-2,
                   center.data = FALSE, 
                   standardize.data = FALSE,
                   standardize.out = FALSE, 
                   standardize = FALSE,
                   mean.observed = TRUE,
                   double = FALSE, 
                   hessian = TRUE),
        qml = list(verbose = FALSE, 
                   optimize = TRUE,
                   nodes = 0, 
                   convergence = 1e-3,
                   center.data = FALSE, 
                   standardize = FALSE,
                   standardize.data = FALSE,
                   standardize.out = FALSE, 
                   mean.observed = TRUE,
                   double = FALSE, 
                   hessian = TRUE)
    )
    
    if  (!method %in% names(settings)) {
        stop2("Unrecognized method")
    }

    args.out <- c(settings[[method]][missingArgs], args[!isMissing])
    
    args.out$standardize.data <- 
      args.out$standardize || args.out$standardize.data
    args.out$standardize.out <- 
      args.out$standardize || args.out$standardize.out
    args.out$mean.observed <- 
      !args.out$standardize && args.out$mean.observed

    args.out
}
