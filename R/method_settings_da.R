getMethodSettingsDA <- function(method, args = NULL) {
    settings <- list(
        lms = list(verbose = FALSE, 
                   optimize = TRUE,
                   nodes = 24, 
                   convergence = 1e-6,
                   center.data = FALSE, 
                   standardize.data = FALSE,
                   standardize.out = FALSE, 
                   standardize = FALSE,
                   mean.observed = TRUE,
                   double = FALSE, 
                   calc.se = TRUE,
                   FIM = "expected",
                   OFIM.hessian = FALSE,
                   EFIM.S = 3e4,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   max.iter = 500, 
                   max.step = 1),
        qml = list(verbose = FALSE, 
                   optimize = TRUE,
                   nodes = 0, 
                   convergence = 1e-6,
                   center.data = FALSE, 
                   standardize = FALSE,
                   standardize.data = FALSE,
                   standardize.out = FALSE, 
                   mean.observed = TRUE,
                   double = FALSE, 
                   calc.se = TRUE,
                   FIM = "observed",
                   OFIM.hessian = TRUE,
                   EFIM.S = 3e4,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   max.iter = 500, 
                   max.step = NULL)
    )

    if (is.null(args)) return(settings[method])

    settingNames <- unique(unlist(lapply(settings, FUN = names))) 
    args <- args[settingNames]
    isMissing <- vapply(args, FUN.VALUE = logical(1L), FUN = is.null)
    missingArgs <- settingNames[isMissing]
    
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
    args.out$OFIM.hessian <- 
      args.out$OFIM.hessian && !args.out$robust.se

    args.out
}



#' default arguments fro LMS and QML approach
#'
#' @param method which method to get the settings for
#' @return modsem_lms or modsem_qml object
#' @export
#' @description
#' This function returns the default settings for the LMS and QML approach.
#' @examples
#' library(modsem)
#' default_settings_da()
default_settings_da <- function(method = c("lms", "qml")) {
  getMethodSettingsDA(method = method, args = NULL)
}
