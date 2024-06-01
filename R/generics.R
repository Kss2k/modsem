#' Extract parameterEstimates from an estimated model
#'
#' @param object An object of class `modsem_pi`, `modsem_lms_qml`, or `modsem_mplus`
#' @param ... Additional arguments passed to other functions
#' @export 
parameter_estimates <- function(object, ...) {
  UseMethod("parameter_estimates")
}
