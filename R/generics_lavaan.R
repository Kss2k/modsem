#' @export
parameter_estimates.lavaan <- function(object, ...) {
  lavaan::parameterEstimates(object, ...)
}


isLavaanObject <- function(x) {
  inherits(x, "lavaan")
}
