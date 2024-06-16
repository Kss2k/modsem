#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param ... arguments passed to lavaan::summary()
#' @rdname summary
#' @export 
summary.modsem_pi <- function(object, ...) {
  cat("modsem: \nMethod =", attributes(object)$method, "\n")
  lavaan::summary(object$lavaan, ...)
}


#' @export
parameter_estimates.modsem_pi <- function(object, ...) {
  lavaan::parameterEstimates(object$lavaan, ...)
}


#' @export 
standardized_estimates.modsem_pi <- function(object, ...) {
  lavaan::standardizedSolution(object$lavaan, ...)
}


#' @export 
modsem_inspect.modsem_pi <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "free"
  lavaan::lavInspect(object$lavaan, what = what, ...)
}


#' @export 
vcov.modsem_pi <- function(object, ...) {
  lavaan::vcov(object$lavaan, ...)
}


#' @export 
coef.modsem_pi <- function(object, ...) {
  lavaan::coef(object$lavaan, ...)
}


#' @export 
coefficients.modsem_pi <- function(object, ...) {
  lavaan::coef(object$lavaan, ...)
}
