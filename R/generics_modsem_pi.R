#' summary.modsem_pi
#'
#' @param object modsem object to summarized
#' @param ... arguments passed to lavaan::summary()
#' @rdname summary
#' @export 
summary.modsem_pi <- function(object, ...) {
  cat("ModSEM: \nMethod =", attributes(object)$method, "\n")
  lavaan::summary(object$lavaan, ...)
}

