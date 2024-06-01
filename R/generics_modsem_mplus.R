#' summary.modsem_mplus
#'
#' @param object modsem object to summarized
#' @param ... arguments passed to other functions
#' @rdname summary
#' @export 
summary.modsem_mplus <- function(object, ...) {
  object
}


#' @export
print.modsem_mplus <- function(x, ...) {
  cat("ModSEM: \nMethod =", attributes(x)$method, "\n")
  print(x$coefParTable)
}


#' @export
parameter_estimates.modsem_mplus <- function(object, ...) {
  object$coefParTable
}
