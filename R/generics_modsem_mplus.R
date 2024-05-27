#' summary.modsem_mplus
#'
#' @param object modsem object to summarized
#' @param ... arguments passed to other functions
#' @rdname summary
#' @export 
summary.modsem_mplus <- function(object, ...) {
  cat("ModSEM: \nMethod =", attributes(object)$method, "\n")
  object$coefParTable
}
