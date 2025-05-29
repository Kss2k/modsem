printModsemPIHeader <- function(approach) {
  cat(paste0("modsem (version ", PKG_INFO$version, ", approach = ",
             approach, "):\n"))
}


#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param ... arguments passed to lavaan::summary()
#' @rdname summary
#' @export
summary.modsem_pi <- function(object, ...) {
  structure(list(lavaan = lavaan::summary(object$lavaan, ...),
                 info   = list(version = PKG_INFO$version,
                               approach = attributes(object)$method)),
            class = c("summary_modsem_pi", "list"))
}


#' @export
print.summary_modsem_pi <- function(x, ...) {
  printModsemPIHeader(x$info$approach)
  print(x$lavaan)
}


#' @export
parameter_estimates.modsem_pi <- function(object, ...) {
  lavaan::parameterEstimates(object$lavaan, ...)
}


#' @export
standardized_estimates.modsem_pi <- function(object, ...) {
  parTable <- lavaan::standardizedSolution(object$lavaan, ...)

  parTable$est <- parTable$est.std
  parTable$est.std <- NULL

  parTable
}


#' @export
modsem_inspect.modsem_pi <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "free"
  lavaan::lavInspect(object$lavaan, what = what, ...)
}


#' @export
#' @importFrom stats vcov
vcov.modsem_pi <- function(object, ...) {
  lavaan::vcov(object$lavaan, ...)
}


#' @export
#' @importFrom stats coef
coef.modsem_pi <- function(object, ...) {
  lavaan::coef(object$lavaan, ...)
}


#' @export
#' @importFrom stats coefficients
coefficients.modsem_pi <- function(object, ...) {
  lavaan::coef(object$lavaan, ...)
}


#' @export
#' @importFrom stats nobs
nobs.modsem_pi <- function(object, ...) {
  lavaan::nobs(object$lavaan, ...)
}


#' @export
print.modsem_pi <- function(x, ...) {
  printModsemPIHeader(attributes(x)$method)
  print(x$lavaan)
}
