#' @export
parameter_estimates.lavaan <- function(object, colon.pi = TRUE, ...) {
  lavaan::parameterEstimates(object, ...)
}


isLavaanObject <- function(x) {
  inherits(x, "lavaan")
}


#' @export
#' @describeIn modsem_inspect Inspect a \code{\link{lavaan}} object
modsem_inspect.lavaan <- function(object, what = "free", ...) {
  lavaan::lavInspect(object, what = what, ...)
}


#' @describeIn centered_estimates Method for \code{lavaan} objects
#'
#' @param monte.carlo Logical. If \code{TRUE}, use Monte Carlo simulation to estimate
#' standard errors; if \code{FALSE}, use the delta method (default).
#' @param mc.reps Number of Monte Carlo repetitions. Default is 10000.
#' @param tolerance.zero Threshold below which standard errors are set to \code{NA}.
#'
#' @export
centered_estimates.lavaan <- function(object, 
                                      correction = FALSE, 
                                      monte.carlo = FALSE,
                                      mc.reps = 10000,
                                      tolerance.zero = 1e-10,
                                      ...) {

  centering <- function(object, grouping = NULL) {

    solution <- centeredSolutionCOEFS(
      object      = object, 
      monte.carlo = monte.carlo, 
      mc.reps     = mc.reps, 
      grouping    = grouping, 
      ...
    ) 

    solution$parTable
  }

  parTable <- parameter_estimates(object, )
 
  hiorder <- isHigherOrderParTable(parTable)
  cluster <- isClustered(object)

  stopif(hiorder, "Centering of higher-order models are not supported!")
  stopif(cluster, "Centering of clustered (multilevel) models is not supported!")

  applyTransformationByGrouping(
    parTable = parTable,
    FUN      = centering,
    object   = object, 
  )
}



