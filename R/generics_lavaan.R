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

  transformEstimatesLavaan(
    object         = object,
    monte.carlo    = monte.carlo,
    mc.reps        = mc.reps,
    tolerance.zero = tolerance.zero,
    transformation = centering
  )
}


#' @describeIn standardized_estimates Method for \code{lavaan} objects
#'
#' @param monte.carlo Logical. If \code{TRUE}, use Monte Carlo simulation to estimate
#' standard errors; if \code{FALSE}, use the delta method (default).
#' @param mc.reps Number of Monte Carlo repetitions. Default is 10000.
#' @param tolerance.zero Threshold below which standard errors are set to \code{NA}.
#'
#' @export
standardized_estimates.lavaan <- function(object, 
                                          correction = FALSE, 
                                          monte.carlo = FALSE,
                                          mc.reps = 10000,
                                          tolerance.zero = 1e-10,
                                          ...) {

  standardization <- function(object, grouping = NULL) {

    solution <- standardizedSolutionCOEFS(
      object      = object, 
      monte.carlo = monte.carlo, 
      mc.reps     = mc.reps, 
      grouping    = grouping, 
      ...
    ) 

    solution$parTable
  }

  transformEstimatesLavaan(
    object         = object,
    monte.carlo    = monte.carlo,
    mc.reps        = mc.reps,
    tolerance.zero = tolerance.zero,
    transformation = standardization 
  )
}


transformEstimatesLavaan <- function(object, 
                                     monte.carlo = FALSE,
                                     mc.reps = 10000,
                                     tolerance.zero = 1e-10,
                                     transformation) {
  parTable <- parameter_estimates(object)
 
  hiorder <- isHigherOrderParTable(parTable)
  cluster <- isClustered(object)

  stopif(hiorder, "Transformations of higher-order models are not supported!")
  stopif(cluster, "Transformations of clustered (multilevel) models is not supported!")

  applyTransformationByGrouping(
    parTable = parTable,
    FUN      = transformation,
    object   = object 
  )
}
