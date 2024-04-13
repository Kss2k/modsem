#' @export
summary.modsemLMS <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}


#' @export 
summary.modsemQML <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}


summaryLmsAndQml <- function(object, ...) {
  # estimates
  est <- object$coefficients

  # standard errors
  stdError <- calcStandardError(object)
  names(stdError) <- names(est)


  tvalue <- est / stdError
  pvalue <- 2 * stats::pnorm(-abs(tvalue))

  coefTable <- data.frame(Parameter = names(est),
                          Estimate = est, `Std. Error` = stdError,
                          `z value` = tvalue, 
                          `Pr(>|z|)` = format.pval(pvalue))
  ans <- list(matricesTheta = fillModel(object$emptyModel, est),
              matricesSE = fillModel(stripModel(object$emptyModel), 
                                     stdError),
              estimates = coefTable,
              iterations = object$iterations,
              finallogLik = object$objective)
  class(ans) <- "summaryModsemLMSQML"
  ans
}


#' @export
print.summaryModsemLMSQML <- function(x, digits = 3, ...) {
  cat("\nSummary for model of class", x$model, "\n")
  cat("\nEstimates:\n")
  est <- lapply(x$estimates, function(col)
    if (is.numeric(col)) round(col, digits) else col) |>
    as.data.frame()
  rownames(est) <- NULL
  print(est)
  cat("\nNumber of iterations:", x$iterations, "\nFinal loglikelihood:",
    round(x$finallogLik, 3), "\n") 
}


calcStandardError <- function(object, ...) {
  UseMethod("calcStandardError")
}


#' @export
calcStandardError.modsemLMS <- function(object, ...) {
  negHessian <- object$negHessian
  stdError <- tryCatch({
      sqrt(diag(solve(negHessian)))
    }, error=function(e) {
      NA
    }, warning=function(w) {
       if (grepl("NaN", conditionMessage(w))) {
         suppressWarnings(sqrt(diag(solve(negHessian))))
      } else{
         sqrt(diag(solve(negHessian)))
      }
    })
    if (all(is.na(stdError))) warning("Standard errors could not be computed, because negative Hessian was either not available or singular.")
    if (any(is.nan(stdError))) warning("Standard errors for some coefficients could not be computed.") 
  stdError
}


#' @export
calcStandardError.modsemQML <- function(object, ...) {
  # not correct yet
  return(calcStandardError.modsemLMS(object, ...))
  H <- object$negHessian 
  invH <- solve(H)
  N <- object$object$info$N
  gradient <- gradientLogLikQml(object$emptyModel, object$coefficients)
  J <- outer(gradient, gradient)
  Jstar <- (1 / N) * (invH %*% J %*% invH)
  sqrt(diag(Jstar))
}
