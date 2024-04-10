#' @param object 
#'
#' @param ... 
#'
#' @export
summary.modsemLMS <- function(object, ...) {

  # estimates
  est <- object$coefficients

  # standard errors
  s.error <- calcStandardError(object$negHessian)
  tvalue <- est / s.error
  pvalue <- 2 * pnorm(-abs(tvalue))

  coefTable <- data.frame(Parameter = names(est),
                          Estimate = est, `Std. Error` = s.error,
                          `z value` = tvalue, 
                          `Pr(>|z|)` = format.pval(pvalue))
  ans <- list(model = object$model.class,
              estimates = coefTable,
              iterations = object$iterations,
              finallogLik = object$objective)

  class(ans) <- "summaryModsemLMS"
  ans
}

#' print lms
#' @param x 
#'
#' @param digits 
#' @param ... 
#'
#' @export
print.summaryModsemLMS <- function(x, digits = 3, ...) {
    
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


calcStandardError <- function(neg.hessian) {
  s.error <- tryCatch({
      sqrt(diag(solve(neg.hessian)))
    }, error=function(e) {
      NA
    }, warning=function(w) {
       if (grepl("NaN", conditionMessage(w))) {
         suppressWarnings(sqrt(diag(solve(neg.hessian))))
      } else{
         sqrt(diag(solve(neg.hessian)))
      }
    })
    if (all(is.na(s.error))) warning("Standard errors could not be computed, because negative Hessian was either not available or singular.")
    if (any(is.nan(s.error))) warning("Standard errors for some coefficients could not be computed.") 
  s.error
}

