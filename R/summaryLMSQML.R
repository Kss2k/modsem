#' @export
summary.modsemLMS <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}


#' @export 
summary.modsemQML <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}



summaryLmsAndQml <- function(object, H0 = TRUE, verbose = TRUE, 
                              r.squared = TRUE, ...) {
  if (inherits(object, "modsemQML")) method <- "qml" 
  else if (inherits(object, "modsemLMS")) method <- "lms" 

  parTable <- object$parTable
  parTable$pvalue <- format.pval(parTable$pvalue)
  names(parTable) <- c("lhs", "op", "rhs", 
                       "est", "std.error", 
                       "t.value", "p.value", #"P(>|z|)", 
                       "ci.lower", "ci.upper")

  out <- list(parTable = parTable,
              data = object$model$data,
              iterations = object$iterations,
              logLik = object$logLik, 
              D = NULL)
  if (H0) {
    estH0 <- estimateNullModel(object$originalParTable, data = out$data, 
                               method = method, verbose = verbose, ...)
    out$D <- calcD(estH0, object)
  } 

  class(out) <- "summaryModsemLMSQML"
  out
}


#' @export
print.summaryModsemLMSQML <- function(x, digits = 3, ...) {
  cat("\nSummary for model of class", x$model, "\n")
  cat("\nEstimates:\n")
  est <- lapply(x$parTable, function(col)
    if (is.numeric(col)) round(col, digits) else col) |>
    as.data.frame()
  print(est)
  cat("\nNumber of iterations:", x$iterations, "\nFinal loglikelihood:",
    round(x$logLik, 3), "\n") 

  if (!is.null(x$D)) {
    cat("Comparative fit to H0 (no interaction effect)\n",
        paste0("D(", x$D$df, ") = ", format(x$D$D, digits = 2), 
               ", p = ", format.pval(x$D$p), "\n"))
  }
}


estimateNullModel <- function(parTable, data, method = "lms", verbose = 
                              FALSE, ...) {
  strippedParTable <- parTable[!grepl(":", parTable$rhs), ]
  syntax <- parTableToSyntax(strippedParTable)
  if (verbose) cat("Estimating null model\n")
  modsem(syntax, data, method, verbose = verbose, ...)
}


calcD <- function(estH0, estH1) {
  df <- length(estH1$theta) - length(estH0$theta) 
  D <- - 2 * (estH1$logLik - estH0$logLik)

  p <- stats::qchisq(p = D, df = df, log = TRUE)
  list(D = D, df = df, p = p) 
}
