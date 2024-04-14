#' @export
summary.modsemLMS <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}


#' @export 
summary.modsemQML <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}


summaryLmsAndQml <- function(object, H0 = TRUE, verbose = TRUE, 
                              r.squared = TRUE, 
                              ...) {
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

  class(out) <- "summaryModsemLMS_QML"
  out
}


#' @export
print.summaryModsemLMS_QML <- function(x, digits = 3, ...) {
  cat("\nModel summary:\n")
  cat("\nNumber of iterations:", x$iterations, "\nFinal loglikelihood:",
    round(x$logLik, 3), "\n") 
  maxWidth <- max(vapply(x$parTable$rhs[x$parTable$rhs != "1"], 
                         FUN.VALUE = numeric(1), FUN = nchar))
  x$parTable$rhs[x$parTable$rhs != "1"] <- 
    stringr::str_pad(x$parTable$rhs[x$parTable$rhs != "1"], 
                     width = maxWidth, side = "left")
  colnames(x$parTable)[[3]] <- "variable"
  if (!is.null(x$D)) {
    cat("\nComparative fit to H0 (no interaction effect)\n",
        paste0("D(", x$D$df, ") = ", format(x$D$D, digits = 2), 
               ", p = ", format.pval(x$D$p), "\n"))
  }

  cat("\nEstimates:\n")
  printLambda(x$parTable)
  printGamma(x$parTable)
  printIntercepts(x$parTable)
  printVariances(x$parTable)

}


#' @export
print.modsemLMS <- function(x, digits = 3, ...) {
  parTable <- x$parTable
  parTable$pvalue <- format.pval(parTable$pvalue)
  names(parTable) <- c("lhs", "op", "rhs", 
                       "est", "std.error", 
                       "t.value", "p.value", #"P(>|z|)", 
                       "ci.lower", "ci.upper")
  est <- lapply(parTable, function(col)
    if (is.numeric(col)) round(col, digits) else col) |>
    as.data.frame()
  print(est)
}


#' @export
print.modsemQML <- function(x, digits = 3, ...) {
  parTable <- x$parTable
  parTable$pvalue <- format.pval(parTable$pvalue)
  names(parTable) <- c("lhs", "op", "rhs", 
                       "est", "std.error", 
                       "t.value", "p.value", #"P(>|z|)", 
                       "ci.lower", "ci.upper")
  est <- lapply(parTable, function(col)
    if (is.numeric(col)) round(col, digits) else col) |>
    as.data.frame()
  print(est)
}


strippedColnames <- c("variable", "est", "std.error", 
                        "t.value", "p.value", "ci.lower", "ci.upper")


printLambda <- function(parTable, digits = 3) {
  # get latent variableiables 
  cat("\nLoadings:\n")
  latentVars <- unique(parTable[parTable$op == "=~", "lhs"])
  for (lv in latentVars) {
    cat("", paste0(" ", lv, ":"), "\n  ")
    out <- capturePrint(parTable[parTable$lhs == lv & 
                                 parTable$op == "=~", 
                                 strippedColnames], 
                        row.names = FALSE,
                        digits = digits)
    cat(stringr::str_replace_all(out, "\n", "\n  "), "\n")
  }
}


printGamma <- function(parTable, digits = 3) {
  # get latent endogenousendogenous  variableiables 
  cat("\nRegressions: \n")
  latentVars <- unique(parTable[parTable$op == "~" &
                                parTable$variable != "1", "lhs"])
  for (lv in latentVars) {
    cat("", paste0(" ", lv, ":"), "\n  ")
    out <- capturePrint(parTable[parTable$lhs == lv &
                                 parTable$op == "~" & 
                                 parTable$variable != "1", 
                                 strippedColnames],
                        row.names = FALSE, 
                        digits = digits)
    cat(stringr::str_replace_all(out, "\n", "\n  "), "\n")
  }
}


printIntercepts <- function(parTable, digits = 3) {
  # get latent endogenousendogenous  variableiables 
  cat("\nIntercepts: \n  ")
  pt <- parTable[parTable$op == "~" &
                 parTable$variable == "1", 
                 !grepl("variable", colnames(parTable))]
  colnames(pt)[[1]] <- "variable"

  out <- capturePrint(pt[ , strippedColnames],
                      row.names = FALSE,
                      digits = digits)
  cat(stringr::str_replace_all(out, "\n", "\n  "), "\n")
}


printVariances <- function(parTable, digits = 3) {
  # get latent endogenousendogenous  variableiables 
  cat("\nVariances:\n  ")
  out <- capturePrint(parTable[parTable$op == "~~", 
                               strippedColnames],
                      row.names = FALSE,
                      digits = digits)
  cat(stringr::str_replace_all(out, "\n", "\n  "), "\n")
}


estimateNullModel <- function(parTable, data, method = "lms", 
                              verbose = FALSE,...) {
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
