#' @export 
parameter_estimates.modsem_lms <- function(object, ...) {
  object$parTable 
}


#' @export 
parameter_estimates.modsem_qml <- function(object, ...) {
  object$parTable 
}


#' summary.modsem_lms
#'
#' @param object modsem object to summarized
#' @param ... additional arguments (see summaryLmsAndQml)
#' @rdname summary
#' @export
summary.modsem_lms <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}


#' summary.modsem_qml
#'
#' @param object modsem object to summarized
#' @param ... additional arguments (see summaryLmsAndQml)
#' @rdname summary
#' @export 
summary.modsem_qml <- function(object, ...) {
  summaryLmsAndQml(object, ...)
}


#' summaryLmsAndQml
#'
#' @param object modsem object estimated with lms or qml
#' @param H0 should a null model be estimated (used for comparison)
#' @param verbose print progress for the estimation of null model
#' @param r.squared calculate R-squared
#' @param digits number of digits to print
#' @param ... additional arguments
#' @rdname summaryLmsAndQml
#' @export 
summaryLmsAndQml <- function(object, H0 = TRUE, verbose = TRUE, 
                             r.squared = TRUE, digits = 3, ...) {
  if (inherits(object, "modsem_qml")) method <- "qml" 
  else if (inherits(object, "modsem_lms")) method <- "lms" 

  parTable <- object$parTable
  parTable$pvalue <- format.pval(parTable$pvalue, digits = digits)
  names(parTable) <- c("lhs", "op", "rhs", 
                       "est", "std.error", 
                       "t.value", "p.value", #"P(>|z|)", 
                       "ci.lower", "ci.upper")

  out <- list(parTable = parTable,
              data = object$data,
              iterations = object$iterations,
              logLik = object$logLik, 
              D = NULL)

  if (H0) {
    estH0 <- estimateNullModel(object$originalParTable, data = out$data, 
                               cov_syntax = object$model$covModel$syntax,
                               method = method, verbose = verbose, ...)
    out$nullModel <- estH0
    out$D <- calcD(estH0, object)
  } else {
    out$D <- NULL
  }

  if (r.squared) {
    out$r.squared <- calcRsquared(parTable)
    if (H0) out$r.squared$H0 <- calcRsquared(estH0$parTable)
  } else {
    out$r.squared <- NULL
  }

  class(out) <- "summary_lms_qml"
  out
}


#' @export
print.summary_lms_qml <- function(x, digits = 3, ...) {
  cat("\n----Model Summary------------------------------------------------------\n")
  cat("\nNumber of iterations:", x$iterations, "\nFinal loglikelihood:",
    round(x$logLik, 3), "\n") 
  # preformating 
  covarianceRows <- x$parTable$op == "~~" & x$parTable$rhs != x$parTable$lhs

  x$parTable$rhs[covarianceRows] <- 
    apply(x$parTable[covarianceRows, c("rhs", "op", "lhs")], MARGIN = 1, 
          FUN = function(x) stringr::str_c(x, collapse = ""))

  maxWidth <- max(vapply(x$parTable$rhs[x$parTable$rhs != "1"], 
                         FUN.VALUE = numeric(1), FUN = nchar))
  
  # using format causes a bug, if the max string-length between lhs and rhs 
  # is different
  x$parTable$rhs[x$parTable$rhs != "1"] <- 
    stringr::str_pad(x$parTable$rhs[x$parTable$rhs != "1"], 
                     width = maxWidth, side = "left")
  x$parTable$lhs[x$parTable$rhs == "1"] <- 
    stringr::str_pad(x$parTable$lhs[x$parTable$rhs == "1"], 
                     width = maxWidth, side = "left")
  x$parTable[!grepl("lhs|op|rhs", colnames(x$parTable))] <- 
    lapplyDf(x$parTable[!grepl("lhs|op|rhs", colnames(x$parTable))], 
             FUN = formatNumeric,
             digits = digits)
  
  colnames(x$parTable)[[3]] <- "variable"
  if (!is.null(x$D)) {
    cat("\nComparative fit to H0 (no interaction effect)\n",
        paste0("Loglikelihood change = ", 
               formatNumeric(x$D$llChange, digits = 2), "\n"),
        paste0("D(", x$D$df, ") = ", formatNumeric(x$D$D, digits = 2), 
               ", p = ", format.pval(x$D$p, digits = digits), "\n"))

  } 

  if (!is.null(x$r.squared)) {
    r.squared <- x$r.squared
    r.squared$Rsqr <- formatNumeric(r.squared$Rsqr, digits = 3)
    maxWidth <- max(vapply(r.squared$Rsqr, 
                           FUN.VALUE = numeric(1), FUN = nchar))
    r.squared$Rsqr <- stringr::str_pad(r.squared$Rsqr, 
                                        width = maxWidth, side = "left")
    cat("\nR-squared:\n")
    for (i in seq_along(r.squared$eta)) {
      cat("  ", r.squared$eta[[i]], "=", r.squared$Rsqr[[i]], "\n")
    }
    if (!is.null(r.squared$H0)) {
      r.squared$H0$Rsqr <- formatNumeric(r.squared$H0$Rsqr, digits = 3)
      maxWidth <- max(vapply(r.squared$H0$Rsqr, 
                             FUN.VALUE = numeric(1), FUN = nchar))
      r.squared$H0$Rsqr <- stringr::str_pad(r.squared$H0$Rsqr, 
                                             width = maxWidth, side = "left")
      cat("R-squared Null-Model (H0):\n")
      for (i in seq_along(r.squared$H0$eta)) {
        cat("  ", r.squared$H0$eta[[i]], "=", 
            formatNumeric(r.squared$H0$Rsqr[[i]], digits = 2), "\n")
      }

      # Calculate Change (using unformatted Rsquared) 
      r.squared$H0$diff <- formatNumeric(x$r.squared$Rsqr - x$r.squared$H0$Rsqr, 
                                  digits = 3)
      cat("R-squared Change:\n")
      for (i in seq_along(r.squared$H0$eta)) {
        cat("  ", x$r.squared$H0$eta[[i]], "=", r.squared$H0$diff[[i]], "\n")
      }
    }
  }

  cat("\n----Estimates----------------------------------------------------------\n")
  printLambda(x$parTable)
  printGamma(x$parTable)
  printIntercepts(x$parTable)
  printVariances(x$parTable)
}


#' @export
print.modsem_lms <- function(x, digits = 3, ...) {
  parTable <- x$parTable
  parTable$pvalue <- format.pval(parTable$pvalue, digits = digits)
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
print.modsem_qml <- function(x, digits = 3, ...) {
  parTable <- x$parTable
  parTable$pvalue <- format.pval(parTable$pvalue, digits = digits)
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


estimateNullModel <- function(parTable, data, method = "lms", cov_syntax = NULL,
                              verbose = FALSE, hessian = FALSE, ...) {
  strippedParTable <- parTable[!grepl(":", parTable$rhs), ]
  syntax <- parTableToSyntax(strippedParTable)
  if (verbose) cat("Estimating null model\n")
  modsem_lms_qml(syntax, data, method, verbose = verbose, cov_syntax = cov_syntax, 
                 hessian = hessian, ...)
}


calcD <- function(estH0, estH1) {
  df <- length(estH1$theta) - length(estH0$theta) 
  D <- - 2 * (estH1$logLik - estH0$logLik)
  if (D > 0) {
    D <- -D 
    warning2("D is positive, this indicates that the null model ",
            "is better than the alternative model. ", 
            "Returning '-D' instead of 'D'")
  }
  p <- stats::qchisq(p = D, df = df, log = TRUE)
  list(D = D, df = df, p = p, llChange = estH1$logLik - estH0$logLik) 
}


calcRsquared <- function(parTable) {
  parTable <- parTable[c("lhs", "op", "rhs", "est")]
  etas <- parTable$lhs[parTable$op == "~" & 
                       parTable$rhs != "1"] |>
    unique()
  # Calculate Variances of Interaction Terms
  intTerms <- parTable[grepl(":", parTable$rhs), ]
  for (i in seq_len(nrow(intTerms))) {
    # interaction term = XY
    XY <- stringr::str_split_fixed(intTerms[i, "rhs"], ":", 2) 
    varX <- calcCovParTable(parTable, XY[[1]], XY[[1]])
    varY <- calcCovParTable(parTable, XY[[2]], XY[[2]])
    covXY <- calcCovParTable(parTable, XY[[1]], XY[[2]])
    newRow <- data.frame(lhs = intTerms[i, "rhs"],
                         op = "~~",
                         rhs = intTerms[i, "rhs"],
                         est = as.character(varX * varY + covXY ^ 2))
    parTable <- rbind(parTable, newRow)
  }
  # Calculate Variances/R squared of Etas
  variances <- residuals <- Rsqr <- vector("numeric", length(etas))
  for (i in seq_along(etas)) {
    variances[[i]] <- calcCovParTable(parTable, etas[[i]], etas[[i]])
    residuals[[i]] <- parTable$est[parTable$lhs == etas[[i]] & 
                                   parTable$op == "~~" &
                                   parTable$rhs == etas[[i]]] |>
      as.numeric()
    Rsqr[[i]] <- 1 - residuals[[i]] / variances[[i]]
  }
  data.frame(eta = etas, variance = variances, 
             residual = residuals, Rsqr = Rsqr) 
}
