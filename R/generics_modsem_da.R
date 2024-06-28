#' @export
parameter_estimates.modsem_lms <- function(object, ...) {
  object$parTable
}


#' @export
parameter_estimates.modsem_qml <- function(object, ...) {
  object$parTable
}


#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param H0 should a null model be estimated (used for comparison)
#' @param verbose print progress for the estimation of null model
#' @param r.squared calculate R-squared
#' @param digits number of digits to print
#' @param scientific print p-values in scientific notation
#' @param ci print confidence intervals
#' @param standardized print standardized estimates
#' @param loadings print loadings
#' @param regressions print regressions
#' @param covariances print covariances
#' @param intercepts print intercepts
#' @param variances print variances
#' @param ... additional arguments
#' @rdname summary
#' @export
#' @examples
#' \dontrun{
#' m1 <- "
#'  # Outer Model
#'  X =~ x1 + x2 + x3
#'  Y =~ y1 + y2 + y3
#'  Z =~ z1 + z2 + z3
#'
#'  # Inner model
#'  Y ~ X + Z + X:Z
#' "
#'
#' est1 <- modsem(m1, oneInt, "lms")
#' summary(est1, ci = TRUE, scientific = TRUE)
#' }
summary.modsem_lms <- function(object,
                               H0 = TRUE,
                               verbose = TRUE,
                               r.squared = TRUE,
                               digits = 3,
                               scientific = FALSE,
                               ci = FALSE,
                               standardized = FALSE,
                               loadings = TRUE,
                               regressions = TRUE,
                               covariances = TRUE,
                               intercepts = TRUE,
                               variances = TRUE,
                               ...) {
  summaryLmsAndQml(object,
    H0 = H0, verbose = verbose,
    r.squared = r.squared, digits = digits,
    scientific = scientific, ci = ci,
    standardized = standardized, 
    loadings = loadings, regressions = regressions,
    covariances = covariances, intercepts = intercepts,
    variances = variances, ...
  )
}


#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param H0 should a null model be estimated (used for comparison)
#' @param verbose print progress for the estimation of null model
#' @param r.squared calculate R-squared
#' @param digits number of digits to print
#' @param scientific print p-values in scientific notation
#' @param ci print confidence intervals
#' @param standardized print standardized estimates
#' @param loadings print loadings
#' @param regressions print regressions
#' @param covariances print covariances
#' @param intercepts print intercepts
#' @param variances print variances
#' @param ... additional arguments
#' @rdname summary
#' @export
#' @examples
#' \dontrun{
#' m1 <- "
#'  # Outer Model
#'  X =~ x1 + x2 + x3
#'  Y =~ y1 + y2 + y3
#'  Z =~ z1 + z2 + z3
#'
#'  # Inner model
#'  Y ~ X + Z + X:Z
#' "
#'
#' est1 <- modsem(m1, oneInt, "qml")
#' summary(est1, ci = TRUE, scientific = TRUE)
#' }
summary.modsem_qml <- function(object,
                               H0 = TRUE,
                               verbose = TRUE,
                               r.squared = TRUE,
                               digits = 3,
                               scientific = FALSE,
                               ci = FALSE,
                               standardized = FALSE, 
                               loadings = TRUE,
                               regressions = TRUE,
                               covariances = TRUE,
                               intercepts = TRUE,
                               variances = TRUE,
                               ...) {
  summaryLmsAndQml(object,
    H0 = H0, verbose = verbose,
    r.squared = r.squared, digits = digits,
    scientific = scientific, ci = ci,
    standardized = standardized, 
    loadings = loadings, regressions = regressions,
    covariances = covariances, intercepts = intercepts,
    variances = variances, ...
  )
}


summaryLmsAndQml <- function(object,
                             H0 = TRUE,
                             verbose = TRUE,
                             r.squared = TRUE,
                             digits = 3,
                             scientific = FALSE,
                             ci = FALSE,
                             standardized = FALSE, 
                             loadings = TRUE,
                             regressions = TRUE,
                             covariances = TRUE,
                             intercepts = TRUE,
                             variances = TRUE,
                             includeVarIntTerm = FALSE, 
                             ...) {
  if (inherits(object, "modsem_qml")) 
    method <- "qml"
  else if (inherits(object, "modsem_lms")) 
    method <- "lms"

  if (standardized) 
    parTable <- standardized_estimates(object)
  else 
    parTable <- parameter_estimates(object)

  if (!includeVarIntTerm) 
    parTable <- removeInteractionVariances(parTable)
  
  args <- object$args
  out <- list(
    parTable = parTable,
    data = object$data,
    iterations = object$iterations,
    logLik = object$logLik,
    AIC = object$AIC,
    D = NULL,
    N = NROW(object$data),
    method = method,
    optimizer = object$optimizer,
    quad = object$info.quad,
    type.se = object$type.se, 
    type.estimates = ifelse(standardized, "standardized", object$type.estimates),
    information = object$information
  )

  if (H0) {
    estH0 <- estimateNullModel(object$originalParTable,
      data = out$data,
      method = method,
      cov.syntax = object$model$covModel$syntax,
      verbose = verbose,
      calc.se = FALSE,
      double = args$double,
      standardize = args$standardize,
      standardize.out = args$standardize.out,
      mean.observed = args$mean.observed
    )

    out$nullModel <- estH0
    if (is.null(estH0)) {
      warning2("Comparative fit to H0 will not be calculated.")
      H0 <- FALSE
      out$D <- NULL 
      out$fitH0 <- NULL

    } else {
      out$D <- compare_fit(estH0, object)
      out$fitH0 <- fit_modsem_da(estH0)
    }
  } else {
    out$D <- NULL
  }

  if (r.squared) {
    out$r.squared <- calcRsquared(parTable)
    if (H0) out$r.squared$H0 <- calcRsquared(estH0$parTable)
  } else {
    out$r.squared <- NULL
  }

  out$format <- list(
    digits = digits,
    scientific = scientific,
    ci = ci,
    loadings = loadings,
    regressions = regressions,
    covariances = covariances,
    intercepts = intercepts,
    variances = variances
  )

  class(out) <- "summary_da"
  out
}


#' @export
print.summary_da <- function(x, digits = 3, ...) {
  width.out <- getWidthPrintedParTable(x$parTable,
                                       scientific = x$format$scientific,
                                       ci = x$format$ci,
                                       digits = x$format$digits,
                                       loadings = x$format$loadings,
                                       regressions = x$format$regressions,
                                       covariances = x$format$covariances,
                                       intercepts = x$format$intercepts,
                                       variances = x$format$variances)
  cat("\nmodsem (version 1.0.1):\n")
  names <- c("Estimator", "Optimization method", "Number of observations", 
             "Number of iterations", "Final loglikelihood", "Akaike (AIC)")
  values <- c(stringr::str_to_upper(c(x$method, x$optimizer)), 
              x$N, x$iterations, round(x$logLik, 3), round(x$AIC, 3))
  cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                   width.out = width.out), "\n")
  
  if (!is.null(x$quad)) {
    cat("Numerical Integration:\n")
    names <- c("Points of integration (per dim)", "Dimensions",
               "Total points of integration")
    values <- c(x$quad$nodes.dim, x$quad$dim, x$quad$nodes.total)
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                     width.out = width.out), "\n")

  }

  if (!is.null(x$D)) {
    cat("Comparative fit to H0 (no interaction effect)\n")
    names <- c("Loglikelihood change", 
               "Bayesian (BIC)",
               "Difference test (D)", 
               "Degrees of freedom (D)", "P-value (D)")
    values <- c(formatNumeric(x$D$llChange, digits = 2), 
                formatNumeric(x$D$BIC, digits = 2),
                formatNumeric(x$D$D, digits = 2),
                x$D$df, 
                format.pval(x$D$p, digits = digits))
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                     width.out = width.out), "\n")
    
    cat("Fit Measures for H0:\n")
    names <- c("Chi-square", "Degrees of Freedom (Chi-square)", 
               "P-value (Chi-square)", "RMSEA") 
    values <- c(formatNumeric(x$fitH0$chisq.value, digits = 2), 
                x$fitH0$chisq.df,
                format.pval(x$fitH0$chisq.pvalue, digits = digits),
                formatNumeric(x$fitH0$RMSEA, digits = 3))
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                     width.out = width.out), "\n")

  }

  if (!is.null(x$r.squared)) {
    r.squared <- x$r.squared
    r.squared$Rsqr <- formatNumeric(r.squared$Rsqr, digits = 3)
    maxWidth <- max(vapply(r.squared$Rsqr, FUN.VALUE = numeric(1), FUN = nchar))
    r.squared$Rsqr <- 
      stringr::str_pad(r.squared$Rsqr, width = maxWidth, side = "left")

    cat("R-Squared:\n")
    names <- r.squared$eta
    values <- character(length(names)) 
    for (i in seq_along(r.squared$eta)) {
      values[[i]] <- r.squared$Rsqr[[i]]
    }
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                     width.out = width.out))

    if (!is.null(r.squared$H0)) {
      r.squared$H0$Rsqr <- formatNumeric(r.squared$H0$Rsqr, digits = 3)
      maxWidth <- 
        max(vapply(r.squared$H0$Rsqr, FUN.VALUE = numeric(1), FUN = nchar))
      r.squared$H0$Rsqr <- 
        stringr::str_pad(r.squared$H0$Rsqr, width = maxWidth, side = "left")

      cat("R-Squared Null-Model (H0):\n")
      names <- r.squared$H0$eta
      for (i in seq_along(names)) {
          values[[i]] <- formatNumeric(r.squared$H0$Rsqr[[i]], digits = 3)
      }
      cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                       width.out = width.out))

      # Calculate Change (using unformatted Rsquared)
      r.squared$H0$diff <- 
        formatNumeric(x$r.squared$Rsqr - x$r.squared$H0$Rsqr, digits = 3)
      cat("R-Squared Change:\n")
      for (i in seq_along(names)) {
        values[[i]] <- r.squared$H0$diff[[i]]
      }
      cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                       width.out = width.out))
    }
  }

  cat("\nParameter Estimates:\n")
  names <- c("Coefficients", "Information", "Standard errors")
  values <- c(x$type.estimates, x$information, x$type.se)
  cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                   width.out = width.out), "\n")

  printParTable(x$parTable,
    scientific = x$format$scientific,
    ci = x$format$ci,
    digits = x$format$digits,
    loadings = x$format$loadings,
    regressions = x$format$regressions,
    covariances = x$format$covariances,
    intercepts = x$format$intercepts,
    variances = x$format$variances)
}


#' @export
print.modsem_lms <- function(x, digits = 3, ...) {
  parTable <- x$parTable
  parTable$p.value <- format.pval(parTable$p.value, digits = digits)
  names(parTable) <- c(
    "lhs", "op", "rhs",
    "est", "std.error",
    "z.value", "p.value", # "P(>|z|)",
    "ci.lower", "ci.upper"
  )
  est <- lapply(parTable, function(col) {
    if (is.numeric(col)) round(col, digits) else col
  }) |>
    as.data.frame()
  print(est)
}


#' @export
print.modsem_qml <- function(x, digits = 3, ...) {
  parTable <- x$parTable
  parTable$p.value <- format.pval(parTable$p.value, digits = digits)
  names(parTable) <- c(
    "lhs", "op", "rhs",
    "est", "std.error",
    "z.value", "p.value", # "P(>|z|)",
    "ci.lower", "ci.upper"
  )
  est <- lapply(parTable, function(col) {
    if (is.numeric(col)) round(col, digits) else col
  }) |>
    as.data.frame()
  print(est)
}


estimateNullModel <- function(parTable,
                              data,
                              method = "lms",
                              cov.syntax = NULL,
                              verbose = FALSE,
                              calc.se = FALSE,
                              standardize = NULL,
                              standardize.out = NULL,
                              mean.observed = NULL,
                              double = NULL,
                              ...) {
  tryCatch({
      strippedParTable <- removeUnknownLabels(parTable[!grepl(":", parTable$rhs), ])
      if (NROW(strippedParTable) == NROW(parTable)) {
        return(NULL)
      }

      syntax <- parTableToSyntax(strippedParTable)
      if (verbose) cat("Estimating null model\n")
      modsem_da(syntax, data, method, 
                verbose = verbose, 
                cov.syntax = cov.syntax, 
                calc.se = calc.se, 
                double = double,
                standardize = standardize, 
                standardize.out = standardize.out,
                mean.observed = mean.observed, ...)
    },
    error = function(e) {
      warning2(
        "Null model could not be estimated. ",
        "Error message: ", e$message
      )
      NULL
    }
  )
}



#' compare model fit for qml and lms models
#'
#' @param estH0 object of class `modsem_lms` or `modsem_qml` representing the
#' null hypothesis model
#' @param estH1 object of class `modsem_lms` or `modsem_qml` representing the
#' @description Compare the fit of two models using the likelihood ratio test.
#' `estH0` representing the null 
#' hypothesis model, and `estH1` the alternative hypothesis model. Importantly, 
#' the function assumes that `estH0` does not have more free parameters 
#' (i.e., degrees of freedom) than `estH1`.
#' alternative hypothesis model
#' @rdname compare_fit
#' @export
#' @examples
#' \dontrun{
#' H0 <- "
#'  # Outer Model
#'  X =~ x1 + x2 + x3
#'  Y =~ y1 + y2 + y3
#'  Z =~ z1 + z2 + z3
#'
#'  # Inner model
#'  Y ~ X + Z
#' "
#'
#' estH0 <- modsem(m1, oneInt, "lms")
#'
#' H1 <- "
#'  # Outer Model
#'  X =~ x1 + x2 + x3
#'  Y =~ y1 + y2 + y3
#'  Z =~ z1 + z2 + z3
#'
#'  # Inner model
#'  Y ~ X + Z + X:Z
#' "
#'
#' estH1 <- modsem(m1, oneInt, "lms")
#' compare_fit(estH0, estH1)
#' }
#' @export
compare_fit <- function(estH0, estH1) {
  if (is.null(estH0) || is.null(estH1)) {
    return(NULL)
  }
  df <- length(coef(estH1)) - length(coef(estH0))
  D <- -2 * (estH0$logLik - estH1$logLik)
  BIC <- 2 * (estH1$logLik - estH0$logLik) + df * log(NROW(estH1$data))
  p <- stats::pchisq(D, df = df, lower.tail = FALSE, log.p = TRUE)
  list(D = D, BIC = BIC, df = df, p = p, llChange = estH1$logLik - estH0$logLik)
}


calcRsquared <- function(parTable) {
  parTable <- var_interactions.data.frame(parTable)
  etas <- unique(parTable$lhs[parTable$op == "~" &
    parTable$rhs != "1"])

  # Calculate Variances/R squared of Etas
  variances <- residuals <- Rsqr <- vector("numeric", length(etas))
  for (i in seq_along(etas)) {
    variances[[i]] <- calcCovParTable(etas[[i]], etas[[i]], parTable)
    residuals[[i]] <- parTable$est[parTable$lhs == etas[[i]] &
      parTable$op == "~~" &
      parTable$rhs == etas[[i]]] |>
      as.numeric()
    Rsqr[[i]] <- 1 - residuals[[i]] / variances[[i]]
  }
  data.frame(
    eta = etas, variance = variances,
    residual = residuals, Rsqr = Rsqr
  )
}


#' @export
var_interactions.modsem_lms <- function(object, ...) {
  var_interactions.data.frame(parameter_estimates(object))
}


#' @export
var_interactions.modsem_qml <- function(object, ...) {
  var_interactions.data.frame(parameter_estimates(object))
}


#' @export 
standardized_estimates.modsem_lms <- function(object, ...) {
  standardized_estimates.data.frame(parameter_estimates(object))
}


#' @export 
standardized_estimates.modsem_qml <- function(object, ...) {
  standardized_estimates.data.frame(parameter_estimates(object))
}


#' @export 
modsem_inspect.modsem_lms <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "default"
  modsem_inspect_da(object, what = what, ...)
}


#' @export 
modsem_inspect.modsem_qml <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "default"
  modsem_inspect_da(object, what = what, ...)
}


#' @export 
#' @importFrom stats vcov
vcov.modsem_qml <- function(object, ...) {
  modsem_inspect_da(object, what = "vcov")[[1]]
}


#' @export 
#' @importFrom stats vcov
vcov.modsem_lms <- function(object, ...) {
  modsem_inspect_da(object, what = "vcov")[[1]]
}


#' @export 
#' @importFrom stats coefficients
coefficients.modsem_qml <- function(object, ...) {
  modsem_inspect_da(object, what = "coefficients")[[1]]
}


#' @export 
#' @importFrom stats coefficients
coefficients.modsem_lms <- function(object, ...) {
  modsem_inspect_da(object, what = "coefficients")[[1]]
}

#' @export 
#' @importFrom stats coef
coef.modsem_qml <- function(object, ...) {
  modsem_inspect_da(object, what = "coefficients")[[1]]
}


#' @export 
#' @importFrom stats coef
coef.modsem_lms <- function(object, ...) {
  modsem_inspect_da(object, what = "coefficients")[[1]]
}
