#' @export
parameter_estimates.modsem_da <- function(object, ...) {
  object$parTable
}


#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param H0 should a null model be estimated (used for comparison)
#' @param verbose print progress for the estimation of null model
#' @param r.squared calculate R-squared
#' @param adjusted.stat should sample size corrected/adjustes AIC and BIC be reported?
#' @param digits number of digits to print
#' @param scientific print p-values in scientific notation
#' @param ci print confidence intervals
#' @param standardized print standardized estimates
#' @param loadings print loadings
#' @param regressions print regressions
#' @param covariances print covariances
#' @param intercepts print intercepts
#' @param variances print variances
#' @param var.interaction if FALSE (default) variances for interaction terms will be removed (if present)
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
summary.modsem_da <- function(object,
                              H0 = TRUE,
                              verbose = TRUE,
                              r.squared = TRUE,
                              adjusted.stat = FALSE,
                              digits = 3,
                              scientific = FALSE,
                              ci = FALSE,
                              standardized = FALSE,
                              loadings = TRUE,
                              regressions = TRUE,
                              covariances = TRUE,
                              intercepts = TRUE,
                              variances = TRUE,
                              var.interaction = FALSE,
                              ...) {
  method <- object$method

  if (standardized) {
    parTable <- standardized_estimates(object)
  } else {
    parTable <- parameter_estimates(object)
  }

  if (!var.interaction) {
    parTable <- removeInteractionVariances(parTable)
  }
  
  args <- object$args
  out <- list(
    parTable       = parTable,
    data           = object$data,
    iterations     = object$iterations,
    logLik         = object$logLik,
    fit            = fit_modsem_da(object, chisq = FALSE),
    D              = NULL,
    N              = NROW(object$data),
    method         = method,
    optimizer      = object$optimizer,
    quad           = object$info.quad,
    type.se        = object$type.se, 
    type.estimates = ifelse(standardized, "standardized", object$type.estimates),
    information    = object$information
  )

  if (H0) {
    estH0 <- estimateNullModel(
        object$originalParTable,
        data            = out$data,
        method          = method,
        cov.syntax      = object$model$covModel$syntax,
        verbose         = verbose,
        calc.se         = FALSE,
        double          = args$double,
        standardize     = args$standardize,
        standardize.out = args$standardize.out,
        mean.observed   = args$mean.observed
    )

    out$nullModel <- estH0
    if (is.null(estH0)) {
      warning2("Comparative fit to H0 will not be calculated.")
      H0        <- FALSE
      out$D     <- NULL
      out$fitH0 <- NULL

    } else {
      out$D     <- compare_fit(estH0, object)
      out$fitH0 <- fit_modsem_da(estH0)
    }
  } else {
    out$D <- NULL
  }

  if (r.squared) {
    out$r.squared <- calcRsquared(parTable)
    if (H0) out$r.squared$H0 <- calcRsquared(estH0$parTable)
  } else out$r.squared <- NULL

  out$format <- list(
    digits        = digits,
    scientific    = scientific,
    adjusted.stat = adjusted.stat,
    ci            = ci,
    loadings      = loadings,
    regressions   = regressions,
    covariances   = covariances,
    intercepts    = intercepts,
    variances     = variances
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
  cat(paste0("\nmodsem (version ", PKG_INFO$version, "):\n"))
  names <- c("Estimator", "Optimization method", "Number of observations", 
             "Number of iterations", "Loglikelihood", 
             "Akaike (AIC)", "Bayesian (BIC)")
  values <- c(stringr::str_to_upper(c(x$method, x$optimizer)), 
              x$N, x$iterations, round(x$logLik, 2), round(x$fit$AIC, 2), 
              round(x$fit$BIC, 2))
  if (x$format$adjusted.stat) {
    names  <- c(names, "Corrected Akaike (AICc)", "Adjusted Bayesian (aBIC)")
    values <- c(values, round(x$fit$AICc, 2), round(x$fit$aBIC, 2))
  }

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
    cat("Fit Measures for H0:\n")
    names <- c("Loglikelihood", "Akaike (AIC)", "Bayesian (BIC)") 
    values <- c(round(x$nullModel$logLik), round(x$fitH0$AIC, 2), round(x$fitH0$BIC, 2))

    if (x$format$adjusted.stat) {
      names <- c(names, "Corrected Akaike (AICc)", "Adjusted Bayesian (aBIC)")
      values <- c(values, round(x$fitH0$AICc, 2), round(x$fitH0$aBIC, 2))
    }

    names <- c(names, "Chi-square", "Degrees of Freedom (Chi-square)", 
               "P-value (Chi-square)", "RMSEA")
    values <- c(values, formatNumeric(x$fitH0$chisq.value, digits = 2), 
                x$fitH0$chisq.df,
                formatPval(x$fitH0$chisq.pvalue, scientific = x$format$scientific),
                formatNumeric(x$fitH0$RMSEA, digits = 3))
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ", 
                     width.out = width.out), "\n")

    cat("Comparative fit to H0 (no interaction effect)\n")
    names <- c("Loglikelihood change", 
               "Difference test (D)", 
               "Degrees of freedom (D)", "P-value (D)")
    values <- c(formatNumeric(x$D$llChange, digits = 2), 
                formatNumeric(x$D$D, digits = 2),
                x$D$df, 
                formatPval(x$D$p, scientific = x$format$scientific))
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
                scientific  = x$format$scientific,
                ci          = x$format$ci,
                digits      = x$format$digits,
                loadings    = x$format$loadings,
                regressions = x$format$regressions,
                covariances = x$format$covariances,
                intercepts  = x$format$intercepts,
                variances   = x$format$variances)
}


#' @export
print.modsem_da <- function(x, digits = 3, ...) {
  parTable         <- x$parTable
  parTable$p.value <- format.pval(parTable$p.value, digits = digits)
  names(parTable)  <- c("lhs", "op", "rhs", "label", "est", "std.error",
                        "z.value", "p.value", "ci.lower", "ci.upper")
  est <- lapply(parTable, FUN = function(col) 
                if (is.numeric(col)) round(col, digits) else col) |>
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
    if (NROW(strippedParTable) == NROW(parTable)) return(NULL)

    syntax <- parTableToSyntax(strippedParTable)
    if (verbose) cat("Estimating null model\n")
    modsem_da(syntax, data = data, 
              method = method, 
              verbose = verbose, 
              cov.syntax = cov.syntax, 
              calc.se = calc.se, 
              double = double,
              standardize = standardize, 
              standardize.out = standardize.out,
              mean.observed = mean.observed, ...)
  },
  error = function(e) {
    warning2("Null model could not be estimated. ",
             "Error message: ", e$message)
    NULL
  })
}



#' compare model fit for qml and lms models
#'
#' @param estH0 object of class `modsem_da` representing the
#' null hypothesis model
#' @param estH1 object of class `modsem_da` representing the
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
  df <- length(coef(estH1, type = "free")) - length(coef(estH0, type = "free"))
  D <- -2 * (estH0$logLik - estH1$logLik)
  p <- stats::pchisq(D, df = df, lower.tail = FALSE, log.p = FALSE)
  list(D = D, df = df, p = p, llChange = estH1$logLik - estH0$logLik)
}


calcRsquared <- function(parTable) {
  parTable <- var_interactions.data.frame(parTable)
  etas     <- unique(parTable$lhs[parTable$op == "~"])

  # Calculate Variances/R squared of Etas
  variances <- residuals <- Rsqr <- vector("numeric", length(etas))
  for (i in seq_along(etas)) {
    variances[[i]] <- calcCovParTable(etas[[i]], etas[[i]], parTable)
    residuals[[i]] <- as.numeric(parTable$est[parTable$lhs == etas[[i]] &
                                              parTable$op == "~~" &
                                              parTable$rhs == etas[[i]]])
    Rsqr[[i]] <- 1 - residuals[[i]] / variances[[i]]
  }

  data.frame(eta = etas, variance = variances,
             residual = residuals, Rsqr = Rsqr)
}


#' @export
var_interactions.modsem_da <- function(object, ...) {
  var_interactions.data.frame(parameter_estimates(object))
}


#' @export 
standardized_estimates.modsem_da <- function(object, ...) {
  standardized_estimates.data.frame(parameter_estimates(object))
}


#' @export 
modsem_inspect.modsem_da <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "default"
  modsem_inspect_da(object, what = what, ...)
}


#' @export 
#' @importFrom stats vcov
vcov.modsem_da <- function(object, ...) {
  modsem_inspect_da(object, what = "vcov")[[1]]
}


#' @export 
#' @importFrom stats coefficients
coefficients.modsem_da <- function(object, type = "all", ...) {
  what <- ifelse(type == "all", yes = "all.coefficients",
                 no = "free.coefficients")
  modsem_inspect_da(object, what = what)[[1]]
}


#' @export 
#' @importFrom stats coef
coef.modsem_da <- function(object, type = "all", ...) {
  coefficients.modsem_da(object, type = type, ...)
}


#' Wrapper for vcov
#'
#' @param object fittet model to inspect
#' @param ... additional arguments
#' @description wrapper for vcov, to be used with modsem::vcov_modsem_da, 
#' since vcov is not in the namespace of modsem, but stats
#' @export
vcov_modsem_da <- function(object, ...) {
  vcov.modsem_da(object, ...)
}


#' Wrapper for coef
#'
#' @param object fittet model to inspect
#' @param ... additional arguments
#' @description wrapper for coef, to be used with modsem::coef_modsem_da, 
#' since coef is not in the namespace of modsem, but stats
#' @export
coef_modsem_da <- function(object, ...) {
  coef.modsem_da(object, ...)
}


#' @export
#' @importFrom stats nobs
nobs.modsem_da <- function(object, ...) {
  modsem_inspect_da(object, what = "N", ...)[[1]]
}
