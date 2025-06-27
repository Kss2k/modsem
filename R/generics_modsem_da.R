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
#' @param monte.carlo should Monte Carlo bootstrapped standard errors be used? Only 
#'   relevant if \code{standardized = TRUE}. If \code{FALSE} delta method is used instead.
#' @param mc.reps number of Monte Carlo repetitions. Only relevant if \code{monte.carlo = TRUE}, 
#'   and \code{standardized = TRUE}.
#' @param loadings print loadings
#' @param regressions print regressions
#' @param covariances print covariances
#' @param intercepts should intercepts be included in the output?
#' If \code{standardized = TRUE} intercepts will by default be excluded. 
#' @param variances print variances
#' @param var.interaction if FALSE variances for interaction terms will be removed 
#' from the output
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
                              verbose = interactive(),
                              r.squared = TRUE,
                              adjusted.stat = FALSE,
                              digits = 3,
                              scientific = FALSE,
                              ci = FALSE,
                              standardized = FALSE,
                              monte.carlo = FALSE,
                              mc.reps = 10000,
                              loadings = TRUE,
                              regressions = TRUE,
                              covariances = TRUE,
                              intercepts = !standardized,
                              variances = TRUE,
                              var.interaction = FALSE,
                              ...) {
  method <- object$method

  if (standardized) {
    parTable <- standardized_estimates(object, intercepts = intercepts, 
                                       monte.carlo = monte.carlo, mc.reps = mc.reps)
    parTableR2 <- parTable

  } else {
    parTable   <- parameter_estimates(object)
    parTableR2 <- var_interactions(centerInteraction(parTable)) # easier to calculate 
                                                                # R-squared when 
                                                                # interactions are mean centeredd
  }

  if (!var.interaction) {
    parTable.out <- removeInteractionVariances(parTable)
  } else parTable.out <- parTable

  args <- object$args
  out <- list(
    parTable       = parTable.out,
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
    if (any(grepl(":", parTable$rhs))) cat("Estimating baseline model (H0)\n")
    est_h0 <- estimate_h0(object, calc.se = FALSE, warn_no_interaction = FALSE)
    
    out$nullModel <- est_h0
    if (is.null(est_h0)) {
      warning2("Comparative fit to H0 will not be calculated.", immediate. = FALSE)
      H0        <- FALSE
      out$D     <- NULL
      out$fitH0 <- NULL

    } else {
      out$D     <- compare_fit(est_h1 = object, est_h0 = est_h0)
      out$fitH0 <- fit_modsem_da(est_h0)
    }
  } else {
    out$D <- NULL
  }

  if (r.squared) {
    out$r.squared <- calcRsquared(parTableR2)
    if (H0) out$r.squared$H0 <- calcRsquared(est_h0$parTable) # no need to center interactions
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
  cat(paste0("\nmodsem (version ", PKG_INFO$version, "):\n\n"))
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
    cat("Fit Measures for Baseline Model (H0):\n")
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

    cat("Comparative Fit to H0 (LRT test):\n")
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

    cat("R-Squared Interaction Model (H1):\n")
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

      cat("R-Squared Baseline Model (H0):\n")
      names <- r.squared$H0$eta
      for (i in seq_along(names)) {
          values[[i]] <- formatNumeric(r.squared$H0$Rsqr[[i]], digits = 3)
      }
      cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                       width.out = width.out))

      # Calculate Change (using unformatted Rsquared)
      r.squared$H0$diff <-
        formatNumeric(x$r.squared$Rsqr - x$r.squared$H0$Rsqr, digits = 3)
      cat("R-Squared Change (H1 - H0):\n")
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


calcRsquared <- function(parTable, recalc.vars = FALSE) {
  if (recalc.vars) parTable <- var_interactions.data.frame(parTable)

  etas     <- unique(parTable$lhs[parTable$op == "~"])

  # Calculate Variances/R squared of Etas
  residuals_df <- parTable[parTable$lhs %in% etas & 
                           parTable$rhs == parTable$lhs, ]
  residuals <- structure(residuals_df$est, names=residuals_df$lhs)[etas]
  variances <- calcVarParTable(etas, parTable)
  Rsqr <- 1 - residuals / variances

  data.frame(eta = etas, variance = variances,
             residual = residuals, Rsqr = Rsqr)
}


#' @export
var_interactions.modsem_da <- function(object, ...) {
  var_interactions.data.frame(parameter_estimates(object), ...)
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


#' @export
#' @importFrom stats nobs
nobs.modsem_da <- function(object, ...) {
  modsem_inspect_da(object, what = "N", ...)[[1]]
}


#' @describeIn standardized_estimates Method for `modsem_da` objects
#'
#' @param monte.carlo Logical. If \code{TRUE}, use Monte Carlo simulation to estimate
#' standard errors; if \code{FALSE}, use the delta method (default).
#' @param mc.reps Number of Monte Carlo repetitions. Default is 10000.
#' @param tolerance.zero Threshold below which standard errors are set to \code{NA}.
#'
#' @export
standardized_estimates.modsem_da <- function(object, 
                                             monte.carlo = FALSE, 
                                             mc.reps = 10000, 
                                             tolerance.zero = 1e-10, ...) {
  stdSolution <- standardizedSolutionCOEFS(
    object, 
    monte.carlo = monte.carlo, 
    mc.reps = mc.reps, 
    tolerance.zero = tolerance.zero, ...
  )
  
  stdSolution$parTable
}


#' @describeIn modsem_predict 
#' Computes (optionally standardised) factor scores via the
#'   regression method using the baseline model unless \code{H0 = FALSE}.
#'
#' @param object \code{\link{modsem_da}} object
#' @param standardized Logical. If \code{TRUE}, return standardized factor scores.
#' @param H0 Logical. If \code{TRUE} (default), use the baseline model to compute factor scores.
#'   If \code{FALSE}, use the model specified in \code{object}. Using \code{H0 = FALSE} is not recommended!
#' @export
modsem_predict.modsem_da <- function(object, standardized = FALSE, H0 = TRUE, ...) {
  modelH1 <- object 

  if (H0) {
    modelH0 <- estimate_h0(modelH1, calc.se = FALSE, warn_no_interaction = FALSE, 
                           verbose = FALSE)
  
    if (is.null(modelH0)) modelH0 <- modelH1
  } else modelH0 <- modelH1

  fitH0      <- fit_modsem_da(modelH0, chisq = TRUE)
  parTableH1 <- parameter_estimates(modelH1)
  parTableH0 <- parameter_estimates(modelH0)

  lVs   <- getLVs(parTableH1)
  sigma <- fitH0$sigma.expected

  sigma.inv <- GINV(sigma)
  lambda    <- getLambdaParTable(parTableH0, rows = colnames(sigma), cols = lVs)
  X         <- apply(as.matrix(modelH0$data), MARGIN = 2, FUN = \(x) x - mean(x, na.rm = TRUE))
  X         <- X[, colnames(sigma), drop = FALSE]

  FSC <- GINV(t(lambda) %*% sigma.inv %*% lambda) %*% (t(lambda) %*% sigma.inv)

  alpha <- matrix(getMeans(lVs, parTable = parTableH1), 
                  nrow = nrow(X), ncol = length(lVs), byrow = TRUE)

  Y <- X %*% t(FSC) + alpha

  if (standardized) {
    mu <- \(x) mean(x, na.rm = TRUE)
    s  <- \(x) stats::sd(x, na.rm = TRUE)
    Y  <- apply(Y, MARGIN = 2, FUN = \(y) (y - mu(y)) / s(y))
  }

  Y
}
