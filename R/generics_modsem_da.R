#' @export
#' @describeIn parameter_estimates Get parameter estimates of a \code{\link{modsem_da}} object
parameter_estimates.modsem_da <- function(object, ...) {
  object$parTable
}


#' summary for modsem objects
#'
#' @param object modsem object to summarized.
#' @param H0 Should a null model be estimated (used for comparison).
#' @param verbose Print progress for the estimation of null model.
#' @param r.squared Calculate R-squared.
#' @param fit Print additional fit measures.
#' @param adjusted.stat Should sample size corrected/adjustes AIC and BIC be reported?
#' @param digits Number of digits to print.
#' @param scientific Print p-values in scientific notation.
#' @param ci Print confidence intervals.
#' @param standardized Print standardized estimates.
#' @param centered Print mean centered estimates.
#' @param monte.carlo Should Monte Carlo bootstrapped standard errors be used? Only
#'   relevant if \code{standardized = TRUE}. If \code{FALSE} delta method is used instead.
#' @param mc.reps Number of Monte Carlo repetitions. Only relevant if \code{monte.carlo = TRUE},
#'   and \code{standardized = TRUE}.
#' @param loadings Print loadings.
#' @param regressions Print regressions.
#' @param covariances Print covariances.
#' @param intercepts Should intercepts be included in the output?
#' If \code{standardized = TRUE} intercepts will by default be excluded.
#' @param variances Print variances.
#' @param var.interaction If FALSE variances for interaction terms will be removed
#' from the output.
#' @param ... Additional arguments.
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
                              fit = FALSE,
                              adjusted.stat = FALSE,
                              digits = 3,
                              scientific = FALSE,
                              ci = FALSE,
                              standardized = FALSE,
                              centered = FALSE,
                              monte.carlo = FALSE,
                              mc.reps = 10000,
                              loadings = TRUE,
                              regressions = TRUE,
                              covariances = TRUE,
                              intercepts = TRUE,
                              variances = TRUE,
                              var.interaction = FALSE,
                              ...) {
  method   <- object$method
  parTable <- parameter_estimates(object)

  extra.cols <- NULL
  if (standardized) {
    std.col     <- "Std.all"
    extra.cols <- c(extra.cols, std.col)

    parTable <- addTransformedEstimatesPT(
      parTable    = parTable,
      values.from = "est",
      values.to   = std.col,
      FUN         = standardized_estimates,
      object      = object,
      monte.carlo = monte.carlo,
      mc.reps     = mc.reps
    )
  }

  if (centered) {
    cnt.col     <- "Cnt.all"
    extra.cols <- c(extra.cols, cnt.col)

    parTable <- addTransformedEstimatesPT(
      parTable    = parTable,
      values.from = "est",
      values.to   = cnt.col,
      FUN         = centered_estimates,
      object      = object,
      monte.carlo = monte.carlo,
      mc.reps     = mc.reps
    )
  }

  if (!var.interaction) {
    parTable.out <- removeInteractionVariances(parTable)
  } else parTable.out <- parTable

  args <- object$args
  out <- list(
    parTable        = parTable.out,
    data            = object$data$data.full,
    iterations      = object$iterations,
    logLik          = object$logLik,
    fit             = fit_modsem_da(object, chisq = FALSE),
    D               = NULL,
    N               = NROW(object$data$data.full),
    method          = method,
    optimizer       = object$optimizer,
    quad            = object$info.quad,
    type.se         = object$type.se,
    type.estimates  = ifelse(standardized, "standardized", object$type.estimates),
    information     = object$information,
    n.fiml.patterns = length(object$data$ids),
    is.fiml         = object$data$is.fiml,
    npar            = length(coef(object, type = "free")),
    convergence.msg = object$convergence.msg
  )

  if (H0) {
    if (any(grepl(":", parTable$rhs)) && verbose)
      cat("Estimating baseline model (H0)\n")

    est_h0 <- estimate_h0(object, calc.se = FALSE, warn_no_interaction = FALSE)

    out$nullModel <- est_h0
    if (is.null(est_h0)) {
      warning2("Comparative fit to H0 will not be calculated.", immediate. = FALSE)
      H0        <- FALSE
      out$D     <- NULL
      out$fitH0 <- NULL

    } else {
      out$D     <- compare_fit(est_h1 = object, est_h0 = est_h0)
      out$fitH0 <- fit_modsem_da(est_h0, lav.fit = TRUE)
    }
  } else {
    out$D <- NULL
  }

  if (r.squared) {
  	out$r.squared <- modsem_inspect(object, "r2.lv")

    if (H0) {
			r.squared.h0 <- modsem_inspect(est_h0, "r2.lv")
			out$r.squared.h0 <- r.squared.h0[names(out$r.squared)] # should't be necessary
			                                                       # but sort just in case...
		} else out$r.squared.h0 <- NULL

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
    variances     = variances,
    extra.cols    = extra.cols,
    extra.fit     = fit,
    scaled.stat   = object$args$robust.se
  )

  class(out) <- "summary_da"
  out
}


#' @export
print.summary_da <- function(x, digits = 3, ...) {
  # We want the width without ci and extra cols
  width.out <- getWidthPrintedParTable(
    parTable    = x$parTable,
    scientific  = x$format$scientific,
    ci          = FALSE,
    digits      = x$format$digits,
    loadings    = x$format$loadings,
    regressions = x$format$regressions,
    covariances = x$format$covariances,
    intercepts  = x$format$intercepts,
    variances   = x$format$variances,
    extra.cols  = NULL
  )

  printf(x$convergence.msg)

  # Convergence and Model Info -------------------------------------------------
  names <- c(
    "Estimator",
    "Optimization method",
    "Number of model parameters",
    "", # blank line
    "Number of observations",
    "Number of missing patterns"
  )

  values <- c(
    stringr::str_to_upper(c(x$method, x$optimizer)),
    x$npar,
    "", # blank line
    x$N,
    x$n.fiml.patterns
  )

  if (!x$is.fiml) {
    fieldFIML <- grepl("Number of missing patterns", names)
    names  <- names[!fieldFIML]
    values <- values[!fieldFIML]
  }

  cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                   width.out = width.out), "\n")

  # Criterion/LogLik -----------------------------------------------------------
  names <- c(
    "Loglikelihood",
    "Akaike (AIC)",
    "Bayesian (BIC)"
  )

  values <- c(
    formatNumeric(x$logLik,  digits = 2),
    formatNumeric(x$fit$AIC, digits = 2),
    formatNumeric(x$fit$BIC, digits = 2)
  )

  if (x$format$adjusted.stat) {
    names  <- c(names, "Corrected Akaike (AICc)", "Adjusted Bayesian (aBIC)")
    values <- c(values, formatNumeric(x$fit$AICc, digits = 2),
                formatNumeric(x$fit$aBIC, digits = 2))
  }

  cat("Loglikelihood and Information Criteria:\n")
  cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                   width.out = width.out), "\n")

  # Intergration ---------------------------------------------------------------
  if (!is.null(x$quad)) {
    cat("Numerical Integration:\n")
    names <- c("Points of integration (per dim)", "Dimensions",
               "Total points of integration")
    values <- c(x$quad$nodes.dim, x$quad$dim, x$quad$nodes.total)
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                     width.out = width.out), "\n")

  }

  # Comparative fit ------------------------------------------------------------
  if (!is.null(x$D)) {
    lav.fit.h0 <- x$fitH0$lav.fit
    fnull      <- \(x) if (is.null(x)) NA else x

    cat("Fit Measures for Baseline Model (H0):\n")

    names <- c("", "Chi-square", "Degrees of Freedom (Chi-square)",
               "P-value (Chi-square)")
    values <- c("Standard",
                formatNumeric(x$fitH0$chisq.value, digits = 2),
                x$fitH0$chisq.df,
                formatPval(x$fitH0$chisq.pvalue, scientific = x$format$scientific))

    if (x$format$scaled.stat) {
      chisq.s <- fnull(lav.fit.h0[["chisq.scaled"]])
      df.s    <- fnull(lav.fit.h0[["df.scaled"]])
      pval.s  <- fnull(lav.fit.h0[["pvalue.scaled"]])
      scale.f <- fnull(lav.fit.h0[["chisq.scaling.factor"]])
      values.scaled <- c("Scaled",
                         formatNumeric(chisq.s, digits = 2),
                         round(df.s),
                         formatPval(pval.s, scientific = x$format$scientific),
                         formatNumeric(scale.f, digits = 3),
                         rep("", 2))
      values <- c(values, rep("", 3))
      names <- c(names, "Scaling correction factor",
                 "  Yuan-Bentler correction (Mplus variant)", "")

    } else values.scaled <- NULL


    names <- c(names, "RMSEA")
    values <- c(values, formatNumeric(x$fitH0$RMSEA, digits = 3))

    if (!is.null(values.scaled)) {
      rmsea.s <- fnull(lav.fit.h0[["rmsea.scaled"]])
      values.scaled <- c(values.scaled, formatNumeric(rmsea.s, digits = 3))
    }

    if (x$format$extra.fit) {
      names  <- c(names, "CFI", "TLI", "SRMR")
      srmr <- x$fitH0$SRMR
      cfi  <- fnull(lav.fit.h0[["cfi"]])
      tli  <- fnull(lav.fit.h0[["tli"]])
      values <- c(values,
                  formatNumeric(cfi, digits = 3),
                  formatNumeric(tli, digits = 3),
                  formatNumeric(srmr, digits = 3))

      if (!is.null(values.scaled)) {
        cfi.s <- fnull(lav.fit.h0[["cfi.scaled"]])
        tli.s <- fnull(lav.fit.h0[["tli.scaled"]])
        values.scaled <- c(values.scaled,
                           formatNumeric(cfi.s, digits = 3),
                           formatNumeric(tli.s, digits = 3),
                           "")
      }
    }

    names <- c(names, "", "Loglikelihood", "Akaike (AIC)", "Bayesian (BIC)")
    values <- c(values, "",
                formatNumeric(x$nullModel$logLik, digits = 2),
                formatNumeric(x$fitH0$AIC, digits = 2),
                formatNumeric(x$fitH0$BIC, digits = 2))

    if (!is.null(values.scaled))
      values.scaled <- c(values.scaled, rep("", 4))

    if (x$format$adjusted.stat) {
      names <- c(names, "Corrected Akaike (AICc)", "Adjusted Bayesian (aBIC)")
      values <- c(values,
                  formatNumeric(x$fitH0$AICc, digits = 2),
                  formatNumeric(x$fitH0$aBIC, digits = 2))

      if (!is.null(values.scaled))
        values.scaled <- c(values.scaled, rep("", 2))
    }

    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                     width.out = width.out, rhs.scaled = values.scaled), "\n")

    cat("Comparative Fit to H0 (LRT test):\n")
    names <- c("Loglikelihood change",
               "Difference test (D)",
               "Degrees of freedom (D)", "P-value (D)")
    values <- c(formatNumeric(x$D$diff.loglik, digits = 2),
                formatNumeric(x$D$D, digits = 2),
                x$D$df,
                formatPval(x$D$p, scientific = x$format$scientific))
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                     width.out = width.out), "\n")

  }

  # R2 -------------------------------------------------------------------------
  if (!is.null(x$r.squared)) {
    r.squared <- formatNumeric(x$r.squared, digits = 3)
    names     <- names(r.squared)

    cat("R-Squared Interaction Model (H1):\n")
    cat(allignLhsRhs(lhs = names, rhs = r.squared,
										 pad = "  ", width.out = width.out))

    if (!is.null(x$r.squared.h0)) {
      r.squared.h0 <- formatNumeric(x$r.squared.h0, digits = 3)
			names.h0     <- names(r.squared.h0)

			cat("R-Squared Baseline Model (H0):\n")
      cat(allignLhsRhs(lhs = names.h0, rhs = r.squared.h0, pad = "  ",
                       width.out = width.out))

      # Calculate Change (using unformatted Rsquared)
     	r.squared.diff <- formatNumeric(x$r.squared - x$r.squared.h0, digits = 3)
			names.diff     <- names(r.squared.diff)
      cat("R-Squared Change (H1 - H0):\n")
      cat(allignLhsRhs(lhs = names.diff, rhs = r.squared.diff,
											 pad = "  ", width.out = width.out))
    }
  }

  # Parameters -----------------------------------------------------------------
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
                variances   = x$format$variances,
                extra.cols  = x$format$extra.cols)
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

  cat(x$convergence.msg)
  print(est)
}


#' @export
var_interactions.modsem_da <- function(object, ...) {
  var_interactions.data.frame(parameter_estimates(object), ...)
}


#' Inspect components of a \code{modsem_da} fit
#'
#' \code{modsem_inspect.modsem_da} Lets you
#' pull matrices, optimiser diagnostics, expected moments, or fit
#' measures from a \code{\link{modsem_da}} object.
#'
#' @param object A fitted object of class \code{"modsem_da"}.
#' @param what   Character scalar selecting what to return (see \emph{Details}).
#'               If \code{NULL} the value \code{"default"} is used.
#' @param ...    Passed straight to \code{modsem_inspect_da()}.
#'
#' @details
#' Below is a list of possible values for the \code{what} argument,
#' organised in several sections.  Keywords are \emph{case-sensitive}.
#'
#' \strong{Presets}
#'
#' \describe{
#'   \item{\code{"default"}}{Everything in \emph{Sample information}, \emph{Optimiser diagnostics}
#'   \emph{Parameter tables}, \emph{Model matrices}, and \emph{Expected-moment matrices} except
#'   the raw \code{data} slot}
#'   \item{\code{"coef"}}{Coefficients and variance-covariance matrix of both free and constrained parameters (same as \code{"coef.all"}).}
#'   \item{\code{"coef.all"}}{Coefficients and variance-covariance matrix of both free and constrained parameters (same as \code{"coef"}).}
#'   \item{\code{"coef.free"}}{Coefficients and variance-covariance matrix of the free parameters.}
#'   \item{\code{"all"}}{All items listed below, including \code{data}.}
#'   \item{\code{"matrices"}}{The model matrices.}
#'   \item{\code{"optim"}}{Only the items under \emph{Optimiser diagnostics}}.
#'   \item{\code{"fit"}}{A list with \code{fit.h0}, \code{fit.h1}, comparative.fit}
#' }
#'
#' \strong{Sample information:}
#'
#' \describe{
#' \item{\code{"N"}}{Number of analysed rows (integer).}
#' }
#'
#' \strong{Parameter estimates and standard errors:}
#'
#' \describe{
#'   \item{\code{"coefficients.free"}}{Free parameter values.}
#'   \item{\code{"coefficients.all"}}{Both free and constrained parameter values.}
#'   \item{\code{"vcov.free"}}{Variance–covariance of free coefficients only.}
#'   \item{\code{"vcov.all"}}{Variance–covariance of both free and constrained coefficients.}
#' }
#'
#' \strong{Optimiser diagnostics:}
#'
#' \describe{
#'   \item{\code{"coefficients.free"}}{Free parameter values.}
#'   \item{\code{"vcov.free"}}{Variance–covariance of free coefficients only.}
#'   \item{\code{"information"}}{Fisher information matrix.}
#'   \item{\code{"loglik"}}{Log-likelihood.}
#'   \item{\code{"iterations"}}{Optimiser iteration count.}
#'   \item{\code{"convergence"}}{\code{TRUE}/\code{FALSE} indicating whether the model converged.}
#' }
#'
#' \strong{Parameter tables:}
#'
#' \describe{
#'   \item{\code{"partable"}}{Parameter table with estimated parameters.}
#'   \item{\code{"partable.input"}}{Parsed model syntax.}
#' }
#'
#' \strong{Model matrices:}
#'
#' \describe{
#'   \item{\code{"lambda"}}{\eqn{\Lambda} – Factor loadings.}
#'   \item{\code{"tau"}}{\eqn{\tau} – Intercepts for indicators.}
#'   \item{\code{"theta"}}{\eqn{\Theta} – Residual (Co-)Variances for indicators.}
#'   \item{\code{"gamma.xi"}}{\eqn{\Gamma_{\xi}} – Structural coefficients between exogenous and endogenous variables.}
#'   \item{\code{"gamma.eta"}}{\eqn{\Gamma_{\eta}} – Structural coefficients between endogenous variables.}
#'   \item{\code{"omega.xi.xi"}}{\eqn{\Omega_{\xi\xi}} – Interaction effects between exogenous variables}
#'   \item{\code{"omega.eta.xi"}}{\eqn{\Omega_{\eta\xi}} – Interaction effects between exogenous and endogenous variables}
#'   \item{\code{"phi"}}{\eqn{\Phi} – (Co-)Variances among exogenous variables.}
#'   \item{\code{"psi"}}{\eqn{\Psi} – Residual (co-)variances among engoenous variables.}
#'   \item{\code{"alpha"}}{\eqn{\alpha} – Intercepts for endogenous variables}
#'   \item{\code{"beta0"}}{\eqn{\beta_0} – Intercepts for exogenous variables}
#' }
#' \strong{Model-implied matrices:}
#'
#' \describe{
#'   \item{\code{"cov.ov"}}{Model-implied covariance of observed variables.}
#'   \item{\code{"cov.lv"}}{Model-implied covariance of latent variables.}
#'   \item{\code{"cov.all"}}{Joint covariance of observed + latent variables.}
#'   \item{\code{"cor.ov"}}{Correlation counterpart of \code{"cov.ov"}.}
#'   \item{\code{"cor.lv"}}{Correlation counterpart of \code{"cov.lv"}.}
#'   \item{\code{"cor.all"}}{Correlation counterpart of \code{"cov.all"}.}
#'   \item{\code{"mean.ov"}}{Expected means of observed variables.}
#'   \item{\code{"mean.lv"}}{Expected means of latent variables.}
#'   \item{\code{"mean.all"}}{Joint mean vector.}
#' }
#'
#' \strong{R-squared and standardized residual variances:}
#'
#' \describe{
#'   \item{\code{"r2.all"}}{R-squared values for both observed (i.e., indicators) and latent endogenous variables.}
#'   \item{\code{"r2.lv"}}{R-squared values for latent endogenous variables.}
#'   \item{\code{"r2.ov"}}{R-squared values for observed (i.e., indicators) variables.}
#'   \item{\code{"res.all"}}{Standardized residuals (i.e., \code{1 - R^2}) for both observed (i.e., indicators) and latent endogenous variables.}
#'   \item{\code{"res.lv"}}{Standardized residuals (i.e., \code{1 - R^2}) for latent endogenous variables.}
#'   \item{\code{"res.ov"}}{Standardized residuals (i.e., \code{1 - R^2}) for observed variables (i.e., indicators).}
#' }
#'
#' \strong{Interaction-specific caveats:}
#'
#' \itemize{
#'   \item If the model contains an \emph{uncentred} latent interaction term it is centred
#'     internally before any \code{cov.*}, \code{cor.*}, or \code{mean.*} matrices are
#'     calculated.
#'   \item These matrices should not be used to compute fit-statistics (e.g.,
#'     chi-square and RMSEA) if there is an interaction term in the model.
#' }
#' @return A named list with the extracted information. If a single piece of information is returned,
#'  it is returned as is; not as a named element in a list.
#'
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
#' est <- modsem(m1, oneInt, "lms")
#'
#' modsem_inspect(est) # everything except "data"
#' modsem_inspect(est, what = "optim")
#' modsem_inspect(est, what = "phi")
#' }
#'
#' @export
#' @describeIn modsem_inspect Inspect a \code{\link{modsem_da}} object
modsem_inspect.modsem_da <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "default"
  modsem_inspect_da(object, what = what, ...)
}


#' @export
#' @importFrom stats vcov
vcov.modsem_da <- function(object, type = c("all", "free"), ...) {
  type <- tolower(type)
  type <- match.arg(type)

  what <- ifelse(type == "all", yes = "vcov.all", no = "vcov.free")
  modsem_inspect_da(object, what = what)
}


#' @export
#' @importFrom stats coefficients
coefficients.modsem_da <- function(object, type = c("all", "free"), ...) {
  type <- tolower(type)
  type <- match.arg(type)

  what <- ifelse(type == "all", yes = "coefficients.all",
                 no = "coefficients.free")
  modsem_inspect_da(object, what = what)
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


#' @describeIn standardized_estimates Method for \code{modsem_da} objects
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
                                             tolerance.zero = 1e-10,
                                             ...) {
  stdSolution <- standardizedSolutionCOEFS(
    object,
    monte.carlo = monte.carlo,
    mc.reps = mc.reps,
    tolerance.zero = tolerance.zero,
    ...
  )

  stdSolution$parTable
}


#' @describeIn centered_estimates Method for \code{modsem_da} objects
#'
#' @param monte.carlo Logical. If \code{TRUE}, use Monte Carlo simulation to estimate
#' standard errors; if \code{FALSE}, use the delta method (default).
#' @param mc.reps Number of Monte Carlo repetitions. Default is 10000.
#' @param tolerance.zero Threshold below which standard errors are set to \code{NA}.
#'
#' @export
centered_estimates.modsem_da <- function(object,
                                         monte.carlo = FALSE,
                                         mc.reps = 10000,
                                         tolerance.zero = 1e-10, ...) {

  stdSolution <- centeredSolutionCOEFS(
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
#' @param newdata Compute factor scores based on a different dataset, than the one used in the model estimation.
#' @param center.data Should data be centered before computing factor scores? Default is \code{TRUE}.
#' @export
modsem_predict.modsem_da <- function(object, standardized = FALSE, H0 = TRUE, newdata = NULL,
                                     center.data = TRUE, ...) {
  modelH1 <- object

  if (H0) {
    modelH0 <- estimate_h0(modelH1, calc.se = FALSE, warn_no_interaction = FALSE,
                           verbose = FALSE)

    if (is.null(modelH0)) modelH0 <- modelH1
  } else modelH0 <- modelH1

  if (!is.null(newdata)) {
    cols <- colnames(modelH0$data$data.full)
    cols.present <- cols %in% colnames(newdata)
    stopif(!all(cols.present), "Missing cols in `newdata`:\n", cols[!cols.present])

    newdata <- as.matrix(newdata)[, cols]

  } else newdata <- modelH0$data$data.full

  transform.x <- if (center.data) \(x) x - mean(x, na.rm = TRUE) else \(x) x

  parTableH1 <- parameter_estimates(modelH1)
  parTableH0 <- parameter_estimates(modelH0)

  lVs   <- getLVs(parTableH1)
  sigma <- modsem_inspect(modelH0, what = "cov.ov")

  sigma.inv <- GINV(sigma)
  lambda    <- getLambdaParTable(parTableH0, rows = colnames(sigma), cols = lVs)

  X <- apply(as.matrix(newdata), MARGIN = 2, FUN = transform.x)
  X <- X[, colnames(sigma), drop = FALSE]

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
