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
                              intercepts = TRUE,
                              variances = TRUE,
                              var.interaction = FALSE,
                              ...) {
  method   <- object$method
  parTable <- parameter_estimates(object)

  if (standardized) {
    std.col  <- "Std.all"
    merge.by <- c("lhs", "op", "rhs")

    parTable.std <- standardized_estimates(
      object, 
      intercepts = intercepts, 
      monte.carlo = monte.carlo, 
      mc.reps = mc.reps
    )

    parTable.std[[std.col]] <- parTable.std$est

    parTable <- leftJoin(x  = parTable, 
                         y  = parTable.std[c(merge.by, std.col)],
                         by = merge.by)
  } else std.col <- NULL

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
      out$fitH0 <- fit_modsem_da(est_h0)
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
    std.col       = std.col
  )

  class(out) <- "summary_da"
  out
}


#' @export
print.summary_da <- function(x, digits = 3, ...) {
  width.out <- getWidthPrintedParTable(x$parTable, # We want the width without ci and std.col
                                       scientific = x$format$scientific,
                                       ci = FALSE, # x$format$ci,
                                       digits = x$format$digits,
                                       loadings = x$format$loadings,
                                       regressions = x$format$regressions,
                                       covariances = x$format$covariances,
                                       intercepts = x$format$intercepts,
                                       variances = x$format$variances,
                                       std.col   = NULL) # x$format$std.col)
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
    values <- c(round(x$nullModel$logLik, 2), round(x$fitH0$AIC, 2), round(x$fitH0$BIC, 2))

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
    values <- c(formatNumeric(x$D$diff.loglik, digits = 2),
                formatNumeric(x$D$D, digits = 2),
                x$D$df,
                formatPval(x$D$p, scientific = x$format$scientific))
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                     width.out = width.out), "\n")

  }

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
                std.col     = x$format$std.col)
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
    cols <- colnames(modelH0$data)
    cols.present <- cols %in% colnames(newdata)
    stopif(!all(cols.present), "Missing cols in `newdata`:\n", cols[!cols.present])
  
    newdata <- as.matrix(newdata)[, cols]

  } else newdata <- modelH0$data
  
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
