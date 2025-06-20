printModsemPIHeader <- function(approach) {
  cat(paste0("modsem (version ", PKG_INFO$version, ", approach = ",
             approach, "):\n"))
}


#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param ... arguments passed to lavaan::summary()
#' @rdname summary
#' @export
summary.modsem_pi <- function(object,
                              H0 = TRUE,
                              r.squared = TRUE,
                              adjusted.stat = FALSE,
                              digits = 3,
                              scientific = FALSE,
                              ci = FALSE,
                              ...) {
  out <- list()
  out$lavaan <- lavaan::summary(extract_lavaan(object), ...)
  out$info <- list(version = PKG_INFO$version, approach = attributes(object)$method)
  out$fit <- lavaan::fitMeasures(extract_lavaan(object))
  out$logLik <- lavaan::fitMeasures(extract_lavaan(object), "logl")
  out$N <- lavaan::nobs(extract_lavaan(object))
  out$format <- list(digits = digits, scientific = scientific, adjusted.stat = adjusted.stat, ci = ci)

  # Check for interaction effect in the model
  parTable <- parameter_estimates(object)
  hasInteraction <- parTableHasInteraction(object$input$parTable)

  out$hasInteraction <- hasInteraction

  # If no interaction effect, skip H0, r.squared, and warn if H0 requested
  if (!hasInteraction) {
    warnif(H0, "Comparative fit to H0 will not be calculated.", immediate. = FALSE)
    class(out) <- c("summary_modsem_pi", "list")
    return(out)
  }

  # Baseline (H0) model
  if (H0) tryCatch({
    est_h0 <- estimate_h0(object, ...)
    out$nullModel <- est_h0

    if (!is.null(est_h0)) {
      # Use compare_fit (modsem generic, which internally uses lavTestLRT/anova for modsem_pi)
      lrt <- compare_fit(object, est_h0)
      out$LRT <- lrt
      out$fitH0 <- lavaan::fitMeasures(extract_lavaan(est_h0))
      out$logLikH0 <- lavaan::fitMeasures(extract_lavaan(est_h0), "logl")
    }
  }, error = \(e) warning2("Baseline model could not be estimated: ", e$message, immediate. = FALSE))

  # R-squared for latent endogenous variables only
  if (r.squared) tryCatch({
    # Get all R2
    r2_all <- modsem_inspect(object, "r2")
    # Get latent variables
    lv_names <- lavaan::lavNames(extract_lavaan(object), type = "lv")
    # Get endogenous (on LHS of ~)
    reg_table <- parameter_estimates(object)
    endo_lhs <- unique(reg_table$lhs[reg_table$op == "~"])
    endo_lv <- intersect(lv_names, endo_lhs)
    out$r.squared <- r2_all[endo_lv]

    if (H0 && !is.null(out$nullModel)) {
      r2_all_H0 <- modsem_inspect(out$nullModel, "r2")
      out$r.squared.H0 <- r2_all_H0[endo_lv]
      out$r.squared.diff <- out$r.squared - out$r.squared.H0
    }
  }, error = \(e) warning2("R-squared could not be computed: ", e$message, immediate. = FALSE))

  class(out) <- c("summary_modsem_pi", "list")
  out
}


#' @export
print.summary_modsem_pi <- function(x, ...) {
  digits <- x$format$digits
  scientific <- x$format$scientific
  adjusted.stat <- x$format$adjusted.stat
  ci <- x$format$ci

  # Compute width for right justification based on lavaan output
  # lavaan always seems to use a width of 54 characters for the main 
  # part of the first header, regardless of the widht of the coefficient 
  # tables... So we don't need to compute it dynamically.
  width.out <- 54
  # lavcat <- utils::capture.output(print(x$lavaan))
  # width.out <- max(nchar(lavcat))

  # Helper for left/right align with indent of 2 spaces for names
  align_lavaan <- function(lhs, rhs, width, indent = 2) {
    lhs <- sprintf("%*s%s", indent, "", lhs)
    pad <- width - nchar(lhs) - nchar(rhs)
    pad <- ifelse(pad < 1, 1, pad)
    paste0(lhs, strrep(" ", pad), rhs)
  }

  printModsemPIHeader(x$info$approach)

  # Only print detailed output if there is an interaction effect
  if (!is.null(x$hasInteraction) && !x$hasInteraction) {
    # The warning is already given in summary, so just return after header and lavaan print.
    print(x$lavaan)
    return(invisible(x))
  }

  # Interaction Model Fit Measures (H1)
  cat("\nInteraction Model Fit Measures (H1):\n")
  fit <- x$fit
  namesH1 <- c("Loglikelihood", "Akaike (AIC)", "Bayesian (BIC)", "Chi-square", "Degrees of Freedom", "P-value (Chi-square)", "RMSEA")
  valuesH1 <- c(
    formatC(x$logLik, digits = 2, format = "f"),
    formatC(fit["aic"], digits = 2, format = "f"),
    formatC(fit["bic"], digits = 2, format = "f"),
    formatC(fit["chisq"], digits = 2, format = "f"),
    as.character(fit["df"]),
    formatC(fit["pvalue"], digits = digits, format = if (scientific) "e" else "f"),
    formatC(fit["rmsea"], digits = 3, format = "f")
  )
  for(i in seq_along(namesH1)) {
    cat(align_lavaan(namesH1[i], valuesH1[i], width.out), "\n")
  }
  cat("\n")

  # Baseline Model (H0) 
  if (!is.null(x$fitH0)) {
    cat("Fit Measures for Baseline Model (H0):\n")
    fitH0 <- x$fitH0
    valuesH0 <- c(
      formatC(x$logLikH0, digits = 2, format = "f"),
      formatC(fitH0["aic"], digits = 2, format = "f"),
      formatC(fitH0["bic"], digits = 2, format = "f"),
      formatC(fitH0["chisq"], digits = 2, format = "f"),
      as.character(fitH0["df"]),
      formatC(fitH0["pvalue"], digits = digits, format = if (scientific) "e" else "f"),
      formatC(fitH0["rmsea"], digits = 3, format = "f")
    )
    for(i in seq_along(namesH1)) {
      cat(align_lavaan(namesH1[i], valuesH0[i], width.out), "\n")
    }
    cat("\n")
  }

  # Comparative Fit (LRT)
  if (!is.null(x$LRT)) {
    cat("Comparative Fit to H0 (LRT test):\n")
    lrt <- x$LRT
    if (nrow(lrt) > 1) {
      lrt_row <- lrt[2, ]
      cat(align_lavaan("Chi-square diff", 
                       formatC(lrt_row[["Chisq diff"]], digits = 3, format = "f"), width.out), "\n")
      cat(align_lavaan("Degrees of freedom diff", 
                       as.character(lrt_row[["Df diff"]]), width.out), "\n")
      cat(align_lavaan("P-value (LRT)", 
                       formatC(lrt_row[["Pr(>Chisq)"]], digits = digits, format = if (scientific) "e" else "f"), width.out), "\n")
    }
    cat("\n")
  }

  # R-Squared
  if (!is.null(x$r.squared)) {
    cat("R-Squared Interaction Model (H1):\n")
    for (i in seq_along(x$r.squared)) {
      cat(align_lavaan(names(x$r.squared)[i], formatC(x$r.squared[i], digits = 3, format = "f"), width.out), "\n")
    }
    if (!is.null(x$r.squared.H0)) {
      cat("R-Squared Baseline Model (H0):\n")
      for (i in seq_along(x$r.squared.H0)) {
        cat(align_lavaan(names(x$r.squared.H0)[i], formatC(x$r.squared.H0[i], digits = 3, format = "f"), width.out), "\n")
      }
      cat("R-Squared Change (H1 - H0):\n")
      for (i in seq_along(x$r.squared.diff)) {
        cat(align_lavaan(names(x$r.squared.diff)[i], formatC(x$r.squared.diff[i], digits = 3, format = "f"), width.out), "\n")
      }
    }
    cat("\n")
  }

  # lavaan
  print(x$lavaan)
}


#' @export
parameter_estimates.modsem_pi <- function(object, colon.pi = FALSE, ...) {
  parTable <- lavaan::parameterEstimates(object$lavaan, ...)

  if (colon.pi) {
    elems <- object$elementsInProdNames

    for (xz in names(elems)) {
      xzColon <- paste0(elems[[xz]], collapse = ":")
      parTable[parTable$rhs == xz, "rhs"] <- xzColon
      parTable[parTable$lhs == xz, "lhs"] <- xzColon # shouldn't be necessary, but just in case
                                                     # the user has done something weird...
    }
  }

  parTable
}


#' @describeIn standardized_estimates Method for \code{modsem_pi} objects
#'
#' @param correction Logical. Whether to apply a correction to the standardized estimates
#' of the interaction term. By default, \code{FALSE}, which standardizes the interaction term
#' such that \code{var(xz) = 1}, consistent with \code{lavaan}'s treatment of latent interactions.
#' If \code{TRUE}, the correction computes the interaction variance using the formula
#' \code{var(xz) = var(x) * var(z) + cov(x, z)^2}, under the assumption of normality and
#' zero-centered variables. This may provide more accurate standardization when \code{x}
#' and \code{z} are correlated.
#'
#' @param std.errors Character string indicating the method used to compute standard errors
#' when \code{correction = TRUE}. Options include:
#' \describe{
#'   \item{"rescale"}{Simply rescales the standard errors. Fastest option.}
#'   \item{"delta"}{Uses the delta method to approximate standard errors.}
#'   \item{"monte.carlo"}{Uses Monte Carlo simulation to estimate standard errors.}
#' }
#' Ignored if \code{correction = FALSE}.
#'
#' @param mc.reps Integer. Number of Monte Carlo replications to use if \code{std.errors = "monte.carlo"}.
#'
#' @param colon.pi Logical. If \code{TRUE}, the interaction terms in the output will be
#' will be formatted using \code{:} to seperate the elements in the interaction term. Default
#' is \code{FALSE}, using the default formatting from \code{lavaan}.
#'
#' @param ... Additional arguments passed on to \code{lavaan::standardizedSolution()}.
#'
#' @export
standardized_estimates.modsem_pi <- function(object, 
                                             correction = FALSE, 
                                             std.errors = c("rescale", "delta", "monte.carlo"),
                                             mc.reps = 10000,
                                             colon.pi = FALSE, 
                                             ...) {
  std.errors <- tolower(std.errors)
  std.errors <- match.arg(std.errors)

  uncorrected <- \(object) rename(lavaan::standardizedSolution(object$lavaan, ...), 
                                  est.std = "est")

  if (correction && std.errors == "rescale") {
    parTable.std <- uncorrected(object)
    parTable.std <- correctStdSolutionPI(object, parTable.std = parTable.std)

  } else if (correction) {
    monte.carlo <- std.errors == "monte.carlo"
    solution <- standardizedSolutionCOEFS(object, monte.carlo = monte.carlo, ...) 
    parTable.std <- solution$parTable
  
  } else {
    parTable.std <- uncorrected(object)
  }

  if (!colon.pi) {
    rm <- \(x) stringr::str_remove_all(x, ":")
    parTable.std$rhs <- rm(parTable.std$rhs)
    parTable.std$lhs <- rm(parTable.std$lhs)
  }

  parTable.std
}


#' @export
modsem_inspect.modsem_pi <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "free"
  lavaan::lavInspect(object$lavaan, what = what, ...)
}


#' @export
#' @importFrom stats vcov
vcov.modsem_pi <- function(object, ...) {
  lavaan::vcov(object$lavaan, ...)
}


#' @export
#' @importFrom stats coef
coef.modsem_pi <- function(object, ...) {
  lavaan::coef(object$lavaan, ...)
}


#' @export
#' @importFrom stats coefficients
coefficients.modsem_pi <- function(object, ...) {
  lavaan::coef(object$lavaan, ...)
}


#' @export
#' @importFrom stats nobs
nobs.modsem_pi <- function(object, ...) {
  lavaan::nobs(object$lavaan, ...)
}


#' @export
print.modsem_pi <- function(x, ...) {
  printModsemPIHeader(attributes(x)$method)
  print(x$lavaan)
}
