printModsemPIHeader <- function(approach) {
  cat(paste0("modsem (version ", PKG_INFO$version, ", approach = ",
             approach, "):\n"))
}


#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param H0 Should the baseline model be estimated, and used to produce
#'  comparative fit?
#' @param digits Number of digits for printed numerical values
#' @param scientific Should scientific format be used for p-values?
#' @param verbose Should messages be printed?
#' @param ... arguments passed to lavaan::summary()
#' @rdname summary
#' @export
summary.modsem_pi <- function(object,
                              H0 = is_interaction_model(object),
                              r.squared = TRUE,
                              adjusted.stat = FALSE,
                              digits = 3,
                              scientific = FALSE,
                              verbose = TRUE,
                              ...) {
  out         <- list()
  out$lavaan  <- lavaan::summary(extract_lavaan(object), ...)
  out$info    <- list(version = PKG_INFO$version, approach = attributes(object)$method)
  out$fit     <- tryCatch(lavaan::fitMeasures(extract_lavaan(object)), error = \(e) NA)
  out$N       <- lavaan::nobs(extract_lavaan(object))
  out$format  <- list(digits = digits, scientific = scientific, adjusted.stat = adjusted.stat)
  out$ngroups <- lavaan::lavInspect(extract_lavaan(object), what = "ngroups")

  # Check for interaction effect in the model
  parTable <- parameter_estimates(object)
  hasInteraction <- parTableHasInteraction(object$input$parTable)

  out$hasInteraction <- hasInteraction

  etas <- getEtas(parTable, isLV = FALSE)

  # If no interaction effect, skip H0, r.squared, and warn if H0 requested
  if (!hasInteraction) {
    warnif(H0, "Comparative fit to H0 will not be calculated.", immediate. = FALSE)
    class(out) <- c("summary_modsem_pi", "list")
    return(out)
  }

  # Baseline (H0) model
  if (H0) tryCatch({
    if (verbose) cat("Estimating baseline model (H0)\n")
    est_h0 <- estimate_h0(object, ...)
    out$nullModel <- est_h0

    if (!is.null(est_h0)) {
      # Use compare_fit (modsem generic, which internally uses lavTestLRT/anova for modsem_pi)
      lrt <- compare_fit(object, est_h0)
      out$LRT <- lrt
      out$fitH0 <- tryCatch(lavaan::fitMeasures(extract_lavaan(est_h0)), error = \(e) NA)
    }
  }, error = \(e) warning2("Baseline model could not be estimated: ", e$message, immediate. = FALSE))

  # R-squared for latent endogenous variables only
  if (r.squared & length(etas)) tryCatch({
    r2All <- getRsqrPI(object)
    out$r.squared <- r2All[etas]

    if (H0 && !is.null(out$nullModel)) {
      r2AllH0 <- getRsqrPI(out$nullModel)

      out$r.squared.H0   <- r2AllH0[etas]
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

  # Compute width for right justification based on lavaan output
  # lavaan always seems to use a width of 54 characters for the main
  # part of the first header, regardless of the widht of the coefficient
  # tables... So we don't need to compute it dynamically??
  # lavcat <- utils::capture.output(print(x$lavaan))
  # ncheck <- 10 # Check only the header
  # width.out <- max(nchar(lavcat[seq_len(ncheck)]))
  width.out <- 54

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
  namesH1 <- c("Loglikelihood",
               "Akaike (AIC)",
               "Bayesian (BIC)",
               "Chi-square",
               "Degrees of Freedom",
               "P-value (Chi-square)",
               "RMSEA",
               "CFI",
               "SRMR")

  valuesH1 <- c(
    tryFormatC(fit["logl"], digits = 2, format = "f"),
    tryFormatC(fit["aic"], digits = 2, format = "f"),
    tryFormatC(fit["bic"], digits = 2, format = "f"),
    tryFormatC(fit["chisq"], digits = 2, format = "f"),
    ifelse(is.na(fit["df"]), yes = "NA", no = as.character(fit["df"])),
    tryFormatC(fit["pvalue"], digits = digits, format = if (scientific) "e" else "f"),
    tryFormatC(fit["rmsea"], digits = 3, format = "f"),
    tryFormatC(fit["cfi"], digits = 3, format = "f"),
    tryFormatC(fit["srmr"], digits = 3, format = "f")
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
      tryFormatC(fitH0["logl"], digits = 2, format = "f"),
      tryFormatC(fitH0["aic"], digits = 2, format = "f"),
      tryFormatC(fitH0["bic"], digits = 2, format = "f"),
      tryFormatC(fitH0["chisq"], digits = 2, format = "f"),
      ifelse(is.na(fitH0["df"]), yes = "NA", no = as.character(fitH0["df"])),
      tryFormatC(fitH0["pvalue"], digits = digits,
              format = if (scientific) "e" else "f"),
      tryFormatC(fitH0["rmsea"], digits = 3, format = "f"),
      tryFormatC(fitH0["cfi"], digits = 3, format = "f"),
      tryFormatC(fitH0["srmr"], digits = 3, format = "f")
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
                       tryFormatC(lrt_row[["Chisq diff"]], digits = 3, format = "f"),
                       width.out), "\n")
      cat(align_lavaan("Degrees of freedom diff",
                       as.character(lrt_row[["Df diff"]]), width.out), "\n")
      cat(align_lavaan("P-value (LRT)",
                       tryFormatC(lrt_row[["Pr(>Chisq)"]], digits = digits,
                               format = if (scientific) "e" else "f"),
                       width.out), "\n")
    }
    cat("\n")
  }

  # R-Squared
  if (!is.null(x$r.squared)) {
    k <- x$ngroups
    if (k > 1) printH <- \(...) printf("Total R-Squared %s (%s):\n", ...)
    else       printH <- \(...) printf("R-Squared %s (%s):\n", ...)

    printH("Interaction Model", "H1")
    for (i in seq_along(x$r.squared)) {
      cat(align_lavaan(names(x$r.squared)[i],
                       tryFormatC(x$r.squared[i], digits = 3, format = "f"),
                       width.out), "\n")
    }

    if (!is.null(x$r.squared.H0)) {
      printH("Baseline Model", "H0")
      for (i in seq_along(x$r.squared.H0)) {
        cat(align_lavaan(names(x$r.squared.H0)[i],
                         tryFormatC(x$r.squared.H0[i], digits = 3, format = "f"),
                         width.out), "\n")
      }

      printH("Change", "H1 - H0")
      for (i in seq_along(x$r.squared.diff)) {
        cat(align_lavaan(names(x$r.squared.diff)[i],
                         tryFormatC(x$r.squared.diff[i], digits = 3, format = "f"),
                         width.out), "\n")
      }
    }
    cat("\n")
  }

  # lavaan
  print(x$lavaan)
}



#' @describeIn parameter_estimates Get parameter estimates of a \code{\link{modsem_pi}} object
#' @param colon.pi Should colons (\code{:}) be added to the interaction terms (E.g., `XZ` -> `X:Z`)?
#' @export
parameter_estimates.modsem_pi <- function(object, colon.pi = FALSE,
                                          label.renamed.prod = FALSE,
                                          high.order.as.measr = NULL, # capture argument
                                          rm.tmp.ov = NULL, # capture argument
                                          is.public = NULL, # capture argument
                                          ...) {
  parTable <- lavaan::parameterEstimates(object$lavaan, ...)

  if (colon.pi) {
    elems <- object$elementsInProdNames

    if (length(elems) && !"label" %in% colnames(parTable))
      parTable$label <- ""

    if (label.renamed.prod)
      origLabels <- getParTableLabels(parTable, labelCol = "label")
    else
      origLabels <- parTable$label

    for (xz in names(elems)) {
      xzColon <- paste0(elems[[xz]], collapse = ":")
      rmatch <- parTable$rhs == xz
      lmatch <- parTable$lhs == xz

      parTable[rmatch | lmatch, "label"] <- origLabels[rmatch | lmatch]

      parTable[rmatch, "rhs"] <- xzColon
      parTable[lmatch, "lhs"] <- xzColon # shouldn't be necessary, but just in case
                                         # the user has done something weird...
    }
  }

  parTable
}


#' @describeIn standardized_estimates Method for \code{modsem_pi} objects
#'
#' @param correction Logical. Whether to apply a correction to the standardized estimates
#' of the interaction term. By default, \code{FALSE}, which standardizes the interaction term
#' such that \eqn{\sigma^2(XZ) = 1}, consistent with \code{lavaan}'s treatment of latent interactions.
#' This is usually wrong, as it does not account for the fact that the interaction term
#' is a product of two variables, which means that the variance of the interaction term
#' of standardized variables (usually) is not equal to 1.
#'
#' Hence the scale of the interaction effect changes, such that
#' the standardized interaction term does not correspond to one (standardized) unit
#' change in the moderating variables. If \code{TRUE}, a correction is applied by
#' computing the interaction term \eqn{b_3 = \frac{B_3 \sigma(X) \sigma(Z)}{\sigma(Y)}} (where
#' \eqn{B_3} is the unstandardized coefficient for the interaction term), and solving for \eqn{\sigma(XZ)}, which
#' is used to correct the variance of the interaction term, and its covariances.
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
#' is \code{FALSE}, using the default formatting from \code{lavaan}. Only relevant if
#' \code{std.errors != "rescale"} and \code{correction = TRUE}.
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
  std.errors  <- tolower(std.errors)
  std.errors  <- match.arg(std.errors)
  monte.carlo <- std.errors == "monte.carlo"

  uncorrected <- function(object) {
    rename(lavaan::standardizedSolution(object$lavaan, ...), est.std = "est")
  }

  cleanOutputParTable <- function(parTable) {
    hasColon <- any(grepl(":", parTable$rhs))

    if (colon.pi && !hasColon) {
      elems <- object$elementsInProdNames

      if (length(elems) && !"label" %in% colnames(parTable))
        parTable$label <- ""

      origLabels <- getParTableLabels(parTable, labelCol = "label")

      for (xz in names(elems)) {
        xzColon <- paste0(elems[[xz]], collapse = ":")
        rmatch <- parTable$rhs == xz
        lmatch <- parTable$lhs == xz

        parTable[rmatch | lmatch, "label"] <- origLabels[rmatch | lmatch]

        parTable[rmatch, "rhs"] <- xzColon
        parTable[lmatch, "lhs"] <- xzColon # shouldn't be necessary, but just in case
                                           # the user has done something weird...
      }

    } else if (!colon.pi) {
      rm <- \(x) stringr::str_remove_all(x, ":")
      parTable$rhs <- rm(parTable$rhs)
      parTable$lhs <- rm(parTable$lhs)

    }

    parTable
  }

  if (correction && std.errors == "rescale") {
    parTable.std.naive <- uncorrected(object)

    correction <- function(object) {
      correctStdSolutionPI(
        object       = object,
        parTable.std = parTable.std.naive
      )
    }

  } else if (correction) {
    correction <- function(object) {

      solution <- standardizedSolutionCOEFS(
        object      = object,
        monte.carlo = monte.carlo,
        mc.reps     = mc.reps,
        center      = FALSE,
        ...
      )

      solution$parTable
    }

  } else return(cleanOutputParTable(uncorrected(object)))

  parTable.ustd <- parameter_estimates(object)

  hiorder <- isHigherOrderParTable(parTable.ustd)
  cluster <- isClustered(object)
  rescale <- std.errors == "rescale"

  stopif(hiorder && !rescale, 'Correction of higher-order models are ',
         'not supported at all with `std.errors = "rescale"`!')
  stopif(cluster, "Correction of clustered (multilevel) models is not supported!")
  warnif(hiorder && rescale, "Correction of higher-order models will likely not work!",
         immediate. = FALSE)

  parTable.std <- correction(object)

  cleanOutputParTable(parTable.std)
}


#' @export
#' @describeIn modsem_inspect Inspect a \code{\link{modsem_pi}} object
modsem_inspect.modsem_pi <- function(object, what = "free", ...) {
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


#' @describeIn modsem_predict Wrapper for \code{lavaan::predict}
#' @export
modsem_predict.modsem_pi <- function(object, ...) {
  lavaan::predict(extract_lavaan(object), ...)
}


#' @export
is_interaction_model.modsem_pi <- function(object) {
  isTRUE(object$has.interaction)
}
