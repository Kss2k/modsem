#' summary for modsem objects
#'
#' @param object modsem object to summarized
#' @param scientific print p-values in scientific notation
#' @param standardized standardize estimates
#' @param ci print confidence intervals
#' @param digits number of digits to print
#' @param loadings print loadings
#' @param regressions print regressions
#' @param covariances print covariances
#' @param intercepts print intercepts
#' @param variances print variances
#' @param ... arguments passed to other functions
#' @rdname summary
#' @export
summary.modsem_mplus <- function(object,
                                 scientific = FALSE,
                                 standardized = FALSE,
                                 ci = FALSE,
                                 digits = 3,
                                 loadings = TRUE,
                                 regressions = TRUE,
                                 covariances = TRUE,
                                 intercepts = TRUE,
                                 variances = TRUE,
                                 ...) {
  extra.cols <- NULL
  if (standardized) {
    object$parTable <- addStdEstimatesMplus(object, ...)
    extra.cols      <- utils::tail(colnames(object$parTable), 1L)
  }

  object$format <- list(digits = digits,
                        scientific = scientific,
                        ci = ci,
                        extra.cols = extra.cols,
                        loadings = loadings,
                        regressions = regressions,
                        covariances = covariances,
                        intercepts = intercepts,
                        variances = variances)
  structure(object, class = "summary_mplus")
}


addStdEstimatesMplus <- function(object, type = "stdyx", ...) {
  # Mplus' STD, STDY and STDYX correspond to lavaan's std.lv, std.nox and std.all
  column <- switch(tolower(type), std = "Std.lv", stdy = "Std.nox", un = "Est", "Std.all")

  parTable <- object$parTable
  std      <- standardized_estimates(object, type = type, ...)

  by  <- intersect(c("lhs", "op", "rhs", "level"), colnames(parTable))
  std <- std[!duplicated(std[by]), c(by, "est"), drop = FALSE]
  colnames(std)[colnames(std) == "est"] <- column

  out <- leftJoin(left = parTable, right = std, by = by)
  out <- out[c(colnames(parTable), column)]
  rownames(out) <- NULL

  out
}


#' @export
print.summary_mplus <- function(x, ...) {
  summaries <- x$model$results$summaries
  getSummary <- \(name) if (name %in% names(summaries)) summaries[[name]][[1L]] else NA

  width.out <- getWidthPrintedParTable(
    parTable   = x$parTable,
    scientific = x$format$scientific,
    ci         = FALSE,
    digits     = x$format$digits,
    extra.cols = x$format$extra.cols
  )

  formatNum <- \(val) formatNumeric(as.numeric(val), digits = 2)

  cat("modsem: \nMethod =", attributes(x)$method, "\n\n")

  info <- c(
    "Mplus version"                = getSummary("Mplus.version"),
    "Estimator"                    = getSummary("Estimator"),
    "Analysis type"                = getSummary("AnalysisType"),
    "Number of observations"       = NROW(x$data),
    "Number of clusters"           = if (!is.null(x$info$cluster))
                                     length(unique(x$data[[x$info$cluster]])),
    "Number of parameters"         = getSummary("Parameters"),
    "Log-likelihood"               = formatNum(getSummary("LL")),
    "Akaike (AIC)"                 = formatNum(getSummary("AIC")),
    "Bayesian (BIC)"               = formatNum(getSummary("BIC")),
    "Deviance (DIC)"               = formatNum(getSummary("DIC")),
    "Effective parameters (pD)"    = formatNum(getSummary("pD")),
    "Posterior predictive p-value" = formatNum(getSummary("PostPred_PValue"))
  )
  info <- info[!is.na(info) & info != "NA"]

  cat(allignLhsRhs(lhs = names(info), rhs = info, pad = "  ",
                   width.out = width.out), "\n", sep = "")

  if (grepl("bayes", getSummary("Estimator"), ignore.case = TRUE)) {
    cat("Std.Error = posterior standard deviation; P-values are one-tailed (Mplus)\n",
        "Confidence intervals are 95% credible intervals\n\n", sep = "")
  }

  printTable <- function(parTable) {
    printParTable(parTable,
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

  if (!"level" %in% names(x$parTable)) {
    printTable(x$parTable)
    return(invisible(x))
  }

  # Defined parameters (level 0) are not part of any level, and are printed last
  levels <- sort(unique(x$parTable$level[x$parTable$level > 0L]))
  for (level in levels) {
    printf("Level %d [%s]:\n\n", level, c("within", "between")[[level]])

    select <- x$parTable$level == level
    if (level == max(levels)) select <- select | x$parTable$level == 0L

    printTable(x$parTable[select, , drop = FALSE])
  }

  invisible(x)
}


#' @export
print.modsem_mplus <- function(x, ...) {
  cat("modsem: \nMethod =", attributes(x)$method, "\n")
  print(x$parTable)
}


#' @export
#' @describeIn parameter_estimates Get parameter estimates of a \code{\link{modsem_mplus}} object
parameter_estimates.modsem_mplus <- function(object,
                                             colon.pi = NULL, # capture
                                             high.order.as.measr = NULL, # capture
                                             rm.tmp.ov = NULL, # capture argument
                                             label.renamed.prod = NULL, # capture argument
                                             is.public = NULL, # capture argument
                                             ...) {
  object$parTable
}


#' @export
var_interactions.modsem_mplus <- function(object, ...) {
  parTable <- parameter_estimates(object)

  if (!"level" %in% names(parTable))
    return(var_interactions.data.frame(parTable, ...))

  # Multilevel models: Interaction terms are computed separately for each level.
  # Defined parameters (level 0) are not part of any level, and are kept as is.
  out <- lapply(sort(unique(parTable$level)), FUN = \(level) {
    parTable.l <- parTable[parTable$level == level, names(parTable) != "level", drop = FALSE]

    if (level > 0L)
      parTable.l <- var_interactions.data.frame(parTable.l, ...)

    parTable.l$level <- level
    parTable.l[names(parTable)]
  })

  modsemParTable(do.call(rbind, out))
}


#' @describeIn standardized_estimates Retrieve standardized estimates from \code{\link{modsem_mplus}} object.
#' @param type Type of standardized estimates to retrieve. Can be one of: \code{"stdyx", "stdy", "std", "un", "modsem"}.
#'   \code{"stdyx"}, \code{"stdy"} and \code{"std"} are retrieved from \code{Mplus}, if available.
#'   Otherwise (or if \code{type = "modsem"}) the standardized estimates are computed by \code{modsem}.
#'   Standardized estimates computed by \code{modsem} are not available for multilevel models.
#' @export
standardized_estimates.modsem_mplus <- function(object, type = "stdyx", mc.reps = 1e6, ...) {
  type       <- tolower(type)
  multilevel <- isTRUE(object$info$spec$multilevel)

  if (type != "modsem") {
    table    <- if (type == "un") "unstandardized" else paste0(type, ".standardized")
    params   <- object$model$results$parameters[[table]]
    parTable <- mplusParTable(params, spec = object$info$spec)

    if (!is.null(parTable))
      return(modsemParTable(parTable))

    mod_stopif(multilevel,
               sprintf("Standardized estimates (`%s`) are not available from `Mplus` for this model!", type))
  }

  mod_stopif(multilevel,
             "Standardized estimates computed by `modsem` are not available for multilevel models!")

  stdSolution <- standardizedSolutionCOEFS(
    object,
    mc.reps = mc.reps,
    ...
  )

  stdSolution$parTable
}


#' @export
#' @importFrom stats nobs
nobs.modsem_mplus <- function(object, ...) {
  NROW(object$data)
}


#' @export
#' @importFrom stats vcov
vcov.modsem_mplus <- function(object, ...) {
  object$vcov
}


#' @export
#' @importFrom stats coef
coef.modsem_mplus <- function(object, ...) {
  object$coefs
}
