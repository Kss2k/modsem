#' Summarize a parameter table from a \code{modsem} model.
#'
#' @param parTable A parameter table, typically obtained from a \code{\link{modsem}} model
#'  using \code{\link{parameter_estimates}} or \code{\link{standardized_estimates}}.
#'
#' @param scientific Logical, whether to print p-values in scientific notation.
#'
#' @param ci Logical, whether to include confidence intervals in the output.
#'
#' @param digits Integer, number of digits to round the estimates to
#'  (default is 3).
#'
#' @param loadings Logical, whether to include factor loadings in the output.
#'
#' @param regressions Logical, whether to include regression coefficients in the output.
#'
#' @param covariances Logical, whether to include covariance estimates in the output.
#'
#' @param intercepts Logical, whether to include intercepts in the output.
#'
#' @param variances Logical, whether to include variance estimates in the output.
#'
#' @return A summary object containing the parameter table and additional information.
#'
#' @examples
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'
#'   # Inner Model
#'   Y ~ X + Z + X:Z
#' '
#' # Double centering approach
#' est_dca <- modsem(m1, oneInt)
#'
#' std <- standardized_estimates(est_dca, correction = TRUE)
#' summarize_partable(std)
#' @export
summarize_partable <- function(parTable,
                               scientific  = FALSE,
                               ci          = FALSE,
                               digits      = 3,
                               loadings    = TRUE,
                               regressions = TRUE,
                               covariances = TRUE,
                               intercepts  = TRUE,
                               variances   = TRUE) {
  if (!"label" %in% colnames(parTable)) parTable$label <- ""

  parTable <- rename(
    .X      = parTable,
    est.std = "est",
    se      = "std.error",
    pvalue  = "p.value",
    p       = "p.value",
    zvalue  = "z.value",
    z       = "z.value"
  )

  standard.cols <- c("lhs", "op", "rhs", "label", "est", 
                     "std.error", "p.value", "z.value",
                     "ci.lower", "ci.upper")
  extra.cols <- setdiff(colnames(parTable), standard.cols)

  width.out <- getWidthPrintedParTable(
    parTable    = parTable,
    scientific  = scientific,
    ci          = ci,
    digits      = digits,
    loadings    = loadings,
    regressions = regressions,
    covariances = covariances,
    intercepts  = intercepts,
    variances   = variances,
    extra.cols  = NULL
  )

  info.names <- c("Number of model parameters",
                  "Number of latent variables",
                  "Number of observed variables")

  info.values <- c(nrow(parTable),
                   length(getLVs(parTable)),
                   length(getOVs(parTable)))

  out <- list(
    parTable    = parTable,
    scientific  = scientific,
    ci          = ci,
    digits      = digits,
    loadings    = loadings,
    regressions = regressions,
    covariances = covariances,
    intercepts  = intercepts,
    variances   = variances,
    width.out   = width.out,
    info.names  = info.names,
    info.values = info.values,
    extra.cols  = extra.cols
  )

  class(out) <- c("list", "modsem_partable_summary")

  out
}


#' @export
print.modsem_partable_summary <- function(x, ...) {
  printf("modsem (version %s)\n\n", PKG_INFO$version)

  pad <- "  "

  cat(allignLhsRhs(
    lhs = x$info.names,
    rhs = x$info.values,
    pad = pad,
    width.out = x$width.out
  ), "\n")

  printParTable(
    parTable    = x$parTable,
    scientific  = x$scientific,
    ci          = x$ci,
    digits      = x$digits,
    loadings    = x$loadings,
    regressions = x$regressions,
    covariances = x$covariances,
    intercepts  = x$intercepts,
    variances   = x$variances,
    extra.cols  = x$extra.cols
  )
}
