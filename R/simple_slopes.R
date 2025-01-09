#' Get the simple slopes of a SEM model
#'
#' This function calculates simple slopes (predicted values of the outcome variable)
#' at user-specified values of the focal predictor (\code{x}) and moderator (\code{z})
#' in a structural equation modeling (SEM) framework. It supports interaction terms
#' (\code{xz}), computes standard errors (SE), and optionally returns confidence or
#' prediction intervals for these predicted values. It also provides p-values for
#' hypothesis testing. This is useful for visualizing and interpreting moderation
#' effects or to see how the slope changes at different values of the moderator.
#'
#' @param x The name of the variable on the x-axis (focal predictor).
#' @param z The name of the moderator variable.
#' @param y The name of the outcome variable.
#' @param xz The name of the interaction term (\code{x:z}). If \code{NULL}, it will
#'   be created by combining \code{x} and \code{z} with a colon (e.g., \code{"x:z"}).
#'   Some backends may remove or alter the colon symbol, so the function tries to
#'   account for that internally.
#' @param model An object of class \code{\link{modsem_pi}}, \code{\link{modsem_da}},
#'   \code{\link{modsem_mplus}}, or a \code{lavaan} object. This should be a fitted SEM
#'   model that includes paths for \code{y ~ x + z + x:z}.
#' @param vals_x Numeric vector of values of \code{x} at which to compute predicted
#'   slopes. Defaults to \code{-3:3}. If \code{rescale = TRUE}, these values are taken
#'   relative to the mean and standard deviation of \code{x}. A higher density of points
#'   (e.g., \code{seq(-3, 3, 0.1)}) will produce smoother curves and confidence bands.
#' @param vals_z Numeric vector of values of the moderator \code{z} at which to compute
#'   predicted slopes. Defaults to \code{-1:1}. If \code{rescale = TRUE}, these values
#'   are taken relative to the mean and standard deviation of \code{z}. Each unique value
#'   of \code{z} generates a separate regression line \code{y ~ x | z}.
#' @param rescale Logical. If \code{TRUE} (default), \code{x} and \code{z} are standardized
#'   according to their means and standard deviations in the model. The values in
#'   \code{vals_x} and \code{vals_z} are interpreted in those standardized units. The
#'   raw (unscaled) values corresponding to these standardized points will be displayed
#'   in the returned table.
#' @param ci_width A numeric value between 0 and 1 indicating the confidence (or
#'   prediction) interval width. The default is 0.95 (i.e., 95\% interval).
#' @param ci_type A string indicating whether to compute a \code{"confidence"} interval
#'   for the predicted mean (\emph{i.e.}, uncertainty in the regression line) or a
#'   \code{"prediction"} interval for individual outcomes. The default is
#'   \code{"confidence"}. If \code{"prediction"}, the residual variance is added to the
#'   variance of the fitted mean, resulting in a wider interval.
#' @param relative_h0 Logical. If \code{TRUE} (default), hypothesis tests for the
#'   predicted values (\code{predicted - h0}) assume \code{h0} is the model-estimated
#'   mean of \code{y}. If \code{FALSE}, the null hypothesis is \code{h0 = 0}.
#' @param ... Additional arguments passed to lower-level functions or other internal
#'   helpers.
#'
#' @details
#' \strong{Computation Steps}  
#' 1. The function extracts parameter estimates (and, if necessary, their covariance
#'    matrix) from the fitted SEM model (\code{model}).  
#' 2. It identifies the coefficients for \code{x}, \code{z}, and \code{x:z} in the model's
#'    parameter table, as well as the variance of \code{x}, \code{z}, and the residual.  
#' 3. If \code{xz} is not provided, it will be constructed by combining \code{x} and
#'    \code{z} with a colon (\code{":"}). In certain SEM software, the colon may be
#'    removed or replaced internally; the function attempts to reconcile that.  
#' 4. A grid of \code{x} and \code{z} values is created from \code{vals_x} and
#'    \code{vals_z}. If \code{rescale = TRUE}, these values are transformed back into raw
#'    metric units for display in the output.  
#' 5. For each point in the grid, a predicted value of \code{y} is computed via
#'    \code{(beta0 + beta_x * x + beta_z * z + beta_xz * x * z)} and, if included, a
#'    mean offset.  
#' 6. The standard error (\code{std.error}) is derived from the covariance matrix of
#'    the relevant parameters, and if \code{ci_type = "prediction"}, adds the residual
#'    variance.  
#' 7. Confidence (or prediction) intervals are formed using \code{ci_width} (defaulting
#'    to 95\%). The result is a table-like data frame with predicted values, CIs,
#'    standard errors, z-values, and p-values.  
#'
#' @return A \code{data.frame} (invisibly inheriting class \code{"simple_slopes"})
#' with columns:
#' \itemize{
#'   \item \code{vals_x}, \code{vals_z}: The requested grid values of \code{x} and \code{z}.
#'   \item \code{predicted}: The predicted value of \code{y} at that combination of
#'         \code{x} and \code{z}.
#'   \item \code{std.error}: The standard error of the predicted value.
#'   \item \code{z.value}, \code{p.value}: The z-statistic and corresponding p-value
#'         for testing the null hypothesis that \code{predicted == h0}.
#'   \item \code{ci.lower}, \code{ci.upper}: Lower and upper bounds of the confidence
#'         (or prediction) interval.
#' }
#' An attribute \code{"variable_names"} (list of \code{x}, \code{z}, \code{y})
#' is attached for convenience. 
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' library(modsem)
#'
#' m1 <- "
#' # Outer Model
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'
#' # Inner model
#'   Y ~ X + Z + X:Z
#' "
#' est1 <- modsem(m1, data = oneInt)
#'
#' # Simple slopes at X in [-3, 3] and Z in [-1, 1], rescaled to the raw metric.
#' simple_slopes(x = "X", z = "Z", y = "Y", model = est1)
#'
#' # If the data or user wants unscaled values, set rescale = FALSE, etc.
#' simple_slopes(x = "X", z = "Z", y = "Y", model = est1, rescale = FALSE)
#' }
simple_slopes <- function(x, 
                          z, 
                          y,
                          xz = NULL,
                          model,
                          vals_x = -3:3,
                          vals_z = -1:1,
                          rescale = TRUE,
                          ci_width = 0.95, 
                          ci_type = "confidence",
                          relative_h0 = TRUE,
                          ...) {
  stopif(!isModsemObject(model) && !isLavaanObject(model), "model must be of class ",
         "'modsem_pi', 'modsem_da', 'modsem_mplus' or 'lavaan'")

  if (is.null(xz))
    xz <- paste(x, z, sep = ":")

  xz <- c(xz, reverseIntTerm(xz))

  if (!inherits(model, c("modsem_da", "modsem_mplus")) &&
      !isLavaanObject(model)) 
    xz <- stringr::str_remove_all(xz, ":")

  parTable <- parameter_estimates(model)
  parTable <- getMissingLabels(parTable)


  if (isLavaanObject(model)) {
    vcov <- lavaan::vcov
    nobs <- lavaan::nobs
    coef <- lavaan::coef
  } 
   
  n     <- nobs(model) 
  VCOV  <- vcov(model)
  coefs <- coef(model)

  if (length(n) > 1) {
    # this won't work for multigroup models
    warning("simple_slopes/plot_interaction does not support multigroup models, be vary of results!")
    n <- n[1]
  }

  # Extract coefficients
  beta_x  <- parTable[parTable$lhs == y & parTable$rhs == x & parTable$op == "~", "est"]
  beta_z  <- parTable[parTable$lhs == y & parTable$rhs == z & parTable$op == "~", "est"]
  beta_xz <- parTable[parTable$lhs == y & parTable$rhs %in% xz & parTable$op == "~", "est"]
  beta0_y <- parTable[parTable$lhs == y & parTable$op == "~1", "est"]
  res_y   <- parTable[parTable$lhs == y & parTable$rhs == y & parTable$op == "~~", "est"]

  var_x   <- calcCovParTable(x, x, parTable)
  var_z   <- calcCovParTable(z, z, parTable)

  stopif(length(var_x) == 0, "Variance for x not found in model")
  stopif(length(var_z) == 0, "Variance for z not found in model")
  stopif(length(beta_x) == 0, "Coefficient for x not found in model")
  stopif(length(beta_z) == 0, "Coefficient for z not found in model")
  stopif(length(beta_xz) == 0, "Coefficient for interaction term not found in model")

  label_beta_x  <- parTable[parTable$lhs == y & parTable$rhs == x & parTable$op == "~", "label"]
  label_beta_z  <- parTable[parTable$lhs == y & parTable$rhs == z & parTable$op == "~", "label"]
  label_beta_xz <- parTable[parTable$lhs == y & parTable$rhs %in% xz & parTable$op == "~", "label"]
  label_beta0_y <- parTable[parTable$lhs == y & parTable$op == "~1", "label"]

  label_beta_x  <- ifelse(length(label_beta_x) == 0, NA, label_beta_x)
  label_beta_z  <- ifelse(length(label_beta_z) == 0, NA, label_beta_z)
  label_beta_xz <- ifelse(length(label_beta_xz) == 0, NA, label_beta_xz)
  label_beta0_y <- ifelse(length(label_beta0_y) == 0, NA, label_beta0_y)

  labels <- c(label_beta0_y, label_beta_x, label_beta_z, label_beta_xz)
  VCOV   <- subsetVCOV(VCOV, labels)

  mean_x <- getMean(x, parTable = parTable)
  mean_z <- getMean(z, parTable = parTable)
  mean_y <- getMean(y, parTable = parTable)
  h0     <- if (relative_h0) mean_y else 0

  if (rescale) {
    vals_x <- vals_x * sqrt(var_x) + mean_x
    vals_z <- vals_z * sqrt(var_z) + mean_z
  }

  alpha  <- 1 - ci_width
  ci.sig <- stats::qnorm(1 - alpha / 2) # two-tailed

  # creating margins
  df           <- structure(expand.grid(x = vals_x, z = vals_z), names= c("vals_x", "vals_z"))
  df$predicted <- mean_y + beta_x * df$vals_x + beta_z * df$vals_z + df$vals_z * df$vals_x * beta_xz
  df$std.error <- calc_se(df, e = res_y, VCOV = VCOV, se_type = ci_type)
  df$z.value   <- (df$predicted - h0) / df$std.error # H0 = mean_y
  df$p.value   <- 2 * stats::pnorm(-abs(df$z.value))
  df$ci.upper  <- df$predicted + ci.sig * df$std.error
  df$ci.lower  <- df$predicted - ci.sig * df$std.error


  variable_names <- c(x = x, z = z, y = y)
  attr(df, "variable_names") <- variable_names
  class(df) <- c("simple_slopes", class(df))
  df
}


calc_se <- function(df, e, VCOV, se_type = "confidence") {
  if (se_type == "prediction") 
    return(sqrt(rep(e, nrow(df))))
  else if (se_type != "confidence") 
    warning("se_type must be 'confidence' or 'prediction', using 'confidence'!")

  vals_x <- df$vals_x
  vals_z <- df$vals_z

  n <- nrow(df)
  i <- rep(1, n)
  X <- matrix(c(i, vals_x, vals_z, vals_x * vals_z), nrow=n)
  V <- calcSESimpleSlopes(X, VCOV)

  sqrt(as.vector(V))
}



printTable <- function(x, header = NULL) {
  if (!NROW(x)) return(NULL)

  for (i in seq_len(nrow(x))) {
    str <- paste(x[i, ], collapse = " | ")
    
    if (i == 1) {
      sep <- paste0(strrep("-", nchar(str)), "\n")
      cat(sep)
      cat(str, "\n")
      cat(sep)
    } else cat(str, "\n")
  }
}


#' @export
print.simple_slopes <- function(x, digits = 2, scientific.p = FALSE, ...) {

  variables  <- attr(x, "variable_names")
  predictors <- variables[1:2]
  outcome    <- variables[3]
  header     <- c(predictors[1], sprintf("Predicted %s", outcome), "Std.Error", "z.value", "p.value", "Conf.Interval")
  ci.lower   <- format(x$ci.lower, digits = digits)
  ci.upper   <- format(x$ci.upper, digits = digits)
  cat_z      <- as.factor(round(x$vals_z, digits))

  X <- data.frame(vals_x    = format(x$vals_x, digits = digits),
                  predicted = format(x$predicted, digits = digits),
                  std.error = format(x$std.error, digits = digits),
                  z.value   = format(x$z.value, digits = digits),
                  p.value   = formatPval(x$p.value, scientific = scientific.p),
                  ci        = paste0("[", ci.lower, ", ", ci.upper, "]"))

  X1 <- matrix(header, nrow = 1)
  X2 <- as.matrix(X)

  X <- apply(rbind(X1, X2), MARGIN = 2, format, digits = digits, justify = "right")
  
  for (z_i in unique(cat_z)) {
    Z1 <- X[1, ]
    Z2 <- X[2:nrow(X),][cat_z == z_i,]
    Z  <- rbind(Z1, Z2)

    printf("\nPredicted %s, given %s = %s:\n", outcome, predictors[2], z_i)
    printTable(Z, header = header)
    cat("\n") 
  }
}


subsetVCOV <- function(VCOV, labels) {
  vlabels <- colnames(VCOV)

  noNa <- labels
  noNa[is.na(noNa) | !noNa %in% vlabels] <- vlabels[1]

  V <- VCOV[noNa, noNa]
  V[is.na(labels), ] <- 0
  V[, is.na(labels)] <- 0

  V
}
