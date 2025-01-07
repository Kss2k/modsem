#' Get the simple slopes of a SEM model
#'
#' @param x The name of the variable on the x-axis
#' @param z The name of the moderator variable
#' @param y The name of the outcome variable
#' @param xz The name of the interaction term. If the interaction term is not specified, it
#' @param model An object of class \code{\link{modsem_pi}}, \code{\link{modsem_da}}, or \code{\link{modsem_mplus}}
#' will be created using \code{x} and \code{z}.
#' @param vals_x The values of the \code{x} variable to plot, the more values the smoother the std.error-area will be.
#' NOTE: \code{vals_x} are measured relative to the mean and standard deviation of \code{x} (overridden by \code{rescale=FALSE}). 
#' The correct values will show up in the table.
#' @param vals_z The values of the moderator variable to plot. A separate regression
#' NOTE: \code{vals_z} are measured relative to the mean and standard deviation of \code{z} (overridden by \code{rescale=FALSE}). 
#' The correct values will show up in the table.
#' line (\code{y ~ x | z}) will be plotted for each value of the moderator variable
#' @param rescale Logical. If \code{TRUE} (default), the values of \code{x} and \code{z} will be rescaled relative to their means and standard deviations.
#' @param ci_width The width of the confidence interval (default is 1.96, corresponding to a 95\% confidence interval)
#' @param relative_h0 Logical. If \code{TRUE} (default), the null hypothesis is that the mean of the outcome variable is equal to the predicted value of the outcome variable. If \code{FALSE}, the null hypothesis is that the outcome variable is equal to zero.
#' @param ... Additional arguments passed to other functions
#' @return A \code{ggplot} object
#' @export
#' @examples
#' library(modsem)
#' \dontrun{
#' m1 <- "
#' # Outer Model
#'   X =~ x1
#'   X =~ x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'
#' # Inner model
#'   Y ~ X + Z + X:Z
#' "
#' est1 <- modsem(m1, data = oneInt)
#' simple_slopes(x = "X", z = "Z", y = "Y", model = est1)
#' }
simple_slopes <- function(x, z, y,
                          xz = NULL,
                          model,
                          vals_x = -3:3,
                          vals_z = -1:1,
                          rescale = TRUE,
                          ci_width = 1.96, 
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
  gamma_x <- parTable[parTable$lhs == x & parTable$op == "~", "est"]

  if (isLavaanObject(model)) {
    # this won't work for multigroup models
    nobs <- unlist(model@Data@nobs)
    warnif(length(nobs) > 1, "plot_interaction is not intended for multigroup models")
    n <- nobs[[1]]

  } else {
    n <- nrow(model$data)
  }

  lVs <- c(x, z, y, xz)
  coefs <- parTable[parTable$op == "~" & parTable$rhs %in% lVs &
                    parTable$lhs == y, ]
  vars <- parTable[parTable$op == "~~" & parTable$rhs %in% lVs &
             parTable$lhs == parTable$rhs, ]
  gamma_x  <- coefs[coefs$rhs == x, "est"]
  var_x    <- calcCovParTable(x, x, parTable)
  gamma_z  <- coefs[coefs$rhs == z, "est"]
  var_z    <- calcCovParTable(z, z, parTable)
  gamma_xz <- coefs[coefs$rhs %in% xz, "est"]
  sd       <- sqrt(vars[vars$rhs == y, "est"]) # residual std.error

  stopif(!length(gamma_x),  "coefficient for x not found in model")
  stopif(!length(var_x),    "variance of x not found in model")
  stopif(!length(gamma_z),  "coefficient for z not found in model")
  stopif(!length(var_z),    "variance of z not found in model")
  stopif(!length(gamma_xz), "coefficient for xz not found in model")
  stopif(!length(sd),       "residual std.error of y not found in model")

  mean_x <- getMean(x, parTable = parTable)
  mean_z <- getMean(z, parTable = parTable)
  mean_y <- getMean(y, parTable = parTable)
  h0     <- if (relative_h0) mean_y else 0

  if (rescale) {
    vals_x <- vals_x * sqrt(var_x) + mean_x
    vals_z <- vals_z * sqrt(var_z) + mean_z
  }

  # creating margins
  df           <- expand.grid(x = vals_x, z = vals_z)
  colnames(df) <- c("vals_x", "vals_z")
  df$predicted <- mean_y + gamma_x * df$vals_x + gamma_z + df$vals_z + df$vals_z * df$vals_x * gamma_xz
  df$std.error <- calc_se(df$vals_x, var = var_x, n = n, s = sd)
  df$z.value   <- (df$predicted - h0) / df$std.error # H0 = mean_y
  df$p.value   <- 2 * stats::pnorm(-abs(df$z.value))
  df$ci.upper  <- df$predicted + ci_width * df$std.error
  df$ci.lower  <- df$predicted - ci_width * df$std.error


  variable_names <- c(x = x, z = z, y = y)
  attr(df, "variable_names") <- variable_names
  class(df) <- c("simple_slopes", class(df))
  df
}


#' @export
print.simple_slopes <- function(x, digits = 2, scientific.p = FALSE, ...) {

  variables  <- attr(x, "variable_names")
  predictors <- variables[1:2]
  outcome    <- variables[3]
  header     <- c(predictors, sprintf("E[%s]", outcome), "Std.Error", "z.value", "p.value", "95% CI")
  ci.lower   <- format(x$ci.lower, digits = digits)
  ci.upper   <- format(x$ci.upper, digits = digits)

  y <- data.frame(vals_x    = format(x$vals_x, digits = digits),
                  vals_z    = format(x$vals_z, digits = digits),
                  predicted = format(x$predicted, digits = digits),
                  std.error = format(x$std.error, digits = digits),
                  z.value   = format(x$z.value, digits = digits),
                  p.value   = formatPval(x$p.value, scientific = scientific.p),
                  ci        = paste0("[", ci.lower, ", ", ci.upper, "]"))

  z1 <- matrix(header, nrow = 1)
  z2 <- as.matrix(y)

  z <- rbind(z1, z2)
  z[, 1:3] <- apply(z[, 1:3], MARGIN = 2, format, digits = digits, justify = "right")
  z[, 4:7] <- apply(z[, 4:7], MARGIN = 2, format, digits = digits, justify = "right")
  
  printf("\nPredicted values of %s\n\n", outcome)

  for (i in seq_len(nrow(z))) {
    str <- paste(z[i, ], collapse = " | ")

    cat(str, "\n")
    if (i == 1) cat(strrep("-", nchar(str)), "\n")
  }
}
