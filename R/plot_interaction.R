#' Plot Interaction Effects
#'
#' @param x The name of the variable on the x-axis
#' @param z The name of the moderator variable
#' @param y The name of the outcome variable
#' @param xz The name of the interaction term. If the interaction term is not specified, it
#' will be created using \code{x} and \code{z}.
#' @param vals_x The values of the \code{x} variable to plot, the more values the smoother the std.error-area will be
#' @param vals_z The values of the moderator variable to plot. A separate regression
#' line (\code{y ~ x | z}) will be plotted for each value of the moderator variable
#' @param model An object of class \code{\link{modsem_pi}}, \code{\link{modsem_da}}, or \code{\link{modsem_mplus}}
#' @param alpha_se The alpha level for the std.error area
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
#' plot_interaction("X", "Z", "Y", "X:Z", -3:3, c(-0.2, 0), est1)
#'
#' tpb <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   # Causal Relationsships
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#'   # BEH ~ ATT:PBC
#'   BEH ~ PBC:INT
#'   # BEH ~ PBC:PBC
#' "
#'
#' est2 <- modsem(tpb, TPB, method = "lms")
#' plot_interaction(x = "INT", z = "PBC", y = "BEH", xz = "PBC:INT",
#'                  vals_z = c(-0.5, 0.5), model = est2)
#' }
plot_interaction <- function(x, z, y, xz = NULL, vals_x = seq(-3, 3, .001) ,
                             vals_z, model, alpha_se = 0.15, ...) {
  if (!isModsemObject(model) && !isLavaanObject(model)) {
    stop2("model must be of class 'modsem_pi', 'modsem_da', 'modsem_mplus' or 'lavaan'")
  }

  if (is.null(xz)) xz <- paste(x, z, sep = ":")
  xz <- c(xz, reverseIntTerm(xz))
  if (!inherits(model, c("modsem_da", "modsem_mplus")) &&
      !isLavaanObject(model)) {
    xz <- stringr::str_remove_all(xz, ":")
  }

  parTable <- parameter_estimates(model)
  gamma_x <- parTable[parTable$lhs == x & parTable$op == "~", "est"]

  if (isLavaanObject(model)) {
    # this won't work for multigroup models
    nobs <- unlist(model@Data@nobs)
    if (length(nobs) > 1) warning2("plot_interaction is not intended for multigroup models")
    n <- nobs[[1]]

  } else {
    n <- nrow(model$data)
  }

  lVs <- c(x, z, y, xz)
  coefs <- parTable[parTable$op == "~" & parTable$rhs %in% lVs &
                    parTable$lhs == y, ]
  vars <- parTable[parTable$op == "~~" & parTable$rhs %in% lVs &
             parTable$lhs == parTable$rhs, ]
  gamma_x <- coefs[coefs$rhs == x, "est"]
  var_x <- calcCovParTable(x, x, parTable)
  gamma_z <- coefs[coefs$rhs == z, "est"]
  var_z <- calcCovParTable(z, z, parTable)
  gamma_xz <- coefs[coefs$rhs %in% xz, "est"]
  sd <- sqrt(vars[vars$rhs == y, "est"]) # residual std.error

  if (length(gamma_x) == 0) stop2("coefficient for x not found in model")
  if (length(var_x) == 0) stop2("variance of x not found in model")
  if (length(gamma_z) == 0) stop2("coefficient for z not found in model")
  if (length(var_z) == 0) stop2("variance of z not found in model")
  if (length(gamma_xz) == 0) stop2("coefficient for xz not found in model")
  if (length(sd) == 0) stop2("residual std.error of y not found in model")

  # creating margins
  df <- expand.grid(x = vals_x, z = vals_z)
  df$se_x <- calc_se(df$x, var = var_x, n = n, s = sd)
  df$proj_y <- gamma_x * df$x + gamma_z + df$z + df$z * df$x * gamma_xz
  df$cat_z <- as.factor(df$z)

  se_x <- df$se_x
  proj_y <- df$proj_y
  cat_z <- df$cat_z
  # plotting margins
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = proj_y, colour = cat_z, group = cat_z,)) +
    ggplot2::geom_smooth(method = "lm", formula = "y ~ x", se = FALSE) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = proj_y - 1.96 * se_x, ymax = proj_y + 1.96 * se_x),
                         alpha = alpha_se, linewidth = 0, linetype = "blank") +
    ggplot2::labs(x = x, y = y, colour = z)
}


# function for calculating std.error of predicted value
calc_se <- function(x, var, n, s) {
  # x = values of x (predictor),
    # this function assumes that 'mean(x) = 0'
  # var = variance of x
  # n = sample size
  # s = residual std.error of model
  SSx <- (n - 1) * var # sum of squares of x
  s * sqrt(1 / n + x ^ 2 / SSx)
}

#' Plot Interaction Effect Using the Johnson-Neyman Technique
#'
#' This function plots the simple slopes of an interaction effect across different values of a moderator variable using the Johnson-Neyman technique. It identifies regions where the effect of the predictor on the outcome is statistically significant.
#'
#' @param x The name of the predictor variable (as a character string).
#' @param z The name of the moderator variable (as a character string).
#' @param y The name of the outcome variable (as a character string).
#' @param xz The name of the interaction term. If not specified, it will be created using \code{x} and \code{z}.
#' @param model A fitted model object of class \code{modsem_da}, \code{modsem_mplus}, \code{modsem_pi}, or \code{lavaan}.
#' @param min_z The minimum value of the moderator variable \code{z} to be used in the plot (default is -3).
#' @param max_z The maximum value of the moderator variable \code{z} to be used in the plot (default is 3).
#' @param alpha The significance level for the confidence intervals (default is 0.05).
#' @param ... Additional arguments (currently not used).
#' @return A \code{ggplot} object showing the interaction plot with regions of significance.
#' @details
#' The function calculates the simple slopes of the predictor variable \code{x} on the outcome variable \code{y} at different levels of the moderator variable \code{z}. It uses the Johnson-Neyman technique to identify the regions of \code{z} where the effect of \code{x} on \code{y} is statistically significant.
#'
#' It extracts the necessary coefficients and variance-covariance information from the fitted model object. The function then computes the critical t-value and solves the quadratic equation derived from the t-statistic equation to find the Johnson-Neyman points.
#'
#' The plot displays:
#' \itemize{
#'   \item The estimated simple slopes across the range of \code{z}.
#'   \item Confidence intervals around the slopes.
#'   \item Regions where the effect is significant (shaded areas).
#'   \item Vertical dashed lines indicating the Johnson-Neyman points.
#'   \item Text annotations providing the exact values of the Johnson-Neyman points.
#' }
#'
#' @examples
#' \dontrun{
#' library(modsem)
#'
#' m1 <-  ' 
#'   visual  =~ x1 + x2 + x3 
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' 
#'   visual ~ speed + textual + speed:textual
#' '
#' 
#' est <- modsem(m1, data = lavaan::HolzingerSwineford1939, method = "ca")
#' plot_jn(x = "speed", z = "textual", y = "visual", model = est, max_z = 6)
#' }
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_vline annotate scale_fill_manual labs theme_minimal
#' @export
plot_jn <- function(x, z, y, xz = NULL, model, min_z = -3, max_z = 3, alpha = 0.05, ...) {
  # Check if model is a valid object
  stopif(!inherits(model, c("modsem_da", "modsem_mplus", "modsem_pi", "lavaan")),
         "model must be of class 'modsem_pi', 'modsem_da', 'modsem_mplus', or 'lavaan'")

  # If interaction term is not specified, create it
  if (is.null(xz)) xz <- paste(x, z, sep = ":")
  xz <- c(xz, reverseIntTerm(xz))
  if (!inherits(model, c("modsem_da", "modsem_mplus")) && !inherits(model, "lavaan")) {
    xz <- stringr::str_remove_all(xz, ":")
  }

  # Extract parameter estimates
  parTable <- parameter_estimates(model)
  parTable <- getMissingLabels(parTable)

  # Extract coefficients
  beta_x <- parTable[parTable$lhs == y & parTable$rhs == x & parTable$op == "~", "est"]
  beta_xz <- parTable[parTable$lhs == y & parTable$rhs %in% xz & parTable$op == "~", "est"]

  stopif(length(beta_x) == 0, "Coefficient for x not found in model")
  stopif(length(beta_xz) == 0, "Coefficient for interaction term not found in model")

  # Extract variance-covariance matrix
  vcov_matrix <- vcov(model)

  label_beta_x  <- parTable[parTable$lhs == y & parTable$rhs == x & parTable$op == "~", "label"]
  label_beta_xz <- parTable[parTable$lhs == y & parTable$rhs %in% xz & parTable$op == "~", "label"]

  # Extract variances and covariances
  var_beta_x <- vcov_matrix[label_beta_x, label_beta_x]
  var_beta_xz <- vcov_matrix[label_beta_xz, label_beta_xz]
  cov_beta_x_beta_xz <- vcov_matrix[label_beta_x, label_beta_xz]

  nobs <- nobs(model)
  npar <- length(coef(model))

  df_resid <- nobs - npar
  
  stopif(df_resid < 1, "Degrees of freedom for residuals must be greater than 0. ",
         "The model may have fewer observations than parameters.")

  # Critical t-value
  t_crit <- stats::qt(1 - alpha / 2, df_resid)

  # Quadratic equation components
  A <- beta_xz^2 - t_crit^2 * var_beta_xz
  B <- 2 * beta_x * beta_xz - 2 * t_crit^2 * cov_beta_x_beta_xz
  C <- beta_x^2 - t_crit^2 * var_beta_x

  discriminant <- B^2 - 4 * A * C

  if (A == 0) {
    # Linear case
    if (B != 0) {
      z_jn <- -C / B
      significant_everywhere <- FALSE
      z_lower <- z_jn
      z_upper <- z_jn
    } else {
      message("No regions where the effect transitions between significant and non-significant.")
      significant_everywhere <- TRUE
    }
  } else if (discriminant < 0) {
    message("No regions where the effect transitions between significant and non-significant.")
    significant_everywhere <- TRUE
  } else if (discriminant == 0) {
    # One real root
    z_jn <- -B / (2 * A)
    significant_everywhere <- FALSE
    z_lower <- z_jn
    z_upper <- z_jn
  } else {
    # Two real roots
    z1 <- (-B + sqrt(discriminant)) / (2 * A)
    z2 <- (-B - sqrt(discriminant)) / (2 * A)
    z_lower <- min(z1, z2)
    z_upper <- max(z1, z2)
    significant_everywhere <- FALSE
  }

  z_range <- seq(min_z, max_z, length.out = 100)

  # Calculate slopes and statistics
  slope <- beta_x + beta_xz * z_range
  SE_slope <- sqrt(var_beta_x + z_range^2 * var_beta_xz + 2 * z_range * cov_beta_x_beta_xz)
  t_value <- slope / SE_slope
  p_value <- 2 * (1 - stats::pt(abs(t_value), df_resid))
  significant <- p_value < alpha

  # Create plotting data frame
  df_plot <- data.frame(z = z_range, slope = slope, SE = SE_slope, t = t_value, p = p_value, significant = significant)

  # Plotting
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = z, y = slope)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = slope - t_crit * SE_slope, ymax = slope + t_crit * SE_slope, fill = significant), alpha = 0.2) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey"), guide = "none") +
    ggplot2::labs(x = z, y = paste("Simple slope of", x, "on", y)) +
    ggplot2::theme_minimal()

  # Add Johnson-Neyman points if applicable
  if (!significant_everywhere) {
    # Ensure JN points are within the range of z
    if (exists("z_jn")) {
      # Single JN point
      if (z_jn >= min_z && z_jn <= max_z) {
        p <- p +
          ggplot2::geom_vline(xintercept = z_jn, linetype = "dashed", color = "red") +
          ggplot2::annotate("text", x = z_jn, y = max(df_plot$slope), label = paste("JN point:", round(z_jn, 2)),
                            hjust = -0.1, vjust = 1, color = "black")
      }
    } else {
      # Two JN points
      if (z_lower >= min_z && z_lower <= max_z && z_upper >= min_z && z_upper <= max_z) {
        # Add light red fill between z_lower and z_upper
        p <- p +
          ggplot2::annotate(geom="rect", xmin = z_lower, xmax = z_upper, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1)
      }
      if (z_lower >= min_z && z_lower <= max_z) {
        p <- p +
          ggplot2::geom_vline(xintercept = z_lower, linetype = "dashed", color = "red") +
          ggplot2::annotate("text", x = z_lower, y = max(df_plot$slope), label = paste("JN point:", round(z_lower, 2)),
                            hjust = -0.1, vjust = 1, color = "black")
      }
      if (z_upper >= min_z && z_upper <= max_z) {
        p <- p +
          ggplot2::geom_vline(xintercept = z_upper, linetype = "dashed", color = "red") +
          ggplot2::annotate("text", x = z_upper, y = max(df_plot$slope), label = paste("JN point:", round(z_upper, 2)),
                            hjust = -0.1, vjust = 1, color = "black")
      }
    }
  }

  p
}

# Helper function to reverse interaction term
reverseIntTerm <- function(xz) {
  if (grepl(":", xz)) {
    vars <- strsplit(xz, ":")[[1]]
    rev_xz <- paste(rev(vars), collapse = ":")
  } else {
    rev_xz <- xz
  }
  rev_xz
}
