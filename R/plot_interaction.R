#' Plot Interaction Effects 
#' 
#' @param x The name of the variable on the x-axis 
#' @param z The name of the moderator variable 
#' @param y The name of the outcome variable 
#' @param xz The name of the interaction term
#' @param vals_x The values of the x variable to plot, the more values the smoother the std.error-area will be
#' @param vals_z The values of the moderator variable to plot. A seperate regression 
#' line ("y ~ x | z") will be plotted for each value of the moderator variable
#' @param model An object of class `modsem_pi`, `modsem_lms_qml`, or `modsem_mplus` 
#' @param alpha_se The alpha level for the std.error area
#' @param ... Additional arguments passed to other functions 
#' @return A ggplot object
#' @export 
plot_interaction <- function(x, z, y, xz, vals_x = seq(-3, 3, .001) , vals_z, model, 
                             alpha_se = 0.15, ...) {
  if (!inherits(model, c("modsem_pi", "modsem_lms", 
                         "modsem_mplus", "modsem_qml"))) {
    stop("model must be of class 'modsem_pi', 'modsem_lms_qml', or 'modsem_mplus'")
  }

  if (!inherits(model, c("modsem_lms", "modsem_qml"))) {
    xz <- stringr::str_remove_all(xz, ":")
  }

  parTable <- parameter_estimates(model)
  gamma_x <- parTable[parTable$lhs == x & parTable$op == "~", "est"] 
  
  data <- model$data
  n <- nrow(data)
  lVs <- c(x, z, y, xz)
  coefs <- parTable[parTable$op == "~" & parTable$rhs %in% lVs, ]
  vars <- parTable[parTable$op == "~~" & parTable$rhs %in% lVs &
             parTable$lhs == parTable$rhs, ]
  gamma_x <- coefs[coefs$rhs == x, "est"]
  var_x <- vars[vars$rhs == x, "est"]
  gamma_z <- coefs[coefs$rhs == z, "est"]
  var_z <- vars[vars$rhs == z, "est"]
  gamma_xz <- coefs[coefs$rhs == xz, "est"]
  sd <- sqrt(vars[vars$rhs == y, "est"]) # residual std.error
  if (length(gamma_x) == 0) stop("coefficient for x not found in model")
  if (length(var_x) == 0) stop("variance of x not found in model")
  if (length(gamma_z) == 0) stop("coefficient for z not found in model")
  if (length(var_z) == 0) stop("variance of z not found in model")
  if (length(gamma_xz) == 0) stop("coefficient for xz not found in model")
  if (length(sd) == 0) stop("residual std.error of y not found in model")

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
                         alpha = alpha_se, linewidth = 0) +
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
