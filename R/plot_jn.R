getJN_PointsConditional <- function(x, z, y, parTable, model, min_z, max_z, sig.level, alpha,
                                    detail, sd.line, standardized, xz, greyscale,
                                    plot.jn.points = TRUE, group = NULL, group.label = NULL,
                                    type = "indirect", mc.reps = 10000, ...) {
  type <- match.arg(type, c("indirect", "total", "direct"))

  mean_z <- getMean(z, parTable = parTable)
  sd_z   <- sqrt(calcCovParTable(x = z, y = z, parTable = parTable))
  stopif(is.na(sd_z), sprintf("Variance for %s not found in model", z))

  z_min <- min_z + mean_z
  z_max <- max_z + mean_z

  labels <- getParTableLabels(parTable, labelCol = "label")
  labels <- stringr::str_replace_all(labels, OP_REPLACEMENTS)
  parTable$label <- labels

  intTerms <- getIntTerms(parTable)

  if (!length(intTerms)) {
    stop2("No interaction terms were found in the model.")
  }

  for (intTerm in intTerms) {
    elems <- stringr::str_split_1(intTerm, pattern = ":")

    if (!z %in% elems)
      next
    else if (length(elems) > 2L)
      stop2("Johnson-Neyman plot is not available for three-way (or higher) interactions!")

    m <- elems[elems != z]

    simpleEffects <- parTable[parTable$rhs == intTerm & parTable$op == "~", , drop = FALSE]

    for (i in seq_len(NROW(simpleEffects))) {
      row   <- simpleEffects[i, ]
      cond1 <- parTable$lhs == row$lhs & parTable$rhs == m
      cond2 <- parTable$lhs == row$lhs & parTable$rhs == intTerm
      beta1 <- parTable[cond1, "label"]
      beta2 <- parTable[cond2, "label"]

      if (length(beta1) && length(beta2)) {
        parTable[cond1, "label"] <- paste0(
          "(", beta1[[1L]], "+", beta2[[1L]], "*", ".NUMERIC_VALUE__Z)"
        )
      }
    }
  }

  parTable.simple <- parTable[!grepl(":", parTable$lhs) &
                              !grepl(":", parTable$rhs), , drop = FALSE]

  paths <- getJN_PathsParTable(x = x, y = y, parTable = parTable.simple)
  stopif(is.null(paths), sprintf("No paths between %s and %s were found!", x, y))

  is.direct <- attr(paths, "is.direct")

  expr_paths <- list(
    total    = paths,
    direct   = paths[is.direct],
    indirect = paths[!is.direct]
  )

  keep <- expr_paths[[type]]
  stopif(!length(keep), sprintf("No %s effect between %s and %s was found!", type, x, y))
  expr <- parse(text = paste0(keep, collapse = "+"))

  V <- tryCatch(vcov(model), error = \(e) NULL)
  coefs <- structure(parTable$est, names = labels)

  if (is.null(V)) {
    k    <- length(coefs)
    pars <- names(coefs)
    V    <- matrix(0, nrow = k, ncol = k, dimnames = list(pars, pars))
  } else {
    dimnames(V) <- list(
      stringr::str_replace_all(rownames(V), OP_REPLACEMENTS),
      stringr::str_replace_all(colnames(V), OP_REPLACEMENTS)
    )
  }

  labels.u <- unique(labels)
  V        <- expandVCOV(V, labels = labels.u)
  coefs    <- coefs[labels.u]

  draws <- mvtnorm::rmvnorm(mc.reps, mean = coefs, sigma = V)
  if (is.null(dim(draws)))
    draws <- matrix(draws, nrow = 1L)

  COEFS    <- as.data.frame(draws)

  getrow <- function(.NUMERIC_VALUE__Z) {
    slopex <- with(COEFS, eval(expr))
    ci <- stats::quantile(slopex, probs = c(sig.level / 2, 1 - sig.level / 2), na.rm = TRUE)
    point <- mean(slopex, na.rm = TRUE)
    c(point, ci, .NUMERIC_VALUE__Z)
  }

  z_range   <- seq(z_min, z_max, length.out = detail)
  jn_points <- do.call(rbind, lapply(z_range, getrow))
  colnames(jn_points) <- c("slope", "lower", "upper", "z")
  jn_points <- as.data.frame(jn_points)
  jn_points$significant <- (jn_points$lower > 0) | (jn_points$upper < 0)

  significant_everywhere <- all(jn_points$significant) || !any(jn_points$significant)

  list(
    grid = jn_points,
    jn_points = approximate_jn_points_from_grid(jn_points),
    significant_everywhere = significant_everywhere,
    mean_z = mean_z,
    sd_z = sd_z,
    min_z = z_min,
    max_z = z_max
  )
}
  

#' Plot Interaction Effect Using the Johnson-Neyman Technique
#'
#' This function plots the simple slopes of an interaction effect across different values of a moderator variable using the Johnson-Neyman technique. It identifies regions where the effect of the predictor on the outcome is statistically significant.
#'
#' @param x The name of the predictor variable (as a character string).
#' @param z The name of the moderator variable (as a character string).
#' @param y The name of the outcome variable (as a character string).
#' @param model A fitted model object of class \code{modsem_da}, \code{modsem_mplus}, \code{modsem_pi}, or \code{lavaan}.
#' @param min_z The minimum value of the moderator variable \code{z} to be used in the plot (default is -3). It is relative to the mean of z.
#' @param max_z The maximum value of the moderator variable \code{z} to be used in the plot (default is 3). It is relative to the mean of z.
#' @param sig.level The alpha-criterion for the confidence intervals (default is 0.05).
#' @param alpha alpha setting used in \code{ggplot} (i.e., the opposite of opacity)
#' @param detail The number of generated data points to use for the plot (default is 1000). You can increase this value for smoother plots.
#' @param sd.line A thick black line showing \code{+/- sd.line * sd(z)}. NOTE: This line will be truncated by \code{min_z} and \code{max_z} if
#' the sd.line falls outside of \code{[min_z, max_z]}.
#' @param standardized Should coefficients be standardized beforehand?
#' @param xz The name of the interaction term. If not specified, it will be created using \code{x} and \code{z}.
#' @param greyscale Logical. If \code{TRUE} the plot is plotted in greyscale.
#' @param plot.jn.points Logical. If \code{TRUE}, omit the numeric annotations for the JN-points from the plot.
#' @param type Which effect to display. One of \code{"direct"}, \code{"indirect"}, or \code{"total"}.
#' @param mc.reps Number of Monte Carlo replicates used to approximate the confidence
#'   bands when \code{type} is \code{"indirect"} or \code{"total"}.
#' @param ... Additional arguments (currently not used).
#'
#' @return A \code{ggplot} object showing the interaction plot with regions of significance.
#' @details
#' The function calculates the simple slopes of the predictor variable \code{x} on the outcome variable \code{y} at different levels of the moderator variable \code{z}. It uses the Johnson-Neyman technique to identify the regions of \code{z} where the effect of \code{x} on \code{y} is statistically significant.
#'
#' When plotting indirect or total effects, the function relies on Monte Carlo
#' draws from the estimated sampling distribution of the parameters to approximate
#' the conditional effect and its confidence interval across the moderator range.
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
#' @export
plot_jn <- function(x, z, y, model, min_z = -3, max_z = 3,
                    sig.level = 0.05, alpha = 0.2, detail = 1000,
                    sd.line = 2, standardized = FALSE, xz = NULL,
                    greyscale = FALSE, plot.jn.points = TRUE,
                    type = c("direct", "indirect", "total"),
                    mc.reps = 10000, ...) {

  stopif(!inherits(model, c("modsem_da", "modsem_mplus", "modsem_pi", "lavaan")),
         "model must be of class 'modsem_pi', 'modsem_da', 'modsem_mplus', or 'lavaan'")

  type <- match.arg(type)

  if (standardized) {
    parTable <- standardized_estimates(model, correction = TRUE)
  } else {
    parTable <- parameter_estimates(
      object = model,
      label.renamed.prod = TRUE
    )
  }

  group.label <- modsem_inspect(model, what = "group.label")
  parTable <- addMissingGroups(getMissingLabels(parTable))

  plots <- list()
  for (g in getGroupsParTable(parTable)) {
    label.g    <- if (length(group.label)) group.label[[g]] else NULL
    parTable.g <- parTable[parTable$group == g, , drop = FALSE]

    plots[[g]] <- plotJN_Group(
      x = x, z = z, y = y, parTable = parTable.g, model = model,
      min_z = min_z, max_z = max_z, sig.level = sig.level,
      alpha = alpha, detail = detail, sd.line = sd.line,
      standardized = standardized, xz = xz, greyscale = greyscale,
      plot.jn.points = plot.jn.points, group = g, group.label = label.g,
      type = type, mc.reps = mc.reps, ...
    )
  }

  if (length(plots) <= 1L)
    return(plots[[1L]])

  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    printf("The `ggpubr` package is needed to arrange Johnson-Neyman plots in multigroup models!\n")
    printf("Do you want to install it? (y/n) ")
    choice <- tolower(substr(readLines(n = 1L), 1L, 1L))

    stopifnot(choice == "y")
    utils::install.packages("ggpubr")
  }

  if (requireNamespace("ggpubr", quietly = TRUE)) { # Make R CMD check happy
    ggpubr::ggarrange(plotlist = plots, labels = group.label)
  } else stop2("The `ggpubr` package is needed to arrange Johnson-Neyman plots in multigroup models!\n")
}


getJN_PointsDirect <- function(x, z, y, parTable, model, min_z, max_z, sig.level, alpha,
                               detail, sd.line, standardized, xz, greyscale,
                               plot.jn.points = TRUE, group = NULL, group.label = NULL, ...) {
  if (is.null(xz))
    xz <- paste(x, z, sep = ":")

  checkInputsSimpleSlopes(x = x, z = z, y = y, xz = xz, parTable = parTable)

  xz <- c(xz, reverseIntTerm(xz))

  if (!inherits(model, c("modsem_da", "modsem_mplus")) && !inherits(model, "lavaan")) {
    xz <- stringr::str_remove_all(xz, ":")
  }

  if (inherits(model, "lavaan")) {
    vcov <- lavaan::vcov
    coef <- lavaan::coef
    nobs <- lavaan::nobs
  }

  # z mean/sd and plotting window
  mean_z <- getMean(z, parTable = parTable)
  sd_z   <- sqrt(calcCovParTable(x = z, y = z, parTable = parTable))
  stopif(is.na(sd_z), sprintf("Variance for %s not found in model", z))
  z_min  <- min_z + mean_z
  z_max  <- max_z + mean_z

  # coefficients
  beta_x  <- parTable[parTable$lhs == y & parTable$rhs == x   & parTable$op == "~", "est"]
  beta_xz <- parTable[parTable$lhs == y & parTable$rhs %in% xz & parTable$op == "~", "est"]
  stopif(length(beta_x)  == 0, "Coefficient for x not found in model")
  stopif(length(beta_xz) == 0, "Coefficient for interaction term not found in model")

  VCOV <- vcov(model)
  label_beta_x  <- parTable[parTable$lhs == y & parTable$rhs == x   & parTable$op == "~", "label"]
  label_beta_xz <- parTable[parTable$lhs == y & parTable$rhs %in% xz & parTable$op == "~", "label"]

  var_beta_x         <- VCOV[label_beta_x,  label_beta_x]
  var_beta_xz        <- VCOV[label_beta_xz, label_beta_xz]
  cov_beta_x_beta_xz <- VCOV[label_beta_x,  label_beta_xz]

  nobs <- nobs(model)
  npar <- length(coef(model))
  df_resid <- nobs - npar
  if (df_resid < 1) {
    warning2("Degrees of freedom for residuals must be greater than 0. ",
             "Using sample size instead of degrees of freedom")
    df_resid <- nobs
  }

  t_crit <- stats::qt(1 - sig.level / 2, df_resid)

  # Johnson–Neyman roots
  A <- beta_xz^2 - t_crit^2 * var_beta_xz
  B <- 2 * beta_x * beta_xz - 2 * t_crit^2 * cov_beta_x_beta_xz
  C <- beta_x^2 - t_crit^2 * var_beta_x
  disc <- B^2 - 4 * A * C

  significant_everywhere <- FALSE
  jn_points <- numeric(0)
  if (A == 0) {
    if (B != 0) {
      z_jn <- -C / B; z_lower <- z_jn; z_upper <- z_jn
      jn_points <- z_jn
    } else {
      message("No regions where the effect transitions between significant and non-significant.")
      significant_everywhere <- TRUE
    }
  } else if (disc < 0) {
    message("No regions where the effect transitions between significant and non-significant.")
    significant_everywhere <- TRUE
  } else if (disc == 0) {
    z_jn <- -B / (2 * A); z_lower <- z_jn; z_upper <- z_jn
    jn_points <- z_jn
  } else {
    z1 <- (-B + sqrt(disc)) / (2 * A)
    z2 <- (-B - sqrt(disc)) / (2 * A)
    z_lower <- min(z1, z2); z_upper <- max(z1, z2)
    jn_points <- c(z_lower, z_upper)
  }

  jn_points <- jn_points[is.finite(jn_points)]

  # grid and simple slopes
  z_range  <- seq(z_min, z_max, length.out = detail)
  slope    <- beta_x + beta_xz * z_range
  SE_slope <- sqrt(var_beta_x + z_range^2 * var_beta_xz + 2 * z_range * cov_beta_x_beta_xz)
  t_value  <- slope / SE_slope
  p_value  <- 2 * (1 - stats::pt(abs(t_value), df_resid))
  significant <- p_value < sig.level

  lower_all <- slope - t_crit * SE_slope
  upper_all <- slope + t_crit * SE_slope

  grid <- data.frame(
    z = z_range,
    slope = as.numeric(slope),
    lower = as.numeric(lower_all),
    upper = as.numeric(upper_all),
    significant = significant
  )

  if (!length(jn_points) && !significant_everywhere) {
    jn_points <- approximate_jn_points_from_grid(grid)
  }

  list(
    grid = grid,
    jn_points = jn_points,
    significant_everywhere = significant_everywhere,
    mean_z = mean_z,
    sd_z = sd_z,
    min_z = z_min,
    max_z = z_max
  )
}


plotJN_Group <- function(x, z, y, parTable, model, min_z, max_z, sig.level, alpha,
                         detail, sd.line, standardized, xz, greyscale,
                         plot.jn.points = TRUE, group = NULL, group.label = NULL,
                         type = c("direct", "indirect", "total"),
                         mc.reps = 10000, ...) {

  type <- match.arg(type)

  result <- if (type == "direct") {
    getJN_PointsDirect(
      x = x, z = z, y = y, parTable = parTable, model = model, min_z = min_z,
      max_z = max_z, sig.level = sig.level, alpha = alpha, detail = detail,
      sd.line = sd.line, standardized = standardized, xz = xz, greyscale = greyscale,
      plot.jn.points = plot.jn.points, group = group, group.label = group.label
    )
  } else {
    getJN_PointsConditional(
      x = x, z = z, y = y, parTable = parTable, model = model, min_z = min_z,
      max_z = max_z, sig.level = sig.level, alpha = alpha, detail = detail,
      sd.line = sd.line, standardized = standardized, xz = xz, greyscale = greyscale,
      plot.jn.points = plot.jn.points, group = group, group.label = group.label,
      type = type, mc.reps = mc.reps
    )
  }

  df_plot   <- result$grid
  mean_z    <- result$mean_z
  sd_z      <- result$sd_z
  z_min     <- result$min_z
  z_max     <- result$max_z
  jn_points <- result$jn_points
  significant_everywhere <- isTRUE(result$significant_everywhere)

  significant <- df_plot$significant
  sig_all <- all(significant)
  sig_any <- any(significant)

  if (sig_all || !sig_any) {
    message("No regions where the effect transitions between significant and non-significant.")
  }

  if (sig_any && !sig_all) {
    format_num <- function(val) formatC(val, format = "f", digits = 2)
    format_sig <- function(val) sub("^0\\.", ".", formatC(val, format = "f", digits = 2))

    sig_points <- if (length(jn_points)) jn_points else df_plot$z[significant]
    interval <- sprintf("[%s, %s]", format_num(min(sig_points)),
                        format_num(max(sig_points)))

    if (!is.null(group.label)) {
      header <- sprintf("Johnson-Neyman Interval (group %s):", group.label)
    } else {
      header <- "Johnson-Neyman Interval:"
    }

    body <- sprintf("When %s is outside the interval %s, the slope of %s is p < %s.",
                    z, interval, x, format_sig(sig.level))
    message(sprintf("%s\n  %s", header, body))
  }


  # split into contiguous runs to avoid polygon bleed
  run_id <- cumsum(c(0, diff(as.integer(significant)) != 0))
  siglabel <- sprintf("p < %s", sig.level)
  significance_chr <- ifelse(significant, siglabel, "n.s.")
  Significance <- factor(significance_chr, levels = c(siglabel, "n.s."))

  x_start <- mean_z - sd.line * sd_z
  x_end   <- mean_z + sd.line * sd_z

  if (x_start < z_min && x_end > z_max) {
    warning2("Truncating SD-range on the right and left!")
  } else if (x_start < z_min) {
    warning2("Truncating SD-range on the left!")
  } else if (x_end > z_max) {
    warning2("Truncating SD-range on the right!")
  }
  x_start <- max(x_start, z_min)
  x_end   <- min(x_end, z_max)
  y_start <- y_end <- 0

  hline_label <- sprintf("+/- %s SDs of %s", sd.line, z)
  data_hline <- data.frame(x_start, x_end, y_start, y_end, hline_label)

  # unified legend for colour/fill
  breaks <- c(siglabel, "n.s.", hline_label)
  values <- structure(c("cyan3", "red", "black"), names = breaks)

  y_range <- range(c(df_plot$lower, df_plot$upper, 0), na.rm = TRUE)
  if (!all(is.finite(y_range))) y_range <- c(-1, 1)

  df_plot$run_id <- run_id
  df_plot$Significance <- Significance

  # plot
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = z, y = slope)) +
    # single ribbon, split into runs; fill follows Significance
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper,
                   fill = Significance, group = run_id),
      alpha = alpha, na.rm = TRUE
    ) +
    # line coloured by Significance, also split into runs for color change
    ggplot2::geom_line(
      ggplot2::aes(color = Significance, group = run_id),
      linewidth = 1, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    suppressWarnings(
      ggplot2::geom_segment(
        mapping = ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end,
                               color = hline_label, fill = hline_label),
        data = data_hline, linewidth = 1.5
      )
    ) +
    ggplot2::ggtitle("Johnson-Neyman Plot") +
    ggplot2::scale_discrete_manual(
      aesthetics = c("colour", "fill"),
      name = "",
      values = values,
      breaks = breaks,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(limits = y_range) +
    ggplot2::labs(x = z, y = paste("Slope of", x, "on", y)) +
    ggplot2::theme_minimal()

  # only show JN lines if there is a transition in the plotted window
  has_transition_in_window <- any(diff(as.integer(df_plot$significant)) != 0, na.rm = TRUE)
  if (!significant_everywhere && has_transition_in_window) {
    jn_points_in_window <- jn_points[jn_points >= z_min & jn_points <= z_max]

    if (!length(jn_points_in_window)) {
      jn_points_in_window <- approximate_jn_points_from_grid(df_plot)
      jn_points_in_window <- jn_points_in_window[jn_points_in_window >= z_min &
                                                 jn_points_in_window <= z_max]
    }

    if (length(jn_points_in_window)) {
      top_y <- suppressWarnings(max(df_plot$slope[is.finite(df_plot$slope)], na.rm = TRUE))
      if (!is.finite(top_y)) top_y <- y_range[2]

      hline_colour <- if (greyscale) "black" else "red"

      for (point in jn_points_in_window) {
        p <- p + ggplot2::geom_vline(xintercept = point, linetype = "dashed", color = hline_colour)

        if (plot.jn.points) {
          p <- p + ggplot2::annotate("text", x = point, y = top_y,
                                     label = paste("JN point:", round(point, 2)),
                                     hjust = -0.1, vjust = 1, color = "black")
        }
      }
    }
  }

  if (greyscale)
    p <- suppressMessages(p + ggplot2::scale_colour_grey() + ggplot2::scale_fill_grey())

  p
}


approximate_jn_points_from_grid <- function(grid) {
  if (is.null(grid) || !NROW(grid)) return(numeric(0))

  sig_vals <- grid$significant
  if (all(sig_vals) || !any(sig_vals)) return(numeric(0))

  transitions <- which(diff(as.integer(sig_vals)) != 0)
  if (!length(transitions)) return(numeric(0))

  approx_zero <- function(z1, z2, v1, v2) {
    if (v1 == v2) return(mean(c(z1, z2)))
    z1 - v1 * (z2 - z1) / (v2 - v1)
  }

  out <- numeric(length(transitions))

  for (i in seq_along(transitions)) {
    idx <- transitions[[i]]
    z1 <- grid$z[idx]
    z2 <- grid$z[idx + 1L]
    lower1 <- grid$lower[idx]
    lower2 <- grid$lower[idx + 1L]
    upper1 <- grid$upper[idx]
    upper2 <- grid$upper[idx + 1L]

    if ((lower1 > 0 && lower2 <= 0) || (lower1 <= 0 && lower2 > 0)) {
      out[[i]] <- approx_zero(z1, z2, lower1, lower2)
    } else if ((upper1 < 0 && upper2 >= 0) || (upper1 >= 0 && upper2 < 0)) {
      out[[i]] <- approx_zero(z1, z2, upper1, upper2)
    } else {
      out[[i]] <- mean(c(z1, z2))
    }
  }

  out[is.finite(out)]
}


getJN_PathsParTable <- function(x, y, parTable) {
  if (x == y) # exit if it's a non-recursive model
    return(NULL)

  cols <- c("lhs", "op", "rhs", "label")
  gamma <- unique(parTable[parTable$lhs == y & parTable$op == "~", cols, drop = FALSE])
  preds <- gamma[, "rhs"]

  if (!length(preds))
    return(NULL)

  paths     <- NULL
  is.direct <- NULL

  for (i in seq_along(preds)) {
    pred <- preds[[i]]

    if (pred == x) {
      paths     <- c(paths, gamma[i, "label"])
      is.direct <- c(is.direct, TRUE)

    } else {
      downstream <- getJN_PathsParTable(x = x, y = pred, parTable = parTable)

      if (!is.null(downstream)) {
        paths <- c(paths, paste0(gamma[i, "label"], "*", downstream))
        is.direct <- c(is.direct, rep(FALSE, length(downstream)))
      }
    }
  }

  attr(paths, "is.direct") <- is.direct
  paths
}
