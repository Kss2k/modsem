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
#' @param var.interaction if FALSE (default) variances for interaction terms will be removed (if present)
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
                              intercepts = !standardized,
                              variances = TRUE,
                              var.interaction = FALSE,
                              ...) {
  method <- object$method

  if (standardized) {
    parTable <- standardized_estimates(object, intercepts = intercepts, 
                                       monte.carlo = monte.carlo, mc.reps = mc.reps)
  } else {
    parTable <- parameter_estimates(object)
  }

  if (!var.interaction) {
    parTable <- removeInteractionVariances(parTable)
  }

  args <- object$args
  out <- list(
    parTable       = parTable,
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
    if (any(grepl(":", parTable$rhs))) cat("Estimating baseline model (H0)\n")
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
    out$r.squared <- calcRsquared(parTable)
    if (H0) out$r.squared$H0 <- calcRsquared(est_h0$parTable)
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
    variances     = variances
  )

  class(out) <- "summary_da"
  out
}


#' @export
print.summary_da <- function(x, digits = 3, ...) {
  width.out <- getWidthPrintedParTable(x$parTable,
                                       scientific = x$format$scientific,
                                       ci = x$format$ci,
                                       digits = x$format$digits,
                                       loadings = x$format$loadings,
                                       regressions = x$format$regressions,
                                       covariances = x$format$covariances,
                                       intercepts = x$format$intercepts,
                                       variances = x$format$variances)
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
    values <- c(round(x$nullModel$logLik), round(x$fitH0$AIC, 2), round(x$fitH0$BIC, 2))

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
    values <- c(formatNumeric(x$D$llChange, digits = 2),
                formatNumeric(x$D$D, digits = 2),
                x$D$df,
                formatPval(x$D$p, scientific = x$format$scientific))
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                     width.out = width.out), "\n")

  }

  if (!is.null(x$r.squared)) {
    r.squared <- x$r.squared
    r.squared$Rsqr <- formatNumeric(r.squared$Rsqr, digits = 3)
    maxWidth <- max(vapply(r.squared$Rsqr, FUN.VALUE = numeric(1), FUN = nchar))
    r.squared$Rsqr <-
      stringr::str_pad(r.squared$Rsqr, width = maxWidth, side = "left")

    cat("R-Squared Interaction Model (H1):\n")
    names <- r.squared$eta
    values <- character(length(names))
    for (i in seq_along(r.squared$eta)) {
      values[[i]] <- r.squared$Rsqr[[i]]
    }
    cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                     width.out = width.out))

    if (!is.null(r.squared$H0)) {
      r.squared$H0$Rsqr <- formatNumeric(r.squared$H0$Rsqr, digits = 3)
      maxWidth <-
        max(vapply(r.squared$H0$Rsqr, FUN.VALUE = numeric(1), FUN = nchar))
      r.squared$H0$Rsqr <-
        stringr::str_pad(r.squared$H0$Rsqr, width = maxWidth, side = "left")

      cat("R-Squared Baseline Model (H0):\n")
      names <- r.squared$H0$eta
      for (i in seq_along(names)) {
          values[[i]] <- formatNumeric(r.squared$H0$Rsqr[[i]], digits = 3)
      }
      cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                       width.out = width.out))

      # Calculate Change (using unformatted Rsquared)
      r.squared$H0$diff <-
        formatNumeric(x$r.squared$Rsqr - x$r.squared$H0$Rsqr, digits = 3)
      cat("R-Squared Change (H1 - H0):\n")
      for (i in seq_along(names)) {
        values[[i]] <- r.squared$H0$diff[[i]]
      }
      cat(allignLhsRhs(lhs = names, rhs = values, pad = "  ",
                       width.out = width.out))
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
                variances   = x$format$variances)
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


calcRsquared <- function(parTable) {
  parTable <- var_interactions.data.frame(parTable)
  etas     <- unique(parTable$lhs[parTable$op == "~"])

  # Calculate Variances/R squared of Etas
  residuals_df <- parTable[parTable$lhs %in% etas & 
                           parTable$rhs == parTable$lhs, ]
  residuals <- structure(residuals_df$est, names=residuals_df$lhs)[etas]
  variances <- calcVarParTable(etas, parTable)
  Rsqr <- 1 - residuals / variances

  data.frame(eta = etas, variance = variances,
             residual = residuals, Rsqr = Rsqr)
}


#' @export
var_interactions.modsem_da <- function(object, ...) {
  var_interactions.data.frame(parameter_estimates(object), ...)
}


#' @export
modsem_inspect.modsem_da <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "default"
  modsem_inspect_da(object, what = what, ...)
}


#' @export
#' @importFrom stats vcov
vcov.modsem_da <- function(object, ...) {
  modsem_inspect_da(object, what = "vcov")[[1]]
}


#' @export
#' @importFrom stats coefficients
coefficients.modsem_da <- function(object, type = "all", ...) {
  what <- ifelse(type == "all", yes = "all.coefficients",
                 no = "free.coefficients")
  modsem_inspect_da(object, what = what)[[1]]
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


#' Get standardized estimates with Monte Carlo bootstrapped standard errors
#'
#' @param object An object of class \code{modsem_da}
#' @param monte.carlo Should standard errors be calculated using Monte Carlo simulation. 
#'   if \code{FALSE} the delta method is used instead. Default is \code{FALSE}.
#' @param mc.reps Number of Monte Carlo repetitions
#' @param tolerance.zero Tolerance for zero values. Standard errors smaller than this 
#'   value will be set to NA.
#' @param ... Additional arguments passed to other functions
#' @details The interaction term is not standardized such that \code{var(xz) = 1}.
#' The interaction term is not an actual variable in the model, meaning that it does not
#' have a variance. It must therefore be calculated from the other parameters in the model.
#' Assuming normality and zero-means, the variance is calculated as
#' \code{var(xz) = var(x) * var(z) + cov(x, z)^2}. Thus setting the variance of the interaction
#' term to 1 would only be 'correct' if the correlation between \code{x} and \code{z} is zero.
#' This means that the standardized estimates for the interaction term will
#' be different from those using \code{lavaan}, since there the interaction term is an
#' actual latent variable in the model, with a standardized variance of 1.
#' @export
standardized_estimates.modsem_da <- function(object, monte.carlo = FALSE, mc.reps = 10000, tolerance.zero = 1e-10, ...) {
  stdSolution <- standardized_solution_COEFS(
    object, monte.carlo = monte.carlo, mc.reps = mc.reps, 
    tolerance.zero = tolerance.zero, ...
  )
  
  stdSolution$parTable
}


standardized_solution_COEFS <- function(object, monte.carlo = FALSE, mc.reps = 10000, tolerance.zero = 1e-10, seed = 123, 
                                        delta.epsilon = 1e-8, ...) {
  set.seed(seed)

  stopif(!inherits(object, "modsem_da"), "The object must be of class 'modsem_da'.")

  parTable <- parameter_estimates(object)
  parTable <- parTable[c("lhs", "op", "rhs", "label", "est", "std.error")]
  parTable <- centerInteraction(parTable) # re-estimate path-coefficients 
                                          # when intercepts are zero
  
  lVs      <- getLVs(parTable)
  intTerms <- getIntTerms(parTable)
  etas     <- getSortedEtas(parTable, isLV = TRUE)
  xis      <- getXis(parTable, etas = etas, isLV = TRUE)
  indsLVs  <- getIndsLVs(parTable, lVs)
  allInds  <- unique(unlist(indsLVs))

  originalLabels <- parTable$label
  labels <- getParTableLabels(parTable, labelCol="label")
  parTable$label <- labels

  V <- vcov(object)
  coefs <- structure(parTable$est, names = labels)

  # parTable$est <- NULL # no longer needed
  parTable$std.error <- NA

  isIntercept <- parTable$op == "~1"
  labels <- unique(labels[!isIntercept]) # remove intercept labels
  parTable <- parTable[!isIntercept, ]

  V <- expandVCOV(V, labels=labels)
  coefs <- coefs[labels]

  legalNames <- stringr::str_replace_all(labels, OP_REPLACEMENTS)

  if (monte.carlo) {
    COEFS <- as.data.frame(mvtnorm::rmvnorm(mc.reps, mean = coefs, sigma = V))
  } else { # delta method
    k <- length(coefs)
    COEFS <- matrix(coefs, nrow=k, ncol=k, byrow=TRUE)
    DeltaEpsilon <- diag(delta.epsilon, k)
    isFixed <- abs(diag(V)) < tolerance.zero
    DeltaEpsilon[isFixed, isFixed] <- 0
    COEFS_p <- as.data.frame(COEFS + DeltaEpsilon)
    COEFS_m <- as.data.frame(COEFS - DeltaEpsilon)
    COEFS <- rbind(COEFS_p, COEFS_m)
  }

  # Get legal parameter names
  names(COEFS) <- legalNames
  names(coefs) <- legalNames
  COEFS <- rbind(as.data.frame(as.list(coefs)), COEFS) # first row is the original values

  colnames(V) <- legalNames
  rownames(V) <- legalNames
  parTable$label <- stringr::str_replace_all(parTable$label, OP_REPLACEMENTS)
  
  parTable <- parTable[c("lhs", "op", "rhs", "label")]
  COEFS <- var_interactions_COEFS(parTable, COEFS) # calculate variances of interaction terms
  # Calculate 

  # get variances
  varianceEquations <- list()
  variances <- list()

  for (x in allInds) {
    # get the variance equation for each variable
    eqVarX <- getCovEqExpr(x=x, y=x, parTable=parTable, measurement.model=TRUE)
    varianceEquations[[x]] <- eqVarX  
    variances[[x]] <- eval(eqVarX, envir = COEFS)
  }
  
  for (x in lVs) {
    # get the variance equation for each variable
    eqVarX <- getCovEqExpr(x=x, y=x, parTable=parTable)
    varianceEquations[[x]] <- eqVarX  
    variances[[x]] <- eval(eqVarX, envir = COEFS)
  }

  for (xz in intTerms) {
    labelVarXZ <- getLabelVarXZ(xz)
    eqVarXZ <- parse(text=labelVarXZ)
    varianceEquations[[xz]] <- eqVarXZ
    variances[[xz]] <- eval(eqVarXZ, envir = COEFS)
  }

  # Factor Loadings
  lambda     <- NULL
  selectRows <- NULL

  for (lV in lVs) {
    for (ind in indsLVs[[lV]]) {
      selectRows  <- parTable$lhs == lV & parTable$op == "=~" & parTable$rhs == ind
      label <- parTable[selectRows, "label"]

      # est in parTable
      scalingCoef <- sqrt(variances[[lV]]) / sqrt(variances[[ind]])
      lambda      <- COEFS[[label]] * scalingCoef

      COEFS[[label]] <- lambda
    }
  }

  # Structural Coefficients
  gamma               <- NULL
  selectStrucExprsEta <- NULL
  structExprsEta      <- NULL
  selectStrucExprs    <- parTable$op == "~" & parTable$lhs %in% etas

  for (eta in etas) {
    selectStrucExprsEta <- selectStrucExprs & parTable$lhs == eta
    structExprsEta      <- parTable[selectStrucExprsEta, ]

    for (xi in structExprsEta$rhs) {
      selectRows  <- selectStrucExprsEta & parTable$rhs == xi
      scalingCoef <- sqrt(variances[[xi]]) / sqrt(variances[[eta]])
      label       <- parTable[selectRows, "label"]
      gamma       <- COEFS[[label]] * scalingCoef
      
      COEFS[[label]] <- gamma
    }
  }

  # (Co-) Variances of xis
  selectCovXis <- parTable$op == "~~" & parTable$lhs %in% xis
  selectRows   <- NULL
  combosXis    <- getUniqueCombos(xis, match = TRUE)

  for (i in seq_len(nrow(combosXis))) {
    xis         <- combosXis[i, , drop = TRUE]
    selectRows  <- selectCovXis & parTable$lhs %in% xis & parTable$rhs %in% xis
    scalingCoef <- sqrt(variances[[xis[[1]]]]) * sqrt(variances[[xis[[2]]]])

    if (xis[[1]] != xis[[2]]) {
      selectRows <- selectRows & parTable$lhs != parTable$rhs
    }

    label <- parTable[selectRows, "label"]
    covs <- COEFS[[label]] / scalingCoef

    COEFS[[label]] <- covs
  }

  # Residual Variances etas
  selectRows <- NULL
  residual   <- NULL

  for (eta in etas) {
    selectRows <- parTable$lhs == eta & parTable$op == "~~" & parTable$rhs == eta
    label <- parTable[selectRows, "label"]
    residual <- COEFS[[label]] / variances[[eta]]

    COEFS[[label]] <- residual
  }

  # residual variances inds
  for (ind in allInds) {
    selectRows <- parTable$lhs == ind & parTable$op == "~~" & parTable$rhs == ind
    label <- parTable[selectRows, "label"]
    residual <- COEFS[[label]] / variances[[ind]]

    COEFS[[label]] <- residual
  }

  # recalculate variance of interaction terms
  # and rescale coefficients for interaction terms
  COEFS <- var_interactions_COEFS(parTable, COEFS)

  for (xz in intTerms) {
    selectRows <- parTable$rhs == xz & parTable$op == "~"
    label <- parTable[selectRows, "label"]
    labelVarXZ <- getLabelVarXZ(xz)
    gamma <- COEFS[[label]] / sqrt(COEFS[[labelVarXZ]])

    COEFS[[label]] <- gamma 
  }

  # recalculate custom parameters
  constrExprs <- sortConstrExprsFinalPt(parTable)
  for (i in seq_len(NROW(constrExprs))) {
    row <- constrExprs[i, , drop=FALSE]
    label <- row$label

    expr      <- parse(text=constrExprs[i, "rhs"])
  
    newVals <- eval(expr, envir = COEFS)
  
    COEFS[[label]] <- newVals
  }

  nzeros <- round(log10(1 / tolerance.zero), 0)
  coefs <- unlist(COEFS[1, ])
  COEFS <- round(COEFS[2:(nrow(COEFS)), ],  nzeros) # skip first row

  if (monte.carlo) {
    vcov <- stats::cov(COEFS)
  } else {
    # delta method
    p <- 1:k
    m <- (k+1):(k + k)

    COEFS_p <- t(COEFS[p, ])
    COEFS_m <- t(COEFS[m, ])

    J <- (COEFS_p - COEFS_m) / (2 * delta.epsilon)
    vcov <- J %*% V %*% t(J)
  }

  # fill parTable
  std.errors <- suppressWarnings(sqrt(diag(vcov)))
  warnFunc <- function(type, row) {
    warning2("Unable to calculate standardized ", type, " for: ", 
             paste0(row$lhs, row$op, row$rhs))
  }

  verboseLabels <- stringr::str_replace_all(parTable$label, OP_REPLACEMENTS)
  for (i in seq_len(nrow(parTable))) {
    row <- parTable[i, , drop = FALSE]
    label <- verboseLabels[[i]]

    if (!label %in% names(coefs)) {
      warnFunc("coefficient", row)
      est <- NA
    } else est <- coefs[[label]]

    if (!label %in% names(std.errors)) {
      warnFunc("std.error", row)
      se <- NA 
    } else se <- std.errors[[label]]

    parTable[i, "est"] <- est
    parTable[i, "std.error"] <- se
  }

  # Fill in NA on zero-standard errors, and create vcov and coefs
  parTable[!is.na(parTable$std.error) & 
           abs(parTable$std.error) < tolerance.zero, "std.error"] <- NA

  isFree <- !is.na(parTable$std.error)
  coefs <- structure(parTable$est[isFree], names = parTable$label[isFree])
  params <- intersect(names(coefs), rownames(vcov))

  coefs <- coefs[params]
  vcov <- vcov[params, params]

  cleanedParamLabels <- stringr::str_replace_all(params, OP_REPLACEMENTS_INV)
  names(coefs) <- cleanedParamLabels 
  colnames(vcov) <- rownames(vcov) <- cleanedParamLabels

  # Finalize parTable
  parTable[!parTable$label %in% originalLabels, "label"] <- "" 
  parTable$z.value  <- parTable$est / parTable$std.error
  parTable$p.value  <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - CI_WIDTH * parTable$std.error
  parTable$ci.upper <- parTable$est + CI_WIDTH * parTable$std.error


  list(parTable = parTable,
       coefs    = coefs, 
       vcov     = vcov,
       COEFS    = COEFS)
}
