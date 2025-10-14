#' Fit measures for QML and LMS models
#'
#' @param model model to be assessed
#' @description Calculates chi-sq test and p-value, as well as RMSEA for
#' the LMS and QML models. Note that the Chi-Square based fit measures should be calculated
#' for the baseline model, i.e., the model without the interaction effect
#'
#' @param model fitted model. Thereafter, you can use 'compare_fit()'
#' to assess the comparative fit of the models. If the interaction effect makes
#' the model better, and e.g., the RMSEA is good for the baseline model,
#' the interaction model likely has a good RMSEA as well.
#'
#' @param chisq should Chi-Square based fit-measures be calculated?
#'
#' @param lav.fit Should fit indices from the \code{lavaan} model used to optimize
#'   the starting parameters be included (if available)? This is usually only approprioate
#'   for linear models (i.e., no interaction effects), where the parameter estimates
#'   for LMS and QML are equivalent to ML estimates from lavaan.
#'
#' @export
fit_modsem_da <- function(model, chisq = TRUE, lav.fit = FALSE) {
  parTable <- model$parTable
  if (isTRUE(chisq) && any(grepl(":", parTable$rhs), na.rm = TRUE)) {
    warning2("Chi-Square based fit-measures for LMS and QML ",
             "should be calculated for baseline model ",
             "i.e., the model without the interaction effect",
             immediate. = FALSE)
  }

  t      <- nFreeInterceptsDA(model)
  mean.s <- model$args$mean.observed || t > 0
  logLik <- model$logLik
  coef   <- coef(model, type = "free")
  k      <- length(coef)

  finalModel <- model$model
  submodels  <- if (!is.null(finalModel$models)) finalModel$models else list()
  if (!length(submodels)) submodels <- list(list(data = model$data))

  n.groups <- max(1L, length(submodels))
  group.labels <- finalModel$info$group.levels
  if (is.null(group.labels) || length(group.labels) != n.groups) {
    group.labels <- paste0("Group", seq_len(n.groups))
  }

  expected.list <- model$expected.matrices
  if (!is.list(expected.list) || length(expected.list) == 0) {
    expected.list <- replicate(n.groups, expected.list, simplify = FALSE)
  }
  if (length(expected.list) < n.groups) {
    expected.list <- c(expected.list,
                       rep(list(NULL), n.groups - length(expected.list)))
  }
  expected.list <- expected.list[seq_len(n.groups)]

  sample_cov_blocks <- vector("list", n.groups)
  expected_cov_blocks <- vector("list", n.groups)
  mu_obs_blocks <- vector("list", n.groups)
  mu_exp_blocks <- vector("list", n.groups)

  N_vec <- numeric(n.groups)
  p_vec <- numeric(n.groups)
  chisq_parts <- rep(NA_real_, n.groups)
  srmr_weight_num <- 0
  srmr_weight_den <- 0

  for (g in seq_len(n.groups)) {
    submodel <- submodels[[g]]
    data_g   <- submodel$data

    data_full <- data_g$data.full
    if (is.null(data_full) || !NROW(data_full) || !NCOL(data_full)) {
      next
    }

    data_full <- as.matrix(data_full)
    N_g <- NROW(data_full)
    p_g <- NCOL(data_full)

    N_vec[g] <- N_g
    p_vec[g] <- p_g

    mu_g_vec <- apply(data_full, 2, mean, na.rm = TRUE)
    mu_g <- matrix(mu_g_vec, ncol = 1,
                   dimnames = list(colnames(data_full), "~1"))
    O_g <- stats::cov(data_full, use = "pairwise.complete.obs")

    expected_g <- expected.list[[g]]
    E_g <- if (!is.null(expected_g) && !is.null(expected_g$sigma.ov)) expected_g$sigma.ov else NULL
    mu_hat_g <- if (!is.null(expected_g) && !is.null(expected_g$mu.ov)) expected_g$mu.ov else NULL

    if (!mean.s || is.null(mu_hat_g)) {
      mu_hat_g <- mu_g
    }

    if (!is.null(E_g)) {
      order_vars <- colnames(data_full)
      E_g <- E_g[order_vars, order_vars, drop = FALSE]
      mu_hat_g <- mu_hat_g[order_vars, , drop = FALSE]
    }

    prefix <- group.labels[[g]]
    sample_cov_blocks[[g]] <- modsemMatrix(O_g, symmetric = TRUE)
    mu_obs_blocks[[g]] <- modsemMatrix(mu_g)

    if (!is.null(E_g)) {
      expected_cov_blocks[[g]] <- modsemMatrix(E_g, symmetric = TRUE)
      mu_exp_blocks[[g]] <- modsemMatrix(mu_hat_g)
    } else {
      expected_cov_blocks[[g]] <- NULL
      mu_exp_blocks[[g]] <- NULL
    }

    if (chisq && !is.null(E_g)) {
      chi_g <- tryCatch(
        calcChiSqr(O = O_g, E = E_g, N = N_g, p = p_g, mu = mu_g, mu.hat = mu_hat_g),
        error = function(e) {
          warning2("Failed to compute chi-square contribution for group '",
                   prefix, "': ", conditionMessage(e), immediate. = FALSE)
          NA_real_
        }
      )
      chisq_parts[g] <- chi_g

      if (!is.na(chi_g)) {
        srmr_g <- tryCatch(
          calcSRMR_Mplus(S = O_g, M = mu_g, Sigma.hat = E_g, Mu.hat = mu_hat_g,
                         mean.structure = mean.s),
          error = function(e) NA_real_
        )

        if (!is.na(srmr_g)) {
          weight <- max(N_g, 1)
          srmr_weight_num <- srmr_weight_num + weight * srmr_g^2
          srmr_weight_den <- srmr_weight_den + weight
        }
      }
    }
  }

  names(sample_cov_blocks) <- group.labels
  names(expected_cov_blocks) <- group.labels
  names(mu_obs_blocks) <- group.labels
  names(mu_exp_blocks) <- group.labels

  N_total <- sum(N_vec)
  df <- getDegreesOfFreedom(p = p_vec[p_vec > 0], coef = coef, mean.structure = mean.s)

  if (chisq) {
    if (all(!is.na(chisq_parts))) {
      chisqValue <- sum(chisq_parts)
      chisqP     <- if (!is.na(df) && df > 0) stats::pchisq(chisqValue, df, lower.tail = FALSE) else NA_real_
      RMSEA      <- if (!is.na(df) && df > 0 && N_total > 0) {
        calcRMSEA(chisqValue, df, N_total)
      } else {
        list(rmsea = NA_real_, lower = NA_real_, upper = NA_real_,
             ci.level = NA_real_, pvalue = NA_real_, close.h0 = NA_real_)
      }
      SRMR       <- if (srmr_weight_den > 0) sqrt(srmr_weight_num / srmr_weight_den) else NA_real_
    } else {
      chisqValue <- NA_real_
      chisqP     <- NA_real_
      RMSEA      <- list(rmsea = NA_real_, lower = NA_real_, upper = NA_real_,
                         ci.level = NA_real_, pvalue = NA_real_, close.h0 = NA_real_)
      SRMR       <- NA_real_
    }
  } else {
    chisqValue <- NULL
    chisqP     <- NULL
    df         <- NULL
    sigma_expected <- NULL
    mu_expected <- NULL
    SRMR       <- NULL
    RMSEA      <- list(
      rmsea = NULL,
      lower = NULL,
      upper = NULL,
      ci.level = NULL,
      pvalue = NULL,
      close.h0 = NULL
    )
  }

  AIC  <- calcAIC(logLik, k = k)
  AICc <- calcAdjAIC(logLik, k = k, N = N_total)
  BIC  <- calcBIC(logLik, k = k, N = N_total)
  aBIC <- calcAdjBIC(logLik, k = k, N = N_total)

  if (lav.fit) {
    lavfit <- tryCatch(
      lavaan::fitMeasures(model$model$lavaan.fit),
      error = function(e) {
        warning2("Unable to retrieve fit measures for lavaan model!\n",
                 "Message: ", conditionMessage(e), immediate. = FALSE)
        NULL
      }
    )
  } else lavfit <- NULL


  list(
    lav.fit = lavfit,

    sigma.observed = sample_cov_blocks,
    sigma.expected = expected_cov_blocks,
    mu.observed    = mu_obs_blocks,
    mu.expected    = mu_exp_blocks,

    chisq.value  = chisqValue,
    chisq.pvalue = chisqP,
    chisq.df     = df,
    chisq.group  = stats::setNames(chisq_parts, group.labels),

    AIC  = AIC,
    AICc = AICc,
    BIC  = BIC,
    aBIC = aBIC,
    SRMR = SRMR,

    RMSEA          = RMSEA$rmsea,
    RMSEA.lower    = RMSEA$lower,
    RMSEA.upper    = RMSEA$upper,
    RMSEA.ci.level = RMSEA$ci.level,
    RMSEA.pvalue   = RMSEA$pvalue,
    RMSEA.close.h0 = RMSEA$close.h0
  )
}


calcChiSqr <- function(O, E, N, p, mu, mu.hat) {
  diff_mu <- mu - mu.hat
  Einv    <- solve(E)
  as.vector(
    (N - 1) * (t(diff_mu) %*% Einv %*% diff_mu +
               tr(O %*% Einv) - log(det(O %*% Einv)) - p)
  )
}


tryCatchUniroot <- function(f, lower, upper, errorVal = NA) {
  tryCatch(
    stats::uniroot(f, lower = lower, upper = upper)$root,
    error = function(e) errorVal
  )
}


calcRMSEA <- function(chi.sq, df, N, ci.level = 0.90, close.h0=0.05) {
  alpha  <- 1 - ci.level

  fLower <- \(lambda) stats::pchisq(chi.sq, df, ncp=lambda) - (1 - alpha/2)
  fUpper <- \(lambda) stats::pchisq(chi.sq, df, ncp=lambda) - (    alpha/2)
  fRMSEA <- \(lambda) sqrt(max(lambda, 0) / (df * (N - 1)))

  point <- chi.sq - df
  lower <- tryCatchUniroot(fLower, lower=0, upper=chi.sq, errorVal=0)
  upper <- tryCatchUniroot(fUpper, lower=0, upper=10*chi.sq, errorVal=df * (N - 1)) # i.e., RMSEA.upper = 1

  rmseaLower  <- fRMSEA(lower)
  rmseaUpper  <- fRMSEA(upper)
  rmseaHat    <- fRMSEA(point)
  rmseaPvalue <- 1 - stats::pchisq(chi.sq, df=df, ncp=df*(N-1)*close.h0^2)

  list(rmsea = rmseaHat, lower = rmseaLower, upper = rmseaUpper,
       ci.level = ci.level, pvalue = rmseaPvalue, close.h0 = close.h0)
}


calcAIC <- function(logLik, k) {
  2 * k - 2 * logLik
}


calcAdjAIC <- function(logLik, k, N) {
  AIC <- calcAIC(logLik, k)
  AIC + (2 * k ^ 2 + 2 * k) / (N - k - 1)
}


calcBIC <- function(logLik, k, N) {
  log(N) * k - 2 * logLik
}


calcAdjBIC <- function(logLik, k, N) {
  log((N + 2) / 24) * k - 2 * logLik
}


calcSRMR_Mplus <- function(S, M, Sigma.hat, Mu.hat, mean.structure = TRUE) {
  # Bollen approach: simply using cov2cor ('correlation residuals')
  S.cor <- cov2cor(S)
  Sigma.cor <- cov2cor(Sigma.hat)
  R.cor <- (S.cor - Sigma.cor)
  nvar  <- NCOL(Sigma.cor)

  # meanstructure
  if (mean.structure) {
    # standardized residual mean vector
    R.cor.mean <- M / sqrt(diag(S)) - Mu.hat / sqrt(diag(Sigma.hat))

    e <- nvar * (nvar + 1) / 2 + nvar
    srmr.mplus <-
      sqrt((sum(R.cor[lower.tri(R.cor, diag = FALSE)]^2) +
            sum(R.cor.mean^2) +
            sum(((diag(S) - diag(Sigma.hat)) / diag(S))^2)) / e)

    e <- nvar * (nvar + 1) / 2
    srmr.mplus.nomean <-
      sqrt((sum(R.cor[lower.tri(R.cor, diag = FALSE)]^2) +
            sum(((diag(S) - diag(Sigma.hat)) / diag(S))^2)) / e)
  } else {
    e <- nvar * (nvar + 1) / 2
    srmr.mplus.nomean <- srmr.mplus <-
      sqrt((sum(R.cor[lower.tri(R.cor, diag = FALSE)]^2) +
            sum(((diag(S) - diag(Sigma.hat)) / diag(S))^2)) / e)
  }

  attr(srmr.mplus, "nomean") <- srmr.mplus.nomean
  srmr.mplus
}
