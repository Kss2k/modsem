#' Fit measures for QML and LMS models
#'
#' @param model model to be assessed
#' @description Calculates chi-sq test and p-value, as well as RMSEA for
#' the LMS and QML models. Note that the Chi-Square based fit measures should be calculated
#' for the baseline model, i.e., the model without the interaction effect
#' @param model fitted model. Thereafter, you can use 'compare_fit()'
#' to assess the comparative fit of the models. If the interaction effect makes
#' the model better, and e.g., the RMSEA is good for the baseline model,
#' the interaction model likely has a good RMSEA as well.
#' @param chisq should Chi-Square based fit-measures be calculated?
#' @export
fit_modsem_da <- function(model, chisq = TRUE) {
  parTable <- model$parTable
  warnif(any(grepl(":", parTable$rhs)) && chisq,
         "Chi-Square based fit-measures for LMS and QML ",
         "should be calculated for baseline model ",
         "i.e., the model without the interaction effect")

  logLik <- model$logLik
  O      <- stats::cov(model$data)
  mu     <- apply(model$data, 2, mean)
  N      <- NROW(model$data)
  p      <- NCOL(model$data)
  coef   <- coef(model, type = "free")
  k      <- length(coef)
  df     <- getDegreesOfFreedom(m = p, coef = coef)

  matrices <- model$model$matrices
  gammaXi  <- matrices$gammaXi
  gammaEta <- matrices$gammaEta
  phi      <- matrices$phi
  psi      <- matrices$psi
  lambdaX  <- matrices$lambdaX
  lambdaY  <- matrices$lambdaY
  thetaY   <- matrices$thetaEpsilon
  thetaX   <- matrices$thetaDelta
  tauX     <- matrices$tauX
  tauY     <- matrices$tauY
  alpha    <- matrices$alpha
  Ieta     <- matrices$Ieta
  Binv     <- solve(Ieta - gammaEta)

  if (chisq) {
    covX  <- lambdaX %*% phi %*% t(lambdaX) + thetaX
    covXY <- lambdaY %*% (Binv %*% gammaXi %*% phi) %*% t(lambdaX)
    covY  <- lambdaY %*%
      (Binv %*% (gammaXi %*% phi %*% t(gammaXi) + psi) %*% t(Binv)) %*%
      t(lambdaY) + thetaY

    E <- rbind(cbind(covX, t(covXY)),
               cbind(covXY, covY))

    if (any(grepl("tau|alpha", names(coef)))) {
      muX   <- tauX
      muY   <- tauY + lambdaY %*% Binv %*% alpha
      muHat <- rbind(muX, muY)
    } else muHat <- mu

    chisqValue <- calcChiSqr(O = O, E = E, N = N, p = p, mu = mu, muHat = muHat)
    chisqP     <- stats::pchisq(chisqValue, df, lower.tail = FALSE)
    RMSEA      <- calcRMSEA(chisqValue, df, N)

  } else {
    E          <- NULL
    chisqValue <- NULL
    chisqP     <- NULL
    df         <- NULL
    muHat      <- NULL
    RMSEA      <- list(
      RMSEA          = NULL, 
      RMSEA.lower    = NULL, 
      RMSEA.upper    = NULL, 
      RMSEA.ci.level = NULL,
      RMSEA.pvalue   = NULL, 
      RMSEA.close.h0 = NULL
    )
  }

  AIC  <- calcAIC(logLik, k = k)
  AICc <- calcAdjAIC(logLik, k = k, N = N)
  BIC  <- calcBIC(logLik, k = k, N = N)
  aBIC <- calcAdjBIC(logLik, k = k, N = N)

  list(
    sigma.observed = O, 
    sigma.expected = E,
    mu.observed    = mu,
    mu.expected    = muHat,
    
    chisq.value  = chisqValue, 
    chisq.pvalue = chisqP,
    chisq.df     = df,

    AIC  = AIC,
    AICc = AICc,
    BIC  = BIC,
    aBIC = aBIC, 

    RMSEA          = RMSEA$rmsea,
    RMSEA.lower    = RMSEA$lower,
    RMSEA.upper    = RMSEA$upper,
    RMSEA.ci.level = RMSEA$ci.level,
    RMSEA.pvalue   = RMSEA$pvalue, 
    RMSEA.close.h0 = RMSEA$close.h0
  )
}


calcChiSqr <- function(O, E, N, p, mu, muHat) {
  diff_mu <- mu - muHat
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
