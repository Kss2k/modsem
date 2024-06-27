#' Fit measures for QML and LMS models
#'
#' @param model model to be assessed
#' @description calculates chi-sq test and p-value, as well as RMSEA for 
#' the LMS and QML models. Note that the fit measures should be calculated
#' for the baseline model, i.e., the model without the interaction effect 
#' @param model fitted model. Thereafter, you can use 'compare_fit()' 
#' to assess the comparative fit of the models. If the interaction effect makes 
#' the model better, and e.g., the RMSEA is good for the baseline model, 
#' the interaction model likely has a good RMSEA as well.
#' @export
fit_modsem_da <- function(model) {
  parTable <- model$parTable
  if (any(grepl(":", parTable$rhs))) {
    warning2("Fit measures for LMS and QML should be calculated for baseline model ", 
             "i.e., the model without the interaction effect")
    parTable <- centerInteraction(parTable)
  }
 
  logLik <- model$logLik
  O <- stats::cov(model$data)
  mu <- apply(model$data, 2, mean)
  N <- NROW(model$data) 
  p <- NCOL(model$data)
  coef <- coef(model)
  df <- getDegreesOfFreedom(m = p, coef = coef)

  matrices <- model$model$matrices
  gammaXi <- matrices$gammaXi
  gammaEta <- matrices$gammaEta
  phi <- matrices$phi
  psi <- matrices$psi
  lambdaX <- matrices$lambdaX
  lambdaY <- matrices$lambdaY
  thetaY <- matrices$thetaEpsilon
  thetaX <- matrices$thetaDelta
  tauX <- matrices$tauX
  tauY <- matrices$tauY 
  alpha <- matrices$alpha
  Ieta <- matrices$Ieta
  Binv <- solve(Ieta - gammaEta)
  
  covX <- lambdaX %*% phi %*% t(lambdaX) + thetaX
  covXY <- lambdaY %*% (Binv %*% gammaXi %*% phi) %*% t(lambdaX) 
  covY <- lambdaY %*% 
    (Binv %*% (gammaXi %*% phi %*% t(gammaXi) + psi) %*% t(Binv)) %*% 
    t(lambdaY) + thetaY

  E <- rbind(cbind(covX, t(covXY)),
              cbind(covXY, covY))

  if (any(grepl("tau|alpha", names(coef)))) {
    muX <- tauX
    muY <- tauY + lambdaY %*% Binv %*% alpha
    muHat <- rbind(muX, muY)
  } else muHat <- mu # diff_mu = 0

  chisqValue <- calcChiSqr(O = O, E = E, N = N, p = p, mu = mu, muHat = muHat)
  chisqP <- stats::pchisq(chisqValue, df, lower.tail = FALSE)
  RMSEA <- calcRMSEA(chisqValue, df, N)
  AIC <- calcAIC(logLik, df)

  list(sigma.observed = O, sigma.expected = E, 
       mu.observed = mu, mu.expected = muHat,
       chisq.value = chisqValue, chisq.pvalue = chisqP, chisq.df = df, 
       AIC = AIC, RMSEA = RMSEA)
}


calcChiSqr <- function(O, E, N, p, mu, muHat) {
  diff_mu <- mu - muHat
  Einv <- solve(E)
  (N - 1) * (t(diff_mu) %*% Einv %*% diff_mu + tr(O %*% Einv) - log(det(O %*% Einv)) - p)
}


calcRMSEA <- function(chi.sq, df, N) {
  ncp <- max(0, chi.sq - df)
  sqrt(ncp / ((N - df) * (N - df)))
}


calcAIC <- function(logLik, df) {
  - 2 * logLik + 2 * df
}
