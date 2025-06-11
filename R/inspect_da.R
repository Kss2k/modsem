inspectDA_Matrices <- c("lambda", "tau", "theta", "gamma_xi",
                        "gamma_eta", "omega_xi_xi",
                        "omega_eta_xi", "phi", "psi", "alpha")

inspectDA_Optim <- c("vcov", "FIM", "data", "coefficients",
                     "loglik", "iterations", "convergence")

modsem_inspect_da <- function(model, what = "default") {
  matrices <- model$model$matrices
  matricesCovModel <- model$model$covModel$matrices
  info <- list(N                 = NROW(model$data),
               vcov              = model$vcov,
               FIM               = model$FIM,
               data              = model$data,
               all.coefficients  = model$coefs,
               free.coefficients = model$theta,
               partable          = model$parTable,
               originalpartable  = model$originalParTable,
               loglik            = model$logLik,
               AIC               = model$AIC,
               iterations        = model$iterations,
               convergence       = model$convergence,

               lambda       = diagPartitionedMat(matrices$lambdaX,
                                                 matrices$lambdaY),
               tau          = diagPartitionedMat(matrices$tauX,
                                                 matrices$tauY),
               theta        = diagPartitionedMat(matrices$thetaDelta,
                                                 matrices$thetaEpsilon),
               gamma_xi     = diagPartitionedMat(matrices$gammaXi,
                                                 matricesCovModel$gammaXi),
               gamma_eta    = diagPartitionedMat(matrices$gammaEta,
                                                 matricesCovModel$gammaEta),
               omega_xi_xi  = diagPartitionedMat(matrices$omegaXiXi,
                                                 matricesCovModel$omegaXiXi),
               omega_eta_xi = diagPartitionedMat(matrices$omegaEtaXi,
                                                 matricesCovModel$omegaEtaXi),

               phi   = diagPartitionedMat(matrices$phi,
                                          matricesCovModel$phi),
               psi   = diagPartitionedMat(matrices$psi,
                                          matricesCovModel$psi),
               alpha = matrices$alpha)

  switch(what,
         default = info[names(info) != "data"],
         all = info,
         matrices = info[inspectDA_Matrices],
         optim = info[inspectDA_Optim],
         fit = {
           h0 <- estimate_h0(model, calc.se = FALSE)
           list(
              fit.h0 = fit_modsem_da(h0, chisq = TRUE),
              fit.h1 = fit_modsem_da(model, chisq = FALSE),
              comparative.fit = compare_fit(est_h1 = model, est_h0 = h0)
           )
         },
         info[what])
}
