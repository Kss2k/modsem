inspectDA_Matrices <- c("lambda", "tau", "theta", "gamma.xi",
                        "gamma.eta", "omega.xi.xi",
                        "omega.eta.xi", "phi", "psi", "alpha", "beta0")


inspectDA_Optim <- c("coefficients.free", "vcov.free", "information", 
                     "loglik", "iterations", "convergence")


modsem_inspect_da <- function(model, what = "default") {
  stopif(!length(what), "`what` is of length zero!")

  matrices          <- model$model$matrices
  matricesCovModel  <- model$model$covModel$matrices
  expected.matrices <- model$expected.matrices

  lambda       <- diagPartitionedMat(matrices$lambdaX,
                                     matrices$lambdaY)
  tau          <- diagPartitionedMat(matrices$tauX,
                                     matrices$tauY)
  theta        <- diagPartitionedMat(matrices$thetaDelta,
                                     matrices$thetaEpsilon)
  gamma.xi     <- diagPartitionedMat(matrices$gammaXi,
                                     matricesCovModel$gammaXi)
  gamma.eta    <- diagPartitionedMat(matrices$gammaEta,
                                     matricesCovModel$gammaEta)
  omega.xi.xi  <- diagPartitionedMat(matrices$omegaXiXi,
                                     matricesCovModel$omegaXiXi)
  omega.eta.xi <- diagPartitionedMat(matrices$omegaEtaXi,
                                     matricesCovModel$omegaEtaXi)
  phi          <- diagPartitionedMat(matrices$phi,
                                     matricesCovModel$phi)
  psi          <- diagPartitionedMat(matrices$psi,
                                     matricesCovModel$psi)

  cov.ov  <- expected.matrices$sigma.ov
  cov.lv  <- expected.matrices$sigma.lv
  cov.all <- expected.matrices$sigma.all

  cor.ov  <- cov2cor(cov.ov)
  cor.lv  <- cov2cor(cov.lv)
  cor.all <- cov2cor(cov.all)

  info <- list(N                 = NROW(model$data),
               vcov.all          = modsemMatrix(model$vcov.all, symmetric = TRUE),
               vcov.free         = modsemMatrix(model$vcov.free, symmetric = TRUE),
               information       = modsemMatrix(model$FIM, symmetric = TRUE),
               data              = model$data,
               coefficients.all  = modsemVector(model$coefs.all),
               coefficients.free = modsemVector(model$coefs.free),
               partable          = modsemParTable(model$parTable),
               partable.input    = model$originalParTable,
               loglik            = model$logLik,
               iterations        = model$iterations,
               convergence       = model$convergence,

               lambda       = modsemMatrix(lambda), 
               tau          = modsemMatrix(tau), 
               theta        = modsemMatrix(theta, symmetric = TRUE),
               gamma.xi     = modsemMatrix(gamma.xi), 
               gamma.eta    = modsemMatrix(gamma.eta), 
               omega.xi.xi  = modsemMatrix(omega.xi.xi), 
               omega.eta.xi = modsemMatrix(omega.eta.xi), 

               phi   = modsemMatrix(phi, symmetric = TRUE), 
               psi   = modsemMatrix(psi, symmetric = TRUE), 

               alpha = modsemMatrix(matrices$alpha),
               beta0 = modsemMatrix(matrices$beta0),

               cov.ov  = modsemMatrix(cov.ov, symmetric = TRUE),
               cov.lv  = modsemMatrix(cov.lv, symmetric = TRUE),
               cov.all = modsemMatrix(cov.all, symmetric = TRUE),

               cor.ov  = modsemMatrix(cor.ov, symmetric = TRUE),
               cor.lv  = modsemMatrix(cor.lv, symmetric = TRUE),
               cor.all = modsemMatrix(cor.all, symmetric = TRUE),

               mean.lv  = modsemMatrix(expected.matrices$mu.lv),
               mean.ov  = modsemMatrix(expected.matrices$mu.ov),
               mean.all = modsemMatrix(expected.matrices$mu.all),

               r2.all  = modsemVector(expected.matrices$r2.all),
               r2.lv   = modsemVector(expected.matrices$r2.lv),
               r2.ov   = modsemVector(expected.matrices$r2.ov),

               res.all = modsemVector(expected.matrices$res.all),
               res.lv  = modsemVector(expected.matrices$res.lv),
               res.ov  = modsemVector(expected.matrices$res.ov)
  )

  FIT <- \() {
    h0 <- estimate_h0(model, calc.se = FALSE)
    list(
      fit.h0 = fit_modsem_da(h0, chisq = TRUE),
      fit.h1 = fit_modsem_da(model, chisq = FALSE),
      comparative.fit = compare_fit(est_h1 = model, est_h0 = h0)
    )
  }

  if (length(what) > 1) {
    fields <- info[what]

  } else {
    fields <- switch(
      EXPR     = what,
      coef     = info[c("vcov.all", "coefficients.all")],
      default  = info[names(info) != "data"],
      all      = info,
      matrices = info[inspectDA_Matrices],
      optim    = info[inspectDA_Optim],
      fit      = FIT(),
      info[[what]]
    )
  }

  nullvalues <- vapply(fields, FUN.VALUE = logical(1L), FUN = is.null)

  warnif(any(nullvalues), "Some fields in `modsem_inspect()` could not be retrieved!",
         immediate. = FALSE)

  fields 
}
