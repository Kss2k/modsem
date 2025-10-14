inspectDA_Matrices <- c("lambda", "tau", "theta", "gamma.xi",
                        "gamma.eta", "omega.xi.xi",
                        "omega.eta.xi", "phi", "psi", "alpha", "beta0")


inspectDA_Optim <- c("coefficients.free", "vcov.free", "information",
                     "loglik", "iterations", "convergence")


modsem_inspect_da <- function(model, what = "default") {
  stopif(!length(what), "`what` is of length zero!")

  finalModel <- model$model
  groupModels <- finalModel$models
  n.groups <- length(groupModels)
  is.multi <- n.groups > 1L

  group.names <- names(groupModels)
  if (!length(group.names) || any(!nzchar(group.names))) {
    level.names <- finalModel$info$group.levels
    if (!is.null(level.names) && length(level.names) == n.groups) {
      group.names <- level.names
    } else {
      group.names <- paste0("Group", seq_len(n.groups))
    }
  }

  expected.raw <- model$expected.matrices
  expected.by.group <- vector("list", length = n.groups)
  names(expected.by.group) <- group.names

  if (is.list(expected.raw) && length(expected.raw)) {
    names.raw <- names(expected.raw)
    for (g in seq_len(n.groups)) {
      idx <- integer()
      if (length(names.raw)) {
        candidates <- unique(c(as.character(g), group.names[[g]]))
        idx <- match(candidates, names.raw, nomatch = NA_integer_)
        idx <- idx[!is.na(idx)]
      }
      if (!length(idx) && g <= length(expected.raw)) {
        idx <- g
      }
      if (length(idx)) expected.by.group[[g]] <- expected.raw[[idx[[1]]]]
    }
  } else if (!is.null(expected.raw)) {
    expected.by.group <- replicate(n.groups, expected.raw, simplify = FALSE)
    names(expected.by.group) <- group.names
  }

  build_expected_payload <- function(expected) {
    if (is.null(expected)) {
      return(list(
        cov.ov  = NULL,
        cov.lv  = NULL,
        cov.all = NULL,
        cor.ov  = NULL,
        cor.lv  = NULL,
        cor.all = NULL,
        mean.lv = NULL,
        mean.ov = NULL,
        mean.all = NULL,
        r2.all  = NULL,
        r2.lv   = NULL,
        r2.ov   = NULL,
        res.all = NULL,
        res.lv  = NULL,
        res.ov  = NULL
      ))
    }

    cov.ov  <- expected$sigma.ov
    cov.lv  <- expected$sigma.lv
    cov.all <- expected$sigma.all

    cor.ov  <- if (!is.null(cov.ov))  cov2cor(cov.ov)   else NULL
    cor.lv  <- if (!is.null(cov.lv))  cov2cor(cov.lv)   else NULL
    cor.all <- if (!is.null(cov.all)) cov2cor(cov.all)  else NULL

    list(
      cov.ov  = modsemMatrix(cov.ov,  symmetric = TRUE),
      cov.lv  = modsemMatrix(cov.lv,  symmetric = TRUE),
      cov.all = modsemMatrix(cov.all, symmetric = TRUE),

      cor.ov  = modsemMatrix(cor.ov,  symmetric = TRUE),
      cor.lv  = modsemMatrix(cor.lv,  symmetric = TRUE),
      cor.all = modsemMatrix(cor.all, symmetric = TRUE),

      mean.lv  = modsemMatrix(expected$mu.lv),
      mean.ov  = modsemMatrix(expected$mu.ov),
      mean.all = modsemMatrix(expected$mu.all),

      r2.all  = modsemVector(expected$r2.all),
      r2.lv   = modsemVector(expected$r2.lv),
      r2.ov   = modsemVector(expected$r2.ov),

      res.all = modsemVector(expected$res.all),
      res.lv  = modsemVector(expected$res.lv),
      res.ov  = modsemVector(expected$res.ov)
    )
  }

  build_group_payload <- function(submodel, expected) {
    matrices         <- submodel$matrices
    matricesCovModel <- submodel$covModel$matrices
    fetchCov <- function(name) {
      if (!is.null(matricesCovModel)) matricesCovModel[[name]] else NULL
    }

    lambda       <- diagPartitionedMat(matrices$lambdaX, matrices$lambdaY)
    theta        <- diagPartitionedMat(matrices$thetaDelta, matrices$thetaEpsilon)
    gamma.xi     <- diagPartitionedMat(matrices$gammaXi,     fetchCov("gammaXi"))
    gamma.eta    <- diagPartitionedMat(matrices$gammaEta,    fetchCov("gammaEta"))
    omega.xi.xi  <- diagPartitionedMat(matrices$omegaXiXi,   fetchCov("omegaXiXi"))
    omega.eta.xi <- diagPartitionedMat(matrices$omegaEtaXi,  fetchCov("omegaEtaXi"))
    phi          <- diagPartitionedMat(matrices$phi,         fetchCov("phi"))
    psi          <- diagPartitionedMat(matrices$psi,         fetchCov("psi"))

    tau   <- rbind(matrices$tauX, matrices$tauY)
    alpha <- matrices$alpha
    beta0 <- matrices$beta0

    if (!is.null(tau))   colnames(tau)   <- "~1"
    if (!is.null(alpha)) colnames(alpha) <- "~1"
    if (!is.null(beta0)) colnames(beta0) <- "~1"

    c(
      list(
        N     = submodel$data$n,
        data  = submodel$data$data.full,
        lambda       = modsemMatrix(lambda),
        tau          = modsemMatrix(tau),
        theta        = modsemMatrix(theta, symmetric = TRUE),
        gamma.xi     = modsemMatrix(gamma.xi),
        gamma.eta    = modsemMatrix(gamma.eta),
        omega.xi.xi  = modsemMatrix(omega.xi.xi),
        omega.eta.xi = modsemMatrix(omega.eta.xi),
        phi          = modsemMatrix(phi, symmetric = TRUE),
        psi          = modsemMatrix(psi, symmetric = TRUE),
        alpha        = modsemMatrix(alpha),
        beta0        = modsemMatrix(beta0)
      ),
      build_expected_payload(expected)
    )
  }

  group.payloads <- Map(build_group_payload, groupModels, expected.by.group)
  names(group.payloads) <- group.names

  collapse_field <- function(field) {
    values <- lapply(group.payloads, `[[`, field)
    if (!is.multi) {
      return(values[[1]])
    }
    names(values) <- group.names
    values
  }

  N.val          <- collapse_field("N")
  data.val       <- collapse_field("data")
  lambda.val     <- collapse_field("lambda")
  tau.val        <- collapse_field("tau")
  theta.val      <- collapse_field("theta")
  gamma.xi.val   <- collapse_field("gamma.xi")
  gamma.eta.val  <- collapse_field("gamma.eta")
  omega.xi.xi.val  <- collapse_field("omega.xi.xi")
  omega.eta.xi.val <- collapse_field("omega.eta.xi")
  phi.val        <- collapse_field("phi")
  psi.val        <- collapse_field("psi")
  alpha.val      <- collapse_field("alpha")
  beta0.val      <- collapse_field("beta0")
  cov.ov.val     <- collapse_field("cov.ov")
  cov.lv.val     <- collapse_field("cov.lv")
  cov.all.val    <- collapse_field("cov.all")
  cor.ov.val     <- collapse_field("cor.ov")
  cor.lv.val     <- collapse_field("cor.lv")
  cor.all.val    <- collapse_field("cor.all")
  mean.lv.val    <- collapse_field("mean.lv")
  mean.ov.val    <- collapse_field("mean.ov")
  mean.all.val   <- collapse_field("mean.all")
  r2.all.val     <- collapse_field("r2.all")
  r2.lv.val      <- collapse_field("r2.lv")
  r2.ov.val      <- collapse_field("r2.ov")
  res.all.val    <- collapse_field("res.all")
  res.lv.val     <- collapse_field("res.lv")
  res.ov.val     <- collapse_field("res.ov")

  info <- list(N                 = N.val,
               vcov.all          = modsemMatrix(model$vcov.all, symmetric = TRUE),
               vcov.free         = modsemMatrix(model$vcov.free, symmetric = TRUE),
               information       = modsemMatrix(model$FIM, symmetric = TRUE),
               data              = data.val,
               coefficients.all  = modsemVector(model$coefs.all),
               coefficients.free = modsemVector(model$coefs.free),
               partable          = modsemParTable(model$parTable),
               partable.input    = model$originalParTable,
               loglik            = model$logLik,
               iterations        = model$iterations,
               convergence       = model$convergence,

               lambda       = lambda.val,
               tau          = tau.val,
               theta        = theta.val,
               gamma.xi     = gamma.xi.val,
               gamma.eta    = gamma.eta.val,
               omega.xi.xi  = omega.xi.xi.val,
               omega.eta.xi = omega.eta.xi.val,

               phi   = phi.val,
               psi   = psi.val,

               alpha = alpha.val,
               beta0 = beta0.val,

               cov.ov  = cov.ov.val,
               cov.lv  = cov.lv.val,
               cov.all = cov.all.val,

               cor.ov  = cor.ov.val,
               cor.lv  = cor.lv.val,
               cor.all = cor.all.val,

               mean.lv  = mean.lv.val,
               mean.ov  = mean.ov.val,
               mean.all = mean.all.val,

               r2.all  = r2.all.val,
               r2.lv   = r2.lv.val,
               r2.ov   = r2.ov.val,

               res.all = res.all.val,
               res.lv  = res.lv.val,
               res.ov  = res.ov.val
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
      EXPR      = what,
      coef.all  = info[c("vcov.all", "coefficients.all")],
      coef      = info[c("vcov.all", "coefficients.all")],
      coef.free = info[c("vcov.free", "coefficients.free")],
      default   = info[names(info) != "data"],
      all       = info,
      matrices  = info[inspectDA_Matrices],
      optim     = info[inspectDA_Optim],
      fit       = FIT(),
      info[[what]]
    )
  }

  nullvalues <- vapply(fields, FUN.VALUE = logical(1L), FUN = is.null)

  warnif(any(nullvalues), "Some fields in `modsem_inspect()` could not be retrieved!",
         immediate. = FALSE)

  fields
}
