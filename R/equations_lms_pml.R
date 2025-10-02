# Observed PML log-likelihood (uses PCS when provided)
obsLogLikLmsPML <- function(theta, model, data, P, sign = 1,
                            pcs_xptr = NULL, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  # Optionally build the PCS pointer if missing (cheap & safe even if some patterns aren't ordinal)
  if (is.null(pcs_xptr)) {
    pcs_xptr <- buildPCS_Xptr(
      dataR         = data$data.split,
      colidxR       = data$colidx0,
      isOrderedEnum = modFilled$info$isOrderedEnum,
      thresholds    = modFilled$matrices$thresholds
    )
  }

  ll <- observedLogLikLmsPMLCpp(
    modelR  = modFilled,
    dataR   = data$data.split,
    colidxR = data$colidx0,
    P       = P,
    pcs_xptr = pcs_xptr,
    n       = data$n.pattern,
    npatterns = data$p,
    ncores  = ThreadEnv$n.threads
  )

  sign * ll
}


# Gradient of observed PML log-likelihood (finite differences under the hood on C++)
gradientObsLogLikLmsPML <- function(theta, model, data, P, sign = 1,
                                    epsilon = 1e-6,
                                    pcs_xptr = NULL, build_pcs_if_null = TRUE) {

  # We'll reuse the same PCS pointer for the closures
  get_pcs <- function(theta_local) {
    if (!is.null(pcs_xptr)) return(pcs_xptr)
    if (!isTRUE(build_pcs_if_null)) return(NULL)
    # build from the model at theta_local to get thresholds shape/isOrderedEnum
    modFilled <- fillModel(model = model, theta = theta_local, method = "lms")
    buildPCS_Xptr(
      dataR         = data$data.split,
      colidxR       = data$colidx0,
      isOrderedEnum = modFilled$info$isOrderedEnum,
      thresholds    = modFilled$matrices$thresholds
    )
  }

  FGRAD <- function(modelR, P, block, row, col, symmetric, colidxR, npatterns,
                    eps, ncores, n, ...) {

    pcs_ptr <- get_pcs(theta)  # build once per gradient eval (fast; can be cached by caller)

    gradObsLogLikLmsPMLCpp(
      modelR = modelR, dataR = data$data.split, colidxR = colidxR, P = P,
      pcs_xptr = pcs_ptr,
      block = block, row = row, col = col, symmetric = symmetric,
      n = n, npatterns = npatterns, eps = eps, ncores = ncores
    )
  }

  FOBJECTIVE <- function(theta, model, colidxR, npatterns, sign, ...) {
    obsLogLikLmsPML(theta = theta, model = model, data = data,
                    P = P, colidxR = colidxR, npatterns = npatterns,
                    sign = sign, pcs_xptr = get_pcs(theta),
                    build_pcs_if_null = FALSE)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, data = data,
                       FGRAD = FGRAD, FOBJECTIVE = FOBJECTIVE)
}


# Hessian of observed PML log-likelihood (FD Hessian on top of the PCS/non-PCS objective)
hessianObsLogLikLmsPML <- function(theta, model, data, P, sign = -1,
                                   .relStep = .Machine$double.eps^(1/5),
                                   pcs_xptr = NULL, build_pcs_if_null = TRUE) {

  get_pcs <- function(theta_local) {
    if (!is.null(pcs_xptr)) return(pcs_xptr)
    if (!isTRUE(build_pcs_if_null)) return(NULL)
    modFilled <- fillModel(model = model, theta = theta_local, method = "lms")
    buildPCS_Xptr(
      dataR         = data$data.split,
      colidxR       = data$colidx0,
      isOrderedEnum = modFilled$info$isOrderedEnum,
      thresholds    = modFilled$matrices$thresholds
    )
  }

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR, n,
                    npatterns, ncores, ...) {
    pcs_ptr <- get_pcs(theta)

    hessObsLogLikLmsPMLCpp(
      modelR = modelR, dataR = data$data.split, colidxR = colidxR, P = P,
      pcs_xptr = pcs_ptr,
      block = block, row = row, col = col, symmetric = symmetric,
      n = n, npatterns = npatterns, relStep = .relStep, minAbs = 0.0, ncores = ncores
    )
  }

  FOBJECTIVE <- function(theta, model, P, sign, data, ...) {
    obsLogLikLmsPML(theta = theta, model = model, data = data,
                    P = P, sign = sign, pcs_xptr = get_pcs(theta),
                    build_pcs_if_null = FALSE)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      data = data, FHESS = FHESS, FOBJECTIVE = FOBJECTIVE,
                      .relStep = .relStep)
}
