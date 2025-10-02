obsLogLikLmsPML <- function(theta, model, data, P, sign = 1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  ll <- observedLogLikLmsPMLCpp(modFilled, dataR = data$data.split, P = P,
                                colidxR = data$colidx0, n = data$n.pattern,
                                npatterns = data$p, ncores = ThreadEnv$n.threads)
  ll * sign
}


gradientObsLogLikLmsPML <- function(theta, model, data, P, sign = 1, epsilon = 1e-6) {
  FGRAD <- function(modelR, P, block, row, col, symmetric, colidxR, npatterns,
                    eps, ncores, n, ...) {
    gradObsLogLikLmsPMLCpp(modelR = modelR, dataR = data$data.split, P = P,
                           block = block, row = row, col = col,
                           symmetric = symmetric, colidxR = colidxR,
                           n = n, npatterns = npatterns, eps = eps,
                           ncores = ncores)
  }

  FOBJECTIVE <- function(theta, model, colidxR, npatterns, sign, ...) {
    obsLogLikLmsPML(theta = theta, model = model, data = data$data.split,
                    P = P, colidxR = colidxR, npatterns = npatterns,
                    sign = sign)
  }

  gradientAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                       epsilon = epsilon, data = data,
                       FGRAD = FGRAD, FOBJECTIVE = FOBJECTIVE)
}


hessianObsLogLikLmsPML <- function(theta, model, data, P, sign = -1,
                                   .relStep = .Machine$double.eps ^ (1/5)) {

  FHESS <- function(modelR, P, block, row, col, symmetric, eps, .relStep, colidxR, n,
                    npatterns, ncores, ...) {
    hessObsLogLikLmsPMLCpp(modelR = modelR, dataR = data$data.split, P = P,
                           block = block, row = row, col = col,
                           symmetric = symmetric, npatterns = npatterns,
                           colidxR = colidxR, n = n, relStep = .relStep,
                           minAbs = 0.0, ncores = ncores)
  }

  FOBJECTIVE <- function(theta, model, P, sign, data, ...) {
    obsLogLikLmsPML(theta = theta, model = model, P = P, data = data, sign = sign)
  }

  hessianAllLogLikLms(theta = theta, model = model, P = P, sign = sign,
                      data = data, FHESS = FHESS, FOBJECTIVE = FOBJECTIVE,
                      .relStep = .relStep)
}
