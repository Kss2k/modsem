calcFIM_da <- function(model,
                       finalModel,
                       theta,
                       data = NULL,
                       method = "lms",
                       calc.se = TRUE,
                       FIM = "observed",
                       robust.se = FALSE,
                       P = NULL,
                       hessian = FALSE,
                       EFIM.parametric = TRUE,
                       NA__ = -999,
                       EFIM.S = 3e4,
                       epsilon = 1e-8,
                       R.max = 1e6,
                       verbose = FALSE) {
  if (!calc.se) return(list(FIM = NULL, vcov = NULL, vcov.sub = NULL, type = "none",
                            raw.labels = names(theta), n.additions = 0))
  if (verbose) printf("Calculating standard errors (%s)\n", FIM)

  I <- switch(method,
     lms =
       switch(FIM,
          observed = calcOFIM_LMS(model, theta = theta, data = data,
                                  epsilon = epsilon, hessian = hessian, P = P),
          expected = calcEFIM_LMS(model, finalModel = finalModel, theta = theta, 
                                  data = data, epsilon = epsilon, S = EFIM.S, 
                                  parametric = EFIM.parametric, verbose = verbose,
                                  R.max = R.max, P = P),
          stop2("FIM must be either expected or observed")),
     qml =
       switch(FIM,
          observed = calcOFIM_QML(model, theta = theta, data = data,
                                  hessian = hessian, epsilon = epsilon),
          expected = calcEFIM_QML(model, finalModel = finalModel, theta = theta, 
                                  data = data, epsilon = epsilon, S = EFIM.S,
                                  parametric = EFIM.parametric, verbose = verbose,
                                  R.max = R.max),
          stop2("FIM must be either expected or observed")),
     stop2("Unrecognized method: ", method)
  )


  if (robust.se) {
    warnif(hessian && FIM == "observed",
           "`robust.se = TRUE` should not be paired with ",
           "`OFIM.hessian = TRUE` and `FIM = \"observed\"`")
    H <- calcHessian(model, theta = theta, data = data, method = method,
                     epsilon = epsilon, P = P)
    invH <- solveFIM(H, NA__ = NA__)

    vcov <- invH %*% I %*% invH

  } else {
    vcov <- solveFIM(I, NA__ = NA__)
  }

  vcov.all <- getVCOV_LabelledParams(vcov = vcov, model = model, theta = theta,
                                     method = method)

  nAdditions   <- ncol(vcov.all) - ncol(vcov)
  lavLabels    <- model$lavLabels
  subLavLabels <- lavLabels[colnames(vcov.all) %in% names(theta)]
  rawLabels    <- colnames(vcov.all)
  dimnames(vcov.all) <- list(lavLabels, lavLabels)
  dimnames(I) <- dimnames(vcov) <- list(subLavLabels, subLavLabels)

  list(FIM = I, vcov = vcov.all, vcov.sub = vcov, type = FIM,
       raw.labels = rawLabels, n.additions = nAdditions)
}


fdHESS <- function(pars, ...) {
  tryCatch(
    nlme::fdHess(pars = pars, ...)$Hessian,
    error = function(e) {
      warning2("Calculation of Hessian matrix failed...\n  ", e$message)
      matrix(NA, nrow = length(pars), ncol = length(pars)) 
    }
  )
}


calcHessian <- function(model, theta, data, method = "lms",
                        epsilon = 1e-8, P = NULL) {
  if (method == "lms") {
    if (is.null(P)) P <- estepLms(model, theta = theta, data = data)
    # negative hessian (sign = -1)
    H <- fdHESS(pars = theta, fun = obsLogLikLms, model = model,
                data = data, P = P, sign = -1,
                .relStep = .Machine$double.eps^(1/5))
    # I assumed this would be faster, but using nlme::fdHess seems faster than
    # using the gradient function...
    # f <- \(theta) gradientObsLogLikLms(theta = theta, model = model, data = data,
    #                                    P = P, sign = -1, eps = epsilon)
    # H <- calcHessFromGradient(gradFun = f, theta = theta, eps = epsilon)

    # I assumed this would be faster, but using nlme::fdHess seems faster than
    # using the gradient function...
    # f <- \(theta) gradientObsLogLikLms(theta = theta, model = model, data = data,
    #                                    P = P, sign = -1, eps = epsilon)
    # H <- calcHessFromGradient(gradFun = f, theta = theta, eps = epsilon)

  } else if (method == "qml") {
    # negative hessian (sign = -1)
    H <- fdHESS(pars = theta, fun = logLikQml, model = model, sign = -1,
                .relStep = .Machine$double.eps^(1/5))
  }

  H
}


solveFIM <- function(H, NA__ = -999) {
  tryCatch(solve(H),
           error = function(e) {
             H[TRUE] <- NA__
             H
           },
           warning = function(w)
             if (grepl("NaN", conditionMessage(w))) suppressWarnings(solve(H)) else solve(H)
  )
}


calcSE_da <- function(calc.se = TRUE, vcov, rawLabels, NA__ = -999) {
  if (!calc.se) return(rep(NA__, length(rawLabels)))
  if (is.null(vcov)) {
    warning2("Fisher Information Matrix (FIM) was not calculated, ",
             "unable to compute standard errors")
    return(rep(NA__, length(rawLabels)))
  }

  se <- suppressWarnings(sqrt(diag(vcov)))

  if (all(is.na(se))) {
    warning2("SE's could not be computed, negative Hessian is singular.")
  } else if (any(is.nan(se))) {
    warning2("SE's for some coefficients could not be computed.")
  }

  if (!is.null(names(se))) names(se) <- rawLabels
  se[is.na(se)] <- NA__
  se
}


calcOFIM_LMS <- function(model, theta, data, hessian = FALSE,
                         epsilon = 1e-6, P = NULL) {
  if (is.null(P)) P <- estepLms(model, theta = theta, data = data)

  N <- nrow(data)
  if (hessian) {
    # negative hessian (sign = -1)
    I <- calcHessian(model, theta = theta, data = data, 
                     method = "lms", epsilon = epsilon, P = P)

    return(I)
  }

  S <- gradientObsLogLikLms_i(theta, model = model, data = data,
                              P = P, sign = 1, epsilon = epsilon)
  crossprod(S)
}


calcEFIM_LMS <- function(model, finalModel = NULL, theta, data,
                         S         = 100,
                         parametric = TRUE,
                         epsilon    = 1e-6,
                         verbose    = FALSE,
                         R.max      = 1e6,
                         P          = NULL) {
  k <- length(theta)                       # number of free parameters
  N <- nrow(data)
  R <- min(R.max, N * S)
  warnif(R.max <= N, "R.max is less than N!")

  if (parametric) {
    stopif(is.null(finalModel), "`finalModel` is needed when parametric = TRUE")
    parTable   <- modelToParTable(finalModel, method = "lms")
    population <- simulateDataParTable(parTable, N = R,
                                       colsOVs = colnames(data))$oV
  } else {
    population <- data[sample(nrow(data), R, replace = TRUE), , drop = FALSE]
  }

  popEstep <- estepLms(model      = model,
                       theta      = theta,
                       data       = population,
                       recalcQuad = TRUE,
                       lastQuad   = if(!is.null(P)) P$quad else NULL)

  J <- gradientObsLogLikLms_i(theta = theta,
                              model = model,
                              data  = population,
                              P     = popEstep,
                              sign  = +1,
                              epsilon = epsilon)      # R × k matrix

  I <- matrix(0, nrow = k, ncol = k)
  for (i in seq_len(S)) {
    if (R == N * S) {
      # non-overlapping split
      idx1 <- (i - 1) * N + 1
      sub  <- idx1:(idx1 + N - 1)
    } else {                    
      sub <- sample(R, N)
    }

    I <- I + crossprod(J[sub, , drop = FALSE])
  }

  if (verbose) cat("\n")

  I / S
}


calcOFIM_QML <- function(model, theta, data, hessian = FALSE,
                         epsilon = 1e-8) {
  N <- nrow(model$data)

  if (hessian) {
    # negative hessian (sign = -1)
    I <- fdHESS(pars = theta, fun = logLikQml, model = model,
                sign = -1, .relStep = .Machine$double.eps^(1/5))
    return(I)
  }

  S <- gradientLogLikQml_i(theta, model = model, sign = 1, epsilon = epsilon)
  crossprod(S)
}


calcEFIM_QML <- function(model, finalModel = NULL, theta, data, S = 100,
                         parametric = TRUE, epsilon = 1e-8, verbose = FALSE,
                         R.max = 1e6) {
  k <- length(theta)                       # number of free parameters
  N <- nrow(data)
  R <- min(R.max, N * S)
  warnif(R.max <= N, "R.max is less than N!")

  if (parametric) {
    if (is.null(finalModel)) stop2("finalModel must be included in calcEFIM_QML")
    parTable <- modelToParTable(finalModel, method = "qml")
    population <- tryCatch(
      simulateDataParTable(parTable, N = R, colsOVs = colnames(data))$oV,
      error = function(e) {
        warning2("Unable to simulate data for EFIM, using stochastic sampling instead")
        calcEFIM_QML(model = model, theta = theta, data = data, S = S, 
                     parametric = FALSE, epsilon = epsilon)
      }
    )

  } else population <- data[sample(R, N, replace = TRUE), ]

  J <- gradientLogLikQml(theta = theta, model = model, sign = +1, 
                         epsilon = epsilon, data = population)

  I <- matrix(0, nrow = k, ncol = k)
  for (i in seq_len(S)) {
    if (R == N * S) {
      # non-overlapping split
      idx1 <- (i - 1) * N + 1
      sub  <- idx1:(idx1 + N - 1)
    } else {                    
      sub <- sample(R, N)
    }

    I <- I + crossprod(J[sub, , drop = FALSE])
  }

  if (verbose) cat("\n")

  I / S
}


getSE_Model <- function(model, se, method, n.additions) {
  model$lenThetaLabel <- model$lenThetaLabel + n.additions
  fillModel(replaceNonNaModelMatrices(model, value = -999),
            theta = se, method = method)
}


calcHessFromGradient <- function(gradFun, theta, eps = 1e-6, ...) {
  p  <- length(theta)
  H  <- matrix(NA_real_, p, p)
  g0 <- gradFun(theta, ...)

  for (j in seq_len(p)) {
    th_j        <- theta
    th_j[j]     <- th_j[j] + eps
    g_fwd       <- gradFun(th_j, ...)
    H[, j]      <- (g_fwd - g0) / eps
  }

  # enforce symmetry against numerical noise
  H <- 0.5 * (H + t(H))

  dimnames(H) <- list(names(theta), names(theta))
  H
}
