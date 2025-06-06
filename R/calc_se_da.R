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
                                  epsilon = epsilon, hessian = hessian),
          expected = calcEFIM_LMS(model, finalModel = finalModel, theta = theta, 
                                  data = data, epsilon = epsilon, S = EFIM.S, 
                                  parametric = EFIM.parametric, verbose = verbose,
                                  R.max = R.max),
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
           "'robust.se = TRUE' should not be paired with ",
           "'EFIM.hessian = TRUE' && 'FIM = \"observed\"'")
    H <- calcHessian(model, theta = theta, data = data, method = method,
                     epsilon = epsilon)
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
                        epsilon = 1e-8) {
  if (method == "lms") {
    P <- estepLms(model, theta = theta, data = data)
    # negative hessian (sign = -1)
    H <- fdHESS(pars = theta, fun = logLikLms, model = model,
                data = data, P = P, sign = -1,
                .relStep = .Machine$double.eps^(1/5))

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
                         epsilon = 1e-6) {
  N <- nrow(data)
  P <- estepLms(model, theta = theta, data = data)
  if (hessian) {
    # negative hessian (sign = -1)
    I <- fdHESS(pars = theta, fun = obsLogLikLms, model = model,
                data = data, P = P, sign = -1,
                .relStep = .Machine$double.eps^(1/5))
    return(I)
  }

  J <- gradientObsLogLikLms_i(theta, model = model, data = data,
                           P = P, sign = 1, epsilon = epsilon)
  I <- matrix(0, nrow = length(theta), ncol = length(theta))

  for (i in seq_len(N)) {
    J_i <- J[i,]
    I <- I + J_i %*% t(J_i)
  }

  I
}


calcEFIM_LMS <- function(model, finalModel = NULL, theta, data, S = 100,
                         parametric = TRUE, epsilon = 1e-6, verbose = FALSE,
                         R.max = 1e6) {
  k <- length(theta)

  if (S <= k) {
    warning2("`EFIM.S` is lower than the number of free parameters! Increasing `EFIM.S` to ", 2 * k)
    S <- 2 * k
  }

  N      <- nrow(model$data)
  R.ceil <- N * S > R.max
  R      <- ifelse(R.ceil, yes = R.max, no = N * S) 
 
  if (R < N) {
    warning2("Population size is smaller than sample size, please increase it using the `R.max` argument!")
    R <- N
  }
  
  if (parametric) {
    if (is.null(finalModel)) stop2("finalModel must be included in calcEFIM_LMS")
    parTable   <- modelToParTable(finalModel, method = "lms")
    population <- tryCatch(
      simulateDataParTable(parTable, N = R, colsOVs = colnames(data))$oV,
      error = function(e) {
        warning2("Unable to simulate data for EFIM, using stochastic sampling instead")
        calcEFIM_LMS(model = model, theta = theta, data = data, S = S, 
                     parametric = FALSE, epsilon = epsilon)
      })

  } else population <- data[sample(R, N, replace = TRUE), ]

  I <- matrix(0, nrow = k, ncol = k)

  for (i in seq_len(S)) {
    if (verbose) printf("\rMonte-Carlo: Iteration = %d/%d", i, S)

    if (!R.ceil) {
      n1  <- (i - 1) * N + 1
      nn  <- n1 + N - 1
      sub <- n1:nn
    } else sub <- sample(R, N)

    P <- estepLms(model = model, theta = theta, data = population[sub, ])
    J <- gradientObsLogLikLms(theta = theta, model = model, data = population[sub, ],
                              P = P, sign = 1, epsilon = epsilon)

    I <- I + J %*% t(J)
  }
  
  if (verbose) {
    cat("\n")
    if (S <= 100) message("Consider increasing the Monte-Carlo Monte-Carloiterations, using the `EFIM.S` argument!")
  }

  I / (S - k) # divide by degrees of freedom
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

  J <- gradientLogLikQml_i(theta, model = model, sign = 1, epsilon = epsilon)
  I <- matrix(0, nrow = length(theta), ncol = length(theta))

  for (i in seq_len(N)) {
    J_i <- J[i,]
    I <- I + J_i %*% t(J_i)
  }

  I
}


calcEFIM_QML <- function(model, finalModel = NULL, theta, data, S = 100,
                         parametric = TRUE, epsilon = 1e-8, verbose = FALSE,
                         R.max = 1e6) {
  k <- length(theta)
  
  if (S <= k) {
    warning2("`EFIM.S` is lower than the number of free parameters! Increasing `EFIM.S` to ", 2 * k)
    S <- 2 * k
  }

  N      <- nrow(model$data)
  R.ceil <- N * S > R.max
  R      <- ifelse(R.ceil, yes = R.max, no = N * S) 
  
  if (R < N) {
    warning2("Population size is smaller than sample size, please increase it using the `R.max` argument!")
    R <- N
  }
  
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


  I <- matrix(0, nrow = k, ncol = k)
  for (i in seq_len(S)) {
    if (verbose) printf("\rMonte-Carlo: Iteration = %d/%d", i, S)
   
    if (!R.ceil) {
      n1  <- (i - 1) * N + 1
      nn  <- n1 + N - 1
      sub <- n1:nn
    } else sub <- sample(R, N)

    J <- gradientLogLikQml(theta = theta, model = model, sign = 1, epsilon = epsilon,
                           data = population[sub, ])
    I <- I + J %*% t(J)
  }

  if (verbose) {
    cat("\n")
    if (S <= 100) message("Consider increasing the Monte-Carlo iterations, using the `EFIM.S` argument!")
  }

  I / (S - k) # divide by degrees of freedom
}


getSE_Model <- function(model, se, method, n.additions) {
  model$lenThetaLabel <- model$lenThetaLabel + n.additions
  fillModel(replaceNonNaModelMatrices(model, value = -999),
            theta = se, method = method)
}
