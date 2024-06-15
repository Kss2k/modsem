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
                       verbose = FALSE) {
  if (!calc.se) return(NULL)
  if (verbose) cat("Calculating standard errors\n")
  
  I <- switch(method, 
     lms = 
       switch(FIM, 
          observed = calcOFIM_LMS(model, theta = theta, data = data, 
                                  hessian = hessian),
          expected = calcEFIM_LMS(model, finalModel = finalModel, 
                                  theta = theta, data = data, 
                                  S = EFIM.S, parametric = EFIM.parametric),
          stop2("FIM must be either expected or observed")),
     qml = 
       switch(FIM, 
          observed = calcOFIM_QML(model, theta = theta, data = data, 
                                  hessian = hessian),
          expected = calcEFIM_QML(model, finalModel = finalModel, 
                                  theta = theta, data = data, 
                                  S = EFIM.S, parametric = EFIM.parametric),
          stop2("FIM must be either expected or observed")),
     stop2("Unrecognized method: ", method) 
  )

  if (robust.se) {
    if (hessian && FIM == "observed")
      warning("'robust.se = TRUE' should not be paired with 'EFIM.hessian = TRUE' && 'FIM = \"observed\"'")
    H <- calcHessian(model, theta = theta, data = data, method = method)
    invH <- solveHessian(H, NA__ = NA__)
    vcov <- invH %*% I %*% invH

  } else {
    vcov <- solveHessian(I, NA__ = NA__)
  }

  list(FIM = I, vcov = vcov)
}

calcHessian <- function(model, theta, data, method = "lms") {
  if (method == "lms") {
    P <- estepLms(model, theta = theta, data = data)
    # negative hessian (sign = -1)
    H <- nlme::fdHess(pars = theta, fun = logLikLms, 
                      model = model, data = data, P = P, sign = -1,
                      .relStep = .Machine$double.eps^(1/5))$Hessian

  } else if (method == "qml") {
    # negative hessian (sign = -1)
    H <- nlme::fdHess(pars = theta, fun = logLikQml, 
                      model = model, sign = -1,
                      .relStep = .Machine$double.eps^(1/5))$Hessian
  }

  H
}


solveHessian <- function(H, NA__ = -999) {
  tryCatch(solve(H), 
           error = function(e) {
             H[TRUE] <- NA__
             H
           }, 
           warning = function(w) 
             if (grepl("NaN", conditionMessage(w))) suppressWarnings(solve(H)) else solve(H)
  )
}


calcSE_da <- function(calc.se = TRUE, vcov, theta, NA__ = -999) {
  if (!calc.se) return(rep(NA__, length(theta)))
  if (is.null(vcov)) {
    warning2("Fisher Information Matrix (FIM) was not calculated, ",
             "unable to compute standard errors")
    return(rep(NA__, length(theta)))
  }

  se <- sqrt(diag(vcov))

  if (all(is.na(se))) 
    warning2("SE's could not be computed, negative Hessian is singular.")
  if (any(is.nan(se))) 
    warning2("SE's for some coefficients could not be computed.") 

  if (!is.null(names(se))) names(se) <- names(theta)
  se[is.na(se)] <- NA__
  se
}


calcOFIM_LMS <- function(model, theta, data, hessian = FALSE) {
  N <- nrow(data)
  P <- estepLms(model, theta = theta, data = data)
  if (hessian) {
    # negative hessian (sign = -1)
    I <- nlme::fdHess(pars = theta, fun = logLikLms, 
                      model = model, data = data, P = P, sign = -1,
                      .relStep = .Machine$double.eps^(1/5))$Hessian
    return(I)
  }
  J <- gradientLogLikLms_i(theta, model = model, data = data, 
                           P = P, sign = 1)  
  I <- matrix(0, nrow = length(theta), ncol = length(theta))
  for (i in seq_len(N)) I <- I + J[i, ] %*% t(J[i, ])

  I
}


calcEFIM_LMS <- function(model, finalModel = NULL, theta, data, S = 3e4, 
                         parametric = TRUE) {
  N <- nrow(data)

  if (parametric) {
    if (is.null(finalModel)) stop2("finalModel must be included in calcEFIM_LMS")
    parTable <- modelToParTable(finalModel, method = "lms")
    population <- simulateDataParTable(parTable, N = S)$oV

  } else {
    population <- data[sample(N, S, replace = TRUE)]
  }

  P <- estepLms(model, theta, data = population)
  J <- gradientLogLikLms_i(theta, model = model, data = population, 
                                     P = P, sign = 1)

  I <- matrix(0, nrow = length(theta), ncol = length(theta))
  for (i in seq_len(S)) I <- I + J[i, ] %*% t(J[i, ])

  I / (S / N)
}


calcOFIM_QML <- function(model, theta, data, hessian = FALSE) {
  N <- nrow(model$data)
  if (hessian) {
    # negative hessian (sign = -1)
    I <- nlme::fdHess(pars = theta, fun = logLikQml, 
                      model = model, sign = -1,
                      .relStep = .Machine$double.eps^(1/5))$Hessian
    return(I)
  }

  J <- gradientLogLikQml_i(theta, model = model, sign = 1)  
  I <- matrix(0, nrow = length(theta), ncol = length(theta))
  for (i in seq_len(N)) I <- I + J[i, ] %*% t(J[i, ])

  I
}


calcEFIM_QML <- function(model, finalModel = NULL, theta, data, S = 3e4, 
                         parametric = TRUE) {
  N <- nrow(model$data)

  if (parametric) {
    if (is.null(finalModel)) stop2("finalModel must be included in calcEFIM_QML")
    parTable <- modelToParTable(finalModel, method = "qml")
    model$data <- simulateDataParTable(parTable, N = S)$oV

  } else {
    model$data <- data[sample(N, S, replace = TRUE)]
  }

  J <- gradientLogLikQml_i(theta, model = model, sign = 1)

  I <- matrix(0, nrow = length(theta), ncol = length(theta))
  for (i in seq_len(S)) I <- I + J[i, ] %*% t(J[i, ])

  I / (S / N)
}
