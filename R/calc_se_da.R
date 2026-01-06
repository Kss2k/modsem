calcFIM_da <- function(model,
                       theta,
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
                       verbose = FALSE,
                       cr1s = TRUE) {
  if (!calc.se) return(list(FIM = NULL, vcov = NULL, vcov.sub = NULL, type = "none",
                            raw.labels = names(theta), n.additions = 0))
  if (verbose) printf("Calculating standard errors (%s)\n", FIM)

  collectCluster <- function(d) {
    if (is.null(d)) return(NULL)

    if (is.list(d) && !is.data.frame(d)) {
      clusters <- lapply(d, function(x) x$cluster)
      clusters <- clusters[lengths(clusters) > 0L]
      if (!length(clusters)) return(NULL)
      do.call(c, clusters)

    } else {
      d$cluster
    }
  }

  DATA <- lapply(model$models, FUN = \(sub) sub$data)
  cluster.vec <- if (robust.se) collectCluster(DATA) else NULL

  I <- switch(method,
     lms =
       switch(FIM,
          observed = calcOFIM_LMS(model, theta = theta, epsilon = epsilon,
                                  hessian = hessian, P = P, robust.se = robust.se,
                                  cluster = cluster.vec, cr1s = cr1s),
          expected = calcEFIM_LMS(model, theta = theta, epsilon = epsilon,
                                  S = EFIM.S, parametric = EFIM.parametric,
                                  verbose = verbose, R.max = R.max, P = P),
          stop2("FIM must be either expected or observed")),
     qml =
       switch(FIM,
          observed = calcOFIM_QML(model, theta = theta, hessian = hessian,
                                  epsilon = epsilon, robust.se = robust.se,
                                  cluster = cluster.vec, cr1s = cr1s),
          expected = calcEFIM_QML(model, theta = theta, epsilon = epsilon,
                                  S = EFIM.S, parametric = EFIM.parametric,
                                  verbose = verbose, R.max = R.max),
          stop2("FIM must be either expected or observed")),
     stop2("Unrecognized method: ", method)
  )


  if (robust.se) {
    warnif(hessian && FIM == "observed",
           "`robust.se = TRUE` should not be paired with ",
           "`OFIM.hessian = TRUE` and `FIM = \"observed\"`")
    H <- calcHessian(model, theta = theta, method = method,
                     epsilon = epsilon, P = P)
    invH <- solveFIM(H, NA__ = NA__)

    vcov <- invH %*% I %*% invH

  } else {
    vcov <- solveFIM(I, NA__ = NA__)
  }

  vcov.all <- getVCOV_LabelledParams(vcov = vcov, model = model, theta = theta,
                                     method = method)

  nAdditions   <- ncol(vcov.all) - ncol(vcov)
  lavLabels    <- model$params$lavLabels
  subLavLabels <- lavLabels[colnames(vcov.all) %in% names(theta)]
  rawLabels    <- colnames(vcov.all)
  dimnames(vcov.all) <- list(lavLabels, lavLabels)
  dimnames(I) <- dimnames(vcov) <- list(subLavLabels, subLavLabels)

  list(FIM = I, vcov.all = vcov.all, vcov.free = vcov, type = FIM,
       raw.labels = rawLabels, n.additions = nAdditions)
}


fdHESS <- function(pars,
                   fun,
                   ...,
                   .relStep = .Machine$double.eps^(1/5),
                   .minAbsPar = 0,
                   .switch.size = 120L,
                   .mem.limit.bytes = 3 * 2 ^ 30) {
  stopifnot(is.numeric(pars))

  npar <- length(pars)
  if (npar == 0)
    return(matrix(0, 0, 0))

  # mirror C++ heuristic: quadratic fit becomes computationally expensive for
  # when there are many free parameters
  m.ls <- 1 + 2 * npar + (npar * (npar - 1)) / 2
  bytes.X <- m.ls * m.ls * 8
  force.full.fd <- npar >= .switch.size || bytes.X > .mem.limit.bytes

  computeFD <- function() {
    tryCatch(
      fdHessFullFD(
        pars = pars,
        fun = fun,
        ...,
        relStep = .relStep,
        minAbsPar = .minAbsPar
      ),
      error = function(e) {
        warning2("Finite-difference Hessian calculation failed...\n  ", e$message)
        matrix(NA, nrow = npar, ncol = npar)
      }
    )
  }

  if (!force.full.fd) {
    tryCatch(
      nlme::fdHess(pars = pars, fun = fun, ..., .relStep = .relStep)$Hessian,
      error = function(e) {
        warning2(
          "Switching to fallback finite-difference Hessian after nlme::fdHess error:\n  ",
          e$message
        )
        computeFD()
      }
    )
  } else computeFD()
}


fdHessFullFD <- function(pars, fun, ..., relStep, minAbsPar) {
  npar <- length(pars)
  par_names <- names(pars)

  relStep <- rep(relStep, length.out = npar)
  minAbsPar <- rep(minAbsPar, length.out = npar)
  base <- as.numeric(pars)
  names(base) <- par_names

  incr <- pmax(abs(base), minAbsPar) * relStep
  incr[incr == 0] <- relStep[incr == 0]

  offset_template <- numeric(npar)
  eval_fun <- function(offset) {
    candidate <- base
    candidate[] <- base + offset
    fun(candidate, ...)
  }

  f0 <- eval_fun(offset_template)
  Hess <- matrix(0, npar, npar)

  for (i in seq_len(npar)) {
    step_plus <- offset_template
    step_minus <- offset_template
    step_plus[i] <- incr[i]
    step_minus[i] <- -incr[i]

    f_ip <- eval_fun(step_plus)
    f_im <- eval_fun(step_minus)

    hi <- incr[i]
    Hess[i, i] <- (f_ip + f_im - 2 * f0) / (hi * hi)
  }

  if (npar > 1) {
    for (i in seq_len(npar - 1)) {
      hi <- incr[i]
      for (j in seq.int(i + 1, npar)) {
        hj <- incr[j]

        step <- offset_template
        step[i] <- hi
        step[j] <- hj
        fpp <- eval_fun(step)

        step[j] <- -hj
        fpm <- eval_fun(step)

        step[i] <- -hi
        fmm <- eval_fun(step)

        step[j] <- hj
        fmp <- eval_fun(step)

        hij <- (fpp - fpm - fmp + fmm) / (4 * hi * hj)
        Hess[i, j] <- hij
        Hess[j, i] <- hij
      }
    }
  }

  if (!is.null(par_names))
    dimnames(Hess) <- list(par_names, par_names)

  Hess
}


calcHessian <- function(model, theta, method = "lms",
                        epsilon = 1e-8, P = NULL) {
  if (method == "lms") {
    if (is.null(P)) P <- estepLms(model, theta = theta)
    # negative hessian (sign = -1)
    fH <- \(model) observedInfoFromLouisLms(model = model, theta = theta,
                                            P = P)$I.obs

    H <- tryCatch(suppressWarnings(fH(model)), error = function(e) {
      warning2("Optimized calculation of Hessian failed, attempting to switch!\n", e)
      model$gradientStruct$hasCovModel <- TRUE

      suppressWarnings(fH(model))
    })

  } else if (method == "qml") {
    # negative hessian (sign = -1)
    suppressWarnings({
      H <- hessianLogLikQml(theta = theta, model = model,
                            sign = -1, .relStep = .Machine$double.eps^(1/5))
    })
  }

  H
}


solveFIM <- function(H, NA__ = -999, use.ginv = FALSE) {
  tryCatch(if (use.ginv) GINV(H) else solve(H),
           error = function(e) {
             if (!use.ginv) return(solveFIM(H, NA__ = NA__, use.ginv = TRUE))

             H[TRUE] <- NA__
             H
           },
           warning = function(w)
             if (grepl("NaN", conditionMessage(w))) suppressWarnings(solve(H)) else solve(H)
  )
}


calcSE_da <- function(calc.se = TRUE, vcov, rawLabels, NA__ = -999) {
  if (!calc.se)
    return(stats::setNames(rep(NA__, length(rawLabels)), nm = rawLabels))

  if (is.null(vcov)) {
    warning2("Fisher Information Matrix (FIM) was not calculated, ",
             "unable to compute standard errors", immediate. = FALSE)
    return(rep(NA__, length(rawLabels)))
  }

  se <- suppressWarnings(sqrt(diag(vcov)))

  if (all(is.na(se))) {
    warning2("Standard errors could not be computed, negative Hessian is singular.",
             immediate. = FALSE)
  } else if (any(is.nan(se))) {
    warning2("Standard errors for some coefficients could not be computed.",
             immediate. = FALSE)
  }

  if (!is.null(names(se))) names(se) <- rawLabels
  se[is.na(se)] <- NA__
  se
}


calcOFIM_LMS <- function(model, theta, hessian = FALSE,
                         epsilon = 1e-6, P = NULL,
                         robust.se = FALSE,
                         cluster   = NULL,
                         cr1s      = TRUE) {
  if (is.null(P))
    P <- estepLms(model, theta = theta)

  if (hessian) {
    # negative hessian (sign = -1)
    I <- calcHessian(model, theta = theta,
                     method = "lms", epsilon = epsilon, P = P)
    return(I)
  }

  # S: N x k matrix of individual score contributions (OPG)
  S <- suppressWarnings(
    gradientObsLogLikLms_i(theta, model = model, P = P, sign = +1,
                           epsilon = epsilon)
  )

  if (!robust.se || is.null(cluster)) {
    # classic OFIM via outer product of gradients (BHHH)
    return(crossprod(S))
  }

  stopif(length(cluster) != nrow(S),
         "Length of 'cluster' must equal the number of rows in the data / scores.")

  f <- as.factor(cluster)
  G <- nlevels(f)
  k <- ncol(S)

  # aggregate scores by cluster: s_g = sum_{i in g} s_i
  Sg <- matrix(0, nrow = G, ncol = k)
  lev <- levels(f)
  for (g in seq_len(G)) {
    idx <- which(f == lev[g])
    Sg[g, ] <- colSums(S[idx, , drop = FALSE])
  }

  B <- crossprod(Sg)  # meat = sum_g s_g s_g'

  # optional CR1S small-sample correction
  if (isTRUE(cr1s)) {
    N <- nrow(S); q <- ncol(S)
    if (G > 1 && N > q) {
      B <- B * (G / (G - 1)) * ((N - 1) / (N - q))
    }
  }

  B
}


calcEFIM_LMS <- function(model,
                         theta,
                         S          = 100,
                         parametric = TRUE,
                         epsilon    = 1e-6,
                         verbose    = FALSE,
                         R.max      = 1e6,
                         P          = NULL) {
  k <- length(theta) # number of free parameters
  N <- sum(vapply(model$models, FUN.VALUE = numeric(1L), FUN = \(sub) sub$data$n))
  G <- model$info$n.groups
  R <- min(R.max, N * S)
  R <- R - R %% G # make R divisble by the number of groups
  R.g <- R / G
  warnif(R.max <= N, "R.max is less than N!")

  ovs <- colnames(model$models[[1L]]$data$data.full)

  if (parametric) {
    # final model (without SEs)
    finalModel <- getFinalModel(model = model, theta = theta, method = "lms")

    for (g in seq_len(model$info$n.groups)) {
      parTable.g <- modelToParTable(finalModel$models[[g]], method = "lms")
      sample.g   <- simulateDataParTable(parTable.g, N = R.g, colsOVs = ovs)$OV[[1L]]

      model$models[[g]]$data <- patternizeMissingDataFIML(sample.g)
    }

  } else for (g in seq_len(model$info$n.groups)) {
    data.g   <- model$models[[g]]$data$data.full
    sample.g <- data.g[sample(NROW(data.g), R.g, replace = TRUE), , drop = FALSE]
    model$models[[g]]$data <- patternizeMissingDataFIML(sample.g)
  }

  popEstep <- estepLms(model      = model,
                       theta      = theta,
                       recalcQuad = TRUE,
                       lastQuad   = if(!is.null(P)) P$quad else NULL)

  suppressWarnings({
    J <- gradientObsLogLikLms_i(theta   = theta,
                                model   = model,
                                P       = popEstep,
                                sign    = +1,
                                epsilon = epsilon)      # R × k matrix
  })

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

  I / S
}


calcEFIM_QML <- function(model, theta, data, S = 100,
                         parametric = TRUE, epsilon = 1e-8, verbose = FALSE,
                         R.max = 1e6) {
  k <- length(theta) # number of free parameters
  N <- sum(vapply(model$models, FUN.VALUE = numeric(1L), FUN = \(sub) sub$data$n))
  G <- model$info$n.groups
  R <- min(R.max, N * S)
  R <- R - R %% G # make R divisble by the number of groups
  R.g <- R / G

  warnif(R.max <= N, "R.max is less than N!")

  ovs <- colnames(model$models[[1L]]$data$data.full)

  if (parametric) {
    # final model (without SEs)
    finalModel <- getFinalModel(model = model, theta = theta, method = "qml")

    for (g in seq_len(model$info$n.groups)) {
      parTable.g <- modelToParTable(finalModel$models[[g]], method = "qml")
      sample.g   <- simulateDataParTable(parTable.g, N = R.g, colsOVs = ovs)$OV[[1L]]

      model$models[[g]]$data <- patternizeMissingDataFIML(sample.g)
    }

  } else for (g in seq_len(model$info$n.groups)) {
    data.g   <- model$models[[g]]$data$data.full
    sample.g <- data.g[sample(NROW(data.g), R.g, replace = TRUE), , drop = FALSE]
    model$models[[g]]$data <- patternizeMissingDataFIML(sample.g)
  }

  for (g in seq_along(model$models)) {
    if (!is.null(model$models[[g]]$matrices$fullU)) {
      fullU     <- model$models[[g]]$matrices$fullU
      fullU_New <- fullU[rep(seq_len(NROW(fullU)), length.out = R.g), , drop = FALSE]

      model$models[[g]]$matrices$fullU <- fullU_New
    }
  }

  suppressWarnings({
    J <- gradientLogLikQml_i(theta = theta, model = model, sign = +1,
                             epsilon = epsilon)
  })

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

  I / S
}


calcOFIM_QML <- function(model, theta,
                         hessian = FALSE,
                         epsilon = 1e-8,
                         robust.se = FALSE,
                         cluster   = NULL,
                         cr1s      = TRUE) {
  N <- sum(vapply(model$models, FUN.VALUE = numeric(1L), FUN = \(sub) NROW(sub$data)))

  if (hessian) {
    # negative hessian (sign = -1)
    I <- calcHessian(model = model, theta = theta,
                     method = "qml", epsilon = epsilon)
    return(I)
  }

  # S: N x k matrix of individual score contributions (sign = +1 => score)
  S <- suppressWarnings(
    gradientLogLikQml_i(theta, model = model, sign = +1, epsilon = epsilon)
  )

  if (!robust.se || is.null(cluster)) {
    # classic OFIM (BHHH / OPG)
    return(crossprod(S))
  }

  stopif(length(cluster) != nrow(S),
         "Length of 'cluster' must equal the number of rows in the data / scores.")

  f <- as.factor(cluster)
  G <- nlevels(f)
  k <- ncol(S)

  # s_g = sum_{i in g} s_i
  Sg <- matrix(0, nrow = G, ncol = k)
  lev <- levels(f)
  for (g in seq_len(G)) {
    idx <- which(f == lev[g])
    Sg[g, ] <- colSums(S[idx, , drop = FALSE])
  }

  B <- crossprod(Sg)  # meat = sum_g s_g s_g'

  # Optional CR1S small-sample correction
  if (isTRUE(cr1s)) {
    q <- ncol(S)
    if (G > 1 && N > q)
      B <- B * (G / (G - 1)) * ((N - 1) / (N - q))
  }

  B
}


getSE_Model <- function(model, se, method, n.additions) {
  params <- model$params

  for (g in seq_len(model$info$n.groups)) {
    SELECT_THETA_LAB  <- params$SELECT_THETA_LAB[[g]]
    SELECT_THETA_COV  <- params$SELECT_THETA_COV[[g]]
    SELECT_THETA_MAIN <- params$SELECT_THETA_MAIN[[g]]

    SELECT_THETA_LAB  <- seq_len(MAX(SELECT_THETA_LAB) + n.additions)
    SELECT_THETA_COV  <- SELECT_THETA_COV  + n.additions
    SELECT_THETA_MAIN <- SELECT_THETA_MAIN + n.additions

    params$SELECT_THETA_LAB[[g]]  <- SELECT_THETA_LAB
    params$SELECT_THETA_COV[[g]]  <- SELECT_THETA_COV
    params$SELECT_THETA_MAIN[[g]] <- SELECT_THETA_MAIN
  }

  model$params <- params

  fillModel(replaceNonNaModelMatrices(model, value = -999),
            theta = se, method = method)
}
