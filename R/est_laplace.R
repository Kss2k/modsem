modsemTMB_MakeADFun <- function(...) {
  msg <- "The `modsemTMB` package is needed to use the LMS-LAPLACE estimator!"

  if (!requireNamespace("modsemTMB", quietly = TRUE)) {
    printf()
    printf("Do you want to install it? (y/n) ")
    choice <- tolower(substr(readLines(n = 1L), 1L, 1L))

    stopifnot(choice == "y")
    remotes::install_github("kss2k/modsemTMB")
  }

  if (requireNamespace("modsemTMB", quietly = TRUE)) { # Make R CMD check happy
    modsemTMB::MakeADFun(...)
  } else mod_msg_stop(msg)
}


boundLaplaceThresholds <- function(model, theta, threshMat_idx, lower, upper,
                                   g = 1L) {
  filled <- fillModel(model, theta, method = "laplace", fillPhi = FALSE)
  threshMat <- filled$models[[g]]$matrices$threshMat
  if (!length(threshMat)) return(list(lower = lower, upper = upper))

  idx <- matrix(threshMat_idx, nrow = NROW(threshMat), ncol = NCOL(threshMat))
  for (i in seq_len(NROW(threshMat))) {
    for (j in seq_len(NCOL(threshMat))) {
      k <- idx[i, j] + 1L
      if (k <= 0L) next

      current <- threshMat[i, j]
      previous <- threshMat[i, j - 1L]
      following <- threshMat[i, j + 1L]

      lower[[k]] <- max(lower[[k]],
                        if (is.finite(previous)) mean(c(previous, current)) else -Inf)
      upper[[k]] <- min(upper[[k]],
                        if (is.finite(following)) mean(c(current, following)) else Inf)
    }
  }

  list(lower = lower, upper = upper)
}


# Build the TMB data list for a single group using a "probe" strategy:
# fill the model with sentinel values (1e6 + 0-based theta index) and read
# back which matrix positions each theta controls.  This is robust to all of
# modsem's naming conventions (position-encoded names, single-element specials,
# equality constraints, etc.).
buildTMBDataGroup <- function(model, g = 1L) {
  submodel <- model$models[[g]]
  M_orig   <- submodel$matrices
  info     <- submodel$data
  data     <- submodel$data

  nXi     <- submodel$info$numXis
  nEta    <- submodel$info$numEtas
  nLatent <- nXi + nEta
  n_theta <- length(model$theta)

  # Probe fillModel with sentinel values
  SENTINEL_BASE <- 1e6
  theta_probe   <- seq(SENTINEL_BASE, SENTINEL_BASE + n_theta - 1L)
  names(theta_probe) <- names(model$theta)

  model_filled <- fillModel(model, theta_probe, method = "laplace", fillPhi = FALSE)
  M_f <- model_filled$models[[g]]$matrices

  # Convert a filled matrix/vector to an integer index array (0-based, -1=fixed)
  asIdx <- function(vals) {
    vals <- as.vector(vals)
    vapply(vals, function(v) {
      if (is.na(v) || !is.finite(v)) return(-1L)
      k <- round(v) - SENTINEL_BASE
      if (k >= 0L && k < n_theta) as.integer(k) else -1L
    }, integer(1L))
  }

  # For symmetric / lower-triangular matrices: only keep lower triangle
  asIdxLTri <- function(mat) {
    n   <- nrow(mat)
    idx <- rep(-1L, length(mat))
    for (j in seq_len(n))
      for (i in j:n) {
        pos <- (j - 1L) * n + i
        v   <- mat[i, j]
        if (!is.na(v) && is.finite(v)) {
          k <- round(v) - SENTINEL_BASE
          if (k >= 0L && k < n_theta) idx[pos] <- as.integer(k)
        }
      }
    idx
  }

  # Fixed skeleton: original matrix with free (NA) positions zeroed out
  skel <- function(mat) { out <- mat; out[is.na(out)] <- 0; out }
  threshSkel <- function(mat) {
    out <- skel(mat)
    out[is.infinite(out) & out < 0] <- -1e30
    out[is.infinite(out) & out > 0] <-  1e30
    out
  }

  # Observation matrix
  nInd <- NROW(M_orig$lambdaX)
  n    <- submodel$data$n

  Y_full <- matrix(NA_real_, nrow = n, ncol = nInd)
  rownames(Y_full) <- seq_len(n)
  colnames(Y_full) <- rownames(M_orig$lambdaX)

  for (p in seq_along(submodel$data$data.split)) {
    dp      <- submodel$data$data.split[[p]]
    obs_col <- which(submodel$data$patterns[p, , drop = TRUE])
    row_idx <- unlist(submodel$data$rowidx[[p]])
    Y_full[row_idx, obs_col] <- dp
  }

  obs_list  <- lapply(seq_len(n), function(i) which(!is.na(Y_full[i, ])) - 1L)
  n_obs_vec <- vapply(obs_list, length, integer(1L))
  max_nObs  <- max(n_obs_vec, 1L)
  obs_idx   <- matrix(-1L, nrow = n, ncol = max_nObs)
  for (i in seq_len(n))
    if (n_obs_vec[i] > 0L)
      obs_idx[i, seq_len(n_obs_vec[i])] <- obs_list[[i]]

  # Ordinal info
  indicator.names <- rownames(M_orig$lambdaX)
  is_ordinal <- as.integer(indicator.names %in% submodel$info$ordered)
  thresh_n <- integer(nInd)
  if (length(submodel$info$ordered)) {
    thresh_n <- vapply(seq_len(nInd), function(j) {
      if (!is_ordinal[[j]]) return(0L)
      as.integer(sum(!is.nan(M_orig$threshMat[j, ])))
    }, integer(1L))
  }

  # Index arrays (via probe)
  Lambda_idx      <- asIdx(M_f$lambdaX)
  Theta_idx       <- asIdxLTri(M_f$thetaDelta)
  tau_idx         <- asIdx(M_f$tauX)
  GammaXi_idx     <- asIdx(M_f$gammaXi)
  GammaEta_idx    <- asIdx(M_f$gammaEta)
  OmegaXiXi_idx   <- asIdx(M_f$omegaXiXi)
  OmegaEtaXi_idx  <- asIdx(M_f$omegaEtaXi)
  alpha_idx       <- asIdx(M_f$alpha)
  beta0_idx       <- asIdx(M_f$beta0)
  A_lower_idx     <- asIdxLTri(M_f$A)
  psi_idx         <- asIdxLTri(M_f$psi)
  covZetaXi_idx   <- asIdx(M_f$covZetaXi)
  threshMat_idx   <- asIdx(M_f$threshMat)

  list(
    nXi  = nXi,
    nEta = nEta,
    nInd = nInd,

    Y       = Y_full,
    obs_idx = obs_idx,
    n_obs   = n_obs_vec,

    is_ordinal = is_ordinal,
    thresh_n   = thresh_n,

    Lambda_f      = skel(M_orig$lambdaX),
    Theta_f       = skel(M_orig$thetaDelta),
    tau_f         = drop(skel(M_orig$tauX)),
    GammaXi_f     = skel(M_orig$gammaXi),
    GammaEta_f    = skel(M_orig$gammaEta),
    OmegaXiXi_f   = skel(M_orig$omegaXiXi),
    OmegaEtaXi_f  = skel(M_orig$omegaEtaXi),
    alpha_f       = drop(skel(M_orig$alpha)),
    beta0_f       = drop(skel(M_orig$beta0)),
    A_lower_f     = skel(M_orig$A),
    psi_f         = skel(M_orig$psi),
    covZetaXi_f   = skel(M_orig$covZetaXi),
    threshMat_f   = threshSkel(M_orig$threshMat),

    Lambda_idx      = Lambda_idx,
    Theta_idx       = Theta_idx,
    tau_idx         = tau_idx,
    GammaXi_idx     = GammaXi_idx,
    GammaEta_idx    = GammaEta_idx,
    OmegaXiXi_idx   = OmegaXiXi_idx,
    OmegaEtaXi_idx  = OmegaEtaXi_idx,
    alpha_idx       = alpha_idx,
    beta0_idx       = beta0_idx,
    A_lower_idx     = A_lower_idx,
    psi_idx         = psi_idx,
    covZetaXi_idx   = covZetaXi_idx,
    threshMat_idx   = threshMat_idx
  )
}


buildLaplaceADFunctions <- function(model, theta, verbose = interactive(), eps = 1e-4) {
  resetOptimizerInfo()

  G <- length(model$models)

  mod_stopif(G > 1,
    "Laplace Approximation is not implemented for multigroup models (yet)!"
  )

  objective <- NULL
  gradient  <- NULL
  hessian   <- NULL
  par       <- NULL

  for (g in seq_len(G)) {
    submodel <- model$models[[g]]
    n        <- submodel$data$n
    nLatent  <- submodel$info$numXis + submodel$info$numEtas

    tmbData <- buildTMBDataGroup(
      model = model,
      g = g
    )

    tmbPar  <- list(
      theta = unname(theta),
      Zeta  = matrix(0.0, nrow = n, ncol = nLatent)
    )

    obj.g <- modsemTMB_MakeADFun(
      data       = tmbData,
      parameters = tmbPar,
      random     = "Zeta"
    )

    objective <- obj.g$fn
    gradient  <- obj.g$gr
    hessian   <- obj.g$hessian
    par       <- obj.g$par
  }

  OBJECTIVE <- function(par) {
    nll <- objective(par) # negative log likelihood
    if (verbose) incrementIterations(-nll)
    nll
  }

  GRADIENT <- function(par) {
    gradient(par)
  }

  HESSIAN <- function(par) {
    # Hessian calculation currently doesn't work with random effects
    H <- tryCatch(obj$he(par), error = function(e) NULL)

    if (is.null(H) || !all(is.finite(H))) {
      p <- length(par)
      H <- matrix(NA_real_, p, p)

      grad0 <- gradient(par)

      for (i in seq_len(p)) {
        pari <- par
        pari[[i]] <- pari[[i]] + eps

        gradi <- gradient(pari)
        H[,i] <- (gradi-grad0)/eps
      }
    }

    H
  }

  list(
    objective = OBJECTIVE,
    gradient  = GRADIENT,
    hessian   = HESSIAN,
    par       = par,
    tmbData   = tmbData,
    tmbPar    = tmbPar
  )
}


# Observed FIM for the Laplace method via the TMB exact Hessian.
# Called by calcFIM_da in calc_se_da.R.
calcOFIM_Laplace <- function(model, theta, ..., eps = 1e-4) {
  AD <- buildLaplaceADFunctions(model = model, theta = theta, eps = eps)

  H <- tryCatch({
    AD$hessian(unname(theta))
  }, error = function(e) {
    mod_msg_warn("Calculation of Hessian failed!")
      p <- length(theta)
      matrix(NA_real_, p, p)
  })

  dimnames(H) <- list(names(theta), names(theta))
  H
}


# Main Laplace estimator.
estLaplace <- function(model,
                       convergence  = 1e-5,
                       verbose      = FALSE,
                       calc.se      = TRUE,
                       epsilon      = 1e-6,
                       optimizer    = "nlminb",
                       inner.control = list(),
                       iter.max     = 1000,
                       ...) {
  mod_stopif(model$info$n.groups > 1L,
    "method=\"laplace\" currently supports single-group models only.")

  theta     <- model$theta
  data      <- lapply(model$models, FUN = \(sub) sub$data)
  submodel  <- model$models[[1L]]
  nXi       <- submodel$info$numXis
  nEta      <- submodel$info$numEtas
  nLatent   <- nXi + nEta
  n         <- submodel$data$n

  AD <- buildLaplaceADFunctions(
    model = model, theta = theta, verbose = verbose
  )

  bounds <- model$params$bounds
  lower  <- unname(bounds$lower)
  upper  <- unname(bounds$upper)
  bounds <- boundLaplaceThresholds(
    model = model, theta = theta,
    threshMat_idx = AD$tmbData$threshMat_idx,
    lower = lower, upper = upper
  )

  lower <- bounds$lower
  upper <- bounds$upper

  if (optimizer == "nlminb") {
    control <- list(
      rel.tol  = convergence,
      iter.max = iter.max,
      eval.max = 5 * iter.max
    )

    suppressWarnings({
      opt <- stats::nlminb(
        start     = AD$par,
        objective = AD$objective,
        gradient  = AD$gradient,
        lower     = lower,
        upper     = upper,
        control   = control
      )
    })

    converged  <- opt$convergence == 0L
    iterations <- opt$iterations

  } else if (optimizer == "L-BFGS-B") {
    opt <- stats::optim(
      par    = AD$par,
      fn     = AD$objective,
      gr     = AD$gradient,
      method = "L-BFGS-B",
      lower  = lower,
      upper  = upper,
      control = list(factr = convergence / .Machine$double.eps)
    )
    converged  <- opt$convergence == 0L
    iterations <- opt$counts[["function"]]
    opt$objective <- opt$value

  } else {
    mod_msg_stop("optimizer must be 'nlminb' or 'L-BFGS-B'.")
  }

  if (verbose) printf("\n")

  theta_hat     <- structure(opt$par, names = names(theta))
  logLik_hat    <- -opt$objective

  finalizeModelEstimatesDA(
    model             = model,
    theta             = theta_hat,
    method            = "laplace",
    data              = data,
    logLik            = logLik_hat,
    iterations        = iterations,
    converged         = converged,
    optimizer         = optimizer,
    calc.se           = calc.se,
    FIM               = "observed",
    OFIM.hessian      = FALSE,
    EFIM.S            = 3e4,
    EFIM.parametric   = TRUE,
    robust.se         = FALSE,
    epsilon           = epsilon,
    cr1s              = FALSE,
    R.max             = 1e6,
    verbose           = verbose,
    includeStartModel = TRUE,
    startModel        = model
  )
}
