# Laplace estimator for modsem DA models (TMB-based).
#
# The TMB template marginalises the per-observation latent variables (Zeta)
# via the Laplace approximation, providing exact AD gradients for the outer
# optimisation over structural parameters and an exact Hessian for SEs.

TMB_DLL_NAME <- "modsem_laplace"
TMB_DLL_ENV  <- new.env(parent = emptyenv())
TMB_DLL_ENV$loaded <- FALSE


# Compile and dyn.load the TMB template (once per session).
loadTMBModel <- function() {
  if (isTRUE(TMB_DLL_ENV$loaded))
    return(invisible(NULL))

  src <- system.file("tmb", paste0(TMB_DLL_NAME, ".cpp"), package = "modsem")
  if (!nzchar(src))
    mod_msg_stop(
      "TMB template not found. Re-install modsem or check inst/tmb/."
    )

  # Compile into a temp directory to avoid conflicts with the main Rcpp DSO
  work_dir <- file.path(tempdir(), "modsem_tmb")
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(src, file.path(work_dir, paste0(TMB_DLL_NAME, ".cpp")),
            overwrite = TRUE)

  old_wd <- setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)

  TMB::compile(paste0(TMB_DLL_NAME, ".cpp"))
  dyn.load(TMB::dynlib(TMB_DLL_NAME))
  TMB_DLL_ENV$loaded <- TRUE

  invisible(NULL)
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

  # ---- Probe fillModel with sentinel values --------------------------------
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

  # ---- Observation matrix ---------------------------------------------------
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

  # ---- Ordinal info ---------------------------------------------------------
  ordinal    <- M_orig$ordinal
  is_ordinal <- as.integer(ordinal$isOrdinal)

  thresh_n <- integer(nInd)
  for (j in seq_len(nInd)) {
    t_j <- ordinal$thresholds[[j]]
    if (!is.null(t_j) && length(t_j) > 0L && is_ordinal[j])
      thresh_n[j] <- length(t_j)
  }
  max_thresh <- max(thresh_n, 1L)
  thresh_mat <- matrix(0.0, nrow = nInd, ncol = max_thresh)
  for (j in seq_len(nInd)) {
    if (thresh_n[j] > 0L) {
      t_j <- ordinal$thresholds[[j]]
      t_j[is.infinite(t_j) & t_j < 0] <- -1e30
      t_j[is.infinite(t_j) & t_j > 0] <-  1e30
      thresh_mat[j, seq_len(thresh_n[j])] <- t_j
    }
  }

  # ---- Index arrays (via probe) --------------------------------------------
  Lambda_idx      <- asIdx(M_f$lambdaX)
  Theta_diag_idx  <- asIdx(diag(M_f$thetaDelta))
  tau_idx         <- asIdx(M_f$tauX)
  GammaXi_idx     <- asIdx(M_f$gammaXi)
  GammaEta_idx    <- asIdx(M_f$gammaEta)
  OmegaXiXi_idx  <- asIdx(M_f$omegaXiXi)
  OmegaEtaXi_idx <- asIdx(M_f$omegaEtaXi)
  alpha_idx       <- asIdx(M_f$alpha)
  beta0_idx       <- asIdx(M_f$beta0)
  A_lower_idx     <- asIdxLTri(M_f$A)
  psi_idx         <- asIdxLTri(M_f$psi)
  covZetaXi_idx  <- asIdx(M_f$covZetaXi)

  list(
    nXi  = nXi,
    nEta = nEta,
    nInd = nInd,

    Y       = Y_full,
    obs_idx = obs_idx,
    n_obs   = n_obs_vec,

    is_ordinal = is_ordinal,
    thresh_n   = thresh_n,
    thresh_mat = thresh_mat,

    Lambda_f      = skel(M_orig$lambdaX),
    Theta_diag_f  = diag(skel(M_orig$thetaDelta)),
    tau_f         = drop(skel(M_orig$tauX)),
    GammaXi_f     = skel(M_orig$gammaXi),
    GammaEta_f    = skel(M_orig$gammaEta),
    OmegaXiXi_f  = skel(M_orig$omegaXiXi),
    OmegaEtaXi_f = skel(M_orig$omegaEtaXi),
    alpha_f       = drop(skel(M_orig$alpha)),
    beta0_f       = drop(skel(M_orig$beta0)),
    A_lower_f     = skel(M_orig$A),
    psi_f         = skel(M_orig$psi),
    covZetaXi_f  = skel(M_orig$covZetaXi),

    Lambda_idx      = Lambda_idx,
    Theta_diag_idx  = Theta_diag_idx,
    tau_idx         = tau_idx,
    GammaXi_idx     = GammaXi_idx,
    GammaEta_idx    = GammaEta_idx,
    OmegaXiXi_idx  = OmegaXiXi_idx,
    OmegaEtaXi_idx = OmegaEtaXi_idx,
    alpha_idx       = alpha_idx,
    beta0_idx       = beta0_idx,
    A_lower_idx     = A_lower_idx,
    psi_idx         = psi_idx,
    covZetaXi_idx  = covZetaXi_idx
  )
}


# Observed FIM for the Laplace method via the TMB exact Hessian.
# Called by calcFIM_da in calc_se_da.R.
calcOFIM_Laplace <- function(model, theta, ..., eps = 1e-4) {
  loadTMBModel()
  submodel <- model$models[[1L]]
  n        <- submodel$data$n
  nLatent  <- submodel$info$numXis + submodel$info$numEtas

  tmb_data <- buildTMBDataGroup(model, g = 1L)
  obj <- TMB::MakeADFun(
    data       = tmb_data,
    parameters = list(theta = unname(theta),
                      Zeta  = matrix(0.0, nrow = n, ncol = nLatent)),
    random     = "Zeta",
    DLL        = TMB_DLL_NAME,
    silent     = TRUE
  )

  H <- tryCatch(obj$he(unname(theta)), error = function(e) NULL)

  if (is.null(H) || !all(is.finite(H))) {
    p <- length(theta)
    H <- matrix(NA_real_, p, p,
                dimnames = list(names(theta), names(theta)))

    grad0 <- obj$gr(unname(theta))
    for (i in seq_len(p)) {
      thetai <- unname(theta)
      thetai[[i]] <- thetai[[i]] + eps
      gradi <- obj$gr(thetai)
      H[,i] <- (gradi-grad0)/eps
    }

    return(H)
  }

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
                       ...) {
  mod_stopif(model$info$n.groups > 1L,
    "method=\"laplace\" currently supports single-group models only.")

  loadTMBModel()

  theta     <- model$theta
  data      <- lapply(model$models, FUN = \(sub) sub$data)
  submodel  <- model$models[[1L]]
  nXi       <- submodel$info$numXis
  nEta      <- submodel$info$numEtas
  nLatent   <- nXi + nEta
  n         <- submodel$data$n

  tmb_data <- buildTMBDataGroup(model, g = 1L)

  tmb_par <- list(
    theta = unname(theta),
    Zeta  = matrix(0.0, nrow = n, ncol = nLatent)
  )

  obj <- TMB::MakeADFun(
    data       = tmb_data,
    parameters = tmb_par,
    random     = "Zeta",
    DLL        = TMB_DLL_NAME,
    silent     = !verbose
  )

  bounds <- model$params$bounds
  lower  <- unname(bounds$lower)
  upper  <- unname(bounds$upper)

  if (verbose) printf("Optimising (method=laplace)\n")

  if (optimizer == "nlminb") {
    control <- list(
      rel.tol  = convergence,
      iter.max = 1000L,
      eval.max = 5000L
    )
    opt <- stats::nlminb(
      start     = obj$par,
      objective = obj$fn,
      gradient  = obj$gr,
      lower     = lower,
      upper     = upper,
      control   = control
    )
    converged  <- opt$convergence == 0L
    iterations <- opt$iterations

  } else if (optimizer == "L-BFGS-B") {
    opt <- stats::optim(
      par    = obj$par,
      fn     = obj$fn,
      gr     = obj$gr,
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
