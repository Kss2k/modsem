getExpectedOrdinalValues <- function(data,
                                     fit,
                                     R = 1e5,
                                     ordered = NULL,
                                     tol = 1e-5,
                                     verbose = TRUE,
                                     num_iter = 2L # Gauss–Seidel sweeps over ordinal items
) {

  # --- 0) Prep: column order, types, and missingness -------------------------
  ovs   <- colnames(modsem_inspect(fit, what = "cov.ov"))
  data  <- as.data.frame(data)[ovs]
  re.ordinalize <- \(x) as.ordered(as.integer(as.ordered(x)))
  if (!is.null(ordered)) data[ordered] <- lapply(data[ordered], re.ordinalize)

  o.cols <- names(data)[vapply(data, is.ordered, logical(1))]
  c.cols <- setdiff(names(data), o.cols)

  missing <- !complete.cases(data)
  if (any(missing)) {
    warning("Missing values are not yet supported; dropping incomplete cases for E-step.")
    data <- data[!missing, , drop = FALSE]
  }

  if (!length(o.cols)) {
    if (isTRUE(verbose)) cat("No ordinal columns; returning numeric copy.\n")
    return(as.matrix(data))
  }

  # --- 1) Model / quadrature (node-wise means & covariances) -----------------
  model <- fit$model
  quad  <- model$quad
  V     <- quad$n            # quadrature nodes (rows)
  w     <- quad$w
  k     <- length(w)

  MU    <- vector("list", length = k)
  SIGMA <- vector("list", length = k)
  for (i in seq_len(k)) {
    MU[[i]]    <- as.numeric(muLmsCpp(model = model, z = V[i, ]))
    SIGMA[[i]] <- as.matrix(sigmaLmsCpp(model = model, z = V[i, ]))
  }

  # --- 2) Keep your simulation-based rescaling EXACTLY as-is -----------------
  inds.xis <- unique(model$info$allIndsXis)
  sim_i    <- simulateDataParTable(parTable = standardized_estimates(fit), N = R)

  rescaled <- rescaleOrderedData(
    data         = data,
    sim.ov       = sim_i$oV,
    cols.ordered = o.cols,
    cols.cont    = c.cols,
    linear.ovs   = inds.xis
  )
  thresholds <- rescaled$thresholds  # named list: each is vector of internal cut-points

  # --- 3) Encode data matrix X; ordinal columns become integer levels --------
  X <- data
  if (length(o.cols)) {
    X[o.cols] <- lapply(X[o.cols], function(x) as.integer(x))
  }
  X <- as.matrix(X)
  rownames(X) <- rownames(data)

  # Output container (will overwrite only ordinal columns)
  Y <- X
  o.idx <- match(o.cols, colnames(X))
  p     <- ncol(X)

  # --- 4) Helpers ------------------------------------------------------------

  # (a) Stable 1D truncated normal mean
  trunc_mean_1d <- function(m, s, a, b) {
    s <- max(s, 1e-12)
    alpha <- (a - m)/s; beta <- (b - m)/s
    Z <- pnorm(beta) - pnorm(alpha)
    if (!is.finite(Z) || Z <= 0) return(m) # fallback
    m + s * (dnorm(alpha) - dnorm(beta)) / Z
  }

  # (b) Conditional Gaussian for scalar j given the rest r, under node t
  mvnorm_cond_scalar <- function(mu, Sigma, idx_j, idx_r, x_r) {
    S_rr <- Sigma[idx_r, idx_r, drop = FALSE]
    S_jr <- Sigma[idx_j, idx_r, drop = FALSE]
    S_rj <- Sigma[idx_r, idx_j, drop = FALSE]
    S_jj <- Sigma[idx_j, idx_j, drop = FALSE]

    # Cholesky solve for stability
    # (add tiny ridge if needed)
    diag(S_rr) <- pmax(diag(S_rr), 1e-10)
    L <- chol(S_rr)
    diff <- x_r - mu[idx_r]
    # A' = S_jr * S_rr^{-1}; compute via solves
    A   <- backsolve(L, t(S_jr), upper.tri = TRUE)
    A   <- backsolve(t(L), A, upper.tri = FALSE)

    m   <- as.numeric(mu[idx_j] + crossprod(A, diff))
    s2  <- as.numeric(S_jj - crossprod(A, S_rj))
    c(m = m, s = sqrt(pmax(1e-12, s2)))
  }

  # (c) Per-row deterministic E-step for all ordinal columns
  expected_ordinals_row <- function(x_row) {
    cur <- x_row

    for (it in seq_len(num_iter)) {
      for (j_local in seq_along(o.idx)) {
        j <- o.idx[j_local]
        colname <- o.cols[j_local]

        # Observed category 1..K and thresholds (add -Inf, +Inf)
        catj <- as.integer(x_row[j])
        tau  <- c(-Inf, thresholds[[colname]], Inf)
        a <- tau[catj]; b <- tau[catj + 1L]
        if (!is.finite(a) && !is.finite(b)) { next } # degenerate, skip

        idx_j <- j
        idx_r <- setdiff(seq_len(p), j)

        # Node posteriors (up to a constant) and node-wise E[c_j | r, y_j, t]
        log_post <- numeric(k)
        E_j_t    <- numeric(k)

        for (t in seq_len(k)) {
          mu_t <- MU[[t]]; S_t <- SIGMA[[t]]

          # Conditional for c_j | r under node t (treat cur[idx_r] as "observed")
          cs  <- mvnorm_cond_scalar(mu_t, S_t, idx_j, idx_r, cur[idx_r])
          m_j <- cs["m"]; s_j <- cs["s"]

          # Univariate interval prob for this ordinal under node t
          pij <- pnorm((b - m_j)/s_j) - pnorm((a - m_j)/s_j)
          pij <- max(pij, 1e-300)

          # Conditional density of r given node t: N(cur_r | mu_r, S_rr)
          S_rr <- S_t[idx_r, idx_r, drop = FALSE]
          diag(S_rr) <- pmax(diag(S_rr), 1e-10)
          ll_r <- mvtnorm::dmvnorm(cur[idx_r], mean = mu_t[idx_r], sigma = S_rr, log = TRUE)

          log_post[t] <- log(w[t]) + ll_r + log(pij)
          E_j_t[t]    <- trunc_mean_1d(m_j, s_j, a, b)
        }

        # Normalize weights and mix
        m  <- max(log_post)
        wt <- exp(log_post - m); wt <- wt / sum(wt)
        cur[j] <- sum(wt * E_j_t)
      }
    }
    cur[o.idx]
  }

  # --- 5) Main loop over rows ------------------------------------------------
  N <- nrow(X)
  if (N == 0L) {
    if (verbose) cat("\n")
    return(Y)
  }

  if (isTRUE(verbose)) {
    cat("E-step (deterministic mixture conditional expectations)…\n")
  }

  tick_every <- max(1L, floor(N / 50L))
  for (n in seq_len(N)) {
    if (verbose && (n %% tick_every == 0L || n == 1L || n == N)) {
      cat(sprintf("\rProgress: %d/%d", n, N)); utils::flush.console()
    }
    Y[n, o.idx] <- expected_ordinals_row(X[n, ])
  }

  if (verbose) cat("\n")
  Y
}



getExpectedOrdinalValuesOld <- function(data,
                                     fit,
                                     R = 1e5,
                                     ordered = NULL,
                                     tol = 1e-5,
                                     verbose = TRUE) {
  # 0) Prep: column order, types, and missingness
  ovs   <- colnames(modsem_inspect(fit, what = "cov.ov"))
  data  <- as.data.frame(data)[ovs]
  re.ordinalize <- \(x) as.ordered(as.integer(as.ordered(x))) # levels [2, 3, 4] -> [1, 2, 3] # important
  if (!is.null(ordered)) data[ordered] <- lapply(data[ordered], re.ordinalize)

  o.cols <- names(data)[vapply(data, is.ordered, logical(1))]
  c.cols <- setdiff(names(data), o.cols)

  missing <- !complete.cases(data)
  if (any(missing)) {
    warning("Missing values are not yet supported, imputing ordinal values (dropping incomplete cases).")
    data <- data[!missing, , drop = FALSE]
  }

  # 1) Model / quadrature
  model <- fit$model
  quad  <- model$quad
  V     <- quad$n            # quadrature nodes (rows)
  w     <- quad$w
  k     <- length(w)

  SIGMA <- vector("list", length = k)
  MU    <- vector("list", length = k)
  for (i in seq_len(k)) {
    MU[[i]]    <- muLmsCpp(model = model, z = V[i, ])
    SIGMA[[i]] <- sigmaLmsCpp(model = model, z = V[i, ])
  }

  # 2) Keep your 100k simulation-based rescaling EXACTLY as-is
  inds.xis <- unique(model$info$allIndsXis)
  sim_i <- simulateDataParTable(parTable = standardized_estimates(fit), N = R)
  rescaled <- rescaleOrderedData(
    data         = data,
    sim.ov       = sim_i$oV,
    cols.ordered = o.cols,
    cols.cont    = c.cols,
    linear.ovs   = inds.xis
  )
  thresholds <- rescaled$thresholds  # named list: each is vector of internal cutpoints

  # 3) Encode data matrix X; ordinal columns become integer levels
  X <- data
  if (length(o.cols)) {
    X[o.cols] <- lapply(X[o.cols], function(x) as.integer(x))
  }
  X <- as.matrix(X)
  Y <- X  # output; will overwrite only ordinal columns

  # Positions of ordinal columns to avoid repeated matching
  o.idx <- match(o.cols, colnames(X))
  p     <- ncol(X)

  # 4) Precompute MVN constants (Cholesky + logdet) for all nodes
  cholL  <- vector("list", k)
  logdet <- numeric(k)
  for (i in seq_len(k)) {
    L <- chol(SIGMA[[i]])
    cholL[[i]] <- L
    logdet[i]  <- 2 * sum(log(diag(L)))
  }
  const <- p * log(2 * pi)
  logw  <- log(w)

  # Small, dependency-free log-sum-exp
  logsumexp <- function(a) {
    m <- max(a)
    if (!is.finite(m)) return(m) # all -Inf
    m + log(sum(exp(a - m)))
  }

  # Fast log mixture density: log sum_i w_i * N(x | MU_i, SIGMA_i)
  llRow <- function(x) {
    # x is full-length row vector, length p
    z <- numeric(k)
    for (i in seq_len(k)) {
      r <- x - MU[[i]]
      y <- backsolve(cholL[[i]], r, transpose = TRUE)
      z[i] <- -0.5 * (sum(y * y) + logdet[i] + const)
    }
    logsumexp(z + logw)
  }

  # 5) Prebuild tau-lists (bounds per ordinal column) once
  # tau_j = (-Inf, thresholds, +Inf)
  tau_list <- lapply(thresholds, function(v) c(-Inf, v, Inf))

  # Helper to compute bounds for a given row's ordinal entries in o.idx order
  compute_bounds_for_row <- function(obsn_ord_levels) {
    # obsn_ord_levels is integer vector of same length as o.idx
    lower <- numeric(length(obsn_ord_levels))
    upper <- numeric(length(obsn_ord_levels))
    for (j in seq_along(obsn_ord_levels)) {
      colj <- o.cols[j]
      tau  <- tau_list[[colj]]
      catj <- obsn_ord_levels[j]
      # Assumes levels are 1..K
      lower[j] <- tau[catj]
      upper[j] <- tau[catj + 1L]
    }

    list(lower = lower, upper = upper)
  }

  # 6) Optimization loop (per-row, box-constrained)
  N <- nrow(X)
  if (N == 0L || length(o.idx) == 0L) {
    if (verbose) cat("\n")
    return(Y)
  }

  # progress printer
  if (isTRUE(verbose)) {
    cat("E-step (per-row bounded optimization)…\n")
  }

  # We’ll use L-BFGS-B for speed; pass full row to llRow (BUGFIX)
  for (n in seq_len(N)) {
    if (verbose && (n %% max(1L, floor(N / 50L)) == 0L || n == 1L || n == N)) {
      # ~50 ticks total
      cat(sprintf("\rProgress: %d/%d", n, N))
      utils::flush.console()
    }

    obsn <- X[n, ]                    # full row (numeric)
    ordv <- obsn[o.idx]               # integer category indices for ord cols

    bnd  <- compute_bounds_for_row(ordv)
    lower <- bnd$lower
    upper <- bnd$upper

    # Start at midpoints of the finite part of the interval
    # For infinite endpoints, use NA to compute rowMeans
    tmp_bounds <- cbind(lower, upper)
    tmp_bounds[is.infinite(tmp_bounds)] <- NA_real_
    start <- rowMeans(tmp_bounds, na.rm = TRUE)

    # If a bound is entirely infinite (rare), fall back to zero start for that dim
    start[!is.finite(start)] <- 0

    objective <- function(vals) {
      obsn[o.idx] <- vals        # update current row's ordinal slots
      -llRow(obsn)               # NEGATIVE log-mixture for minimizer
    }

    opt <- try(
      stats::optim(
        par     = start,
        fn      = objective,
        method  = "L-BFGS-B",
        lower   = lower,
        upper   = upper,
        control = list(factr = tol / .Machine$double.eps)
      ),
      silent = TRUE
    )

    if (inherits(opt, "try-error") || !is.finite(opt$value)) {
      browser()
      # Fallback to nlminb if L-BFGS-B had issues (very rarec)
      opt <- stats::nlminb(
        start     = start,
        objective = objective,
        lower     = lower,
        upper     = upper,
        control   = list(rel.tol = tol)
      )
      Y[n, o.idx] <- opt$par
    } else {
      Y[n, o.idx] <- opt$par
    }
  }

  if (verbose) cat("\n")
  Y
}


ordered_lms <- function(model, data, ordered = NULL, tol = 1e-9,
                        # capture args
                        method = NA,
                        MAX_ITER = 100,
                        ...) {
  ovs <- getOVs(modsemify(model))
  data <- data[, ovs]

  isOrdered <- sapply(data, FUN = is.ordered)
  ordered <- union(colnames(data)[isOrdered], ordered)

  if (!length(ordered)) return(modsem(model, data, method = "lms", ...))

  std1 <- \(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

  data.int <- as.data.frame(data)
  data.int[ordered] <- lapply(data.int, FUN = as.integer)
  
  Y <- data.int
  Y[ordered] <- lapply(Y, FUN = std1)
 
  theta0 <- NULL

  for (i in seq_len(MAX_ITER)) {
    printf("Iteration %d...\n", i)
    fit_sam <- lavaan::sam(model, data = Y)
    centered_estimates(fit_sam)
    fit <- modsem(model, data = Y, method = "lms", ...)
    Y <- getExpectedOrdinalValues(data = data, fit = fit, ordered = ordered, tol = tol)
    Y[, ordered] <- apply(Y[, ordered], MARGIN = 2, FUN = std1)

    theta1 <- standardized_estimates(fit)$est

    if (!is.null(theta0)) {
      diff <- mean(abs(theta1 - theta0))
      if (diff < 1e-5) break
    }
  }

  fit.out <- fit
  sim <- simulateDataParTable(parTable = parameter_estiamtes(fit),
                              N        = N)
  rescaled <- rescaleOrderedData(data = data.int,
                                 sim.ov = sim$oV,
                                 cols.ordered = ordered,
                                 cols.cont = setdiff(colnames(data), ordered),
                                 linear.ovs = fit$model$info$allIndsXis,
                                 thresholds.vcov = TRUE)

  vcov.t <- NULL
  coef.t <- NULL
  for (col in cols.ordered) {
    vcov.t <- diagPartitionedMat(vcov.t, rescaled$vcov[[col]])
    coef.t <- c(coef.t, rescaled$thresholds[[col]])
  }

  parTable <- fit.out$parTable
  cols.t <- c("lhs", "op", "rhs", "label", "est", "std.error")
  cols.m <- setdiff(colnames(parTable), cols.t)

  lhs.t = stringr::str_split_i(names(coef.t), pattern = "\\|", i = 1L)
  rhs.t = stringr::str_split_i(names(coef.t), pattern = "\\|", i = 2L)

  parTable.t <- data.frame(lhs = lhs.t, op = "|", rhs = rhs.t, label = "",
                           est = coef.t, std.error = sqrt(diag(vcov.t)))
  parTable.t[cols.m] <- NA
  parTable.t <- addZStatsParTable(parTable.t)

  parTable <- sortParTableDA(rbind(parTable, parTable.t), model = fit.out$model)

  fit.out$coef.all <- c(fit.out$coef.all, coef.t)
  fit.out$vcov.all <- diagPartitionedMat(fit.out$vcov.all, vcov.t)
  fit.out$parTable <- parTable

  fit.out
}
