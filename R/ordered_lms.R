getExpectedOrdinalValues <- function(data,
                                     fit,
                                     R = 1e5,
                                     ordered = NULL,
                                     tol = 1e-5,
                                     verbose = TRUE,
                                     num_iter = 2L,         # Gauss–Seidel sweeps
                                     parallel = TRUE,        # set future::plan() outside
                                     workers = 8) {
  platform <- tolower(.Platform$OS.type)

  if (platform == "unix") {
    future::plan(future::multicore, workers = workers)
  } else if (platform == "windows") {
    future::plan(future::multisession, workers = workers)
  }

  # --- 0) Prep ---------------------------------------------------------------
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

  # --- 1) Quadrature: MU/SIGMA/w & precompute node facts --------------------
  model <- fit$model
  quad  <- model$quad
  V     <- quad$n
  w     <- quad$w
  k     <- length(w)

  MU    <- vector("list", k)
  SIGMA <- vector("list", k)
  CHOL  <- vector("list", k)    # Cholesky of Sigma (for ll_full)
  LOGDET<- numeric(k)
  QINV  <- vector("list", k)    # Precision matrix Sigma^{-1}

  for (i in seq_len(k)) {
    MU[[i]]    <- as.numeric(muLmsCpp(model = model, z = V[i, ]))
    S          <- as.matrix(sigmaLmsCpp(model = model, z = V[i, ]))
    # tiny ridge for safety
    diag(S)    <- pmax(diag(S), 0)
    S          <- S + diag(1e-8, nrow(S))
    SIGMA[[i]] <- S
    L          <- chol(S)
    CHOL[[i]]  <- L
    LOGDET[i]  <- 2 * sum(log(diag(L)))
    # precision via Cholesky: Sigma^{-1} = (L^{-T} L^{-1})
    QINV[[i]]  <- chol2inv(L)
  }
  const <- ncol(data) * log(2 * pi)
  logw  <- log(w)

  # --- 2) Simulation-based thresholds (your original path) -------------------
  inds.xis <- unique(model$info$allIndsXis)
  sim_i    <- simulateDataParTable(parTable = standardized_estimates(fit), N = R)
  rescaled <- rescaleOrderedData(
    data         = data,
    sim.ov       = sim_i$oV,
    cols.ordered = o.cols,
    cols.cont    = c.cols,
    linear.ovs   = inds.xis
  )
  thresholds <- rescaled$thresholds
  tau_list   <- lapply(thresholds, function(v) c(-Inf, v, Inf))  # cache with infinities

  # --- 3) Encode numeric matrix X -------------------------------------------
  X <- data
  X[o.cols] <- lapply(X[o.cols], function(x) as.integer(x))
  X <- as.matrix(X)
  rownames(X) <- rownames(data)

  Y     <- X                            # output
  p     <- ncol(X)
  o.idx <- match(o.cols, colnames(X))

  # --- 4) Helpers (vectorized/fast) -----------------------------------------
  # Truncated-normal mean (robust)
  trunc_mean_1d <- function(m, s, a, b) {
    s <- max(s, 1e-12)
    alpha <- (a - m)/s; beta <- (b - m)/s
    alpha <- pmax(alpha, -1e2); beta <- pmin(beta, 1e2)
    Z <- pnorm(beta) - pnorm(alpha)
    if (!is.finite(Z) || Z <= 0) return(pmin(pmax(m, a), b))
    m + s * (dnorm(alpha) - dnorm(beta)) / Z
  }

  # Conditional via precision: for scalar j | rest, under node with precision Q = Sigma^{-1}
  # m = mu_j - (Q_j,-j / Q_jj) %*% (x_-j - mu_-j),   s2 = 1 / Q_jj
  cond_from_precision <- function(mu, Q, j, x) {
    Qjj <- Q[j, j]
    s2  <- 1 / max(Qjj, 1e-12)
    qj  <- Q[j, -j, drop = FALSE]
    m   <- as.numeric(mu[j] - (qj %*% (x[-j] - mu[-j])) / Qjj)
    c(m = m, s = sqrt(s2))
  }

  # Node-wise full joint log-density for one x (vector) using precomputed CHOL/LOGDET
  ll_full_per_node <- function(x) {
    # returns numeric length-k
    out <- numeric(k)
    for (t in seq_len(k)) {
      r  <- x - MU[[t]]
      y  <- backsolve(CHOL[[t]], r, transpose = TRUE)
      out[t] <- -0.5 * (sum(y * y) + LOGDET[t] + const) + logw[t]
    }
    out
  }

  # One row, all ordinals (Gauss–Seidel sweeps), using precision and ll_full reuse
  expected_ordinals_row <- function(x_row) {
    cur <- x_row

    # init ordinals at interval midpoints
    for (j_local in seq_along(o.idx)) {
      j <- o.idx[j_local]; colname <- o.cols[j_local]
      catj <- as.integer(x_row[j])
      tau  <- tau_list[[colname]]
      a <- tau[catj]; b <- tau[catj + 1L]
      if (is.finite(a) && is.finite(b)) {
        cur[j] <- 0.5 * (a + b)
      } else if (is.infinite(a) && is.finite(b)) {
        cur[j] <- b - 1.0
      } else if (is.finite(a) && is.infinite(b)) {
        cur[j] <- a + 1.0
      } else cur[j] <- 0
    }

    for (it in seq_len(num_iter)) {
      # precompute full-join log-density per node for the current cur (reused across j)
      log_joint <- ll_full_per_node(cur)   # length-k

      for (j_local in seq_along(o.idx)) {
        j <- o.idx[j_local]; colname <- o.cols[j_local]
        catj <- as.integer(x_row[j])
        tau  <- tau_list[[colname]]
        a <- tau[catj]; b <- tau[catj + 1L]

        # node-wise conditional mean/var via precision; also conditional log-density for j
        log_post   <- numeric(k)
        E_j_given  <- numeric(k)

        for (t in seq_len(k)) {
          cs  <- cond_from_precision(MU[[t]], QINV[[t]], j, cur)
          mj  <- cs["m"]; sj <- cs["s"]

          # soft cap to keep numbers sane prior to truncation
          if (is.finite(a) && mj < a) mj <- max(mj, a - 8 * sj)
          if (is.finite(b) && mj > b) mj <- min(mj, b + 8 * sj)

          # univariate interval probability under node t
          pij <- pnorm((b - mj)/sj) - pnorm((a - mj)/sj)
          pij <- max(pij, 1e-300)

          # p(x_-j|t) = p(x|t) / p(x_j | x_-j, t)
          # conditional density of current value cur[j]
          ll_cond_j <- dnorm(cur[j], mean = mj, sd = sj, log = TRUE)

          log_post[t]  <- (log_joint[t] - ll_cond_j) + log(pij) # ∝ log w_t + log p(x_-j|t) + log pij
          E_j_given[t] <- trunc_mean_1d(mj, sj, a, b)
        }

        # normalize weights and mix
        m  <- max(log_post)
        wt <- exp(log_post - m); wt <- wt / sum(wt)
        cur[j] <- sum(wt * E_j_given)
      }
    }
    cur[o.idx]
  }

  # --- 5) Parallel over rows -------------------------------------------------
  N <- nrow(X)
  if (N == 0L) {
    if (verbose) cat("\n")
    return(Y)
  }
  if (isTRUE(verbose)) {
    cat("E-step (parallel mixture conditional expectations)…\n")
  }

  row_fun <- function(n) expected_ordinals_row(X[n, ])

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    # Use existing future::plan() set by the caller
    chunks <- future.apply::future_lapply(seq_len(N), function(n) row_fun(n))
  } else if (parallel && .Platform$OS.type != "windows") {
    # Fallback to mclapply on Unix-alikes
    chunks <- parallel::mclapply(seq_len(N), row_fun, mc.cores = max(1L, parallel::detectCores() - 1L))
  } else {
    # Serial
    if (isTRUE(verbose)) {
      tick_every <- max(1L, floor(N / 50L))
    }
    chunks <- vector("list", N)
    for (n in seq_len(N)) {
      if (verbose && (n %% tick_every == 0L || n == 1L || n == N)) {
        cat(sprintf("\rProgress: %d/%d", n, N)); utils::flush.console()
      }
      chunks[[n]] <- row_fun(n)
    }
    if (verbose) cat("\n")
  }

  # stitch back
  ord_mat <- do.call(rbind, chunks)
  Y[, o.idx] <- ord_mat
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
                        method = NA,
                        MAX_ITER = 100,
                        verbose = TRUE,
                        calc.se = TRUE, # capture
                        ...) {

  # ---- helpers ----
  printf <- function(...) cat(sprintf(...))
  std1   <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  scale_once <- function(df) {
    m <- vapply(df, function(x) mean(x, na.rm = TRUE), 0.0)
    s <- vapply(df, function(x) sd(x,   na.rm = TRUE), 0.0); s[s == 0] <- 1
    list(
      transform = function(df) as.data.frame(mapply(function(x, mu, sd) (x - mu) / sd, df, m, s, SIMPLIFY = FALSE)),
      inverse   = function(z)  as.data.frame(mapply(function(x, mu, sd)  x * sd + mu, z,  m, s, SIMPLIFY = FALSE)),
      mean = m, sd = s
    )
  }

  # observed-data log-likelihood under LMS node mixture
  loglik_obs <- function(X, fit) {
    # Build MU/SIGMA/w from the fitted model's quadrature
    mdl  <- fit$model
    quad <- mdl$quad
    V    <- quad$n
    w    <- quad$w
    k    <- length(w)
    MU    <- vector("list", k)
    SIGMA <- vector("list", k)
    for (i in seq_len(k)) {
      MU[[i]]    <- as.numeric(muLmsCpp(model = mdl, z = V[i, ]))
      SIGMA[[i]] <- as.matrix(sigmaLmsCpp(model = mdl, z = V[i, ]))
    }
    # Cholesky precompute
    cholL  <- lapply(SIGMA, function(S) {
      # tiny ridge for safety
      diag(S) <- pmax(diag(S), 0); S <- S + diag(1e-8, nrow(S))
      chol(S)
    })
    logdet <- vapply(cholL, function(L) 2 * sum(log(diag(L))), 0.0)
    const  <- ncol(X) * log(2 * pi)
    lw     <- log(w)
    lse <- function(a) { m <- max(a); m + log(sum(exp(a - m))) }
    s <- 0.0
    for (i in seq_len(nrow(X))) {
      x <- as.numeric(X[i, ])
      z <- numeric(length(w))
      for (q in seq_along(w)) {
        r <- x - MU[[q]]
        y <- backsolve(cholL[[q]], r, transpose = TRUE)
        z[q] <- -0.5 * (sum(y * y) + logdet[q] + const)
      }
      s <- s + lse(z + lw)
    }
    s
  }

  # ---- prep data/cols ----
  ovs  <- getOVs(modsemify(model))
  data <- as.data.frame(data)[ovs]
  N    <- nrow(data)

  isOrdered <- vapply(data, is.ordered, logical(1))
  ordered   <- union(names(data)[isOrdered], ordered)
  ordered   <- intersect(ordered, colnames(data))  # sanitize

  if (!length(ordered)) {
    if (verbose) printf("No ordinal variables detected; fitting LMS directly.\n")
    return(modsem(model, data, method = "lms", ...))
  }

  # Make an integer-coded copy for ordinals (1..K). Keep others numeric as-is.
  data.int <- data
  data.int[ordered] <- lapply(data.int[ordered], function(x) as.integer(as.ordered(x)))

  # --- Standardize ONCE (optional but matches your original intent) ----
  # If you prefer to stay on the raw scale, set data.sc <- data.int and skip scaling.
  sc      <- scale_once(data.int)
  Y       <- sc$transform(data.int)
  Y       <- as.data.frame(Y)  # LMS expects a data.frame

  theta0  <- NULL
  ell0    <- NA_real_

  if (verbose) {
    printf("==================================================================\n")
    printf("Starting EM (max %d iterations). Ordinals: %s\n",
           MAX_ITER, paste(ordered, collapse = ", "))
  }

  for (iter in seq_len(MAX_ITER)) {
    if (verbose) {
      printf("------------------------------------------------------------------\n")
      printf("Iteration %d\n", iter)
    }

    # --- M-step: fit LMS on current Y
    fit <- modsem(model, data = Y, method = "lms", calc.se = FALSE, ...)
    fit_std <- standardized_estimates(fit)
    # --- Report current observed-data log-likelihood under node mixture
    ell <- tryCatch(loglik_obs(as.matrix(Y), fit),
                    error = function(e) { if (verbose) message("LL calc error: ", e$message); NA_real_ })
    if (verbose) printf("Log-likelihood (observed-data, mixture): %.6f\n", ell)

    if (verbose) {
      printf("Standardized estimates:\n")
      print(fit_std)
    }

    # --- E-step: deterministic mixture expectation for ordinals
    # IMPORTANT: pass the **current-scale** data used to fit the model
    Y <- getExpectedOrdinalValues(
      data    = Y,
      fit     = fit,
      ordered = ordered,
      tol     = tol,
      verbose = TRUE
    )
    Y <- as.data.frame(Y)

    # ---- Convergence check
    theta1 <- fit_std$est
    diff_theta <- if (is.null(theta0)) NA_real_ else mean(abs(theta1 - theta0))
    diff_ell   <- if (is.na(ell0) || is.na(ell)) NA_real_ else abs(ell - ell0)

    if (verbose) {
      printf("Δθ (mean |Δ|): %s   Δℓ: %s\n",
             ifelse(is.na(diff_theta), "NA", format(diff_theta, digits = 5)),
             ifelse(is.na(diff_ell),   "NA", format(diff_ell,   digits = 5)))
    }

    if (!is.null(theta0) && !is.na(diff_theta) && diff_theta < 1e-5) {
      if (verbose) printf("Converged by parameter change < 1e-5.\n")
      break
    }
    theta0 <- theta1
    ell0   <- ell
  }

  # ---- Post-processing: attach threshold estimates & vcov from rescaling ----
  fit.out <- modsem(model, data = Y, method = "lms", calc.se = calc.se, ...)

  # Simulate on the **final model** to build threshold VCOVs
  sim <- simulateDataParTable(
    parTable = parameter_estimates(fit.out),  # fixed typo
    N        = N
  )

  rescaled <- rescaleOrderedData(
    data           = data.int,
    sim.ov         = sim$oV,
    cols.ordered   = ordered,
    cols.cont      = setdiff(colnames(data.int), ordered),
    linear.ovs     = fit.out$model$info$allIndsXis,
    thresholds.vcov= TRUE
  )

  # Pack thresholds into parTable with SEs
  vcov.t <- NULL
  coef.t <- NULL
  for (col in ordered) {
    vcov.t <- diagPartitionedMat(vcov.t, rescaled$vcov[[col]])
    coef.t <- c(coef.t, rescaled$thresholds[[col]])
  }

  parTable <- fit.out$parTable
  cols.t   <- c("lhs", "op", "rhs", "label", "est", "std.error")
  cols.m   <- setdiff(colnames(parTable), cols.t)

  lhs.t <- stringr::str_split_i(names(coef.t), pattern = "\\|", i = 1L)
  rhs.t <- stringr::str_split_i(names(coef.t), pattern = "\\|", i = 2L)

  parTable.t <- data.frame(
    lhs = lhs.t, op = "|", rhs = rhs.t, label = "",
    est = coef.t, std.error = sqrt(diag(vcov.t))
  )
  parTable.t[cols.m] <- NA
  parTable.t <- addZStatsParTable(parTable.t)

  parTable <- sortParTableDA(rbind(parTable, parTable.t), model = fit.out$model)

  fit.out$coef.all <- c(fit.out$coef.all, coef.t)
  fit.out$vcov.all <- diagPartitionedMat(fit.out$vcov.all, vcov.t)
  fit.out$parTable <- parTable

  if (verbose) {
    printf("==================================================================\n")
    printf("Finished. Final log-likelihood: %.6f\n", ell0)
  }

  fit.out
}
