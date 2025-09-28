getExpectedOrdinalValues <- function(
  data, fit,
  R = 1e5,
  ordered = NULL,
  tol = 1e-5,
  verbose = TRUE,
  num_iter = 2L,                 # Gauss–Seidel sweeps
  parallel = TRUE,               # set future::plan() outside to control backend, or let this choose
  impute = c("mean","stochastic_exact","stochastic_gaussian"),
  return_variances = TRUE,
  seed = 123,
  workers = 8
) {
  impute <- match.arg(impute)
  if (!is.null(seed)) set.seed(seed)

  # choose a parallel plan only if none set
  if (parallel && length(future::plan("list")) == 0) {
    if (tolower(.Platform$OS.type) == "unix") {
      future::plan(future::multicore, workers = workers)
    } else {
      future::plan(future::multisession, workers = workers)
    }
  }

  # --- 0) Prep ----------------------------------------------------------------
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
    ans <- as.matrix(data)
    return(if (return_variances) list(values = ans, variances = matrix(0, nrow(ans), 0)) else ans)
  }

  # --- 1) Quadrature: MU/SIGMA/w with precomputations -------------------------
  model <- fit$model
  quad  <- model$quad
  V     <- quad$n
  w     <- quad$w
  k     <- length(w)

  MU    <- vector("list", k)
  SIGMA <- vector("list", k)
  CHOL  <- vector("list", k)
  LOGDET<- numeric(k)
  QINV  <- vector("list", k)

  for (i in seq_len(k)) {
    MU[[i]]    <- as.numeric(muLmsCpp(model = model, z = V[i, ]))
    S          <- as.matrix(sigmaLmsCpp(model = model, z = V[i, ]))
    diag(S)    <- pmax(diag(S), 0)
    S          <- S + diag(1e-8, nrow(S))
    SIGMA[[i]] <- S
    L          <- chol(S)
    CHOL[[i]]  <- L
    LOGDET[i]  <- 2 * sum(log(diag(L)))
    QINV[[i]]  <- chol2inv(L)  # precision
  }
  p     <- ncol(SIGMA[[1L]])
  const <- p * log(2 * pi)
  logw  <- log(w)

  # --- 2) Thresholds via your rescaling (unchanged) ---------------------------
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
  tau_list   <- lapply(thresholds, function(v) c(-Inf, v, Inf))

  # --- 3) Numeric matrices ----------------------------------------------------
  X <- data
  X[o.cols] <- lapply(X[o.cols], function(x) as.integer(as.ordered(x)))  # ensures 1..K
  X <- as.matrix(X)
  rownames(X) <- rownames(data)

  Y_vals <- X                       # will overwrite ordinal cols
  Y_var  <- matrix(0, nrow(X), length(o.cols),
                   dimnames = list(rownames(X), o.cols))
  o.idx  <- match(o.cols, colnames(X))

  # --- 4) Helpers (NA/underflow-proof) ---------------------------------------
  eps <- 1e-12

  softmax_safe <- function(lp) {
    if (anyNA(lp)) lp[is.na(lp)] <- -Inf
    m <- max(lp)
    if (!is.finite(m)) return(rep(1/length(lp), length(lp)))
    w <- exp(lp - m); s <- sum(w)
    if (!is.finite(s) || s <= 0) rep(1/length(lp), length(lp)) else w / s
  }

  trunc_mean_var_1d <- function(m, s, a, b) {
    if (!is.finite(m) || !is.finite(s)) return(c(mean = m, var = 0))
    s <- max(s, 1e-12)
    if (!is.finite(a)) a <- -Inf
    if (!is.finite(b)) b <-  Inf
    alpha <- (a - m)/s; beta <- (b - m)/s
    alpha <- pmax(alpha, -1e2); beta <- pmin(beta,  1e2)
    Za <- pnorm(alpha); Zb <- pnorm(beta); Z <- Zb - Za
    if (!is.finite(Z) || Z <= 0) return(c(mean = pmin(pmax(m, a), b), var = 0))
    kappa <- (dnorm(alpha) - dnorm(beta)) / Z
    m_out <- m + s * kappa
    v_out <- s^2 * (1 + (alpha*dnorm(alpha) - beta*dnorm(beta))/Z - kappa^2)
    c(mean = m_out, var = max(v_out, 0))
  }

  # robust truncated-normal sampler (inverse-CDF), CDF clipped into (eps, 1-eps)
  rtruncnorm1 <- function(m, s, a, b) {
    s <- max(s, 1e-12)
    if (!is.finite(a)) a <- -Inf
    if (!is.finite(b)) b <-  Inf
    alpha <- (a - m)/s; beta <- (b - m)/s
    Fa <- pnorm(alpha); Fb <- pnorm(beta)
    # if interval collapsed numerically, return mean
    if (!is.finite(Fa) || !is.finite(Fb) || Fb <= Fa) return(m)
    u <- runif(1)
    p <- Fa + u * (Fb - Fa)
    p <- min(max(p, eps), 1 - eps)     # <- critical: avoid qnorm(0/1)
    m + s * qnorm(p)
  }

  # conditional via precision, Σ-fallback, then marginal fallback
  cond_from_precision <- function(mu, Q, j, x, Sigma_fallback = NULL) {
    Qjj <- Q[j, j]
    if (is.finite(Qjj) && Qjj > 0 && all(is.finite(Q[j, -j]))) {
      s2 <- 1 / Qjj
      m  <- as.numeric(mu[j] - (Q[j, -j, drop=FALSE] %*% (x[-j] - mu[-j])) / Qjj)
      if (is.finite(m) && is.finite(s2) && s2 > 0) return(c(m = m, s = sqrt(s2)))
    }
    if (!is.null(Sigma_fallback)) {
      idx_r <- setdiff(seq_along(mu), j)
      S <- Sigma_fallback
      S_rr <- S[idx_r, idx_r, drop = FALSE]
      S_jr <- S[j, idx_r, drop = FALSE]
      S_rj <- S[idx_r, j, drop = FALSE]
      diag(S_rr) <- pmax(diag(S_rr), 0); S_rr <- S_rr + diag(1e-8, nrow(S_rr))
      cs <- try({
        L <- chol(S_rr)
        diff <- x[idx_r] - mu[idx_r]
        v <- backsolve(L, forwardsolve(t(L), diff))
        u <- backsolve(L, forwardsolve(t(L), S_rj))
        m <- as.numeric(mu[j] + S_jr %*% v)
        s2 <- as.numeric(S[j, j] - S_jr %*% u)
        c(m = m, s = sqrt(max(s2, 1e-12)))
      }, silent = TRUE)
      if (!inherits(cs, "try-error") && all(is.finite(cs))) return(cs)
    }
    # last resort
    if (!is.null(Sigma_fallback)) {
      return(c(m = mu[j], s = sqrt(max(Sigma_fallback[j, j], 1e-12))))
    } else {
      return(c(m = mu[j], s = 1))
    }
  }

  # log p(x | node t)
  ll_full_per_node <- function(x) {
    out <- numeric(k)
    for (t in seq_len(k)) {
      r  <- x - MU[[t]]
      y  <- backsolve(CHOL[[t]], r, transpose = TRUE)
      out[t] <- -0.5 * (sum(y * y) + LOGDET[t] + const) + logw[t]
    }
    out
  }

  # Row kernel
  expected_ordinals_row <- function(x_row) {
    cur <- x_row
    # init: interval midpoints (NA-safe)
    for (j_local in seq_along(o.idx)) {
      j <- o.idx[j_local]; col <- o.cols[j_local]
      tau <- tau_list[[col]]; K <- length(tau) - 1L
      catj <- as.integer(x_row[j])
      if (is.na(catj) || catj < 1L || catj > K) { a <- -Inf; b <- Inf } else { a <- tau[catj]; b <- tau[catj + 1L] }
      cur[j] <- if (is.finite(a) && is.finite(b)) 0.5*(a+b) else if (is.infinite(a) && is.finite(b)) b-1 else if (is.finite(a) && is.infinite(b)) a+1 else 0
    }

    mean_out <- setNames(numeric(length(o.idx)), o.cols)
    var_out  <- setNames(numeric(length(o.idx)), o.cols)

    for (it in seq_len(num_iter)) {
      log_joint <- ll_full_per_node(cur)

      for (j_local in seq_along(o.idx)) {
        j <- o.idx[j_local]; col <- o.cols[j_local]
        tau <- tau_list[[col]]; K <- length(tau) - 1L
        catj <- as.integer(x_row[j])
        if (is.na(catj) || catj < 1L || catj > K) { a <- -Inf; b <- Inf } else { a <- tau[catj]; b <- tau[catj + 1L] }

        log_post <- numeric(k)
        e_t      <- numeric(k)
        v_t      <- numeric(k)

        for (t in seq_len(k)) {
          cs <- cond_from_precision(MU[[t]], QINV[[t]], j, cur, SIGMA[[t]])
          mj <- unname(cs["m"]); sj <- unname(cs["s"])

          # NA-safe soft caps (pre-truncation)
          if (isTRUE(is.finite(a) && is.finite(mj) && mj < a)) mj <- max(mj, a - 8 * sj)
          if (isTRUE(is.finite(b) && is.finite(mj) && mj > b)) mj <- min(mj, b + 8 * sj)

          tv <- trunc_mean_var_1d(mj, sj, a, b)
          e_t[t] <- if (is.finite(tv["mean"])) tv["mean"] else mj
          v_t[t] <- if (is.finite(tv["var"]))  tv["var"]  else 0

          # node posterior up to const: log w_t + log p(x_-j|t) + log P(a<b | t)
          ll_cond_j <- dnorm(cur[j], mean = mj, sd = sj, log = TRUE)
          pij <- pnorm((b - mj)/sj) - pnorm((a - mj)/sj)
          if (!is.finite(pij)) pij <- 0
          log_post[t] <- (log_joint[t] - ll_cond_j) + log(pmax(pij, 1e-300))
        }

        wt <- softmax_safe(log_post)

        mbar <- sum(wt * e_t); if (!is.finite(mbar)) mbar <- mean(e_t[is.finite(e_t)])
        vbar <- sum(wt * (v_t + (e_t - mbar)^2)); if (!is.finite(vbar)) vbar <- 0
        mean_out[j_local] <- mbar
        var_out[j_local]  <- vbar

        # overwrite cur[j] for next sweep
        if (it < num_iter) {
          if (impute == "stochastic_exact") {
            tsel <- sample.int(k, size = 1L, prob = wt)
            cs   <- cond_from_precision(MU[[tsel]], QINV[[tsel]], j, cur, SIGMA[[tsel]])
            cur[j] <- rtruncnorm1(m = unname(cs["m"]), s = unname(cs["s"]), a = a, b = b)
          } else if (impute == "stochastic_gaussian") {
            cur[j] <- mbar + rnorm(1, sd = sqrt(max(vbar, 0)))
          } else {
            cur[j] <- mbar
          }
        } else {
          cur[j] <- mbar
        }
      }
    }

    # optional final sampling
    if (impute != "mean") {
      for (j_local in seq_along(o.idx)) {
        j <- o.idx[j_local]; col <- o.cols[j_local]
        tau <- tau_list[[col]]; K <- length(tau) - 1L
        catj <- as.integer(x_row[j])
        if (is.na(catj) || catj < 1L || catj > K) { a <- -Inf; b <- Inf } else { a <- tau[catj]; b <- tau[catj + 1L] }

        log_joint <- ll_full_per_node(cur)
        log_post <- numeric(k); e_t <- numeric(k); v_t <- numeric(k)
        for (t in seq_len(k)) {
          cs <- cond_from_precision(MU[[t]], QINV[[t]], j, cur, SIGMA[[t]])
          mj <- unname(cs["m"]); sj <- unname(cs["s"])
          tv <- trunc_mean_var_1d(mj, sj, a, b)
          e_t[t] <- if (is.finite(tv["mean"])) tv["mean"] else mj
          v_t[t] <- if (is.finite(tv["var"]))  tv["var"]  else 0
          ll_cond_j <- dnorm(cur[j], mean = mj, sd = sj, log = TRUE)
          pij <- pnorm((b - mj)/sj) - pnorm((a - mj)/sj); if (!is.finite(pij)) pij <- 0
          log_post[t] <- (log_joint[t] - ll_cond_j) + log(pmax(pij, 1e-300))
        }
        wt <- softmax_safe(log_post)

        if (impute == "stochastic_exact") {
          tsel <- sample.int(k, size = 1L, prob = wt)
          cs   <- cond_from_precision(MU[[tsel]], QINV[[tsel]], j, cur, SIGMA[[tsel]])
          mean_out[j_local] <- rtruncnorm1(m = unname(cs["m"]), s = unname(cs["s"]), a = a, b = b)
        } else {
          mean_out[j_local] <- mean_out[j_local] + rnorm(1, sd = sqrt(max(var_out[j_local], 0)))
        }
      }
    }

    list(mean = mean_out, var = var_out)
  }

  # --- 5) Parallel map over rows ---------------------------------------------
  N <- nrow(X)
  if (N == 0L) {
    ans <- Y_vals
    return(if (return_variances) list(values = ans, variances = Y_var) else ans)
  }
  if (isTRUE(verbose)) cat("E-step (", impute, ", mixture expectations)…\n", sep = "")

  row_fun <- function(n) expected_ordinals_row(X[n, ])

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    chunks <- future.apply::future_lapply(seq_len(N), function(n) row_fun(n),
                                          future.seed = TRUE)
  } else if (parallel && .Platform$OS.type != "windows") {
    chunks <- parallel::mclapply(seq_len(N), row_fun, mc.cores = max(1L, parallel::detectCores()-1L))
  } else {
    tick <- max(1L, floor(N/50L))
    chunks <- vector("list", N)
    for (n in seq_len(N)) {
      if (verbose && (n %% tick == 0L || n == 1L || n == N)) {
        cat(sprintf("\rProgress: %d/%d", n, N)); utils::flush.console()
      }
      chunks[[n]] <- row_fun(n)
    }
    if (verbose) cat("\n")
  }

  # stitch outputs
  Mmat <- do.call(rbind, lapply(chunks, `[[`, "mean"))
  Vmat <- do.call(rbind, lapply(chunks, `[[`, "var"))

  Y_vals[, o.idx] <- Mmat
  if (return_variances) {
    Y_var[,] <- Vmat
    return(list(values = Y_vals, variances = Y_var))
  } else {
    return(Y_vals)
  }
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
    Expectation <- getExpectedOrdinalValues(
      data    = Y,
      fit     = fit,
      ordered = ordered,
      tol     = tol,
      verbose = TRUE,
      impute  = "stochastic_exact"
    )
    Y <- as.data.frame(Expectation$values)

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
