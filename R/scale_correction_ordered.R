modsemOrderedScaleCorrection <- function(model.syntax,
                                         data,
                                         method = "lms",
                                         ordered = NULL,
                                         calc.se = TRUE,
                                         # convergence controls
                                         R = 50,
                                         tol = 1e-4,                 # stopping tolerance
                                         tol_kind = c("theta","data"),
                                         patience = 3,               # need 'patience' consecutive passes under tol
                                         lambda = 0.9,               # damping in (0,1]
                                         # simulation controls
                                         N = 1e6,
                                         smooth_eps = 0,             # e.g. 0.5 for Laplace smoothing of category probs
                                         verbose = interactive(),
                                         optimize = TRUE,
                                         start = NULL,
                                         ...) {
  message("Scale correcting ordinal variables. ",
          "This is an experimental feature!\n",
          "See `help(modsem_da)` for more information.")

  if (is.null(verbose)) verbose <- TRUE
  tol_kind <- match.arg(tol_kind)

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]
  if (length(cols.ordered) == 0L) {
    warning("No ordered variables found; returning a standard fit.")
    return(modsem(model.syntax, data, method, calc.se, verbose, optimize, start, ...))
  }

  # Empirical re-scaling from simulated y to category-conditional means
  # Quantile-cutpoint style with left-closed/right-open bins (last bin inclusive)
  # + exact centering to kill residual mean drift.
  rescaleOrderedVariable <- function(name, data, sim.ov,
                                     smooth_eps = 0,
                                     eps_expand = 1e-12) {
    x <- as.integer(as.ordered(data[[name]]))
    y <- sim.ov[, name]

    # standardize for scale stability (no distributional assumption)
    y_mean <- mean(y, na.rm = TRUE)
    y_sd   <- stats::sd(y,  na.rm = TRUE)
    if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
    y <- (y - y_mean) / y_sd

    # observed category proportions (optional Laplace smoothing)
    tab <- as.numeric(table(x))
    K   <- length(tab)
    n   <- sum(tab)
    if (smooth_eps > 0) {
      p_obs <- (tab + smooth_eps) / (n + K * smooth_eps)
    } else {
      p_obs <- tab / n
    }

    # empirical CDF from data; quantile cutpoints on simulated y
    cdf_obs <- c(0, cumsum(p_obs))
    q <- stats::quantile(y, probs = cdf_obs, names = FALSE, type = 7)

    # compute conditional means in each [q[i], q[i+1]) (last bin inclusive)
    mu <- numeric(K)
    for (i in seq_len(K)) {
      lo <- q[i]
      hi <- q[i + 1]

      # left-closed/right-open, last bin inclusive
      if (i < K) {
        idx <- (y >= lo & y < hi)
      } else {
        idx <- (y >= lo & y <= hi)
      }

      # guard for degenerate bins (ties / finite-N quantile artifacts)
      if (!any(idx)) {
        # expand a hair around the cutpoints to catch ties
        if (i < K) {
          idx <- (y >= (lo - eps_expand) & y < (hi + eps_expand))
        } else {
          idx <- (y >= (lo - eps_expand) & y <= (hi + eps_expand))
        }
        # still empty? fallback to midpoint
        if (!any(idx)) mu[i] <- 0.5 * (lo + hi)
      }

      if (any(idx)) mu[i] <- mean(y[idx])
    }

    # exact centering to remove residual drift (weighted by observed category probs)
    mu <- mu - sum(p_obs * mu, na.rm = TRUE)

    # map each observed category to its conditional mean
    out <- rep(NA_real_, length(x))
    for (i in seq_len(K)) out[x == i] <- mu[i]
    out
  }

  rescaleOrderedData <- function(data, sim.ov, smooth_eps = 0) {
    data.y <- data
    for (col in cols.ordered) {
      data.y[[col]] <- rescaleOrderedVariable(col, data, sim.ov, smooth_eps)
    }
    data.y
  }

  dist_theta <- function(a, b) mean(abs(a - b))
  dist_data  <- function(A, B, cols) {
    s <- 0; m <- 0
    for (c in cols) {
      aa <- A[[c]]; bb <- B[[c]]
      ok <- is.finite(aa) & is.finite(bb)
      if (any(ok)) { s <- s + mean(abs(aa[ok] - bb[ok])); m <- m + 1 }
    }
    if (m == 0) NA_real_ else s / m
  }

  # initialization

  data.x <- data
  # start with a one-shot simulation-based rescale using a naive fit
  printedLines <- utils::capture.output(split = TRUE, {
    if (verbose)
      printf("Initializing scale-correction with a naive fit...\n")

    fit0 <- modsem(
      model.syntax = model.syntax,
      data         = data.x,
      method       = method,
      calc.se      = FALSE,
      verbose      = verbose,
      optimize     = optimize,
      start        = start,
      ...
    )
  })

  sim0 <- simulateDataParTable(
    parTable = parameter_estimates(fit0),
    N        = N
  )

  data.y <- rescaleOrderedData(data = data.x, sim.ov = sim0$oV, smooth_eps = smooth_eps)

  prev.theta  <- NULL
  prev.data.y <- data.y
  best.theta  <- fit0$theta
  best.fit    <- fit0
  passes <- 0

  for (iter in seq_len(R)) {
    nprinted <- length(printedLines)
    eraseConsoleLines(nprinted)

    printedLines <- utils::capture.output(split = TRUE, {
      if (verbose)
        printf("Scale-correction iter %d/%d...\n", iter, R)

      fit <- modsem(
        model.syntax = model.syntax,
        data         = data.y,
        method       = method,
        calc.se      = FALSE,
        verbose      = verbose,
        optimize     = if (is.null(prev.theta)) optimize else FALSE,
        start        = if (is.null(prev.theta)) start else prev.theta,
        ...
      )
    })

    theta_new <- fit$theta

    sim <- simulateDataParTable(parTable = parameter_estimates(fit), N = N)

    data.y_new <- rescaleOrderedData(data = data.x, sim.ov = sim$oV,
                                     smooth_eps = smooth_eps)

    # damping to promote contraction
    for (c in cols.ordered) {
      data.y[[c]] <- (1 - lambda) * data.y[[c]] + lambda * data.y_new[[c]]
    }

    # convergence metric
    met <- NA_real_
    if (!is.null(prev.theta))
      met <- dist_theta(theta_new, prev.theta)

    if (is.finite(met)) {
      if (met < tol) {
        passes <- passes + 1
      } else {
        passes <- 0
      }

      if (passes >= patience) {
        nprinted <- length(printedLines)
        eraseConsoleLines(nprinted)

        if (verbose)
          printf("Converged (metric %.3g) for %d consecutive checks.\n", met, passes)

        best.theta <- theta_new; best.fit <- fit
        break
      }
    }

    prev.theta  <- theta_new
    prev.data.y <- data.y
    best.theta  <- theta_new
    best.fit    <- fit
  }

  # final fit on frozen scaling
  modsem(
    model.syntax = model.syntax,
    method       = method,
    data         = data.y,
    calc.se      = calc.se,
    verbose      = verbose,
    start        = best.theta,   # warm-start
    optimize     = FALSE,
    ...
  )
}
