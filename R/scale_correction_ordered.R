modsemOrderedScaleCorrection <- function(model.syntax,
                                         data,
                                         method = "lms",
                                         ordered = NULL,
                                         calc.se = TRUE,
                                         iter = 75L,
                                         warmup = 25L,
                                         N = max(NROW(data), 1e5), # 10,000 initally
                                         se = "simple",
                                         diff.theta.tol = 1e-10,
                                         ordered.mean.observed = FALSE,
                                         # Capture args
                                         verbose = interactive(),
                                         optimize = TRUE,
                                         start = NULL,
                                         mean.observed = NULL,
                                         scaling.factor.int = NULL,
                                         ...) {
  message("Correcting scale of ordered variables...\n",
          "This is an experimental feature, ",
          "see `help(modsem_da)` for more information!")

  standardize <- \(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  if (is.null(verbose))
    verbose <- TRUE # default

  parTable.in <- modsemify(model.syntax)
  xis <- getXis(parTable.in, checkAny = FALSE)
  inds.xis <- unique(unlist(getIndsLVs(parTable.in, lVs = xis)))

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]
  cols.cont    <- cols[!(cols %in% cols.ordered) & vapply(data, is.numeric, TRUE)]

  ERROR <- \(e) {warning2(e, immediate. = FALSE); NULL}

  data.x <- data
  data.y <- data
  data.y[cols.ordered] <- lapply(data.y[cols.ordered], FUN = as.integer)

  stopif(iter <= warmup, "`ordered.iter` must be larger than `ordered.warmup`!")

  iter.keep <- max(1L, floor(iter - warmup))

  data_i <- data
  data_i[cols.ordered] <- lapply(data_i[cols.ordered],
                                 FUN = \(x) standardize(as.integer(x)))
  thresholds <- NULL
  theta_i    <- NULL

  for (i in seq_len(iter)) {
    printedLines <- utils::capture.output(split = TRUE, {
      j <- i - warmup
      mode <- if (j <= 0) "warmup" else "sampling"

      if (verbose)
        printf("Iterations %d/%d [%s]...\n", i, iter, mode)

      if (is.null(theta_i)) {
        optimize_i <- TRUE
        start_i    <- NULL
      } else {
        optimize_i <- FALSE
        start_i   <- theta_i
      }

      fit_i <- tryCatch(
        modsem_da(
          model.syntax     = model.syntax,
          data             = data_i,
          method           = method,
          start            = start_i,
          verbose          = verbose,
          optimize         = optimize_i,
          calc.se          = FALSE,
          mean.observed    = ordered.mean.observed,
          ...
        ), error = \(e) {cat("\n"); print(e); NULL}
      )

      stopif(is.null(fit_i), "Model estimation failed!")

      theta_k <- theta_i
      theta_i <- fit_i$theta

      pars <- parameter_estimates(fit_i)

      if (!is.null(scaling.factor.int)) {
        # bias represents the degree to which the coefficient is
        # over/under estimated. we use it to get a scaling factor
        is.int <- grepl(":", pars$rhs)
        pars[is.int, "est"] <- scaling.factor.int * pars[is.int, "est"]
      }

      sim_i <- simulateDataParTable(
        parTable = pars,
        N        = N
      )

      rescaled <- rescaleOrderedData(data = data.x,
                                     sim.ov = sim_i$oV,
                                     cols.ordered = cols.ordered,
                                     cols.cont = cols.cont,
                                     linear.ovs = inds.xis)

      if (j >= 1L && !is.null(thresholds) && !is.null(data_i)) {
        data_j   <- rescaled$data
        lambda_j <- 1 / j
        lambda_i <- 1 - lambda_j

        for (col in cols.ordered)
          data_i[[col]] <- lambda_i * data_i[[col]] + lambda_j * data_j[[col]]

      } else {
        data_i     <- rescaled$data
        thresholds <- rescaled$thresholds
      }
    })

    nprinted <- length(printedLines)
    if (i < iter) eraseConsoleLines(nprinted)

    epsilon    <- theta_i - theta_k
    diff.theta <- mean(abs(epsilon))

    if (diff.theta <= diff.theta.tol && j > 1L) {
      printf("Solution converged (iter = %d, warmup = %d, diff = %.2g)\n",
             i, warmup, diff.theta)

      break
    }
  }

  fit.out <- modsem_da(model.syntax     = model.syntax,
                       data             = data_i,
                       method           = method,
                       start            = theta_i,
                       verbose          = verbose,
                       optimize         = FALSE,
                       calc.se          = calc.se,
                       mean.observed    = ordered.mean.observed,
                       ...)

  rescaled <- rescaleOrderedData(data = data.x,
                                 sim.ov = sim_i$oV,
                                 cols.ordered = cols.ordered,
                                 cols.cont = cols.cont,
                                 linear.ovs = inds.xis,
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
  fit.out$args$optimize <- optimize
  fit.out$args$start    <- start

  fit.out
}


rescaleOrderedVariableMonteCarlo <- function(name,
                                             data,
                                             sim.ov,
                                             smooth_eps = 0,
                                             eps_expand = 1e-12,
                                             standardize = TRUE,
                                             thresholds.vcov = FALSE,
                                             vcov.boot = 250) {
  x   <- as.ordered(data[[name]])
  x.i <- as.integer(x)
  y   <- sim.ov[, name]

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
  cdf_obs <- cumsum(p_obs[-length(p_obs)]) # drop last, we know that it sums to 1
  q <- stats::quantile(y, probs = cdf_obs, names = FALSE, type = 7)
  q <- c(-Inf, q, Inf) # [0, ..., last(p_obs)]

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

  if (thresholds.vcov) {
    Q <- matrix(NA_real_, nrow = vcov.boot, ncol = length(cdf_obs),
                dimnames = list(NULL, paste0(name, "|t", seq_along(cdf_obs))))

    for (i in seq_len(vcov.boot)) {
      N <- length(y)
      n <- length(x)

      replace  <- N <= 5L * n
      y.sample <- y[sample(N, n, replace = replace)]
      Q[i, ] <- stats::quantile(y.sample, probs = cdf_obs, names = FALSE, type = 7)
    }

    vcov <- stats::cov(Q)
  } else vcov <- NULL

  # exact centering to remove residual drift (weighted by observed category probs)
  mu <- mu - sum(p_obs * mu, na.rm = TRUE)

  # map each observed category to its conditional mean
  x.out <- rep(NA_real_, length(x))
  for (i in seq_len(K))
    x.out[x.i == i] <- mu[i]

  if (standardize)
    x.out <- std1(x.out)

  labels.t <- paste0(name, "|t", seq_len(sum(is.finite(q))))
  thresholds <- stats::setNames(q[is.finite(q)], labels.t)

  list(values = x.out, thresholds = thresholds, vcov = vcov)
}


rescaleOrderedVariableAnalytic <- function(name,
                                           data,
                                           sim.ov = NULL, # ignored
                                           smooth_eps = 0,
                                           standardize = TRUE,
                                           thresholds.vcov = FALSE,
                                           vcov.method = c("delta", "bootstrap"),
                                           vcov.boot = 250) {
  vcov.method <- match.arg(vcov.method)

  x   <- as.ordered(data[[name]])
  x.i <- as.integer(x)

  # observed category proportions (optional Laplace smoothing)
  tab <- as.numeric(table(x))
  K   <- length(tab)
  n   <- sum(tab)
  if (K < 2) stop("Need at least 2 ordered categories for '", name, "'.")

  if (smooth_eps > 0) {
    p_hat <- (tab + smooth_eps) / (n + K * smooth_eps)
  } else {
    p_hat <- tab / n
  }

  # cumulative probs for interior thresholds (K-1 of them)
  C <- cumsum(p_hat)[-K]  # drop last; sums to 1
  # numeric guards
  eps <- .Machine$double.eps^0.5
  C <- pmin(pmax(C, eps), 1 - eps)

  # thresholds on standard-normal scale
  t_int <- stats::qnorm(C)                 # length K-1
  q     <- c(-Inf, t_int, Inf)             # length K+1 for interval bounds

  # truncated-normal means for each category interval (a_i, b_i]
  # mu_i = (phi(a_i) - phi(b_i)) / (Phi(b_i) - Phi(a_i))
  a <- q[1:K]
  b <- q[2:(K+1)]
  Phi_a <- stats::pnorm(a)
  Phi_b <- stats::pnorm(b)
  phi_a <- stats::dnorm(a)
  phi_b <- stats::dnorm(b)
  denom <- pmax(Phi_b - Phi_a, eps)
  mu    <- (phi_a - phi_b) / denom

  # exact centering to remove small drift from smoothing/finite-n
  mu <- mu - sum(p_hat * mu)

  # map each observed category to its conditional mean
  x.out <- rep(NA_real_, length(x.i))
  for (i in seq_len(K)) x.out[x.i == i] <- mu[i]

  # optional standardization of the mapped scores (sample standardization)
  if (standardize)
    x.out <- std1(x.out)

  # labels for interior thresholds
  labels.t   <- paste0(name, "|t", seq_len(K - 1))
  thresholds <- stats::setNames(t_int, labels.t)

  vcov <- NULL
  if (thresholds.vcov && vcov.method == "delta") {
    # Delta method for t = qnorm(C), where C are cumulative probs from multinomial
    # Multinomial covariance of cell proportions: Var(p) = (diag(p) - pp^T)/n
    # Let C_j = sum_{i<=j} p_i, j=1..K-1.
    # J_p->C is lower-triangular with ones:
    J_pc <- matrix(0, nrow = K - 1, ncol = K)
    for (j in 1:(K - 1)) J_pc[j, 1:j] <- 1

    V_p  <- (diag(p_hat) - tcrossprod(p_hat, p_hat)) / n  # K x K
    V_C  <- J_pc %*% V_p %*% t(J_pc)                      # (K-1) x (K-1)

    # t = qnorm(C); dt/dC = 1 / dnorm(t) evaluated at t
    D_tC <- diag(1 / pmax(stats::dnorm(t_int), eps), nrow = K - 1)
    V_t  <- D_tC %*% V_C %*% D_tC

    colnames(V_t) <- rownames(V_t) <- labels.t
    vcov <- V_t

  } else if (thresholds.vcov && vcov.method == "bootstrap") {
    Q <- matrix(NA_real_, nrow = vcov.boot, ncol = K - 1,
                dimnames = list(NULL, labels.t))
    for (b_ix in seq_len(vcov.boot)) {
      # parametric bootstrap from Multinomial(n, p_hat)
      counts_b <- as.numeric(stats::rmultinom(1, size = n, prob = p_hat))
      p_b      <- counts_b / n
      C_b      <- cumsum(p_b)[-K]
      C_b      <- pmin(pmax(C_b, eps), 1 - eps)
      Q[b_ix,] <- stats::qnorm(C_b)
    }
    vcov <- stats::cov(Q)
  }

  list(values = x.out, thresholds = thresholds, vcov = vcov)
}


rescaleOrderedData <- function(data, sim.ov, cols.ordered, cols.cont,
                               thresholds.vcov = FALSE, linear.ovs = NULL,
                               ...) {
  if (!length(cols.ordered)) # nothing to do
    return(cols.ordered)

  # ensure column order consistent with data
  sim.mat <- as.data.frame(sim.ov)[, colnames(data), drop = FALSE]

  # make a full output copy right away, so continuous vars are preserved
  out <- as.data.frame(data, stringsAsFactors = FALSE)

  data.cat <- as.data.frame(data)
  data.cat[cols.ordered] <- lapply(data.cat[cols.ordered], as.integer)

  # thresholds + univariate conditional means for each ordered variable
  NamedList  <- \(nm) stats::setNames(vector("list", length(nm)), nm = nm)
  thresholds     <- NamedList(cols.ordered)
  univarEst      <- NamedList(cols.ordered)
  thresholdsVcov <- NamedList(cols.ordered)

  # standardized versions of the ordered cols in sim (these are the values we copy)
  sim.std.ord <- lapplyNamed(cols.ordered, \(col) std1(sim.mat[[col]]),
                             names = cols.ordered)
  sim.std.ord <- as.data.frame(sim.std.ord)

  for (col in cols.ordered) {
    if (col %in% linear.ovs) .f <- rescaleOrderedVariableAnalytic
    else                     .f <- rescaleOrderedVariableMonteCarlo

    scaling <- .f(name = col, data = data.cat, sim.ov = sim.mat,
                  thresholds.vcov = thresholds.vcov, ...)
    univarEst[[col]]      <- scaling$values
    thresholds[[col]]     <- scaling$thresholds
    thresholdsVcov[[col]] <- scaling$vcov
  }

  out[cols.ordered] <- as.data.frame(univarEst)[cols.ordered]
  # preserve original column order & types for continuous vars
  out <- out[, colnames(data), drop = FALSE]

  list(data = out, thresholds = thresholds, vcov = thresholdsVcov)
}


std1 <- function(v) {
  mu    <- mean(v, na.rm = TRUE)
  sigma <- stats::sd(v, na.rm = TRUE)

  if (!is.finite(sigma) || sigma == 0)
    sigma <- 1

  (v - mu) / sigma
}
