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


modsemOrderedScaleCorrectionV2 <- function(model.syntax,
                                         data,
                                         method = "lms",
                                         ordered = NULL,
                                         calc.se = TRUE,
                                         iter = 75L,
                                         warmup = 25L,
                                         N = max(NROW(data), 1e5),
                                         se = "simple",
                                         diff.theta.tol = 1e-10,
                                         ordered.mean.observed = FALSE,
                                         verbose = interactive(),
                                         optimize = TRUE,
                                         start = NULL,
                                         mean.observed = NULL,
                                         scaling.factor.int = NULL,
                                         ...) {

  message("Correcting scale of ordered variables via ordered-probit CE for endogenous indicators...")

  standardize <- \(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  if (is.null(verbose)) verbose <- TRUE

  parTable.in <- modsemify(model.syntax)
  xis <- getXis(parTable.in, checkAny = FALSE)

  inds <- getInds(parTable.in)
  parTable.cfa <- parTable.in[parTable.in$op == "=~" |
                              (parTable.in$lhs %in% inds &
                               parTable.in$rhs %in% inds), ,
                              drop = FALSE]
  syntax.cfa <- parTableToSyntax(parTable.cfa)
  fit.cfa    <- lavaan::cfa(syntax.cfa, data = data, ordered = ordered)
  fscores <- as.data.frame(lavaan::lavPredict(fit.cfa))

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]
  cols.cont    <- cols[!(cols %in% cols.ordered) & vapply(data, is.numeric, TRUE)]

  stopif(iter <= warmup, "`ordered.iter` must be larger than `ordered.warmup`!")
  iter.keep <- max(1L, floor(iter - warmup))

  data.x <- data
  data_i <- data
  data_i[cols.ordered] <- lapply(data_i[cols.ordered], \(x) standardize(as.integer(x)))

  thresholds <- NULL
  theta_i    <- NULL
  last_rescale <- NULL

  for (i in seq_len(iter)) {
    printedLines <- utils::capture.output(split = TRUE, {
      j <- i - warmup
      mode <- if (j <= 0) "warmup" else "sampling"

      if (verbose) printf("Iterations %d/%d [%s]...\n", i, iter, mode)

      optimize_i <- is.null(theta_i)
      start_i    <- if (is.null(theta_i)) NULL else theta_i

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
        is.int <- grepl(":", pars$rhs)
        pars[is.int, "est"] <- scaling.factor.int * pars[is.int, "est"]
      }

      # -------- Alternative A: ordered-probit CE rescaling ----------

      rescaled <- rescaleOrderedData_OP(
        data          = data.x,
        cols.ordered  = cols.ordered,
        parTable.in   = parTable.in,
        fscores    = fscores,
        thresholds.vcov = FALSE
      )

      last_rescale <- rescaled  # remember thresholds from last pass

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
      # -------------------------------------------------------------

      if (!is.null(theta_k)) {
        epsilon    <- theta_i - theta_k
        diff.theta <- mean(abs(epsilon))
        if (diff.theta <= diff.theta.tol && j > 1L) {
          printf("Solution converged (iter = %d, warmup = %d, diff = %.2g)\n",
                 i, warmup, diff.theta)
          break
          # break naturally after capture
        }
      }
    })

    nprinted <- length(printedLines)
    if (i < iter) eraseConsoleLines(nprinted)
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

  # Use the last pass' thresholds (probit for endogenous; analytic for others)
  rescaled <- last_rescale
  vcov.t <- NULL
  coef.t <- NULL
  for (col in cols.ordered) {
    if (!is.null(rescaled$vcov[[col]]) && length(rescaled$vcov[[col]])) {
      vcov.t <- diagPartitionedMat(vcov.t, rescaled$vcov[[col]])
    }
    if (!is.null(rescaled$thresholds[[col]]) && length(rescaled$thresholds[[col]])) {
      coef.t <- c(coef.t, rescaled$thresholds[[col]])
    }
  }

  parTable <- fit.out$parTable
  cols.t <- c("lhs", "op", "rhs", "label", "est", "std.error")
  cols.m <- setdiff(colnames(parTable), cols.t)

  if (length(coef.t)) {
    lhs.t <- stringr::str_split_i(names(coef.t), pattern = "\\|", i = 1L)
    rhs.t <- stringr::str_split_i(names(coef.t), pattern = "\\|", i = 2L)
    parTable.t <- data.frame(lhs = lhs.t, op = "|", rhs = rhs.t, label = "",
                             est = coef.t,
                             std.error = if (!is.null(vcov.t)) sqrt(diag(vcov.t)) else NA_real_)
    parTable.t[cols.m] <- NA
    parTable.t <- addZStatsParTable(parTable.t)
    parTable <- sortParTableDA(rbind(parTable, parTable.t), model = fit.out$model)
    fit.out$coef.all <- c(fit.out$coef.all, coef.t)
    if (!is.null(vcov.t)) fit.out$vcov.all <- diagPartitionedMat(fit.out$vcov.all, vcov.t)
  }

  fit.out$parTable <- parTable
  fit.out$args$optimize <- optimize
  fit.out$args$start    <- start
  fit.out
}


# Closed-form CE for ordered-probit (vectorized)
# y: integer categories 1..K, eta: linear predictor, tau: interior thresholds (length K-1)
ce_oprobit <- function(y, eta, tau, eps = 1e-12) {
  K <- max(y, na.rm = TRUE)
  a <- c(-Inf, tau)[y]
  b <- c(tau,  Inf)[y]
  al <- a - eta
  bl <- b - eta
  num <- stats::dnorm(al) - stats::dnorm(bl)
  den <- pmax(stats::pnorm(bl) - stats::pnorm(al), eps)
  eta + num / den
}

impute_with_oprobit <- function(y_ord, X, standardize_out = TRUE, weights = NULL) {
  # y_ord: integer/ordered factor of categories (1..K)
  # X    : data.frame/matrix of predictors (already includes interactions if desired)
  # weights: optional nonnegative case weights for the CE aggregation within categories
  
  # 1) Fit ordered probit (no explicit intercept; thresholds capture location)
  y_fac <- if (is.factor(y_ord) || is.ordered(y_ord)) y_ord else factor(y_ord, ordered = TRUE)
  df    <- as.data.frame(X)
  fit   <- suppressWarnings(MASS::polr(y_fac ~ ., data = df, method = "probit", Hess = FALSE))
  
  # Extract pieces
  beta <- stats::coef(fit)            # slope coefficients
  tau  <- as.numeric(fit$zeta)        # interior thresholds (K-1)
  K    <- length(levels(y_fac))
  yint <- as.integer(y_fac)
  
  # 2) Linear predictor eta_i = x_i' beta
  Xmm  <- model.matrix(~ . , data = df)[, names(beta), drop = FALSE]
  eta  <- as.numeric(Xmm %*% beta)
  sigma <- 1.0                         # homoskedastic probit here
  
  # 3) Compute conditional latent means m_ij for all i,j
  # Bounds a_j, b_j with a_1=-Inf, b_K=+Inf
  a <- c(-Inf, tau)
  b <- c(tau,  Inf)
  n <- length(eta)
  m <- matrix(NA_real_, nrow = n, ncol = K)
  
  # truncated-normal mean:
  # E[y* | a<y*<=b, eta, sigma] = eta + sigma * {phi((a-eta)/sigma) - phi((b-eta)/sigma)} /
  #                                           {Phi((b-eta)/sigma) - Phi((a-eta)/sigma)}
  eps <- 1e-12
  for (j in seq_len(K)) {
    alpha <- (a[j] - eta) / sigma
    betaJ <- (b[j] - eta) / sigma
    numer <- dnorm(alpha) - dnorm(betaJ)
    denom <- pnorm(betaJ) - pnorm(alpha)
    denom <- pmax(denom, eps)
    m[, j] <- eta + sigma * (numer / denom)
  }
  colnames(m) <- levels(y_fac)
  
  # 4) Aggregate to ONE CE per category: average over *observed* members of that category
  if (is.null(weights)) {
    w <- rep(1, n)
  } else {
    if (length(weights) != n) stop("weights must have length nrow(X).")
    w <- as.numeric(weights)
  }
  
  mu <- numeric(K)
  for (j in seq_len(K)) {
    idx <- which(yint == j)
    if (length(idx) == 0L) {
      # Fallback: if a category is unobserved, use population-weighted average
      pj <- pnorm((b[j] - eta)/sigma) - pnorm((a[j] - eta)/sigma)
      mu[j] <- sum(m[, j] * pj * w) / max(sum(pj * w), eps)
    } else {
      wj <- w[idx]
      mu[j] <- sum(m[idx, j] * wj) / sum(wj)
    }
  }
  names(mu) <- levels(y_fac)
  
  # 5) Impute: replace each observed category with its CE constant
  y_ce <- mu[yint]
  
  # 6) Optional standardization of the imputed series to mean 0, sd 1
  if (isTRUE(standardize_out)) {
    s <- stats::sd(y_ce)
    if (is.finite(s) && s > 0) {
      y_ce <- (y_ce - mean(y_ce)) / s
    } else {
      warning("Standardization skipped: zero or non-finite sd.")
    }
  }

  list(
    values      = y_ce,   # imputed vector: one constant per observed category
    mu_per_cat  = mu,     # the category constants (on probit scale)
    tau         = tau,    # probit thresholds
    eta         = eta,    # linear predictor (for reference)
    fit         = fit     # the polr fit (in case you want to inspect)
  )
}

# Build X for a given endogenous latent's indicator:
# includes the latent's composite itself + all structural RHS terms for that latent, with ':' expanded as products
.build_X_for_indicator <- function(latent, rhs_terms, fscores) {
  X <- list()  # don't include latent self
  if (length(rhs_terms)) {
    for (trm in rhs_terms) {
      if (grepl(":", trm, fixed = TRUE)) {
        parts <- strsplit(trm, ":", fixed = TRUE)[[1]]
        stopifnot(all(parts %in% colnames(fscores)))
        X[[trm]] <- Reduce(`*`, lapply(parts, function(p) fscores[[p]]))
      } else {
        stopifnot(trm %in% colnames(fscores))
        X[[trm]] <- fscores[[trm]]
      }
    }
  }

  X <- as.data.frame(X)
  # drop duplicate columns if any (name collisions)
  X[, !duplicated(names(X)), drop = FALSE]
}

# Hybrid rescaler:
# - For ordered indicators of ENDOGENOUS latents: ordered-probit CE using fscores + interactions
# - For other ordered indicators: analytic (normal) mapping
rescaleOrderedData_OP <- function(data, cols.ordered, parTable.in, fscores,
                                  thresholds.vcov = FALSE, ...) {
  if (!length(cols.ordered)) return(list(data = data, thresholds = list(), vcov = list()))

  out <- as.data.frame(data, stringsAsFactors = FALSE)
  data.cat <- as.data.frame(data)
  data.cat[cols.ordered] <- lapply(data.cat[cols.ordered], as.integer)

  # maps: indicator -> latent (from loadings)
  meas <- subset(parTable.in, op == "=~")
  ind2lat <- stats::setNames(as.character(meas$lhs), meas$rhs)

  # structural terms per endogenous latent
  struc <- subset(parTable.in, op == "~")
  rhs_by_lat <- tapply(struc$rhs, struc$lhs, c, simplify = FALSE)

  etas <- getEtas(parTable.in, checkAny = FALSE)
  ordered_endog <- intersect(cols.ordered, names(ind2lat)[ind2lat %in% etas])

  thresholds <- vector("list", length(cols.ordered)); names(thresholds) <- cols.ordered
  vcovs      <- vector("list", length(cols.ordered)); names(vcovs)      <- cols.ordered

  for (col in cols.ordered) {
    if (col %in% ordered_endog) {
      # indicator belongs to an endogenous latent -> ordered-probit CE
      latent <- ind2lat[[col]]
      rhs_terms <- rhs_by_lat[[latent]]; if (is.null(rhs_terms)) rhs_terms <- character(0)
      X <- .build_X_for_indicator(latent, rhs_terms, fscores)
      imp <- impute_with_oprobit(y_ord = data[[col]], X = X, standardize_out = TRUE)
      out[[col]] <- imp$values
      # give names to thresholds like "col|t1", "col|t2", ...
      if (length(imp$tau)) {
        tn <- paste0(col, "|t", seq_along(imp$tau))
        thresholds[[col]] <- stats::setNames(imp$tau, tn)
      } else thresholds[[col]] <- numeric(0)
      vcovs[[col]] <- NULL  # lightweight: no vcov for probit thresholds here
    } else {
      # exogenous or non-endogenous ordered indicator -> analytic normal mapping
      ana <- rescaleOrderedVariableAnalytic(name = col, data = data.cat,
                                            thresholds.vcov = thresholds.vcov, ...)
      out[[col]]      <- ana$values
      thresholds[[col]] <- ana$thresholds
      vcovs[[col]]      <- ana$vcov
    }
  }

  list(data = out, thresholds = thresholds, vcov = vcovs)
}
