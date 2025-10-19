modsemOrderedScaleCorrection <- function(..., ordered.v2 = TRUE) {
  message("Transforming ordered variables to interval scale...\n",
          "This is an experimental feature, see `help(modsem_da)` for more information!")

  # v2 seems to be the best. Left here for testing purposes
  # The user can pass the `ordered.v2` argument to `modsem_da`
  # without it being publicly exposed as an option
  if (ordered.v2) FUN <- modsemOrderedScaleCorrectionV2
  else            FUN <- modsemOrderedScaleCorrectionV1

  FUN(...)
}


modsemOrderedScaleCorrectionV1 <- function(model.syntax,
                                           data,
                                           method = "lms",
                                           ordered = NULL,
                                           calc.se = TRUE,
                                           group = NULL,
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
                                           probit.correction = FALSE,
                                           ...) {
  standardize <- \(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  stopif(!is.null(group), "Multigroup models do not work with: ",
         "`modsemOrderedScaleCorrectionV1()`!")

  if (is.null(verbose))
    verbose <- TRUE # default

  parTable.in <- modsemify(model.syntax)
  ovs <- getOVs(parTable.in)
  xis <- getXis(parTable.in, checkAny = FALSE)
  inds.xis <- unique(unlist(getIndsLVs(parTable.in, lVs = xis)))

  data <- as.data.frame(data)[ovs]
  cols <- ovs
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

  pcorr <- polychor(vars = cols, data = data,
                    thresholds = NULL, # rescaled$thresholds,
                    ordered    = cols.ordered)

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
          group            = group,
          ...
        ), error = \(e) {cat("\n"); print(e); NULL}
      )

      stopif(is.null(fit_i), "Model estimation failed!")

      theta_k <- theta_i
      theta_i <- fit_i$theta

      pars <- parameter_estimates(fit_i)

      sim_i <- simulateDataParTable(
        parTable = pars,
        N        = N
      )

      rescaled <- rescaleOrderedData(data = data.x,
                                     sim.ov = sim_i$oV,
                                     cols.ordered = cols.ordered,
                                     cols.cont = cols.cont,
                                     linear.ovs = inds.xis)

      if (probit.correction) data_j <- fitX2Cov(X = rescaled$data, S = pcorr)
      else                   data_j <- rescaled$data

      if (j >= 1L && !is.null(thresholds) && !is.null(data_i)) {
        lambda_j <- 1 / j
        lambda_i <- 1 - lambda_j

        for (col in cols.ordered)
          data_i[[col]] <- lambda_i * data_i[[col]] + lambda_j * data_j[[col]]

      } else {
        data_i     <- data_j
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
                                           group = NULL,
                                           N = max(NROW(data), 1e5),
                                           se = "simple",
                                           diff.theta.tol = 1e-10,
                                           ordered.mean.observed = FALSE,
                                           verbose = interactive(),
                                           mean.observed = NULL,
                                           probit.correction = FALSE,
                                           ...) {
  standardize <- \(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  if (is.null(verbose)) verbose <- TRUE

  parTable.in <- modsemify(model.syntax)
  xis  <- getXis(parTable.in, checkAny = FALSE)
  inds <- getInds(parTable.in)

  parTable.cfa <- parTable.in[parTable.in$op == "=~" |
                              (parTable.in$lhs %in% inds &
                               parTable.in$rhs %in% inds), ,
                              drop = FALSE]

  data <- as.data.frame(data)
  cols <- getOVs(parTable.in)
  cols.ordered <- cols[cols %in% ordered | sapply(data[cols], is.ordered)]
  cols.cont    <- cols[!(cols %in% cols.ordered) & vapply(data[cols], is.numeric, TRUE)]
  ncat <- vapply(cols.ordered, FUN.VALUE = numeric(1L), FUN = \(x) length(unique(data[[x]])))

  warnif(any(ncat <= 2L) && !probit.correction,
         "Some ordered variables only have two categories.\n",
         "Consider passing `ordered.probit.correction=TRUE`")

  ordinalize <- function(x) {
    if   (is.ordered(x)) return(as.ordered(as.integer(x))) # re-order levels, e.g., [2, 3, 4] -> [1, 2, 3]
    else                 return(as.ordered(x))
  }

  data.x <- data
  data.x[cols.ordered] <- lapply(data.x[cols.ordered], FUN = ordinalize)

  if (verbose) printf("Estimating factor scores...\n")
  syntax.cfa <- parTableToSyntax(parTable.cfa)

  if (is.null(group)) group.vals <- rep(1L, NROW(data.x))
  else                group.vals <- data[[group]]

  data.s   <- NULL
  RESCALED <- list()

  for (g in unique(group.vals)) {
    data.xg <- data.x[group.vals == g, , drop = FALSE]

    fit.cfa <- lavaan::cfa(syntax.cfa, data = data.xg, ordered = cols.ordered, se = "none")
    fscores <- as.data.frame(lavaan::lavPredict(fit.cfa))

    rescaled.g <- rescaleOrderedData_OP(
      data            = data.xg,
      cols.ordered    = cols.ordered,
      parTable.in     = parTable.in,
      fscores         = fscores,
      thresholds.vcov = TRUE
    )

    RESCALED[[g]] <- rescaled.g
    data.sg <- rescaled.g$data[cols]

    pcorr <- polychor(vars = cols, data = data.xg,
                      thresholds = NULL, ordered = cols.ordered)

    if (probit.correction)
      data.sg <- fitX2Cov(X = data.sg, S = pcorr)

    data.sg <- as.data.frame(data.sg)

    if (!is.null(group))
      data.sg[[group]] <- g

    data.s <- rbind(data.s, data.sg)
  }

  fit.out <- modsem_da(model.syntax     = model.syntax,
                       data             = data.s,
                       method           = method,
                       verbose          = verbose,
                       calc.se          = calc.se,
                       mean.observed    = ordered.mean.observed,
                       group            = group,
                       ...)

  for (g in seq_along(RESCALED)) {
    rescaled.g <- RESCALED[[g]]

    vcov.t <- NULL
    coef.t <- NULL

    for (col in cols.ordered) {
      if (!is.null(rescaled.g$vcov[[col]]) && length(rescaled.g$vcov[[col]]))
        vcov.t <- diagPartitionedMat(vcov.t, rescaled.g$vcov[[col]])

      if (!is.null(rescaled.g$thresholds[[col]]) && length(rescaled.g$thresholds[[col]]))
        coef.t <- c(coef.t, rescaled.g$thresholds[[col]])
    }

    parTable <- fit.out$parTable
    cols.t <- c("lhs", "op", "rhs", "group", "label", "est", "std.error")
    cols.m <- setdiff(colnames(parTable), cols.t)

    if (length(coef.t)) {
      lhs.t <- stringr::str_split_i(names(coef.t), pattern = "\\|", i = 1L)
      rhs.t <- stringr::str_split_i(names(coef.t), pattern = "\\|", i = 2L)

      parTable.t <- data.frame(lhs = lhs.t, op = "|", rhs = rhs.t, label = "",
                               group = g, est = coef.t,
                               std.error = if (!is.null(vcov.t)) sqrt(diag(vcov.t)) else NA_real_)
      parTable.t[cols.m] <- NA

      parTable.t <- addZStatsParTable(parTable.t)
      parTable   <- sortParTableDA(rbind(parTable, parTable.t), model = fit.out$model)
      fit.out$coef.all <- c(fit.out$coef.all, coef.t)

      if (!is.null(vcov.t))
        fit.out$vcov.all <- diagPartitionedMat(fit.out$vcov.all, vcov.t)
    }
  }

  fit.out$parTable <- parTable
  fit.out
}


# Closed-form CE for ordered-probit (vectorized)
# y: integer categories 1..K, eta: linear predictor, tau: interior thresholds (length K-1)
calcCondExpOrdProbit <- function(y, eta, tau, eps = 1e-12) {
  K <- max(y, na.rm = TRUE)
  a <- c(-Inf, tau)[y]
  b <- c(tau,  Inf)[y]
  al <- a - eta
  bl <- b - eta
  num <- stats::dnorm(al) - stats::dnorm(bl)
  den <- pmax(stats::pnorm(bl) - stats::pnorm(al), eps)
  eta + num / den
}


imputeWithOrdProbit <- function(y.ord, X, standardize = TRUE, weights = NULL,
                                std.thresholds = TRUE,
                                vcov.method = c("delta", "simple", "none")) {
  vcov.method <- match.arg(vcov.method)

  # 1) Fit ordered probit with Hessian
  y.fac <- if (is.factor(y.ord) || is.ordered(y.ord)) y.ord else factor(y.ord, ordered = TRUE)
  df    <- as.data.frame(X)
  fit   <- suppressWarnings(MASS::polr(y.fac ~ ., data = df, weights = weights,
                                       method = "probit", Hess = TRUE))

  # 2) Extract coefficients and rebuild design as used by polr
  beta <- stats::coef(fit)                 # length P
  tau_probit <- as.numeric(fit$zeta)       # length Q = K-1
  K    <- length(levels(y.fac)); Q <- K - 1L
  yint <- as.integer(y.fac)

  Xmm_all <- stats::model.matrix(~ ., data = df)
  Xmm     <- if (length(beta)) Xmm_all[, names(beta), drop = FALSE] else
             matrix(0, nrow = nrow(df), ncol = 0) # degenerate
  eta     <- as.numeric(Xmm %*% beta)
  n       <- length(eta)

  # weights handling
  if (is.null(weights)) {
    w <- rep(1, n); wsum <- n
  } else {
    if (length(weights) != n) stop("weights must have length nrow(X).")
    w <- as.numeric(weights); wsum <- sum(w)
  }

  # 3) CE (conditional expectations) on probit scale; used only for imputed values
  a <- c(-Inf, tau_probit); b <- c(tau_probit,  Inf)
  eps <- 1e-12
  m <- matrix(NA_real_, nrow = n, ncol = K,
              dimnames = list(NULL, levels(y.fac)))
  for (j in seq_len(K)) {
    alpha <- (a[j] - eta)
    betaJ <- (b[j] - eta)
    numer <- stats::dnorm(alpha) - stats::dnorm(betaJ)
    denom <- stats::pnorm(betaJ) - stats::pnorm(alpha); denom <- pmax(denom, eps)
    m[, j] <- eta + numer / denom
  }

  # Aggregate to a single CE constant per observed category
  mu_cat <- numeric(K)
  for (j in seq_len(K)) {
    idx <- which(yint == j)
    if (!length(idx)) {
      pj <- stats::pnorm(b[j] - eta) - stats::pnorm(a[j] - eta)
      mu_cat[j] <- sum(m[, j] * pj * w) / max(sum(pj * w), eps)
    } else {
      wj <- w[idx]
      mu_cat[j] <- sum(m[idx, j] * wj) / sum(wj)
    }
  }
  names(mu_cat) <- levels(y.fac)
  y_ce <- mu_cat[yint]
  if (standardize) y_ce <- std1(y_ce)

  # 4) Thresholds: rescale from probit (Var eps =1) to SEM unit-Var(y*)
  #    y* variance = Var(eta) + 1 ; mean = E[eta]
  mu_eta <- sum(w * eta) / wsum
  # weighted variance (population form)
  var_eta <- sum(w * (eta - mu_eta)^2) / wsum
  s <- sqrt(var_eta + 1)

  tau.ystar <- if (std.thresholds) (tau_probit - mu_eta) / s else tau_probit

  # 5) Vcov for thresholds
  vcov.full <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  vcov.tau.probit <- vcov.tau.ystar <- NULL
  if (!is.null(vcov.full) && Q > 0) {
    P <- length(beta)
    vcov.tau.probit <- as.matrix(vcov.full[P + seq_len(Q), P + seq_len(Q), drop = FALSE])

    if (vcov.method == "none") {
      vcov.tau.ystar <- NULL

    } else if (!std.thresholds || s <= 0) {
      # no rescale -> same vcov; if s invalid keep NULL
      vcov.tau.ystar <- vcov.tau.probit

    } else if (vcov.method == "simple") {
      # only scale by 1/s^2, ignore uncertainty in mu_eta and s
      vcov.tau.ystar <- vcov.tau.probit / (s^2)

    } else if (vcov.method == "delta") {
      # Full delta: tau' = (tau - mu_eta)/s, with mu_eta = E[X]^T beta,
      # s = sqrt(beta^T Sx beta + 1). Build J and compute J V J^T.
      EX <- if (P) colSums(Xmm * w) / wsum else numeric(0)
      # weighted covariance of X columns
      if (P) {
        Xc <- sweep(Xmm, 2, EX, FUN = "-")
        Sx <- (t(Xc * w) %*% Xc) / wsum   # P x P
        Sb <- as.numeric(Sx %*% beta)     # P-vector: Sx * beta
      } else {
        Sx <- matrix(0, 0, 0); Sb <- numeric(0)
      }

      J <- matrix(0, nrow = Q, ncol = P + Q)  # rows: thresholds; cols: [beta, tau]
      if (P) {
        # d tau'_k / d beta = -(EX)/s - (tau_k - mu_eta) * (Sx beta) / s^3
        for (k in seq_len(Q)) {
          J[k, 1:P] <- -(EX) / s - (tau_probit[k] - mu_eta) * (Sb) / (s^3)
        }
      }
      # d tau'_k / d tau_j = (1/s) if j==k else 0
      for (k in seq_len(Q)) J[k, P + k] <- 1 / s

      V <- as.matrix(vcov.full)
      vcov.tau.ystar <- J %*% V %*% t(J)
      # ensure symmetry
      vcov.tau.ystar <- 0.5 * (vcov.tau.ystar + t(vcov.tau.ystar))
    }
  }

  # Named outputs for thresholds (we’ll re-prefix with indicator name upstream)
  names(tau_probit) <- paste0("t", seq_len(Q))
  names(tau.ystar)  <- paste0("t", seq_len(Q))
  if (!is.null(vcov.tau.probit)) dimnames(vcov.tau.probit) <- list(names(tau_probit), names(tau_probit))
  if (!is.null(vcov.tau.ystar))  dimnames(vcov.tau.ystar)  <- list(names(tau.ystar),  names(tau.ystar))

  list(
    values          = y_ce,              # imputed y* CE values (you standardize again downstream if desired)
    mu_per_cat      = mu_cat,
    tau_probit      = tau_probit,        # thresholds on probit-error scale
    vcov.tau.probit = vcov.tau.probit,
    tau.ystar       = tau.ystar,         # thresholds on SEM y* unit-variance scale
    vcov.tau.ystar  = vcov.tau.ystar,    # vcov for tau.ystar (delta/simple/none)
    eta             = eta,
    mu_eta          = mu_eta,
    s_ystar         = s,                 # scaling used: s = sqrt(Var(eta)+1)
    fit             = fit
  )
}


.buildX_ForIndicator <- function(latent, rhsTerms, fscores) {
  X <- list()  # don't include latent self
  stopif(!length(rhsTerms), sprintf("expected predictors for eta (%s)!", latent))

  for (trm in rhsTerms) {
    if (grepl(":", trm, fixed = TRUE)) {
      parts <- strsplit(trm, ":", fixed = TRUE)[[1]]
      stopifnot(all(parts %in% colnames(fscores)))
      X[[trm]] <- Reduce(`*`, lapply(parts, function(p) fscores[[p]]))

    } else {
      stopifnot(trm %in% colnames(fscores))
      X[[trm]] <- fscores[[trm]]
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
  meas    <- parTable.in[parTable.in$op == "=~", , drop = FALSE]
  ind2lat <- stats::setNames(as.character(meas$lhs), meas$rhs)

  # structural terms per endogenous latent
  struc       <- parTable.in[parTable.in$op == "~", , drop = FALSE]
  rhsByLatent <- tapply(struc$rhs, struc$lhs, c, simplify = FALSE)

  etas <- getEtas(parTable.in, checkAny = FALSE)
  ordered_endog <- intersect(cols.ordered, names(ind2lat)[ind2lat %in% etas])

  thresholds <- stats::setNames(vector("list", length(cols.ordered)), nm = cols.ordered)
  vcovs      <- stats::setNames(vector("list", length(cols.ordered)), nm = cols.ordered)

  for (col in cols.ordered) {
    if (col %in% ordered_endog) {
      # indicator belongs to an endogenous latent -> ordered-probit CE
      latent   <- ind2lat[[col]]
      rhsTerms <- rhsByLatent[[latent]]

      if (is.null(rhsTerms)) rhsTerms <- character(0)

      X   <- .buildX_ForIndicator(latent, rhsTerms, fscores)
      imp <- imputeWithOrdProbit(y.ord = data[[col]], X = X, standardize = TRUE)

      out[[col]] <- imp$values
      if (length(imp$tau.ystar)) {
        tn <- paste0(col, "|t", seq_along(imp$tau.ystar))
        thresholds[[col]] <- stats::setNames(imp$tau.ystar, tn)

        if (!is.null(imp$vcov.tau.ystar)) {
          dimnames(imp$vcov.tau.ystar) <- list(tn, tn)
          vcovs[[col]] <- imp$vcov.tau.ystar

        } else vcovs[[col]] <- NULL

      } else {
        thresholds[[col]] <- numeric(0); vcovs[[col]] <- NULL
      }

    } else {
      # exogenous or non-endogenous ordered indicator -> analytic normal mapping
      ana <- rescaleOrderedVariableAnalytic(name = col, data = data.cat,
                                            thresholds.vcov = thresholds.vcov, ...)
      out[[col]]        <- ana$values
      thresholds[[col]] <- ana$thresholds
      vcovs[[col]]      <- ana$vcov
    }
  }

  list(data = out, thresholds = thresholds, vcov = vcovs)
}


fitX2Cov <- function(X, S) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  stopif(!all(dim(S) == c(p, p)), "S must be a square matrix with dimensions equal to the number of columns in X")

  # Compute sample mean and covariance of X
  mu_X <- colMeans(X)
  X_centered <- sweep(X, 2, mu_X, FUN = "-")

  Sigma_X <- stats::cov(X)

  # Check positive definiteness
  if (any(eigen(Sigma_X, symmetric = TRUE)$values <= 0)) {
    warning("Sample covariance matrix of X is not positive-definite. Consider regularization.")
    return(X)
  }
  if (any(eigen(S, symmetric = TRUE)$values <= 0)) {
    warning("Target covariance matrix S must be positive-definite.")
    return(X)
  }

  # Matrix square root function using eigen-decomposition
  matrix_sqrt <- function(M) {
    eig <- eigen(M, symmetric = TRUE)
    eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  }

  # Compute inverse square root of Sigma_X
  Sigma_X_inv_sqrt <- solve(matrix_sqrt(Sigma_X))

  # Compute square root of S
  S_sqrt <- matrix_sqrt(S)

  # Transformation matrix A
  A <- Sigma_X_inv_sqrt %*% S_sqrt

  # Transform the centered data
  Y_centered <- X_centered %*% A

  # Add target mean back
  Y <- sweep(Y_centered, 2, mu_X, FUN = "+")

  colnames(Y) <- colnames(X)
  as.data.frame(Y)
}


polychor <- function(vars, data, thresholds = NULL, ordered = NULL) {
  combos <- getUniqueCombos(vars, match = FALSE)
  syntax <- paste0(sprintf("%s~~%s", combos$V1, combos$V2), collapse = "\n")

  for (var in names(thresholds)) {
    tau  <- thresholds[[var]]
    ftau <- format(tau[is.finite(tau)], scientific = FALSE)

    constr <- sprintf("%s|%s*t%d", var, ftau, seq_along(tau))
    syntax <- paste(syntax, paste0(constr, collapse = "\n"), sep = "\n")
  }

  lavaan::lavInspect(lavaan::sem(syntax, data, ordered = ordered, se = "none"),
                     what = "cov.ov")
}
