modsemOrderedScaleCorrection <- function(..., type = c("simple", "rubin")) {
  type <- tolower(type)
  type <- match.arg(type)

  message("Scale correcting ordinal variables. ",
          "This is an experimental feature!\n",
          "See `help(modsem_da)` for more information.")

  switch(type,
         simple = modsemOrderedScaleCorrectionSimple(...),
         rubin  = modsemOrderedScaleCorrectionRubin(...),
         stop2("Unrecognized type: ", type))
}


modsemOrderedScaleCorrectionSimple <- function(model.syntax,
                                               data,
                                               method = "lms",
                                               ordered = NULL,
                                               calc.se = TRUE,
                                               R = 50L,
                                               warmup = floor(R / 2L),
                                               N = 1e5, # 100,000
                                               verbose = interactive(),
                                               optimize = TRUE,
                                               start = NULL,
                                               lambda = 1,
                                               ordered.tol = 1e-6,
                                               ...) {
  if (is.null(verbose))
    verbose <- TRUE # default

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]

  rescaleOrderedData <- function(data, sim.ov) {
    data.y <- data

    thresholds <- stats::setNames(vector("list", length(cols.ordered)),
                                  nm = cols.ordered)
    for (col in cols.ordered) {
      scaling <- rescaleOrderedVariable(name = col, data = data, sim.ov = sim.ov)
      data.y[[col]]     <- scaling$values
      thresholds[[col]] <- scaling$thresholds
    }

    list(data = data.y, thresholds = thresholds)
  }

  ERROR <- \(e) {warning2(e, immediate. = FALSE); NULL}

  data.x <- data
  data.y <- data
  data.y[cols.ordered] <- lapply(data.y[cols.ordered], FUN = as.integer)

  thresholds <- stats::setNames(lapply(cols.ordered,  FUN = \(x) 0),
                                nm = cols.ordered)

  prev.theta <- NULL
  for (r in seq_len(R)) {
    printedLines <- utils::capture.output(split = TRUE, {
      if (verbose) printf("Bootstrapping scale-correction %d/%d...\n", r, R)

      fit.naive <- modsem(
        model.syntax = model.syntax,
        data         = data.y,
        method       = method,
        calc.se      = FALSE,
        verbose      = verbose,
        optimize     = if (is.null(prev.theta) || r <= 3L) optimize else FALSE,
        start        = if (is.null(prev.theta) || r <= 3L) start else prev.theta,
        ...
      )

      new.theta <- fit.naive$theta

      sim <- simulateDataParTable(
        parTable = parameter_estimates(fit.naive),
        N        = N
      )

      rescaled <- rescaleOrderedData(data = data.x, sim.ov = sim$oV)
      data.y.new     <- rescaled$data
      thresholds.new <- rescaled$thresholds

      # damping to promote contraction
      lambda.i <- if (r > warmup) lambda/(r - warmup) else 1 # no mixing before warmup
      for (c in cols.ordered) {
        data.y[[c]]     <- (1 - lambda.i) * data.y[[c]] +
                                lambda.i  * data.y.new[[c]]
        thresholds[[c]] <- (1 - lambda.i) * thresholds.new[[c]] +
                                lambda.i  * thresholds.new[[c]]
      }
    })

    nprinted <- length(printedLines)
    eraseConsoleLines(nprinted)

    if (!is.null(prev.theta)) {
      eucdist   <- \(x) sqrt(c(t(x) %*% x))
      normalize <- \(x) x / eucdist(x)

      epsilon <- normalize(prev.theta) - normalize(new.theta)

      delta <- eucdist(epsilon)
      if (r > warmup & delta <= ordered.tol) {
        printf("Scale correction converged (iter = %d, warmup = %d, tau = %.2g)\n",
               r, warmup, delta)
        break
      }
    }

    prev.theta <- new.theta
  }

  out <- modsem(
    model.syntax = model.syntax,
    method       = method,
    data         = data.y,
    calc.se      = calc.se,
    verbose      = verbose,
    start        = new.theta,
    optimize     = FALSE,
    ...
  )

  parTable <- parameter_estimates(out)
  na.cols <- setdiff(colnames(parTable), c("lhs", "op", "rhs", "est", "label"))

  for (col in cols.ordered) {
    t <- thresholds[[col]]
    t <- t[is.finite(t)]

    newRows <- data.frame(lhs = col,
                          op  = "|",
                          rhs = paste0("t", seq_along(t)),
                          est = t, label = "")
    newRows[na.cols] <- NA

    parTable <- rbind(parTable, newRows)
  }

  out$parTable <- sortParTableDA(parTable, model = out$model)
  out$args$optimize <- optimize # restore input arguments
  out$args$start    <- start

  out
}


modsemOrderedScaleCorrectionRubin <- function(model.syntax,
                                              data,
                                              method = "lms",
                                              ordered = NULL,
                                              calc.se = TRUE,
                                              R = 50,
                                              warmup = floor(R / 2L),
                                              verbose = interactive(),
                                              optimize = TRUE,
                                              start = NULL,
                                              lambda = 1,
                                              ordered.tol = 1e-6,
                                              ...) {
  if (is.null(verbose))
    verbose <- TRUE # default

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]

  rescaleOrderedData <- function(data, sim.ov) {
    data.y <- data

    thresholds <- stats::setNames(vector("list", length(cols.ordered)),
                                  nm = cols.ordered)
    for (col in cols.ordered) {
      scaling <- rescaleOrderedVariable(name = col, data = data, sim.ov = sim.ov)
      data.y[[col]]     <- scaling$values
      thresholds[[col]] <- scaling$thresholds
    }

    list(data = data.y, thresholds = thresholds)
  }

  ERROR <- \(e) {warning2(e, immediate. = FALSE); NULL}

  data.x <- data
  data.y <- data
  data.y[cols.ordered] <- lapply(data.y[cols.ordered], FUN = as.integer)

  prev.theta <- NULL
  COEF.ALL  <- vector("list", length = warmup - 1L)
  COEF.FREE <- vector("list", length = warmup - 1L)

  thresholds <- stats::setNames(lapply(cols.ordered,  FUN = \(x) 0),
                                nm = cols.ordered)
  for (r in seq_len(R)) {
    printedLines <- utils::capture.output(split = TRUE, {
      if (verbose) printf("Bootstrapping scale-correction %d/%d...\n", r, R)

      fit.naive <- modsem(
        model.syntax = model.syntax,
        data         = data.y,
        method       = method,
        calc.se      = FALSE,
        verbose      = verbose,
        optimize     = if (is.null(prev.theta) || r <= 3L) optimize else FALSE,
        start        = if (is.null(prev.theta) || r <= 3L) start else prev.theta,
        ...
      )

      new.theta <- fit.naive$theta

      sim <- simulateDataParTable(
        parTable = parameter_estimates(fit.naive),
        N        = NROW(data.y)
      )

      rescaled <- rescaleOrderedData(data = data.x, sim.ov = sim$oV)
      data.y.new     <- rescaled$data
      thresholds.new <- rescaled$thresholds

      # damping to promote contraction
      lambda.i <- if (r > warmup) lambda/(r - warmup) else 1 # no mixing before warmup
      for (c in cols.ordered) {
        data.y[[c]]     <- (1 - lambda.i) * data.y[[c]] +
                                lambda.i  * data.y.new[[c]]
        thresholds[[c]] <- (1 - lambda.i) * thresholds.new[[c]] +
                                lambda.i  * thresholds.new[[c]]
      }

      if (!is.null(prev.theta)) {
        eucdist   <- \(x) sqrt(c(t(x) %*% x))
        normalize <- \(x) x / eucdist(x)

        epsilon <- normalize(prev.theta) - normalize(new.theta)
        # normalize

        delta <- eucdist(epsilon)
        if (r > warmup & delta <= ordered.tol) {
          break
        }
      }

      if (r > 1L && r <= warmup) { # accumulate to correct std-errors
        COEF.ALL[[r-1L]]  <- coef(fit.naive, type = "all")
        COEF.FREE[[r-1L]] <- coef(fit.naive, type = "free")
      }

      prev.theta <- new.theta
    })

    nprinted <- length(printedLines)
    eraseConsoleLines(nprinted)
  }

  fit.final <- modsem(
    model.syntax = model.syntax,
    method       = method,
    data         = data.y,
    calc.se      = calc.se,
    verbose      = verbose,
    start        = new.theta,
    optimize     = FALSE,
    ...
  )

  vcov.naive.all  <- vcov(fit.final, type = "all") # final fit, but naive vcov
  vcov.naive.free <- vcov(fit.final, type = "free")

  VCOV.ALL  <- vector("list", length = warmup - 1L)
  VCOV.FREE <- vector("list", length = warmup - 1L)

  for (i in seq_len(warmup - 1L)) {
    VCOV.ALL[[i]]  <- vcov.naive.all
    VCOV.FREE[[i]] <- vcov.naive.free
  }

  pool.all  <- rubinPool(COEF.ALL,  VCOV.ALL)
  pool.free <- rubinPool(COEF.FREE, VCOV.FREE)

  vcov.all  <- pool.all$Tvcov # naive vcov adjusted for variance in coefs
  vcov.free <- pool.free$Tvcov

  fit.final$vcov.all  <- vcov.all
  fit.final$vcov.free <- vcov.free
  fit.final$FIM       <- solve(vcov.free)

  parTable1   <- parameter_estimates(fit.final)
  orig.labels <- parTable1$label
  parTable1   <- getMissingLabels(parTable1)
  parTableT   <- data.frame(label = colnames(vcov.all),
                            std.error.t = sqrt(diag(vcov.all)))

  parTable <- leftJoin(left = parTable1, right = parTableT, by = "label")
  match    <- !is.na(parTable$std.error.t)

  parTable$std.error[match] <- parTable$std.error.t[match]
  parTable$std.error.t <- NULL
  parTable$label[!parTable$label %in% orig.labels] <- ""

  na.cols <- setdiff(colnames(parTable), c("lhs", "op", "rhs", "est", "label"))

  for (col in cols.ordered) {
    t <- thresholds[[col]]
    t <- t[is.finite(t)]

    newRows <- data.frame(lhs = col,
                          op  = "|",
                          rhs = paste0("t", seq_along(t)),
                          est = t, label = "")
    newRows[na.cols] <- NA

    parTable <- rbind(parTable, newRows)
  }

  fit.final$parTable <- sortParTableDA(parTable, model = fit.final$model)
  fit.final$args$optimize <- optimize # restore input arguments
  fit.final$args$start    <- start

  fit.final
}


rescaleOrderedVariable <- function(name, data, sim.ov,
                                   smooth_eps = 0,
                                   eps_expand = 1e-12) {
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

  # exact centering to remove residual drift (weighted by observed category probs)
  mu <- mu - sum(p_obs * mu, na.rm = TRUE)

  # map each observed category to its conditional mean
  x.out <- rep(NA_real_, length(x))
  for (i in seq_len(K))
    x.out[x.i == i] <- mu[i]

  list(values = x.out, thresholds = q)
}
