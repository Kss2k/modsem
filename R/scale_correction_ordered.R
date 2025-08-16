modsemOrderedScaleCorrection <- function(model.syntax,
                                         data,
                                         method = "lms",
                                         ordered = NULL,
                                         calc.se = TRUE,
                                         iter = 75L,
                                         warmup = 25L,
                                         N.start = max(NROW(data), 1e4), # 10,000 initally
                                         N.max   = max(NROW(data), 1e6), # 1,000,000
                                         target.coverage = 0.7,
                                         lambda = 1,
                                         se = "simple",
                                         ordered.standardize.data = TRUE,
                                         ordered.mean.observed = FALSE,
                                         # Capture args
                                         verbose = interactive(),
                                         optimize = TRUE,
                                         start = NULL,
                                         standardize.data = NULL,
                                         mean.observed    = NULL,
                                         ...) {
  message("Sampling values for ordered variables, ",
          "this is an experimental feature!\n",
          "Consider increasing the `ordered.iter` and `ordered.warmup` arguments!")

  standardize <- \(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  if (is.null(verbose))
    verbose <- TRUE # default

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]
  cols.cont    <- cols[!(cols %in% cols.ordered) & vapply(data, is.numeric, TRUE)]

  ERROR <- \(e) {warning2(e, immediate. = FALSE); NULL}

  data.x <- data
  data.y <- data
  data.y[cols.ordered] <- lapply(data.y[cols.ordered], FUN = as.integer)

  thresholds <- stats::setNames(lapply(cols.ordered,  FUN = \(x) 0),
                                nm = cols.ordered)

  stopif(iter <= warmup, "`ordered.iter` must be larger than `ordered.warmup`!")

  iter.keep <- max(1L, floor(iter - warmup))
  COEF.ALL  <- vector("list", length = iter.keep)
  COEF.FREE <- vector("list", length = iter.keep)
  VCOV.ALL  <- vector("list", length = iter.keep)
  VCOV.FREE <- vector("list", length = iter.keep)
  fits      <- vector("list", length = iter.keep)
  imputed   <- vector("list", length = iter.keep)

  data_i <- data
  data_i[cols.ordered] <- lapply(data_i[cols.ordered],
                                 FUN = \(x) standardize(as.integer(x)))

  THETA     <- NULL
  PARS      <- NULL
  N         <- N.start
  W         <- 0       # cumulative weight
  maxLL_run <- -Inf

  for (i in seq_len(iter)) {
    printedLines <- utils::capture.output(split = TRUE, {
      j <- i - warmup
      mode <- if (j <= 0) "warmup" else "sampling"

      if (verbose)
        printf("Iterations %d/%d [%s]...\n", i, iter, mode)

      if (is.null(calc.se))
        calc.se_i <- ifelse(j < 1L || (se == "simple" && j > 1L), yes = FALSE, no = TRUE)
      else
        calc.se_i <- calc.se

      if (is.null(THETA)) {
        optimize <- TRUE
        start    <- NULL
      } else {
        optimize <- FALSE
        start    <- apply(THETA, MARGIN = 2, FUN = mean, na.rm = TRUE) # na.rm should't be necessary, ever...
      }

      if (j > 0L) imputed[[j]] <- data_i

      fit_i <- tryCatch(
        modsem_da(
          model.syntax     = model.syntax,
          data             = data_i,
          method           = method,
          start            = start,
          verbose          = verbose,
          optimize         = optimize,
          calc.se          = calc.se_i,
          standardize.data = ordered.standardize.data,
          mean.observed    = ordered.mean.observed,
          ...
        ), error = \(e) {cat("\n"); print(e); NULL}
      )

      if (is.null(fit_i)) next

      pars_i <- parameter_estimates(fit_i)

      if (is.null(PARS) || mode == "warmup") {
        PARS <- pars_i

      } else {
        ll <- as.numeric(fit_i$logLik)

        if (!is.finite(ll)) {
          # skip weighting when logLik is NA/Inf
          w <- 0

        } else {
          # if we have a new maximum, rescale cumulative weight to the new baseline
          if (ll > maxLL_run) {
            if (is.finite(maxLL_run)) {
              scale <- exp((maxLL_run - ll) / lambda)  # < 1
              W <- W * scale
            }
            maxLL_run <- ll
          }

          # softmax-style weight in (0, 1], peaked at current best
          w <- exp((ll - maxLL_run) / lambda)
        }

        # online likelihood-weighted average update
        if (w > 0) {
          W <- W + w
          PARS$est <- PARS$est + (w / W) * (pars_i$est - PARS$est)
        }
      }

      sim_i <- simulateDataParTable(
        parTable = PARS,
        N        = N
      )

      rescaled   <- rescaleOrderedData(data = data.x, sim.ov = sim_i$oV,
                                       cols.ordered = cols.ordered,
                                       cols.cont = cols.cont)
      coverage   <- rescaled$coverage
      n.scalef   <- target.coverage / coverage
      N          <- min(N.max, floor(n.scalef * N))
      data_i     <- rescaled$data
      thresholds <- rescaled$thresholds

      THETA <- rbind(
        THETA,
        matrix(fit_i$theta, nrow = 1, dimnames = list(NULL, names(fit_i$theta)))
      )

      if (j >= 1L) {
        fits[[j]]      <- fit_i
        COEF.FREE[[j]] <- coef(fit_i, type = "free")
        COEF.ALL[[j]]  <- addThresholdsCoef(coef(fit_i, type = "all"),
                                            thresholds = thresholds)

        if (j > 1L && se == "simple") {
          VCOV.ALL[[j]]  <- VCOV.ALL[[1]]
          VCOV.FREE[[j]] <- VCOV.FREE[[1]]

        } else if (calc.se_i) {
          VCOV.FREE[[j]]  <- vcov(fit_i, type = "free")
          VCOV.ALL[[j]] <- addThresholdsVcov(vcov(fit_i, type = "all"),
                                              thresholds = thresholds)

        } else {
          k.all  <- length(COEF.ALL[[j]]) + length(unlist(thresholds))
          k.free <- length(COEF.FREE[[j]])
          d.all  <- c(names(COEF.ALL[[j]]), getLabelsThresholds(thresholds))
          d.free <- names(COEF.FREE[[j]])

          VCOV.ALL[[j]]  <- matrix(0, nrow = k.all, ncol = k.all,
                                   dimnames = list(d.all, d.all))
          VCOV.FREE[[j]] <- matrix(0, nrow = k.free, ncol = k.free,
                                   dimnames = list(d.free, d.free))
        }
      }
    })

    nprinted <- length(printedLines)
    if (i < iter) eraseConsoleLines(nprinted)
  }

  failed <- vapply(fits, FUN.VALUE = logical(1L), FUN = is.null)

  if (any(failed)) {
    warning2(sprintf("Model estimation failed in %d out of %d impuations!",
                     sum(failed), R), immediate. = FALSE)

    fits      <- fits[!failed]
    COEF.ALL  <- COEF.ALL[!failed]
    COEF.FREE <- COEF.FREE[!failed]
    VCOV.ALL  <- VCOV.ALL[!failed]
    VCOV.FREE <- VCOV.FREE[!failed]
    R         <- sum(!failed)
  }

  pool.all  <- rubinPool(COEF.ALL,  VCOV.ALL)
  pool.free <- rubinPool(COEF.FREE, VCOV.FREE)

  coef.all  <- pool.all$theta.bar
  vcov.all  <- pool.all$Tvcov

  coef.free <- pool.free$theta.bar
  vcov.free <- pool.free$Tvcov

  # Re-do parameter estimates
  parTable1   <- addThresholdsParTable(parameter_estimates(fits[[1]]),
                                       thresholds = thresholds)
  orig.labels <- parTable1$label
  parTable1   <- getMissingLabels(parTable1)
  parTableT   <- data.frame(label = names(coef.all),
                            est.t = coef.all,
                            std.error.t = sqrt(diag(vcov.all)))

  parTable    <- leftJoin(left = parTable1, right = parTableT, by = "label")
  match       <- !is.na(parTable$est.t)

  parTable$est[match]       <- parTable$est.t[match]
  parTable$std.error[match] <- parTable$std.error.t[match]
  parTable$est.t       <- NULL
  parTable$std.error.t <- NULL
  parTable$label[!parTable$label %in% orig.labels] <- ""
  parTable <- parTable[c("lhs", "op", "rhs", "label", "est", "std.error")] # remove z-statistics
  parTable <- addZStatsParTable(parTable)


  parTable <- modsemParTable(sortParTableDA(parTable, model = fits[[1L]]$model))

  matrices    <- aggregateMatrices(fits, type = "main")
  covMatrices <- aggregateMatrices(fits, type = "cov")
  expected.matrices <- aggregateMatrices(fits, type = "expected")

  getScalarFit <- function(fit, field, dtype = numeric)
    vapply(fits, FUN.VALUE = dtype(1L), \(fit) fit[[field]])

  fit.out <- fits[[1]]
  fit.out$coefs.all       <- coef.all
  fit.out$coefs.free      <- coef.free
  fit.out$vcov.all        <- vcov.all
  fit.out$vcov.free       <- vcov.free
  fit.out$parTable        <- parTable
  fit.out$information     <- sprintf("Rubin-corrected (m=%d)", iter.keep)
  fit.out$FIM             <- solve(vcov.free)
  fit.out$theta           <- apply(THETA, MARGIN = 2, FUN = mean, na.rm = TRUE)
  fit.out$iterations      <- sum(getScalarFit(fits, field = "iterations"))
  fit.out$logLik          <- mean(getScalarFit(fits, field = "logLik"))
  fit.out$convergence     <- all(getScalarFit(fits, field = "convergence",
                                              dtype = logical))
  fit.out$convergence.msg <- getConvergenceMessage(fit.out$convergence,
                                                   fit.out$iterations)
  fit.out$model$matrices          <- matrices
  fit.out$model$covModel$matrices <- covMatrices
  fit.out$expected.matrices       <- expected.matrices

  fit.out$imputations <- list(fitted = fits, data = imputed)

  # restore passed arguments
  fit.out$args$optimize <- optimize
  fit.out$args$start    <- NULL

  fit.out
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
  y.sample <- y[sample(NROW(sim.ov), NROW(data))] # sample y to get appropriate sampling
                                                  # variation of thresholds
  q <- stats::quantile(y.sample,
                       probs = cdf_obs, names = FALSE, type = 7)
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


addThresholdsParTable <- function(parTable, thresholds) {
  na.cols <- setdiff(colnames(parTable), c("lhs", "op", "rhs", "est", "label"))

  for (col in names(thresholds)) {
    t <- thresholds[[col]]
    t <- t[is.finite(t)]

    newRows <- data.frame(lhs = col,
                          op  = "|",
                          rhs = paste0("t", seq_along(t)),
                          est = t, label = "")
    newRows[na.cols] <- NA

    parTable <- rbind(parTable, newRows)
  }

  parTable
}


addThresholdsCoef <- function(coef, thresholds) {
  for (col in names(thresholds)) {
    t <- thresholds[[col]]
    t <- t[is.finite(t)]

    names(t) <- paste0(col, "|t", seq_along(t))
    coef <- c(coef, t)
  }

  coef
}


getLabelsThresholds <- function(thresholds) {
  labels <- NULL
  for (col in names(thresholds)) {
    t <- thresholds[[col]]
    t <- t[is.finite(t)]

    labels <- c(labels, paste0(col, "|t", seq_along(t)))
  }

  labels
}


addThresholdsVcov <- function(vcov, thresholds) {
  for (col in names(thresholds)) {
    t <- thresholds[[col]]
    t <- t[is.finite(t)]
    k <- length(t)

    lab    <- paste0(col, "|t", seq_along(t))
    vcov.t <- matrix(0, nrow = k, ncol = k, dimnames = list(lab, lab))
    vcov   <- diagPartitionedMat(vcov, vcov.t)
  }

  vcov
}


rescaleOrderedData <- function(data, sim.ov, cols.ordered, cols.cont) {
  if (!length(cols.ordered)) # nothing to do
    return(cols.ordered)

  if (length(cols.cont)) {
    obsContScaled <- scale(as.matrix(data[, cols.cont, drop = FALSE]))
    contCenter <- attr(obsContScaled, "scaled:center")
    contScale  <- attr(obsContScaled, "scaled:scale")
    use.kNN <- TRUE
    kNN_k   <- 25L
  } else {
    obsContScaled <- NULL
    contCenter <- contScale <- numeric(0L)
    use.kNN <- FALSE
    kNN_k   <- 0L
  }

  # ensure column order consistent with data
  sim.mat <- as.data.frame(sim.ov)[, colnames(data), drop = FALSE]

  # make a full output copy right away, so continuous vars are preserved
  out <- as.data.frame(data, stringsAsFactors = FALSE)

  data.cat <- as.data.frame(data)
  data.cat[cols.ordered] <- lapply(data.cat[cols.ordered], as.integer)

  # thresholds + univariate conditional means for each ordered variable
  NamedList  <- \(nm) stats::setNames(vector("list", length(nm)), nm = nm)
  thresholds <- NamedList(cols.ordered)
  univarEst  <- NamedList(cols.ordered)

  std1 <- function(v) {
    mu    <- mean(v, na.rm = TRUE)
    sigma <- stats::sd(v, na.rm = TRUE)

    if (!is.finite(sigma) || sigma == 0)
      sigma <- 1

    (v - mu) / sigma
  }

  # standardized versions of the ordered cols in sim (these are the values we copy)
  sim.std.ord <- lapplyNamed(cols.ordered, \(col) std1(sim.mat[[col]]),
                             names = cols.ordered)
  sim.std.ord <- as.data.frame(sim.std.ord)

  sim.cat <- sim.mat
  for (col in cols.ordered) {
    scaling <- rescaleOrderedVariable(name = col, data = data.cat, sim.ov = sim.mat)
    univarEst[[col]]  <- scaling$values
    thresholds[[col]] <- scaling$thresholds

    # categories via findInterval for speed (on standardized scale)
    sim.cat[[col]] <- findInterval(sim.std.ord[[col]],
                                   vec = scaling$thresholds,
                                   left.open = TRUE, rightmost.closed = TRUE)
  }

  # start with the univariate estimates for ordered variables only
  multivarEst <- as.data.frame(univarEst)

  # ordered-combo keys
  obs.code <- interaction(data.cat[cols.ordered], drop = TRUE, lex.order = TRUE)
  sim.code <- interaction(sim.cat[cols.ordered],  drop = TRUE, lex.order = TRUE)
  combos.both <- intersect(levels(obs.code), levels(sim.code))
  combos.miss <- setdiff(levels(obs.code), levels(sim.code))

  denominator <- length(levels(obs.code))
  numerator   <- denominator - length(combos.miss)
  coverage <- numerator / denominator

  # pre-split indices
  obs.splits <- split(seq_len(nrow(data.cat)), obs.code, drop = TRUE)
  sim.splits <- split(seq_len(nrow(sim.cat)),  sim.code, drop = TRUE)

  # precompute standardized continuous block for sim (to match observed)
  if (length(cols.cont)) {
    X <- as.matrix(sim.mat[, cols.cont, drop = FALSE])

    simContScaled <- as.matrix(
      sweep(sweep(X, 2, contCenter, "-"), 2, contScale, "/")
    )
  }

  for (combo in combos.both) {
    oi <- obs.splits[[combo]]
    si <- sim.splits[[combo]]
    if (!length(oi) || !length(si)) next

    # kNN within combo using continuous variables (if available)
    # I.e., if there are continuous variables available we also want
    # to incorporate their information, when sampling continuous values
    # for the ordinal variables
    if (use.kNN && length(cols.cont)) {
      Oc <- obsContScaled[oi, , drop = FALSE]
      Sc <- simContScaled[si, , drop = FALSE]

      if (ncol(Oc) && nrow(Sc)) {
        k_eff <- max(1L, min(kNN_k, nrow(Sc)))
        picked <- integer(length(oi))

        nn <- RANN::nn2(Sc, Oc, k = k_eff, searchtype = "standard")
        for (r in seq_along(oi)) {
          nei <- nn$nn.idx[r, seq_len(k_eff)]
          d   <- nn$nn.dists[r, seq_len(k_eff)]
          w   <- 1 / pmax(d, 1e-8)  # inverse-distance weights
          picked[r] <- si[ nei[ sample.int(k_eff, 1L, prob = w) ] ]
        }

        # write replacements for ordered columns only
        multivarEst[oi, cols.ordered] <- sim.std.ord[picked, cols.ordered, drop = FALSE]
        next
      }
    }

    # fallback: uniform within combo
    k <- length(oi)
    m <- length(si)
    M <- si[sample.int(m, k, replace = (m < k))]
    multivarEst[oi, cols.ordered] <- sim.std.ord[M, cols.ordered, drop = FALSE]
  }

  # Put the multivariate ordered replacements back onto a FULL data frame.
  out[cols.ordered] <- multivarEst[cols.ordered]

  # preserve original column order & types for continuous vars
  out <- out[, colnames(data), drop = FALSE]

  list(data = out, thresholds = thresholds, coverage = coverage)
}
