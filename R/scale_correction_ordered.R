modsemOrderedScaleCorrection <- function(model.syntax,
                                         data,
                                         method = "lms",
                                         ordered = NULL,
                                         calc.se = TRUE,
                                         iter = 50L,
                                         warmup = floor(iter / 2L),
                                         N = 1e5, # 500,000
                                         verbose = interactive(),
                                         optimize = TRUE,
                                         start = NULL,
                                         lambda = 1,
                                         ordered.tol = 1e-6,
                                         se = "simple",
                                         standardize.data = NULL, # override
                                         ...) {
  message("Bootstrapping continuous values for ordinal variables...\n",
          "This is an experimental feature!\n",
          "See `help(modsem_da)` for more information.")

  standardize <- \(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

  if (is.null(verbose))
    verbose <- TRUE # default

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]

  rescaleOrderedData <- function(data, sim.ov) {

    sim.cat  <- as.data.frame(sim.ov)
    data.cat <- as.data.frame(data)
    data.cat[cols.ordered] <- lapply(data.cat[cols.ordered], as.integer)
    data.cont <- data.cat
    data.cont[cols.ordered] <- NA_real_

    sim.cat <- sim.cat[colnames(data.cat)] # ensure correct order
    sim.std <- lapplyDf(sim.cat, FUN = standardize)

    thresholds <- stats::setNames(vector("list", length(cols.ordered)),
                                  nm = cols.ordered)
    univarEst <- stats::setNames(vector("list", length(cols.ordered)),
                                   nm = cols.ordered)
    for (col in cols.ordered) {
      scaling <- rescaleOrderedVariable(name = col, data = data.cat, sim.ov = sim.ov)

      univarEst[[col]]  <- scaling$values # univariate mean-estimates
      thresholds[[col]] <- scaling$thresholds
      sim.cat[[col]]    <- cut(sim.std[[col]], breaks = scaling$thresholds,
                               ordered_result = TRUE, labels = FALSE)
    }

    univarEst <- as.data.frame(univarEst)

    # right now this will work best if all variables are categorical
    # Optimally we would want to use information from numerical variables
    # as well...
    combostring <- \(...) paste(..., sep = "-")
    combos.obs  <- do.call(combostring, data.cat[cols.ordered])
    combos.sim  <- do.call(combostring, sim.cat[cols.ordered])
    combos.both <- intersect(combos.obs, combos.sim)

    multivarEst <- univarEst # replace univariate estimates with multi-variate ones
                             # where possible

    # sim.agg <- as.data.frame(sim.std)[cols.ordered] |>
    #   dplyr::mutate(.combo__ = combos.sim) |> # pick a weird name
    #   dplyr::filter(.combo__ %in% combos.both) |>
    #   dplyr::group_by_at(".combo__") |>
    #   dplyr::summarize_at(.vars = cols.ordered, .funs = mean)
    #
    # sim.agg.n <- as.data.frame(sim.std)[cols.ordered] |>
    #   dplyr::mutate(.combo__ = combos.sim) |> # pick a weird name
    #   dplyr::filter(.combo__ %in% combos.both) |>
    #   dplyr::group_by_at(".combo__") |>
    #   dplyr::summarize(n = dplyr::n())

    for (combo in combos.both) {
      obs.mask <- combos.obs == combo
      sim.mask <- combos.sim == combo

      # mu <- sim.agg[sim.agg$.combo__ == combo, cols.ordered, drop = FALSE]
      k <- sum(obs.mask)
      m <- sum(sim.mask)
      M <- sample(m, k, replace = TRUE)

      multivarEst[obs.mask, cols.ordered] <-
        sim.std[sim.mask, cols.ordered, drop = FALSE][M, , drop = FALSE]
    }

    list(data = multivarEst, thresholds = thresholds)
  }

  ERROR <- \(e) {warning2(e, immediate. = FALSE); NULL}

  data.x <- data
  data.y <- data
  data.y[cols.ordered] <- lapply(data.y[cols.ordered], FUN = as.integer)

  thresholds <- stats::setNames(lapply(cols.ordered,  FUN = \(x) 0),
                                nm = cols.ordered)

  stopif(iter <= warmup, "`ordered.boot` must be larger than `ordered.warmup`!")

  iter.keep <- max(1L, floor(iter - warmup))
  COEF.ALL  <- vector("list", length = iter.keep)
  COEF.FREE <- vector("list", length = iter.keep)
  VCOV.ALL  <- vector("list", length = iter.keep)
  VCOV.FREE <- vector("list", length = iter.keep)
  fits      <- vector("list", length = iter.keep)
  imputed   <- vector("list", length = iter.keep)

  THETA <- NULL

  data_i <- data
  data_i[cols.ordered] <- lapply(data_i[cols.ordered],
                                 FUN = \(x) standardize(as.integer(x)))

  PARS <- NULL
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
          model.syntax = model.syntax,
          data         = data_i,
          method       = method,
          start        = start,
          verbose      = verbose,
          optimize     = optimize,
          calc.se      = calc.se_i,
          standardize.data = TRUE,
          ...
        ), error = \(e) {cat("\n");print(e); NULL}
      )

      if (is.null(fit_i)) next

      # pars_i <- parameter_estimates(fit_i)

      # if (is.null(PARS) || mode == "warmup") {
      #   PARS <- pars_i
      # } else { # get more an more stable estimates of the true parameter values
      #   lambda <- (1 / j)
      #   PARS$est <- (1 - lambda) * PARS$est + lambda * pars_i$est
      # }

      # sim_i <- simulateDataParTable(
      #   parTable = PARS,
      #   N        = N
      # )

      sim_i <- simulateDataParTable(
        parTable = parameter_estimates(fit_i),
        N        = N
      )

      rescaled   <- rescaleOrderedData(data = data.x, sim.ov = sim_i$oV)
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
