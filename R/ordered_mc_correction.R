modsemOrderedMCCorrection <- function(model.syntax,
                                      data,
                                      method = "lms",
                                      ordered = NULL,
                                      calc.se = TRUE,
                                      group = NULL,
                                      verbose = interactive(),
                                      optimize = TRUE,
                                      start = NULL,
                                      mc.reps = NULL,
                                      ordered.mc.reps = mc.reps,
                                      ordered.min.iter = 20L,
                                      ordered.max.iter = 100L,
                                      ordered.tol = 1e-4,
                                      ordered.rng.seed = NULL,
                                      ordered.fixed.seed = FALSE,
                                      ordered.polyak.juditsky = TRUE,
                                      ordered.pj.extrapolate = TRUE,
                                      ordered.fn.args = list(),
                                      ordered.se = c("penalty", "naive", "delta"),
                                      ordered.se.penalty = 0.5,
                                      ordered.delta.reps = ordered.mc.reps,
                                      ordered.delta.epsilon = 1e-2,
                                      ordered.standardize = TRUE,
                                      standardize = NULL,
                                      standardize.out = NULL,
                                      cluster = NULL,
                                      sampling.weights = NULL,
                                      sampling.weights.normalization = NULL,
                                      mean.observed = NULL, # capture
                                      ...) {
  method <- tolower(method)
  ordered.se <- match.arg(ordered.se)
  mod_stopif(
    !method %in% c("lms", "qml"),
    "MC ordered correction is only available for LMS and QML."
  )

  if (is.null(verbose))
    verbose <- interactive()
  if (is.null(calc.se))
    calc.se <- TRUE
  if (is.null(ordered.mc.reps))
    ordered.mc.reps <- max(NROW(data), 2000L)
  if (is.null(ordered.delta.reps))
    ordered.delta.reps <- max(ordered.mc.reps, 10000L)

  ordered.mc.reps <- as.integer(ordered.mc.reps)
  ordered.delta.reps <- as.integer(ordered.delta.reps)

  mod_stopif(ordered.mc.reps < 2L, "`ordered.mc.reps` must be at least 2.")
  mod_stopif(ordered.min.iter < 1L, "`ordered.min.iter` must be at least 1.")
  mod_stopif(ordered.max.iter < ordered.min.iter,
         "`ordered.max.iter` must be larger than `ordered.min.iter`.")
  mod_stopif(!is.list(ordered.fn.args), "`ordered.fn.args` must be a list.")

  if (ordered.fixed.seed && is.null(ordered.rng.seed)) {
    ordered.rng.seed <- floor(stats::runif(1L, min = 1, max = 9999999))
    if (verbose)
      mod_msg_note(sprintf("Using fixed MC ordered seed %i.", ordered.rng.seed))
  }

  parTable.in <- modsemify(model.syntax)
  ovs <- getOVs(parTable.in)

  data <- as.data.frame(data)
  cols.ordered <- mcOrderedCols(data = data, ovs = ovs, ordered = ordered)
  mod_stopif(!length(cols.ordered),
    "No ordered variables were found for MC ordered correction."
  )

  if (verbose) {
    mod_msg_note(
      sprintf(
        "Estimating Monte-Carlo corrected estimates (method=%s, R=%d).\n",
        toupper(method), ordered.mc.reps
      ),
      "This is an experimental feature!",
      "See `help(\"modsem_da\")` for more information.",
      "Consider using the MC-OrdPLSc estimator in the `plssem` package."
    )
  }

  observed <- mcPrepareOrderedData(
    data = data,
    ovs = ovs,
    ordered = cols.ordered,
    group = group,
    standardize = ordered.standardize
  )

  PROBS <- mcOrderedProbs(
    data = observed$raw, ordered = cols.ordered,
    group = group
  )

  fit0 <- mcFitModsemOrdered(
    model.syntax = model.syntax,
    data = observed$data,
    method = method,
    group = group,
    calc.se = isTRUE(calc.se),
    verbose = verbose,
    optimize = optimize,
    start = start,
    cluster = cluster,
    sampling.weights = sampling.weights,
    sampling.weights.normalization = sampling.weights.normalization,
    standardize = FALSE,
    standardize.out = FALSE,
    mean.observed = FALSE,
    ...
  )

  std.info <- mcStandardizedStateInfo(fit0)
  theta0 <- std.info$state
  p <- theta0
  if (!is.null(start) && length(start) == length(theta0))
    p[] <- start

  f <- function(p.i, reps = ordered.mc.reps, seed = ordered.rng.seed,
                verbose.fit = FALSE) {

    fit.sim <- mcFitSimulatedOrdinalStandardized(
      fit0 = fit0,
      state = p.i,
      std.info = std.info,
      model.syntax = model.syntax,
      method = method,
      group = group,
      ovs = ovs,
      ordered = cols.ordered,
      PROBS = PROBS,
      mc.reps = reps,
      seed = seed,
      standardize = ordered.standardize,
      sampling.weights = sampling.weights,
      sampling.weights.normalization = sampling.weights.normalization,
      mean.observed = FALSE,
      verbose = verbose.fit,
      ...
    )

    fit.sim$state - theta0
  }

  if (isTRUE(ordered.polyak.juditsky) && !isTRUE(ordered.pj.extrapolate)) {
    if (verbose) mod_msg_note("Warming up MC ordered solver.")

    warmup <- mcRobbinsMonro(
      p = p,
      f = f,
      tol = 10 * ordered.tol,
      min.iter = min(ordered.min.iter, 5L),
      max.iter = min(max(ordered.min.iter, 5L), 20L),
      verbose = verbose,
      polyak.juditsky = FALSE,
      fn.args = ordered.fn.args
    )

    p <- as.vector(warmup$root)
    names(p) <- names(theta0)
  }

  mcfit <- mcRobbinsMonro(
    p = p,
    f = f,
    tol = ordered.tol,
    min.iter = ordered.min.iter,
    max.iter = ordered.max.iter,
    verbose = verbose,
    polyak.juditsky = ordered.polyak.juditsky,
    pj.extrapolate = ordered.pj.extrapolate,
    fn.args = ordered.fn.args
  )

  theta.mc <- as.vector(mcfit$root)
  names(theta.mc) <- names(theta0)

  vcov.free <- NULL
  type.se <- if (isTRUE(calc.se)) fit0$type.se else "none"

  naive.se.label   <- "naive"
  penalty.se.label <- "penalized"
  delta.se.label   <- "mc-delta"

  if (isTRUE(calc.se)) {
    vcov.free <- std.info$vcov.free

    if (!is.null(vcov.free) &&
        all(dim(vcov.free) == c(length(theta.mc), length(theta.mc)))) {

      if (identical(ordered.se, "delta")) {
        if (verbose)
          mod_msg_note("Calculating MC ordered delta-method standard errors.")

        vcov.free <- tryCatch({
          delta.seed <- ordered.rng.seed

          if (is.null(delta.seed)) {
            delta.seed <- floor(stats::runif(1L, min = 1, max = 9999999))

            if (verbose)
              mod_msg_note(sprintf("Using fixed MC ordered delta seed %i.", delta.seed))
          }

          vcov.delta <- mcAlignVcov(vcov.free, labels = names(theta.mc))

          H <- mcDeltaJacobian(
            p = theta.mc,
            f = f,
            reps = ordered.delta.reps,
            seed = delta.seed,
            epsilon = ordered.delta.epsilon,
            verbose = verbose
          )

          H.inv <- mcStableInverse(H)
          V <- H.inv %*% vcov.delta %*% t(H.inv)
          V <- mcStableVcov(V)

          dimnames(V) <- list(names(theta.mc), names(theta.mc))
          attr(V, "delta.jacobian") <- H
          attr(V, "delta.seed") <- delta.seed
          attr(V, "delta.inverse.method") <- attr(H.inv, "method")
          attr(V, "delta.condition") <- attr(H.inv, "condition")
          attr(V, "delta.psd.adjusted") <- attr(V, "psd.adjusted")

          type.se <<- delta.se.label

          V

        }, error = function(e) {
          mod_msg_warn(paste0(
            "Delta-method MC correction failed; using ordered MC penalty ",
            "standard errors instead. Message: ", conditionMessage(e)
          ))

          penalized <- mcOrderedPenaltyVcov(
            theta.mc = theta.mc,
            theta0 = theta0,
            vcov.free = std.info$vcov.free,
            lambda = ordered.se.penalty
          )

          type.se <<- penalty.se.label
          penalized
        })

      } else if (identical(ordered.se, "naive")) {
        if (verbose)
          mod_msg_note("Using naive standard errors.")

        vcov.free <- mcOrderedNaiveScalingVcov(
          theta.mc = theta.mc,
          theta0 = theta0,
          vcov.free = vcov.free
        )

        type.se <- naive.se.label

      } else if (identical(ordered.se, "penalty")) {
        if (verbose)
          mod_msg_note("Using penalized standard errors.")

        vcov.free <- mcOrderedPenaltyVcov(
          theta.mc = theta.mc,
          theta0 = theta0,
          vcov.free = vcov.free,
          lambda = ordered.se.penalty
        )

        type.se <- penalty.se.label
      }
    }
  }

  fit.out <- mcOverwriteFitWithStdState(
    fit = fit0,
    state = theta.mc,
    vcov.free = vcov.free,
    std.info = std.info,
    calc.se = calc.se,
    type.se = type.se
  )

  thresholds <- mcOrderedThresholdsDelta(
    data = observed$raw,
    ordered = cols.ordered,
    group = group,
    calc.se = calc.se
  )

  fit.out <- mcAppendThresholds(
    fit = fit.out,
    thresholds = thresholds
  )

  fit.out$expected.matrices <- mcOrderedExpectedMatrices(
    fit = fit.out,
    std.info = std.info,
    state = theta.mc,
    observed = observed$data,
    ovs = ovs,
    ordered = cols.ordered,
    PROBS = PROBS,
    group = group,
    mc.reps = ordered.mc.reps,
    seed = ordered.rng.seed,
    standardize = ordered.standardize
  )

  fit.out$iterations <- mcfit$iter
  fit.out$convergence <- isTRUE(mcfit$converged)
  fit.out$convergence.msg <- getConvergenceMessage(
    fit.out$convergence, fit.out$iterations
  )

  fit.out$args$ordered <- cols.ordered
  fit.out$args$ordered.mc.reps <- ordered.mc.reps
  fit.out$args$ordered.min.iter <- ordered.min.iter
  fit.out$args$ordered.max.iter <- ordered.max.iter
  fit.out$args$ordered.tol <- ordered.tol
  fit.out$args$ordered.rng.seed <- ordered.rng.seed
  fit.out$args$ordered.polyak.juditsky <- ordered.polyak.juditsky
  fit.out$args$ordered.pj.extrapolate <- ordered.pj.extrapolate
  fit.out$args$ordered.fn.args <- ordered.fn.args
  fit.out$args$ordered.se <- ordered.se
  fit.out$args$ordered.se.penalty <- ordered.se.penalty
  fit.out$args$ordered.delta.reps <- ordered.delta.reps
  fit.out$args$ordered.delta.epsilon <- ordered.delta.epsilon
  fit.out$args$optimize <- optimize
  fit.out$args$start <- start
  fit.out$ordered.mc <- mcfit

  fit.out$args$standardize <- TRUE
  fit.out$args$standardize.out <- TRUE

  fit.out
}


mcOrderedCols <- function(data, ovs, ordered = NULL) {
  ordered <- ordered %||% character(0L)
  cols <- ovs[ovs %in% colnames(data)]
  cols[cols %in% ordered | vapply(data[cols], is.ordered, logical(1L))]
}


mcPrepareOrderedData <- function(data, ovs, ordered, group = NULL,
                                 standardize = TRUE) {
  out <- as.data.frame(data)
  raw <- as.data.frame(data)

  if (is.null(group)) {
    out[ordered] <- lapply(out[ordered], mcOrderedInteger)
    if (standardize)
      out[ordered] <- lapply(out[ordered], std1)
  } else {
    group.vals <- out[[group]]
    for (col in ordered) {
      x.out <- numeric(NROW(out))
      for (g in unique(group.vals)) {
        idx <- which(group.vals == g)
        x <- mcOrderedInteger(out[[col]][idx])
        if (standardize) x <- std1(x)
        x.out[idx] <- x
      }
      out[[col]] <- x.out
    }
  }

  list(data = out, raw = raw[unique(c(ovs, group))])
}


mcOrderedInteger <- function(x) {
  if (is.ordered(x) || is.factor(x)) as.integer(as.ordered(x))
  else as.integer(as.ordered(x))
}


mcOrderedProbs <- function(data, ordered, group = NULL) {
  if (is.null(group)) {
    out <- list("1" = mcOrderedProbsGroup(data, ordered))
  } else {
    vals <- unique(data[[group]])
    out <- stats::setNames(vector("list", length(vals)), as.character(vals))
    for (g in vals)
      out[[as.character(g)]] <- mcOrderedProbsGroup(data[data[[group]] == g, ],
                                                    ordered)
  }

  out
}


mcOrderedProbsGroup <- function(data, ordered) {
  stats::setNames(lapply(ordered, function(col) {
    x <- as.ordered(data[[col]])
    tab <- as.numeric(table(x))
    p <- tab / sum(tab)
    p[p > 0]
  }), ordered)
}


mcFitModsemOrdered <- function(model.syntax, data, method, group, calc.se,
                               verbose, optimize, start, cluster = NULL,
                               sampling.weights = NULL,
                               sampling.weights.normalization = NULL,
                               standardize = FALSE,
                               standardize.out = FALSE, ...) {
  modsem_da(
    model.syntax = model.syntax,
    data = data,
    method = method,
    group = group,
    cluster = cluster,
    sampling.weights = sampling.weights,
    sampling.weights.normalization = sampling.weights.normalization,
    standardize = standardize,
    standardize.out = standardize.out,
    ordered = NULL,
    calc.se = calc.se,
    verbose = verbose,
    optimize = optimize,
    start = start,
    ...
  )
}


mcFitSimulatedOrdinalStandardized <- function(fit0,
                                              state,
                                              std.info,
                                              model.syntax,
                                              method,
                                              group,
                                              ovs,
                                              ordered,
                                              PROBS,
                                              mc.reps,
                                              seed,
                                              standardize = TRUE,
                                              sampling.weights = NULL,
                                              sampling.weights.normalization = NULL,
                                              verbose = FALSE,
                                              ...) {
  parTable <- mcStdParTableFromState(std.info = std.info, state = state)
  sim <- simulateDataParTableStandardized(
    parTable = parTable,
    N = mc.reps,
    colsOVs = ovs,
    seed = seed
  )

  data.sim <- mcOrdinalizeSimulatedData(
    sim = sim,
    ovs = ovs,
    ordered = ordered,
    PROBS = PROBS,
    group = group,
    standardize = standardize
  )

  if (!is.null(sampling.weights))
    data.sim[[sampling.weights]] <- 1

  fit.sim <- mcFitModsemOrdered(
    model.syntax = model.syntax,
    data = data.sim,
    method = method,
    group = group,
    calc.se = FALSE,
    verbose = verbose,
    optimize = FALSE,
    start = fit0$theta,
    cluster = NULL,
    sampling.weights = sampling.weights,
    sampling.weights.normalization = sampling.weights.normalization,
    ...
  )

  list(
    fit = fit.sim,
    state = mcExtractStandardizedState(fit.sim, std.info = std.info)
  )
}


mcOrdinalizeSimulatedData <- function(sim, ovs, ordered, PROBS, group = NULL,
                                      standardize = TRUE) {
  groups <- seq_along(sim$OV)
  out <- NULL

  for (g in groups) {
    data.g <- as.data.frame(sim$OV[[g]])[ovs]
    prob.g <- if (length(PROBS) == 1L) PROBS[[1L]] else PROBS[[g]]

    for (col in ordered) {
      data.g[[col]] <- mcOrdinalize(data.g[[col]], probs = prob.g[[col]])
      if (standardize)
        data.g[[col]] <- std1(data.g[[col]])
    }

    if (!is.null(group)) {
      group.name <- names(PROBS)[[g]]
      data.g[[group]] <- utils::type.convert(group.name, as.is = TRUE)
    }

    out <- rbind(out, data.g)
  }

  out
}


mcOrdinalize <- function(x, probs) {
  probs <- probs / sum(probs)
  if (length(probs) < 2L)
    return(rep(1L, length(x)))

  breaks <- stats::quantile(x, probs = cumsum(probs[-length(probs)]),
                            names = FALSE, type = 7, na.rm = TRUE)
  as.integer(findInterval(x, vec = breaks) + 1L)
}


mcParTableWithTheta <- function(fit, theta, method) {
  coefs <- getLavCoefs(model = fit$start.model, theta = theta, method = method)$all
  parTable <- as.data.frame(fit$parTable)
  key <- mcParTableKeys(parTable)
  idx <- match(key, names(coefs))
  hit <- !is.na(idx)
  parTable$est[hit] <- coefs[idx[hit]]
  parTable[parTable$op != "|", , drop = FALSE]
}


mcParTableKeys <- function(parTable) {
  key <- paste0(parTable$lhs, parTable$op, parTable$rhs)
  if ("group" %in% colnames(parTable)) {
    g <- parTable$group
    key <- ifelse(!is.na(g) & g > 1L, paste0(key, ".g", g), key)
  }

  key
}


mcRobbinsMonro <- function(p, f, tol, min.iter, max.iter, verbose = FALSE,
                           polyak.juditsky = FALSE, pj.extrapolate = TRUE,
                           fn.args = list(), k = 3L) {
  if (max.iter < min.iter) max.iter <- min.iter

  fn.a <- fn.args$fn.a
  if (is.null(fn.a))
    fn.a <- function(iter, a = 1, b = 1/2, c = 0, ...) a / (iter + c)^b
  fn.args$fn.a <- NULL

  history <- rbind(p, matrix(NA_real_, nrow = max.iter, ncol = length(p)))
  colnames(history) <- names(p)
  k.succ <- 0L
  pbar.last <- pbar <- p

  for (i in seq_len(max.iter)) {
    p.old <- p
    fp <- f(p)
    fp[!is.finite(fp)] <- 0
    a <- do.call(fn.a, c(list(iter = i), fn.args))
    p <- p - a * fp
    names(p) <- names(p.old)

    history[i + 1L, ] <- p
    change <- max(abs(p.old - p), na.rm = TRUE)

    if (isTRUE(polyak.juditsky)) {
      pbar.last <- pbar
      pbar <- mcPKAverage(history)
      names(pbar) <- names(p)
      change <- max(abs(pbar.last - pbar), na.rm = TRUE)
    }

    if (isTRUE(verbose)) {
      if (isTRUE(polyak.juditsky)) msg <- "\rMC iter %d; \u0394E(\u03B8) = %.4g"
      else msg <- "\rMC iter %d; \u0394\u03B8 = %.4g"

      printf(msg, i, change)
    }

    if (i >= min.iter && is.finite(change) && change < tol) {
      k.succ <- k.succ + 1L
      if (k.succ >= k) break
    } else {
      k.succ <- 0L
    }
  }

  history <- history[seq_len(i + 1L), , drop = FALSE]
  root <- if (isTRUE(polyak.juditsky)) pbar else p
  if (isTRUE(polyak.juditsky) && isTRUE(pj.extrapolate))
    root <- mcGetConvergencePoints(history)

  names(root) <- colnames(history)

  if (isTRUE(verbose))
    printf("\n")

  list(root = root, iter = i, converged = i < max.iter, history = history,
       polyak.juditsky = polyak.juditsky,
       pj.extrapolate = pj.extrapolate)
}


mcPKAverage <- function(history) {
  colMeans(history, na.rm = TRUE)
}


mcGetConvergencePoint <- function(y, t = seq_along(y)) {
  c.aitken <- mcAitkenAccelerate(y)
  span <- diff(range(y))

  if (length(y) < 4L || !is.finite(span) ||
      span < .Machine$double.eps^0.5)
    return(c.aitken)

  tryCatch({
    fit <- stats::nls(
      y ~ c + a * exp(-k * t),
      start = list(
        c = mean(utils::tail(y, 3L)),
        a = y[[1L]] - mean(utils::tail(y, 3L)),
        k = 0.1
      ),
      algorithm = "port",
      lower = c(c = -Inf, a = -Inf, k = 0)
    )

    c.nls <- stats::coef(fit)[["c"]]
    k.fit <- stats::coef(fit)[["k"]]

    bad.range <- c.nls < min(y) - span || c.nls > max(y) + span
    bad.k <- k.fit < sqrt(.Machine$double.eps)
    bad.agree <- abs(c.nls - c.aitken) > 0.5 * span

    if (bad.range || bad.k || bad.agree) c.aitken else c.nls
  }, error = function(e) c.aitken)
}


mcAitkenAccelerate <- function(y) {
  n <- length(y)
  if (n < 3L) return(mean(y, na.rm = TRUE))

  ests <- vapply(seq_len(n - 2L), FUN.VALUE = numeric(1L), FUN = function(i) {
    p0 <- y[[i]]
    p1 <- y[[i + 1L]]
    p2 <- y[[i + 2L]]

    denom <- p2 - 2 * p1 + p0
    if (abs(denom) < .Machine$double.eps^0.5)
      return(NA_real_)

    p0 - (p1 - p0)^2 / denom
  })

  valid <- ests[is.finite(ests)]
  if (!length(valid)) return(mean(y, na.rm = TRUE))
  stats::median(valid)
}


mcGetConvergencePoints <- function(history) {
  history <- history[stats::complete.cases(history), , drop = FALSE]
  out <- apply(history, MARGIN = 2L, FUN = mcGetConvergencePoint)
  names(out) <- colnames(history)
  out
}


mcDeltaJacobian <- function(p, f, reps, seed, epsilon = 1e-3, verbose = interactive()) {
  k <- length(p)
  J <- matrix(NA_real_, nrow = k, ncol = k,
              dimnames = list(names(p), names(p)))

  for (j in seq_len(k)) {
    if (verbose) printf("\rCalculating Jacobian %d/%d...", j, k)

    step <- epsilon * max(1, abs(p[[j]]))
    p.plus <- p
    p.minus <- p
    p.plus[[j]] <- p.plus[[j]] + step
    p.minus[[j]] <- p.minus[[j]] - step

    f.plus <- f(p.plus, reps = reps, seed = seed, verbose.fit = FALSE)
    f.minus <- f(p.minus, reps = reps, seed = seed, verbose.fit = FALSE)

    if (!all(is.finite(f.plus)) || !all(is.finite(f.minus)))
      stop("Non-finite finite-difference function value.")

    J[, j] <- (f.plus - f.minus) / (2 * step)
  }

  if (verbose)
    printf("\n")

  if (!all(is.finite(J)))
    stop("Non-finite finite-difference Jacobian.")

  J
}


mcAlignVcov <- function(vcov, labels) {
  if (is.null(vcov))
    stop("Missing variance-covariance matrix.")

  dn <- dimnames(vcov)
  if (is.null(dn) || is.null(dn[[1L]]) || is.null(dn[[2L]]))
    stop("Variance-covariance matrix must have dimnames.")

  missing <- setdiff(labels, rownames(vcov))
  if (length(missing))
    stop("Variance-covariance matrix is missing parameters: ",
         paste(missing, collapse = ", "))

  V <- vcov[labels, labels, drop = FALSE]
  V <- 0.5 * (V + t(V))
  dimnames(V) <- list(labels, labels)
  V
}


mcStableInverse <- function(X, rcond.tol = 1e-8) {
  if (NROW(X) != NCOL(X))
    stop("Jacobian must be square.")
  if (!all(is.finite(X)))
    stop("Jacobian contains non-finite values.")

  sv <- svd(X)
  d <- sv$d
  d.max <- max(d)
  d.min <- min(d)
  if (!is.finite(d.max) || d.max <= 0)
    stop("Jacobian is numerically zero.")

  rcond <- d.min / d.max
  if (is.finite(rcond) && rcond > rcond.tol) {
    out <- solve(X)
    method <- "solve"
    condition <- 1 / rcond
  } else {
    d.reg <- pmax(d, d.max * rcond.tol)
    out <- sv$v %*% (diag(1 / d.reg, nrow = length(d.reg)) %*% t(sv$u))
    method <- "svd-ridge"
    condition <- 1 / max(rcond, rcond.tol, na.rm = TRUE)
  }

  dimnames(out) <- rev(dimnames(X))
  attr(out, "method") <- method
  attr(out, "condition") <- condition
  out
}


mcStableVcov <- function(V, eigen.tol = .Machine$double.eps^0.5) {
  nms <- dimnames(V)
  V <- 0.5 * (V + t(V))
  if (!all(is.finite(V)))
    stop("Delta-method variance-covariance matrix contains non-finite values.")

  eig <- eigen(V, symmetric = TRUE)
  scale <- max(abs(eig$values), 1)
  lower <- -eigen.tol * scale
  adjusted <- any(eig$values < 0)
  materially.adjusted <- any(eig$values < lower)
  values <- pmax(eig$values, 0)

  if (adjusted) {
    V <- eig$vectors %*% (diag(values, nrow = length(values)) %*% t(eig$vectors))
    V <- 0.5 * (V + t(V))
  }

  dimnames(V) <- nms
  attr(V, "psd.adjusted") <- adjusted
  attr(V, "psd.materially.adjusted") <- materially.adjusted
  V
}


mcOrderedNaiveScalingVcov <- function(theta.mc, theta0, vcov.free,
                                      eps = .Machine$double.eps^0.5) {
  if (is.null(vcov.free)) return(NULL)

  ids <- intersect(names(theta.mc), intersect(names(theta0), rownames(vcov.free)))
  if (!length(ids)) return(vcov.free)

  scale <- rep(1, length(ids))
  names(scale) <- ids

  denom <- theta0[ids]
  numer <- theta.mc[ids]
  ok <- is.finite(denom) & is.finite(numer) & abs(denom) > eps
  scale[ok] <- numer[ok] / denom[ok]

  D <- diag(scale, nrow = length(scale))
  dimnames(D) <- list(ids, ids)
  V <- vcov.free[ids, ids, drop = FALSE]
  V <- D %*% V %*% D
  V <- 0.5 * (V + t(V))
  dimnames(V) <- list(ids, ids)
  expandVCOV(V, labels = rownames(vcov.free))
}


mcOrderedPenaltyVcov <- function(theta.mc, theta0, vcov.free, lambda = 1) {
  if (is.null(vcov.free)) return(NULL)

  ids <- intersect(names(theta.mc), intersect(names(theta0), rownames(vcov.free)))
  if (!length(ids)) return(vcov.free)

  d <- theta.mc[ids] - theta0[ids]
  V <- vcov.free[ids, ids, drop = FALSE]
  V <- V + lambda * diag(d^2, nrow = length(d))
  V <- 0.5 * (V + t(V))
  dimnames(V) <- list(ids, ids)
  expandVCOV(V, labels = rownames(vcov.free))
}


simulateDataParTableStandardized <- function(parTable, N, colsOVs = NULL,
                                             colsLVs = NULL, seed = NULL) {
  simulateDataParTable(
    parTable = parTable,
    N = N,
    colsOVs = colsOVs,
    colsLVs = colsLVs,
    seed = seed
  )
}


mcOrderedExpectedMatrices <- function(fit,
                                      std.info,
                                      state,
                                      observed,
                                      ovs,
                                      ordered,
                                      PROBS,
                                      group = NULL,
                                      mc.reps,
                                      seed = NULL,
                                      standardize = TRUE) {

  parTable <- mcStdParTableFromState(
    std.info = std.info, state = state
  )

  sim <- simulateDataParTableStandardized(
    parTable = parTable,
    N = mc.reps,
    colsOVs = ovs,
    seed = seed
  )

  sim.ord <- mcOrdinalizeSimulatedData(
    sim = sim,
    ovs = ovs,
    ordered = ordered,
    PROBS = PROBS,
    group = group,
    standardize = standardize
  )

  observed.blocks <- mcSplitOrderedBlocks(
    data = observed,
    ovs = ovs,
    group = group
  )

  expected.blocks <- mcSplitOrderedBlocks(
    data = sim.ord,
    ovs = ovs,
    group = group
  )

  out <- vector("list", length(observed.blocks))
  for (i in seq_along(observed.blocks)) {
    obs.g <- observed.blocks[[i]]
    exp.g <- expected.blocks[[min(i, length(expected.blocks))]]

    sigma.obs <- stats::cov(obs.g, use = "pairwise.complete.obs")
    sigma.exp <- stats::cov(exp.g, use = "pairwise.complete.obs")

    mu.obs <- matrix(
      colMeans(obs.g, na.rm = TRUE), ncol = 1,
      dimnames = list(colnames(obs.g), "~1")
    )

    mu.exp <- matrix(
      colMeans(exp.g, na.rm = TRUE), ncol = 1,
      dimnames = list(colnames(exp.g), "~1")
    )

    out[[i]] <- list(
      sigma.all = sigma.exp,
      sigma.lv = NULL,
      sigma.ov = sigma.exp,
      mu.all = mu.exp,
      mu.lv = NULL,
      mu.ov = mu.exp,
      r2.all = NULL,
      r2.lv = NULL,
      r2.ov = NULL,
      res.all = NULL,
      res.lv = NULL,
      res.ov = NULL,
      lambda = NULL,
      gammaXi = NULL,
      gammaEta = NULL,
      psi = NULL,
      phi = NULL,
      theta = NULL,
      sigma.ord.observed = sigma.obs,
      sigma.ord.expected = sigma.exp,
      mu.ord.observed = mu.obs,
      mu.ord.expected = mu.exp,
      sim.ov.ord = exp.g
    )
  }

  out
}


mcSplitOrderedBlocks <- function(data, ovs, group = NULL) {
  if (is.null(group)) {
    return(list(as.data.frame(data)[, ovs, drop = FALSE]))
  }

  vals <- unique(data[[group]])
  stats::setNames(
    lapply(vals, function(g) {
      as.data.frame(data[data[[group]] == g, ovs, drop = FALSE])
    }),
    as.character(vals)
  )
}


mcStandardizedStateInfo <- function(fit) {
  solution <- standardizedSolutionCOEFS(
    fit,
    monte.carlo = FALSE
  )
  parTable <- as.data.frame(solution$parTable)
  parTable$.state.id <- getParTableLabels(parTable, labelCol = "label")
  parTable$.state.redundant <- FALSE

  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  has.residual.cov <- mcHasResidualCovariances(parTable)
  keep <- op %in% c("=~", "<~", "~", "~~")

  if (!has.residual.cov) {
    redundant <- (
      keep & (
        (op == "~~" & lhs == rhs) |
        (op == "~~" & (grepl(":", lhs) | grepl(":", rhs)))
      )
    )

    parTable$.state.redundant[redundant] <- TRUE
    keep <- keep & !redundant
  }

  ids <- unique(parTable$.state.id[keep])
  est <- stats::setNames(parTable$est, parTable$.state.id)
  est <- est[ids]

  vcov.free <- solution$vcov
  if (!is.null(vcov.free) && length(vcov.free)) {
    vcov.free <- expandVCOV(vcov.free, labels = ids)
  } else {
    vcov.free <- NULL
  }

  list(
    template = parTable,
    state = est,
    vcov.free = vcov.free,
    reduced = !has.residual.cov,
    has.residual.cov = has.residual.cov
  )
}


mcHasResidualCovariances <- function(parTable) {
  is.cov <- parTable$op == "~~" & parTable$lhs != parTable$rhs
  if (!any(is.cov)) return(FALSE)

  lVs <- getLVs(parTable)
  indsLVs <- getIndsLVs(parTable, lVs = lVs)
  inds <- unique(unlist(indsLVs))
  endogenous <- unique(parTable[parTable$op == "~", "lhs"])

  covs <- parTable[is.cov, , drop = FALSE]
  any(
    (covs$lhs %in% inds & covs$rhs %in% inds) |
      covs$lhs %in% endogenous |
      covs$rhs %in% endogenous
  )
}


mcExtractStandardizedState <- function(fit, std.info) {
  solution <- standardizedSolutionCOEFS(fit, monte.carlo = FALSE)
  parTable <- as.data.frame(solution$parTable)

  parTable$.state.id <- getParTableLabels(parTable, labelCol = "label")
  est <- stats::setNames(parTable$est, parTable$.state.id)
  out <- est[names(std.info$state)]
  names(out) <- names(std.info$state)

  out
}


mcStdParTableFromState <- function(std.info, state) {
  parTable <- std.info$template
  idx <- match(parTable$.state.id, names(state))
  hit <- !is.na(idx)
  parTable$est[hit] <- state[idx[hit]]

  if (isTRUE(std.info$reduced))
    parTable <- mcStdRefreshRedundant(parTable)

  parTable
}


mcStdRefreshRedundant <- function(parTable) {
  # parTable <- var_interactions(parTable)

  groups <- getGroupsParTable(addMissingGroups(parTable))
  out <- parTable

  for (g in groups) {
    mask.g <- out$group == g
    parTable.g <- out[mask.g, , drop = FALSE]
    if (!NROW(parTable.g)) next

    lVs.g <- getLVs(parTable.g)
    intTerms.g <- getIntTerms(parTable.g)
    etas.g <- getSortedEtas(parTable.g, isLV = FALSE)
    xis.g <- getXis(parTable.g, etas = etas.g, isLV = FALSE)
    indsLVs.g <- getIndsLVs(parTable.g, lVs = lVs.g)
    allInds.g <- unique(unlist(indsLVs.g))

    exogenous <- unique(c(xis.g, intTerms.g))
    for (x in exogenous) {
      sel <- parTable.g$lhs == x & parTable.g$rhs == x & parTable.g$op == "~~"
      if (any(sel)) parTable.g[sel, "est"] <- 1
    }

    for (eta in etas.g) {
      sel <- parTable.g$lhs == eta & parTable.g$rhs == eta & parTable.g$op == "~~"
      if (!any(sel)) next
      total <- calcVarParTable(eta, parTable = parTable.g)
      parTable.g[sel, "est"] <- pmax(parTable.g[sel, "est"] + (1 - total), 1e-8)
    }

    for (ind in allInds.g) {
      sel <- parTable.g$lhs == ind & parTable.g$rhs == ind & parTable.g$op == "~~"
      if (!any(sel)) next
      total <- calcVarParTable(ind, parTable = parTable.g, measurement.model = TRUE)
      parTable.g[sel, "est"] <- pmax(parTable.g[sel, "est"] + (1 - total), 1e-8)
    }

    out[mask.g, ] <- parTable.g
  }

  out
}


mcOverwriteFitWithStdState <- function(fit, state, vcov.free, std.info, calc.se,
                                       type.se = fit$type.se) {
  parTable <- mcStdParTableFromState(std.info = std.info, state = state)

  se <- stats::setNames(rep(NA_real_, length(state)), names(state))
  if (isTRUE(calc.se) && !is.null(vcov.free)) {
    ids.v <- intersect(names(state), rownames(vcov.free))
    if (length(ids.v)) {
      se[ids.v] <- sqrt(pmax(diag(vcov.free[ids.v, ids.v, drop = FALSE]), 0))
    }
  }

  parTable$std.error <- NA_real_
  idx.se <- match(parTable$.state.id, names(se))
  hit.se <- !is.na(idx.se)
  parTable$std.error[hit.se] <- se[idx.se[hit.se]]

  parTable <- addZStatsParTable(parTable)
  parTable$.state.id <- NULL
  parTable$.state.redundant <- NULL

  fit$parTable <- modsemParTable(sortParTableDA(parTable, model = fit$model))
  fit$coefs.all <- stats::setNames(parTable$est, getParTableLabels(parTable, labelCol = "label", replace.dup = TRUE))
  fit$coefs.free <- state
  fit$vcov.all <- if (!is.null(vcov.free))
    expandVCOV(vcov.free, labels = names(fit$coefs.all))
  else NULL
  fit$vcov.free <- vcov.free
  fit$type.se <- if (isTRUE(calc.se)) type.se else "none"
  fit$type.estimates <- "standardized"
  fit
}


mcReplaceFreeVcovAll <- function(fit, vcov.free) {
  if (is.null(fit$vcov.all) || is.null(dimnames(fit$vcov.all)))
    return(fit)

  lav <- mcFreeLavLabels(fit = fit, theta = fit$theta)
  keep <- lav %in% rownames(fit$vcov.all)
  if (!any(keep)) return(fit)

  V <- fit$vcov.all
  V[lav[keep], lav[keep]] <- vcov.free[keep, keep, drop = FALSE]
  fit$vcov.all <- V
  fit
}


mcSEFromVcov <- function(fit, theta, vcov.free, zero.tol = 1e-10) {
  if (is.null(vcov.free))
    return(stats::setNames(rep(NA_real_, length(fit$coefs.all)),
                           names(fit$coefs.all)))

  lav <- mcFreeLavLabels(fit = fit, theta = theta)

  se <- sqrt(pmax(diag(vcov.free), 0))
  se[se < zero.tol] <- NA_real_

  stats::setNames(se, lav)
}


mcFreeLavLabels <- function(fit, theta) {
  lav <- fit$start.model$params$lavLabels
  if (!is.null(names(lav))) {
    out <- lav[names(theta)]
  } else {
    out <- lav[seq_along(theta)]
  }
  unname(out)
}


mcOrderedThresholdsDelta <- function(data, ordered, group = NULL,
                                     calc.se = TRUE) {
  if (is.null(group)) {
    groups <- list("1" = data)
  } else {
    vals <- unique(data[[group]])
    groups <- stats::setNames(lapply(vals, function(g) data[data[[group]] == g, ]),
                              as.character(vals))
  }

  out <- NULL
  V.out <- NULL
  for (g in seq_along(groups)) {
    data.g <- groups[[g]]
    for (col in ordered) {
      x <- as.ordered(data.g[[col]])
      tab <- as.numeric(table(x))
      K <- length(tab)
      n <- sum(tab)
      if (K < 2L) next

      p <- tab / n
      C <- cumsum(p)[-K]
      eps <- .Machine$double.eps^0.5
      C <- pmin(pmax(C, eps), 1 - eps)
      tau <- stats::qnorm(C)

      labs <- paste0(col, "|t", seq_len(K - 1L))
      se.tau <- rep(NA_real_, length(tau))

      out <- rbind(out, data.frame(
        lhs = col,
        op = "|",
        rhs = paste0("t", seq_len(K - 1L)),
        label = "",
        group = g,
        est = tau,
        std.error = se.tau,
        stringsAsFactors = FALSE,
        row.names = NULL
      ))
    }
  }

  if (!is.null(out)) {
    out <- addZStatsParTable(out)
    attr(out, "vcov.thresholds") <- V.out
  }
  out
}


mcAppendThresholds <- function(fit, thresholds) {
  if (is.null(thresholds) || !NROW(thresholds))
    return(fit)

  V.tau <- attr(thresholds, "vcov.thresholds")
  parTable <- as.data.frame(fit$parTable)
  missing.cols <- setdiff(colnames(parTable), colnames(thresholds))
  thresholds[missing.cols] <- NA
  thresholds <- thresholds[colnames(parTable)]

  fit$parTable <- modsemParTable(sortParTableDA(rbind(parTable, thresholds),
                                                model = fit$model))
  coef.t <- stats::setNames(thresholds$est,
                            paste0(thresholds$lhs, thresholds$op,
                                   thresholds$rhs))
  fit$coefs.all <- c(fit$coefs.all, coef.t)

  if (!is.null(V.tau))
    fit$vcov.all <- diagPartitionedMat(fit$vcov.all, V.tau)

  fit
}


`%||%` <- function(x, y) if (is.null(x)) y else x
