# Categorical LMS -----------------------------------------------------------
#
# The categorical estimator deliberately has its own likelihood backend.  The
# ordinary LMS backend integrates only a conditioning subset and then uses a
# multivariate-normal density. The ordinal backend instead evaluates normal-ogive
# category probabilities conditional on the latent state.

lmsCatControl <- function(control = list(), nodes = NULL) {
  defaults <- list(
    order = nodes %||% 9L,
    sparse.level = 4L,
    qmc.points = 4096L,
    qmc.replicates = 4L,
    qmc.max.points = 65536L,
    qmc.loglik.tol = 0.05,
    qmc.adapt = TRUE,
    gradient = TRUE,
    gradient.analytical = TRUE,
    gradient.eps = 1e-6,
    gradient.check = TRUE,
    gradient.check.tol = 1e-3,
    algorithm = "NLMINB",
    ema.max.step = 0.5,
    ema.min.step = 1e-6,
    ema.max.extrapolation = 10,
    max.nodes = 50000L,
    adapt = TRUE,
    adaptive.frequency = 3L,
    rel.tol = 1e-6,
    seed = 19790417L,
    start.max.iter = 50L,
    progress.every = 10L,
    progress.seconds = 1
  )
  utils::modifyList(defaults, control)
}


lmsCatOrderedCols <- function(data, model.syntax, ordered = NULL) {
  ovs <- getOVs(modsemify(model.syntax))
  auto <- names(data)[vapply(data, is.ordered, logical(1L))]
  cols <- unique(c(ordered, auto))
  cols <- intersect(cols, ovs)
  mod_stopif(!length(cols),
             "`lms-cat` requires at least one binary or ordered indicator.")
  missing <- setdiff(cols, names(data))
  mod_stopif(length(missing), paste("Ordered variables not found:",
                                    paste(missing, collapse = ", ")))
  cols
}


lmsCatScoreData <- function(data, ordered) {
  out <- data
  for (v in ordered) {
    x <- data[[v]]
    if (is.factor(x)) x <- as.integer(x)
    else {
      lev <- sort(unique(x[!is.na(x)]))
      x <- match(x, lev)
    }
    z <- as.numeric(scale(x))
    z[is.na(x)] <- NA_real_
    out[[v]] <- z
  }
  out
}


lmsCatThresholdSpec <- function(data, ordered, group = NULL) {
  if (is.null(group)) groups <- rep(1L, NROW(data))
  else groups <- match(data[[group]], unique(data[[group]]))

  specs <- list()
  delta <- numeric()
  for (g in sort(unique(groups))) {
    for (v in ordered) {
      x <- data[[v]][groups == g]
      if (is.factor(x)) code <- as.integer(x)
      else code <- match(x, sort(unique(x[!is.na(x)])))
      K <- max(code, na.rm = TRUE)
      mod_stopif(!is.finite(K) || K < 2L,
                 sprintf("Ordered indicator `%s` has fewer than two categories in group %d.", v, g))
      probs <- vapply(seq_len(K - 1L), function(k) mean(code <= k, na.rm = TRUE), numeric(1L))
      probs <- pmin(pmax(probs, 1e-5), 1 - 1e-5)
      tau <- stats::qnorm(probs)
      gaps <- diff(tau)
      d <- c(tau[1L], log(expm1(pmax(gaps, 1e-6))))
      labs <- paste0(v, "|t", seq_len(K - 1L), if (g > 1L) paste0(".g", g) else "")
      names(d) <- paste0("thresholdDelta:", labs)
      from <- length(delta) + seq_along(d)
      delta <- c(delta, d)
      specs[[paste(g, v, sep = "::")]] <- list(
        group = g, variable = v, K = K, code = code, index = from, labels = labs
      )
    }
  }
  list(specs = specs, delta = delta, groups = groups)
}


lmsCatDeltaToThresholds <- function(delta, spec) {
  d <- delta[spec$index]
  if (length(d) == 1L) return(d)
  c(d[1L], d[1L] + cumsum(log1p(exp(-abs(d[-1L]))) + pmax(d[-1L], 0)))
}


lmsCatThresholdJacobian <- function(delta, threshold.info) {
  n <- length(delta)
  J <- matrix(0, n, n)
  for (spec in threshold.info$specs) {
    idx <- spec$index
    J[idx[1L], idx[1L]] <- 1
    if (length(idx) > 1L) {
      for (r in seq_along(idx)[-1L]) {
        J[idx[r], idx[1L]] <- 1
        J[idx[r], idx[2:r]] <- stats::plogis(delta[idx[2:r]])
      }
    }
  }
  J
}


# Attach ordinal parameters to the ordinary DA model-matrix layout. The
# rectangular matrices use NaN for structurally absent thresholds so that a
# shorter item does not accidentally acquire parameters belonging to an item
# with more response categories.
lmsCatAttachThresholdMatrices <- function(model, threshold.info, delta,
                                           vcov.delta = NULL) {
  for (g in seq_along(model$models)) {
    submodel <- model$models[[g]]
    indicators <- rownames(submodel$matrices$lambdaX)
    specs <- Filter(function(x) identical(x$group, g), threshold.info$specs)
    max.thresholds <- max(c(0L, vapply(specs, function(x) x$K - 1L, integer(1L))))
    columns <- paste0("t", seq_len(max.thresholds))

    thresholds <- matrix(NaN, length(indicators), max.thresholds,
                         dimnames = list(indicators, columns))
    threshold.delta <- thresholds
    threshold.labels <- matrix("", length(indicators), max.thresholds,
                               dimnames = dimnames(thresholds))
    threshold.se <- thresholds

    for (spec in specs) {
      count <- spec$K - 1L
      cols <- seq_len(count)
      thresholds[spec$variable, cols] <- lmsCatDeltaToThresholds(delta, spec)
      threshold.delta[spec$variable, cols] <- delta[spec$index]
      threshold.labels[spec$variable, cols] <- spec$labels
      if (!is.null(vcov.delta)) {
        J <- lmsCatThresholdJacobian(delta, threshold.info)
        threshold.se[spec$variable, cols] <- sqrt(pmax(diag(
          J[spec$index, , drop = FALSE] %*% vcov.delta %*%
            t(J[spec$index, , drop = FALSE])
        ), 0))
      }
    }

    submodel$matrices$thresholds <- thresholds
    submodel$matrices$thresholdDelta <- threshold.delta
    submodel$labelMatrices$thresholds <- threshold.labels
    submodel$labelMatrices$thresholdDelta <- threshold.labels
    if (!is.null(submodel$matricesNA)) {
      skeleton <- thresholds
      skeleton[is.finite(skeleton)] <- NA_real_
      submodel$matricesNA$thresholds <- skeleton
      submodel$matricesNA$thresholdDelta <- skeleton
    }
    if (!is.null(submodel$matricesSE)) {
      submodel$matricesSE$thresholds <- threshold.se
      submodel$matricesSE$thresholdDelta <- threshold.se
    }
    submodel$info$ordered <- names(specs) |>
      vapply(function(key) threshold.info$specs[[key]]$variable, character(1L)) |>
      unname()
    model$models[[g]] <- submodel
  }
  model$info$ordered <- unique(unlist(lapply(model$models, function(x) x$info$ordered)))
  model
}


lmsCatRegisterThresholdParameters <- function(model, threshold.info) {
  old.theta <- model$theta
  model <- lmsCatAttachThresholdMatrices(
    model, threshold.info, threshold.info$delta
  )

  threshold.starts <- numeric()
  for (g in seq_along(model$models)) {
    submodel <- model$models[[g]]
    free <- is.finite(submodel$matrices$thresholdDelta)
    starts.g <- submodel$matrices$thresholdDelta[free]
    names(starts.g) <- getParamNamesMatrix(
      submodel$matrices$thresholdDelta, "thresholdDelta"
    )[free]
    if (g > 1L) names(starts.g) <- sprintf("%s.g%d", names(starts.g), g)
    threshold.starts <- c(threshold.starts, starts.g)
    submodel$matrices$thresholdDelta[free] <- NA_real_
    submodel$labelMatrices$thresholdDelta[] <- ""
    submodel$labelMatrices$thresholds[] <- ""
    model$models[[g]] <- submodel
  }

  params <- createTheta(model, parTable.in = model$parTable)
  common <- intersect(names(old.theta), names(params$theta))
  params$theta[common] <- old.theta[common]
  threshold.index <- grep("^thresholdDelta", names(params$theta))
  mod_stopif(length(threshold.index) != length(threshold.info$delta),
             sprintf(paste0("Threshold parameter registration produced an invalid ",
                            "parameter count (%d registered, %d expected)."),
                     length(threshold.index), length(threshold.info$delta)))
  params$theta[names(threshold.starts)] <- threshold.starts
  model$params[names(params)] <- params
  model$theta <- params$theta
  model$params$bounds <- getParamBounds(model)
  model$params$gradientStruct <- getGradientStruct(
    model, theta = model$theta, method = "lms"
  )
  model
}


lmsCatThresholdThetaIndex <- function(model, threshold.info) {
  index <- integer()
  for (spec in threshold.info$specs) {
    delta.matrix <- model$models[[spec$group]]$matrices$thresholdDelta
    parameter.names <- matrix(
      getParamNamesMatrix(delta.matrix, "thresholdDelta"),
      nrow = NROW(delta.matrix), ncol = NCOL(delta.matrix),
      dimnames = dimnames(delta.matrix)
    )
    names.spec <- parameter.names[spec$variable, seq_len(spec$K - 1L)]
    if (spec$group > 1L)
      names.spec <- sprintf("%s.g%d", names.spec, spec$group)
    index <- c(index, match(names.spec, names(model$theta)))
  }
  mod_stopif(anyNA(index), "Unable to map thresholds into the DA parameter vector.")
  index
}


lmsCatParameterJacobian <- function(theta, model) {
  structure <- model$params$gradientStruct
  J <- refreshCovModelJacobian(theta, model, structure$Jacobian)$J
  if (length(structure$nlinDerivs)) {
    eval.theta <- structure$evalTheta
    parameter.full <- stringr::str_split_i(colnames(J), pattern = "#", i = 1L)
    parameter.part <- rownames(J)
    theta.environment <- list2env(as.list(eval.theta(theta)))
    for (dependent in names(structure$nlinDerivs)) {
      for (independent in names(structure$nlinDerivs[[dependent]])) {
        derivative <- eval(structure$nlinDerivs[[dependent]][[independent]],
                           envir = theta.environment)
        J[parameter.part == independent, parameter.full == dependent] <- derivative
      }
    }
  }
  J
}


lmsCatTensorRule <- function(order, dimension) {
  q <- quadrature(m = as.integer(order), k = dimension, adaptive = FALSE)
  list(nodes = q$n, weights = q$w, log.adjust = rep(0, NROW(q$n)),
       method = "gh", signed = FALSE,
       dimension = dimension, nodes.total = NROW(q$n))
}


lmsCatAdaptiveTarget <- function(u, log.kernel) {
  val <- log.kernel(matrix(u, nrow = 1L))
  if (length(val) != 1L || !is.finite(val)) return(.Machine$double.xmax / 100)
  -(val - .5 * sum(u^2))
}


lmsCatAdaptiveRule <- function(rule, log.kernel, start = NULL,
                               mode.tol = 1e-8, maxit = 100L) {
  # `rule` integrates against a standard-normal density.  We locate the mode
  # and curvature of log{L(u) phi(u)}, then express the same integral as an
  # expectation under N(mode, curvature^-1).  `log.adjust` is the exact
  # phi(u)/q(u) importance correction; consequently adaptation changes only
  # numerical efficiency, never the estimand.
  d <- NCOL(rule$nodes)
  if (!d) return(rule)
  if (is.null(start)) start <- rep(0, d)
  mode <- stats::optim(start, lmsCatAdaptiveTarget, log.kernel = log.kernel,
                       method = "BFGS",
                       control = list(reltol = mode.tol, maxit = maxit))
  H <- tryCatch(stats::optimHess(mode$par, lmsCatAdaptiveTarget,
                                 log.kernel = log.kernel),
                error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) return(rule)
  H <- (H + t(H)) / 2
  R <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R)) return(rule)
  L <- backsolve(R, diag(d))
  # L L' = H^-1 for the upper-triangular R returned by chol(H).
  z <- rule$nodes
  u <- sweep(z %*% t(L), 2L, mode$par, "+")
  log.det <- sum(log(diag(L)))
  log.adjust <- -.5 * rowSums(u^2) + .5 * rowSums(z^2) + log.det
  rule$nodes <- u
  rule$log.adjust <- log.adjust
  rule$method <- "aghq"
  rule$adaptation <- list(mode = mode$par, covariance = L %*% t(L),
                          converged = mode$convergence == 0L)
  rule
}


lmsCatGenzKeisterRule <- function(level) {
  # Nested GK22 family, normalized from exp(-x^2) to N(0, 1).  Levels 1--4
  # have orders 1, 3, 9, and 19 and polynomial precision 1, 5, 15, and 29.
  level <- as.integer(level)
  mod_stopif(level < 1L || level > 4L,
             "Genz-Keister levels must currently be between 1 and 4.")
  nodes <- list(
    0,
    c(0, 1.2247448713915889),
    c(0, 0.52403354748695763, 1.2247448713915889,
      2.0232301911005157, 2.9592107790638380),
    c(0, 0.52403354748695763, 0.87004089535290285,
      1.2247448713915889, 1.8357079751751868, 2.0232301911005157,
      2.2665132620567876, 2.9592107790638380, 3.6677742159463378,
      4.4995993983103881)
  )[[level]]
  weights <- list(
    1.7724538509055159,
    c(1.1816359006036772, 0.29540897515091930),
    c(0.45014700975378197, 0.47869428549114124,
      0.16811892894767771, 0.014173117873979098,
      1.6708826306882348e-4),
    c(0.53788160700510168, 0.36924643368920851,
      0.10838861955003017, 0.11360729895748269,
      0.032055243099445879, -0.011232438489069229,
      0.0051133174390883855, 1.0656589772852267e-4,
      1.0802767206624762e-6, 1.5295717705322357e-9)
  )[[level]]
  if (length(nodes) > 1L) {
    nodes <- c(-rev(nodes[-1L]), nodes)
    weights <- c(rev(weights[-1L]), weights)
  }
  list(n = matrix(sqrt(2) * nodes, ncol = 1L),
       w = weights / sqrt(pi), level = level)
}


lmsCatSparseRule <- function(level, dimension, max.nodes = 50000L) {
  # Smolyak construction over the nested Genz-Keister sequence.  Identical
  # abscissas are merged after combination, giving substantial node reuse.
  level <- min(as.integer(level), 4L)
  q <- as.integer(level) + dimension - 1L
  grids <- list(); weights <- list(); at <- 0L
  recurse <- function(prefix, left, slots) {
    if (slots == 1L) return(list(c(prefix, left)))
    out <- list()
    for (i in seq_len(left - slots + 1L))
      out <- c(out, recurse(c(prefix, i), left - i, slots - 1L))
    out
  }
  for (s in seq.int(max(dimension, q - dimension + 1L), q)) {
    coef <- (-1)^(q - s) * choose(dimension - 1L, q - s)
    if (coef == 0) next
    for (idx in recurse(integer(), s, dimension)) {
      rules <- lapply(idx, function(i) {
        lmsCatGenzKeisterRule(i)
      })
      n <- do.call(expand.grid, lapply(rules, `[[`, "n")) |> as.matrix()
      w <- Reduce(kronecker, rev(lapply(rules, `[[`, "w"))) * coef
      at <- at + 1L; grids[[at]] <- n; weights[[at]] <- w
    }
  }
  nodes <- do.call(rbind, grids); w <- unlist(weights, use.names = FALSE)
  key <- apply(round(nodes, 13L), 1L, paste, collapse = ":")
  sums <- rowsum(w, key, reorder = FALSE)[, 1L]
  nodes <- nodes[match(rownames(rowsum(w, key, reorder = FALSE)), key), , drop = FALSE]
  w <- sums
  nz <- abs(w) > .Machine$double.eps
  nodes <- nodes[nz, , drop = FALSE]; w <- w[nz]
  mod_stopif(NROW(nodes) > max.nodes,
             sprintf("Sparse grid requires %d nodes (maximum %d).", NROW(nodes), max.nodes))
  list(nodes = nodes, weights = w, log.adjust = rep(0, NROW(nodes)),
       method = "sparse", family = "genz-keister", level = level,
       signed = any(w < 0), dimension = dimension, nodes.total = NROW(nodes))
}


lmsCatQmcRule <- function(points, dimension, replicates = 4L, seed = 19790417L) {
  points <- as.integer(points); replicates <- as.integer(replicates)
  mod_stopif(points < 1L, "`qmc.points` must be at least 1.")
  mod_stopif(replicates < 1L, "`qmc.replicates` must be at least 1.")
  nodes <- vector("list", replicates)
  for (r in seq_len(replicates)) {
    u <- qrng::sobol(n = points, d = dimension, randomize = "digital.shift",
                     seed = as.integer(seed + r - 1L))
    u <- matrix(u, ncol = dimension)
    u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    nodes[[r]] <- stats::qnorm(u)
  }
  list(nodes = do.call(rbind, nodes),
       weights = rep(1 / (points * replicates), points * replicates),
       replicate = rep(seq_len(replicates), each = points), method = "qmc",
       signed = FALSE, dimension = dimension, nodes.total = points * replicates,
       points = points, replicates = replicates, seed = seed)
}


lmsCatCompressGroup <- function(raw, indicators, ordered, weights = NULL) {
  if (!all(indicators %in% ordered))
    return(list(data = raw, weights = weights, compressed = FALSE))
  encoded <- lapply(raw[, indicators, drop = FALSE], function(x) {
    if (is.factor(x)) as.integer(x) else match(x, sort(unique(x[!is.na(x)])))
  })
  encoded <- as.data.frame(encoded)
  encoded[] <- lapply(encoded, function(x) replace(x, is.na(x), 0L))
  key <- do.call(paste, c(encoded, sep = ":"))
  keep <- !duplicated(key)
  if (is.null(weights)) weights <- rep(1, NROW(raw))
  summed <- rowsum(weights, key, reorder = FALSE)[, 1L]
  idx <- match(rownames(rowsum(weights, key, reorder = FALSE)), key)
  list(data = raw[idx, , drop = FALSE], weights = summed,
       compressed = TRUE, n.original = NROW(raw), n.patterns = sum(keep))
}


lmsCatIntegrationRule <- function(method, dimension, control) {
  requested <- method
  reason <- sprintf("explicit `%s` selection", method)
  if (method == "auto") {
    dense.nodes <- control$order^dimension
    if (dimension <= 3L && dense.nodes <= control$max.nodes) {
      method <- "aghq"
      reason <- sprintf("product rule is feasible (%d dimensions, %d nodes)",
                        dimension, dense.nodes)
    } else if (dimension <= 6L) {
      method <- "sparse"
      reason <- sprintf("product rule is too large; using a nested sparse grid in %d dimensions",
                        dimension)
    } else {
      method <- "qmc"
      reason <- sprintf("integration dimension %d exceeds the sparse-grid cutoff", dimension)
    }
  }
  if (method == "aghq") {
    mod_stopif(control$order^dimension > control$max.nodes,
               "Product Gauss-Hermite rule exceeds `max.nodes`; use sparse or QMC integration.")
    ans <- lmsCatTensorRule(control$order, dimension)
    ans$requested <- requested; ans$selection.reason <- reason
    return(ans)
  }
  if (method == "sparse") {
    ans <- tryCatch(lmsCatSparseRule(control$sparse.level, dimension, control$max.nodes),
                    error = function(e) NULL)
    if (!is.null(ans)) {
      ans$requested <- requested; ans$selection.reason <- reason
      return(ans)
    }
    method <- "qmc"
    reason <- paste0(reason, "; sparse grid exceeded `max.nodes`, so QMC was used")
  }
  ans <- lmsCatQmcRule(control$qmc.points, dimension, control$qmc.replicates, control$seed)
  ans$requested <- requested; ans$selection.reason <- reason
  ans
}


lmsCatOrdinalNuisanceTheta <- function(model, ordered) {
  loc <- model$params$gradientStruct$locations
  idx <- as.integer(sub(".*#", "", loc$param))
  nuisance <- (loc$block == DA_BLOCKS$tauX &
                 rownames(model$models[[1L]]$matrices$lambdaX)[loc$row + 1L] %in% ordered) |
    (loc$block == DA_BLOCKS$thetaDelta &
       (rownames(model$models[[1L]]$matrices$lambdaX)[loc$row + 1L] %in% ordered |
        rownames(model$models[[1L]]$matrices$lambdaX)[loc$col + 1L] %in% ordered))
  unique(idx[nuisance & !idx %in% idx[!nuisance]])
}


lmsCatActiveTheta <- function(model, ordered) {
  setdiff(seq_along(model$theta), lmsCatOrdinalNuisanceTheta(model, ordered))
}


lmsCatStateOne <- function(M, node) {
  p <- NROW(M$A); r <- NROW(M$psi)
  ux <- matrix(node[seq_len(p)], nrow = 1L)
  xi <- c(M$beta0) + c(ux %*% t(M$A))
  if (!r) return(xi)
  zproj <- solve(M$A, t(M$covZetaXi)) |> t()
  psi.orth <- M$psi - zproj %*% t(zproj)
  psi.orth <- (psi.orth + t(psi.orth)) / 2
  L <- tryCatch(t(chol(psi.orth)), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  ue <- matrix(node[p + seq_len(r)], nrow = 1L)
  zeta <- c(ux %*% t(zproj) + ue %*% t(L))
  Ie <- diag(r)
  x <- matrix(xi, ncol = 1L)
  kz <- kronecker(Ie, x)
  B <- Ie - M$gammaEta - t(kz) %*% M$omegaEtaXi
  rhs <- M$alpha + M$gammaXi %*% x + t(kz) %*% M$omegaXiXi %*% x +
    matrix(zeta, ncol = 1L)
  c(xi, solve(B, rhs))
}


lmsCatLatentStates <- function(M, nodes) {
  ans <- t(vapply(seq_len(NROW(nodes)), function(i) lmsCatStateOne(M, nodes[i, ]),
                  numeric(NROW(M$A) + NROW(M$psi))))
  if (any(!is.finite(ans))) NULL else ans
}


lmsCatReduction <- function(submodel, ordered, tol = 1e-8) {
  M <- submodel$matrices
  factors <- c(submodel$info$xis, submodel$info$etas)
  cat.factors <- names(Filter(function(x) any(x %in% ordered),
                             c(submodel$info$indsXis, submodel$info$indsEtas)))
  cat.rows <- match(cat.factors, factors)
  d <- NROW(M$A) + NROW(M$psi)
  explicit <- integer()
  probes <- rbind(rep(0.37, d), seq(-.41, .43, length.out = d))

  # Include every independent source that can alter a categorical factor.  A
  # nonzero probe avoids missing pure interaction dependencies at the origin.
  for (j in seq_len(d)) {
    depends <- FALSE
    for (b in seq_len(NROW(probes))) {
      u0 <- probes[b, ]; u1 <- u0; u1[j] <- u1[j] + .61
      g0 <- lmsCatStateOne(M, u0); g1 <- lmsCatStateOne(M, u1)
      if (any(abs(g1[cat.rows] - g0[cat.rows]) > tol)) depends <- TRUE
    }
    if (depends) explicit <- c(explicit, j)
  }

  # Conditional on `explicit`, the remainder must be affine.  If a quadratic
  # or cross-product remains, condition one member of that term and repeat.
  repeat {
    analytic <- setdiff(seq_len(d), explicit)
    if (!length(analytic)) break
    base <- numeric(d); g0 <- lmsCatStateOne(M, base)
    offending <- integer()
    for (j in analytic) {
      up <- base; um <- base; up[j] <- 1; um[j] <- -1
      if (max(abs(lmsCatStateOne(M, up) + lmsCatStateOne(M, um) - 2 * g0)) > tol)
        offending <- c(offending, j)
    }
    if (!length(offending) && length(analytic) > 1L) {
      pairs <- utils::combn(analytic, 2L)
      for (k in seq_len(NCOL(pairs))) {
        j <- pairs[1L, k]; h <- pairs[2L, k]
        uj <- base; uh <- base; ujh <- base
        uj[j] <- 1; uh[h] <- 1; ujh[c(j, h)] <- 1
        cross <- lmsCatStateOne(M, ujh) - lmsCatStateOne(M, uj) -
          lmsCatStateOne(M, uh) + g0
        if (max(abs(cross)) > tol) { offending <- j; break }
      }
    }
    if (!length(offending)) break
    explicit <- sort(unique(c(explicit, offending[1L])))
  }
  list(explicit = sort(unique(explicit)),
       analytic = setdiff(seq_len(d), explicit), full.dimension = d,
       factor.names = factors)
}


lmsCatReducedMoments <- function(M, nodes, reduction) {
  Q <- NROW(nodes); p <- NROW(M$A) + NROW(M$psi)
  means <- matrix(NA_real_, Q, p)
  covs <- vector("list", Q)
  for (q in seq_len(Q)) {
    u <- numeric(reduction$full.dimension)
    u[reduction$explicit] <- nodes[q, ]
    g0 <- lmsCatStateOne(M, u)
    means[q, ] <- g0
    if (!length(reduction$analytic)) covs[[q]] <- matrix(0, p, p)
    else {
      J <- vapply(reduction$analytic, function(j) {
        uj <- u; uj[j] <- 1
        lmsCatStateOne(M, uj) - g0
      }, numeric(p))
      if (is.null(dim(J))) J <- matrix(J, ncol = 1L)
      covs[[q]] <- J %*% t(J)
    }
  }
  list(mean = means, cov = covs)
}


lmsCatLogLikGroup <- function(submodel, raw, ordered, threshold.info, delta,
                              rule, reduction, group.index = 1L, weights = NULL,
                              compiled = TRUE, score = FALSE, details = FALSE,
                              posterior.details = TRUE) {
  M <- submodel$matrices
  indicators <- rownames(M$lambdaX)
  continuous <- setdiff(indicators, ordered)
  moments <- if (compiled && !length(continuous)) {
    means <- tryCatch(
      lmsCatLatentMeansCpp(M, rule$nodes, reduction$explicit,
                           reduction$full.dimension),
      error = function(e) NULL
    )
    if (is.null(means)) NULL else list(mean = means, cov = NULL)
  } else {
    tryCatch(lmsCatReducedMomentsCpp(M, rule$nodes, reduction$explicit,
                                     reduction$analytic,
                                     reduction$full.dimension),
             error = function(e) NULL)
  }
  if (is.null(moments) || any(!is.finite(moments$mean))) return(-Inf)
  raw <- raw[, indicators, drop = FALSE]
  mu <- sweep(moments$mean %*% t(M$lambdaX), 2L, c(M$tauX), "+")
  N <- NROW(raw); Q <- NROW(rule$nodes)

  if (compiled && !length(continuous)) {
    codes <- do.call(cbind, lapply(indicators, function(v) {
      x <- raw[[v]]
      if (is.factor(x)) as.integer(x)
      else match(x, sort(unique(x[!is.na(x)])))
    }))
    codes <- matrix(codes, nrow = N, ncol = length(indicators))
    thresholds <- lapply(indicators, function(v) {
      spec <- threshold.info$specs[[paste(group.index, v, sep = "::")]]
      lmsCatDeltaToThresholds(delta, spec)
    })
    if (is.null(weights)) weights <- rep(1, N)
    ans <- ordinalLogLikCpp(
      codes, mu, thresholds, rule$weights,
      rule$log.adjust %||% rep(0, Q), weights, gradient = score,
      details = details, posteriorDetails = posterior.details
    )
    if (score) {
      ans$mu <- mu
      return(ans)
    }
    if (details) return(ans)
    return(ans$logLik)
  }

  if (compiled && length(continuous)) {
    categorical <- indicators[indicators %in% ordered]
    categorical.idx <- match(categorical, indicators)
    continuous.idx <- match(continuous, indicators)
    codes <- do.call(cbind, lapply(categorical, function(v) {
      x <- raw[[v]]
      if (is.factor(x)) as.integer(x)
      else match(x, sort(unique(x[!is.na(x)])))
    }))
    codes <- matrix(codes, nrow = N, ncol = length(categorical))
    thresholds <- lapply(categorical, function(v) {
      spec <- threshold.info$specs[[paste(group.index, v, sep = "::")]]
      lmsCatDeltaToThresholds(delta, spec)
    })
    sigma.nodes <- lapply(seq_len(Q), function(q) {
      full <- M$lambdaX %*% moments$cov[[q]] %*% t(M$lambdaX) + M$thetaDelta
      full[continuous.idx, continuous.idx, drop = FALSE]
    })
    continuous.data <- as.matrix(raw[, continuous, drop = FALSE])
    storage.mode(continuous.data) <- "double"
    if (is.null(weights)) weights <- rep(1, N)
    ans <- mixedLogLikCpp(
      codes, mu[, categorical.idx, drop = FALSE], thresholds,
      continuous.data, mu[, continuous.idx, drop = FALSE], sigma.nodes,
      rule$weights, rule$log.adjust %||% rep(0, Q), weights,
      gradient = score, details = details,
      posteriorDetails = posterior.details
    )
    if (score) {
      if (!is.finite(ans$logLik)) return(ans)
      predictor.score <- matrix(0, Q, length(indicators),
                                dimnames = list(NULL, indicators))
      predictor.score[, categorical.idx] <- ans$categorical.score
      predictor.score[, continuous.idx] <- ans$continuous.score
      ans$predictor.score <- predictor.score
      ans$mu <- mu
      ans$continuous <- continuous
      ans$continuous.index <- continuous.idx
      ans$covariance <- sigma.nodes
      return(ans)
    }
    if (details) return(ans)
    return(ans$logLik)
  }
  logc <- matrix(0, N, Q)

  for (v in intersect(ordered, indicators)) {
    spec <- threshold.info$specs[[paste(group.index, v, sep = "::")]]
    tau <- lmsCatDeltaToThresholds(delta, spec)
    x <- raw[[v]]
    code <- if (is.factor(x)) as.integer(x) else match(x, sort(unique(x[!is.na(x)])))
    j <- match(v, indicators)
    for (i in which(!is.na(code))) {
      lo <- if (code[i] == 1L) -Inf else tau[code[i] - 1L]
      hi <- if (code[i] == spec$K) Inf else tau[code[i]]
      pr <- stats::plogis(hi - mu[, j]) - stats::plogis(lo - mu[, j])
      logc[i, ] <- logc[i, ] + log(pmax(pr, .Machine$double.xmin))
    }
  }

  if (length(continuous)) {
    idx <- match(continuous, indicators)
    sigma.nodes <- lapply(seq_len(Q), function(q) {
      M$lambdaX %*% moments$cov[[q]] %*% t(M$lambdaX) + M$thetaDelta
    })
    for (i in seq_len(N)) {
      obs <- which(!is.na(raw[i, continuous, drop = TRUE]))
      if (!length(obs)) next
      ii <- idx[obs]
      y <- as.numeric(raw[i, continuous[obs], drop = TRUE])
      for (q in seq_len(Q)) {
        S <- sigma.nodes[[q]][ii, ii, drop = FALSE]
        Li <- tryCatch(chol(S), error = function(e) NULL)
        if (is.null(Li)) return(-Inf)
        d <- mu[q, ii] - y
        z <- backsolve(Li, d, transpose = TRUE)
        logc[i, q] <- logc[i, q] - .5 *
          (length(obs) * log(2 * pi) + 2 * sum(log(diag(Li))) + sum(z^2))
      }
    }
  }

  marginal <- numeric(N)
  if (!rule$signed) {
    lw <- log(rule$weights) + (rule$log.adjust %||% 0)
    for (i in seq_len(N)) {
      a <- logc[i, ] + lw; m <- max(a)
      marginal[i] <- m + log(sum(exp(a - m)))
    }
  } else {
    adjusted.weights <- rule$weights * exp(rule$log.adjust %||% 0)
    pos <- adjusted.weights > 0; neg <- adjusted.weights < 0
    for (i in seq_len(N)) {
      m <- max(logc[i, ])
      val <- sum(adjusted.weights[pos] * exp(logc[i, pos] - m)) -
        sum(-adjusted.weights[neg] * exp(logc[i, neg] - m))
      if (!is.finite(val) || val <= 0) return(-Inf)
      marginal[i] <- m + log(val)
    }
  }
  if (is.null(weights)) weights <- rep(1, N)
  sum(weights * marginal)
}


lmsCatAdaptiveLogLikGroup <- function(submodel, raw, ordered, threshold.info,
                                      delta, rules, reduction,
                                      group.index = 1L, weights = NULL) {
  M <- submodel$matrices
  indicators <- rownames(M$lambdaX)
  continuous <- setdiff(indicators, ordered)
  categorical <- intersect(indicators, ordered)
  N <- NROW(raw)
  mod_stopif(length(rules) != N,
             "One adaptive quadrature rule is required per response pattern.")
  nodes.per.pattern <- vapply(rules, function(x) NROW(x$nodes), integer(1L))
  mod_stopif(length(unique(nodes.per.pattern)) != 1L,
             "Adaptive rules must currently have equal node counts.")
  Q <- nodes.per.pattern[1L]
  nodes <- do.call(rbind, lapply(rules, `[[`, "nodes"))
  quad.weights <- unlist(lapply(rules, `[[`, "weights"), use.names = FALSE)
  log.adjust <- unlist(lapply(rules, function(x) x$log.adjust %||% rep(0, Q)),
                       use.names = FALSE)
  if (is.null(weights)) weights <- rep(1, N)
  raw <- raw[, indicators, drop = FALSE]
  moments <- if (!length(continuous)) {
    means <- tryCatch(lmsCatLatentMeansCpp(
      M, nodes, reduction$explicit, reduction$full.dimension
    ), error = function(e) NULL)
    if (is.null(means)) NULL else list(mean = means, cov = NULL)
  } else {
    tryCatch(lmsCatReducedMomentsCpp(M, nodes, reduction$explicit,
                                     reduction$analytic,
                                     reduction$full.dimension),
             error = function(e) NULL)
  }
  if (is.null(moments) || any(!is.finite(moments$mean))) return(-Inf)
  mu <- sweep(moments$mean %*% t(M$lambdaX), 2L, c(M$tauX), "+")
  codes <- do.call(cbind, lapply(categorical, function(v) {
    x <- raw[[v]]
    if (is.factor(x)) as.integer(x)
    else match(x, sort(unique(x[!is.na(x)])))
  }))
  codes <- matrix(codes, nrow = N, ncol = length(categorical))
  thresholds <- lapply(categorical, function(v) {
    spec <- threshold.info$specs[[paste(group.index, v, sep = "::")]]
    lmsCatDeltaToThresholds(delta, spec)
  })
  if (!length(continuous)) {
    return(adaptiveOrdinalLogLikCpp(
      codes, mu, thresholds, quad.weights, log.adjust, weights, Q
    ))
  }
  categorical.idx <- match(categorical, indicators)
  continuous.idx <- match(continuous, indicators)
  sigma.nodes <- lapply(seq_len(NROW(nodes)), function(q) {
    full <- M$lambdaX %*% moments$cov[[q]] %*% t(M$lambdaX) + M$thetaDelta
    full[continuous.idx, continuous.idx, drop = FALSE]
  })
  continuous.data <- as.matrix(raw[, continuous, drop = FALSE])
  storage.mode(continuous.data) <- "double"
  adaptiveMixedLogLikCpp(
    codes, mu[, categorical.idx, drop = FALSE], thresholds,
    continuous.data, mu[, continuous.idx, drop = FALSE], sigma.nodes,
    quad.weights, log.adjust, weights, Q
  )
}


modsemLmsCat <- function(model.syntax, data, ordered = NULL, group = NULL,
                         sampling.weights = NULL,
                         sampling.weights.normalization = "total",
                         missing = NULL, nodes = NULL,
                         integration = "auto", integration.control = list(),
                         calc.se = TRUE, robust.se = FALSE, cluster = NULL,
                         start = NULL, optimizer = "nlminb", max.iter = 500L,
                         convergence.abs = 1e-5, verbose = interactive(),
                         adaptive.frequency = NULL,
                         requested.method = "lms-cat", cov.syntax = NULL,
                         auto.fix.first = NULL, auto.fix.single = NULL,
                         n.threads = NULL, ...) {
  data <- as.data.frame(data)
  ordered <- lmsCatOrderedCols(data, model.syntax, ordered)
  ctl <- lmsCatControl(integration.control, nodes)
  ctl$adaptive.frequency <- adaptive.frequency %||% ctl$adaptive.frequency
  calc.se <- isTRUE(calc.se %||% TRUE)
  verbose <- isTRUE(verbose %||% interactive())
  if (isTRUE(robust.se) || !is.null(cluster)) {
    mod_msg_warn("Robust and cluster-robust standard errors are not available yet for `lms-cat`; using the observed Hessian.")
  }

  progress <- new.env(parent = emptyenv())
  progress$evaluations <- 0L
  progress$best <- -Inf
  progress$started <- proc.time()[["elapsed"]]
  progress$last <- progress$started
  progress$stage <- "initialization"
  progress$score.evaluations <- 0L
  progress$light.evaluations <- 0L
  beginStage <- function(stage, message) {
    progress$stage <- stage
    if (verbose) {
      clearConsoleLine()
      printf("lms-jv: %s\n", message)
    }
  }

  beginStage("initialization", "fitting the continuous-model initializer")
  scored <- lmsCatScoreData(data, ordered)
  start.fit <- modsem_da(
    model.syntax = model.syntax, data = scored, method = "lms", ordered = NULL,
    group = group, cov.syntax = cov.syntax, sampling.weights = sampling.weights,
    sampling.weights.normalization = sampling.weights.normalization,
    missing = missing, nodes = min(12L, ctl$order), calc.se = FALSE,
    max.iter = ctl$start.max.iter, convergence.abs = 1e-3,
    auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
    n.threads = n.threads, verbose = FALSE
  )
  model <- start.fit$start.model
  theta0 <- start.fit$theta
  # A cumulative-link item has thresholds but no separately identified
  # intercept or Gaussian residual variance.  The continuous initializer
  # estimates both, so remove them before constructing threshold starts.  In
  # particular, leaving the initializer's item intercept fixed would shift all
  # reported thresholds away from the usual zero-intercept parameterization.
  ordinal.nuisance <- lmsCatOrdinalNuisanceTheta(model, ordered)
  theta0[ordinal.nuisance] <- 0
  model$theta <- theta0
  threshold.info <- lmsCatThresholdSpec(data, ordered, group)
  model <- lmsCatRegisterThresholdParameters(model, threshold.info)
  theta0 <- model$theta
  threshold.active <- lmsCatThresholdThetaIndex(model, threshold.info)
  active <- setdiff(lmsCatActiveTheta(model, ordered), threshold.active)

  # The full latent innovation representation is the correctness oracle.  A
  # subsequent reduction pass may replace this set without changing the kernel.
  filled0 <- fillModel(model, theta0, fillPhi = TRUE, method = "lms")
  reductions <- lapply(filled0$models, lmsCatReduction, ordered = ordered)
  dims <- vapply(reductions, function(x) length(x$explicit), integer(1L))
  mod_stopif(length(unique(dims)) != 1L,
             "All groups must currently have the same categorical LMS integration dimension.")
  rule <- lmsCatIntegrationRule(integration, dims[1L], ctl)
  if (verbose) {
    printf("lms-jv: integration=%s, explicit dimensions=%d, analytic dimensions=%d, nodes=%d\n",
           rule$method, dims[1L], length(reductions[[1L]]$analytic), rule$nodes.total)
  }

  p0 <- c(theta0[active], theta0[threshold.active])
  lower <- c(model$params$bounds$lower[active],
             model$params$bounds$lower[threshold.active])
  upper <- c(model$params$bounds$upper[active],
             model$params$bounds$upper[threshold.active])
  ng <- length(model$models)
  raw.groups <- if (is.null(group)) list(data) else {
    lev <- unique(data[[group]])
    lapply(lev, function(x) data[data[[group]] == x, , drop = FALSE])
  }
  compressed <- vector("list", ng)
  for (g in seq_len(ng)) {
    indicators <- rownames(model$models[[g]]$matrices$lambdaX)
    wg <- if (!is.null(sampling.weights)) raw.groups[[g]][[sampling.weights]] else NULL
    compressed[[g]] <- lmsCatCompressGroup(raw.groups[[g]], indicators, ordered, wg)
  }

  adaptive.rules <- NULL
  baseLogLik <- function(p, evaluation.rule = rule) {
    theta <- theta0; theta[active] <- p[seq_along(active)]
    delta <- p[length(active) + seq_along(threshold.active)]
    theta[threshold.active] <- delta
    filled <- tryCatch(fillModel(model, theta, fillPhi = TRUE, method = "lms"),
                       error = function(e) NULL)
    if (is.null(filled)) return(-Inf)
    values <- vapply(seq_len(ng), function(g) {
      lmsCatLogLikGroup(filled$models[[g]], compressed[[g]]$data, ordered,
                        threshold.info, delta, evaluation.rule, reductions[[g]], g,
                        compressed[[g]]$weights)
    }, numeric(1L))
    if (any(!is.finite(values))) -Inf else sum(values)
  }
  qmcDiagnostics <- function(p, qmc.rule) {
    reps <- sort(unique(qmc.rule$replicate))
    estimates <- vapply(reps, function(r) {
      use <- qmc.rule$replicate == r
      rr <- qmc.rule
      rr$nodes <- qmc.rule$nodes[use, , drop = FALSE]
      rr$weights <- rep(1 / sum(use), sum(use))
      rr$replicate <- rep(1L, sum(use))
      rr$nodes.total <- sum(use)
      rr$log.adjust <- rep(0, sum(use))
      baseLogLik(p, rr)
    }, numeric(1L))
    finite <- is.finite(estimates)
    standard.error <- if (sum(finite) > 1L)
      stats::sd(estimates[finite]) / sqrt(sum(finite)) else NA_real_
    list(replicate.logLik = estimates,
         mean.logLik = if (any(finite)) mean(estimates[finite]) else NA_real_,
         standard.error = standard.error,
         range = if (any(finite)) range(estimates[finite]) else c(NA_real_, NA_real_),
         points = qmc.rule$points, replicates = length(reps),
         tolerance = ctl$qmc.loglik.tol,
         reliable = is.finite(standard.error) && standard.error <= ctl$qmc.loglik.tol)
  }
  if (identical(integration, "auto") && identical(rule$method, "sparse") &&
      !is.finite(baseLogLik(p0, rule))) {
    reason <- paste0(rule$selection.reason,
                     "; the sparse likelihood was invalid at the initializer, so QMC was used")
    rule <- lmsCatQmcRule(ctl$qmc.points, dims[1L], ctl$qmc.replicates, ctl$seed)
    rule$requested <- "auto"
    rule$selection.reason <- reason
    if (verbose) printf("lms-jv: sparse-grid check failed; switching to QMC (%d nodes)\n",
                        rule$nodes.total)
  }
  objective <- function(p) {
    progress$evaluations <- progress$evaluations + 1L
    theta <- theta0; theta[active] <- p[seq_along(active)]
    delta <- p[length(active) + seq_along(threshold.active)]
    theta[threshold.active] <- delta
    filled <- tryCatch(fillModel(model, theta, fillPhi = TRUE, method = "lms"),
                       error = function(e) NULL)
    if (is.null(filled)) return(.Machine$double.xmax / 100)
    ll <- 0
    for (g in seq_len(ng)) {
      if (is.null(adaptive.rules)) {
        value <- lmsCatLogLikGroup(filled$models[[g]], compressed[[g]]$data, ordered,
                                   threshold.info, delta, rule, reductions[[g]], g,
                                   compressed[[g]]$weights)
      } else {
        value <- lmsCatAdaptiveLogLikGroup(
          filled$models[[g]], compressed[[g]]$data, ordered,
          threshold.info, delta, adaptive.rules[[g]], reductions[[g]], g,
          compressed[[g]]$weights
        )
      }
      if (!is.finite(value)) return(.Machine$double.xmax / 100)
      ll <- ll + value
    }
    progress$best <- max(progress$best, ll)
    -ll
  }

  all.categorical <- all(vapply(seq_len(ng), function(g) {
    indicators <- rownames(model$models[[g]]$matrices$lambdaX)
    all(indicators %in% ordered)
  }, logical(1L)))
  use.gradient <- isTRUE(ctl$gradient)
  analytical.predictor <- isTRUE(ctl$gradient.analytical)
  ema.cache <- new.env(parent = emptyenv())
  ema.cache$par <- NULL
  ema.cache$state <- NULL
  gradientFunction <- if (use.gradient) function(p) {
    cached <- !is.null(ema.cache$par) && identical(c(p), c(ema.cache$par))
    if (cached) {
      theta <- ema.cache$state$theta
      delta <- ema.cache$state$delta
      filled <- ema.cache$state$filled
      scores <- ema.cache$state$scores
    } else {
      theta <- theta0; theta[active] <- p[seq_along(active)]
      delta <- p[length(active) + seq_along(threshold.active)]
      theta[threshold.active] <- delta
      filled <- tryCatch(fillModel(model, theta, fillPhi = TRUE, method = "lms"),
                         error = function(e) NULL)
      if (is.null(filled)) return(rep(NA_real_, length(p)))
      scores <- vector("list", ng)
      for (g in seq_len(ng)) {
        scores[[g]] <- lmsCatLogLikGroup(
          filled$models[[g]], compressed[[g]]$data, ordered,
          threshold.info, delta, rule, reductions[[g]], g,
          compressed[[g]]$weights, compiled = TRUE, score = TRUE
        )
        if (!is.finite(scores[[g]]$logLik)) return(rep(NA_real_, length(p)))
      }
    }
    tau.score <- numeric(length(delta))
    for (g in seq_len(ng)) {
      indicators <- rownames(filled$models[[g]]$matrices$lambdaX)
      indicators <- indicators[indicators %in% ordered]
      offset <- 0L
      for (v in indicators) {
        spec <- threshold.info$specs[[paste(g, v, sep = "::")]]
        count <- length(spec$index)
        tau.score[spec$index] <- scores[[g]]$threshold.score[offset + seq_len(count)]
        offset <- offset + count
      }
    }

    if (analytical.predictor) {
      locations <- model$params$gradientStruct$locations
      location.score <- numeric(NROW(locations))
      for (g in seq_len(ng)) {
        use <- which(locations$group == g)
        loc <- locations[use, , drop = FALSE]
        location.score[use] <- lmsCatMatrixScoreCpp(
          filled$models[[g]]$matrices, rule$nodes, reductions[[g]]$explicit,
          reductions[[g]]$full.dimension, scores[[g]]$predictor.score,
          scores[[g]]$covariance.score %||% list(),
          scores[[g]]$continuous.index %||% integer(),
          reductions[[g]]$analytic,
          loc$block, loc$row, loc$col, loc$symmetric
        )
      }
      parameter.score <- drop(lmsCatParameterJacobian(theta, model) %*% location.score)
      theta.score <- parameter.score[active]
    } else {
      theta.score <- numeric(length(active))
      for (k in seq_along(active)) {
        index <- active[k]
        h <- ctl$gradient.eps * max(1, abs(theta[index]))
        theta.plus <- theta.minus <- theta
        theta.plus[index] <- theta.plus[index] + h
        theta.minus[index] <- theta.minus[index] - h
        plus <- tryCatch(fillModel(model, theta.plus, fillPhi = TRUE, method = "lms"),
                         error = function(e) NULL)
        minus <- tryCatch(fillModel(model, theta.minus, fillPhi = TRUE, method = "lms"),
                          error = function(e) NULL)
        if (is.null(plus) || is.null(minus)) return(rep(NA_real_, length(p)))
        for (g in seq_len(ng)) {
          momentsFor <- function(submodel) {
            M <- submodel$matrices
            indicators <- rownames(M$lambdaX)
            continuous <- setdiff(indicators, ordered)
            if (!length(continuous)) {
              latent <- lmsCatLatentMeansCpp(
                M, rule$nodes, reductions[[g]]$explicit,
                reductions[[g]]$full.dimension
              )
              return(list(
                mu = sweep(latent %*% t(M$lambdaX), 2L, c(M$tauX), "+"),
                covariance = NULL
              ))
            }
            moments <- lmsCatReducedMomentsCpp(
              M, rule$nodes, reductions[[g]]$explicit,
              reductions[[g]]$analytic, reductions[[g]]$full.dimension
            )
            mu <- sweep(moments$mean %*% t(M$lambdaX), 2L, c(M$tauX), "+")
            idx <- match(continuous, indicators)
            covariance <- lapply(seq_len(NROW(rule$nodes)), function(q) {
              full <- M$lambdaX %*% moments$cov[[q]] %*% t(M$lambdaX) +
                M$thetaDelta
              full[idx, idx, drop = FALSE]
            })
            list(mu = mu, covariance = covariance)
          }
          plus.moments <- momentsFor(plus$models[[g]])
          minus.moments <- momentsFor(minus$models[[g]])
          derivative <- (plus.moments$mu - minus.moments$mu) / (2 * h)
          theta.score[k] <- theta.score[k] +
            sum(scores[[g]]$predictor.score * derivative)
          if (length(scores[[g]]$covariance.score)) {
            for (q in seq_along(scores[[g]]$covariance.score)) {
              covariance.derivative <-
                (plus.moments$covariance[[q]] - minus.moments$covariance[[q]]) /
                (2 * h)
              theta.score[k] <- theta.score[k] +
                sum(scores[[g]]$covariance.score[[q]] * covariance.derivative)
            }
          }
        }
      }
    }
    delta.score <- drop(t(lmsCatThresholdJacobian(delta, threshold.info)) %*% tau.score)
    -c(theta.score, delta.score)
  } else NULL

  adaptiveGradientFunction <- if (use.gradient && analytical.predictor) function(p) {
    progress$score.evaluations <- progress$score.evaluations + 1L
    theta <- theta0
    theta[active] <- p[seq_along(active)]
    delta <- p[length(active) + seq_along(threshold.active)]
    theta[threshold.active] <- delta
    filled <- tryCatch(fillModel(model, theta, fillPhi = TRUE, method = "lms"),
                       error = function(e) NULL)
    if (is.null(filled)) return(rep(NA_real_, length(p)))
    locations <- model$params$gradientStruct$locations
    location.score <- numeric(NROW(locations))
    tau.score <- numeric(length(delta))
    for (g in seq_len(ng)) {
      use.locations <- which(locations$group == g)
      loc <- locations[use.locations, , drop = FALSE]
      indicators <- rownames(filled$models[[g]]$matrices$lambdaX)
      categorical <- indicators[indicators %in% ordered]
      wi <- compressed[[g]]$weights %||% rep(1, NROW(compressed[[g]]$data))
      for (i in seq_len(NROW(compressed[[g]]$data))) {
        score <- lmsCatLogLikGroup(
          filled$models[[g]], compressed[[g]]$data[i, , drop = FALSE], ordered,
          threshold.info, delta, adaptive.rules[[g]][[i]], reductions[[g]], g,
          wi[i], compiled = TRUE, score = TRUE
        )
        if (!is.finite(score$logLik)) return(rep(NA_real_, length(p)))
        location.score[use.locations] <- location.score[use.locations] +
          lmsCatMatrixScoreCpp(
            filled$models[[g]]$matrices, adaptive.rules[[g]][[i]]$nodes,
            reductions[[g]]$explicit, reductions[[g]]$full.dimension,
            score$predictor.score, score$covariance.score %||% list(),
            score$continuous.index %||% integer(), reductions[[g]]$analytic,
            loc$block, loc$row, loc$col, loc$symmetric
          )
        offset <- 0L
        for (v in categorical) {
          spec <- threshold.info$specs[[paste(g, v, sep = "::")]]
          count <- length(spec$index)
          tau.score[spec$index] <- tau.score[spec$index] +
            score$threshold.score[offset + seq_len(count)]
          offset <- offset + count
        }
      }
    }
    parameter.score <- drop(lmsCatParameterJacobian(theta, model) %*% location.score)
    theta.score <- parameter.score[active]
    delta.score <- drop(t(lmsCatThresholdJacobian(delta, threshold.info)) %*% tau.score)
    -c(theta.score, delta.score)
  } else NULL

  gradient.check <- NULL
  if (use.gradient && isTRUE(ctl$gradient.check)) {
    direction <- cos(seq_along(p0))
    direction <- direction / sqrt(sum(direction^2))
    h <- sqrt(ctl$gradient.eps)
    analytical <- tryCatch(sum(gradientFunction(p0) * direction),
                           error = function(e) NA_real_)
    numerical <- tryCatch(
      (objective(p0 + h * direction) - objective(p0 - h * direction)) / (2 * h),
      error = function(e) NA_real_
    )
    relative.error <- abs(analytical - numerical) / max(1, abs(numerical))
    passed <- is.finite(relative.error) && relative.error <= ctl$gradient.check.tol
    gradient.check <- list(
      type = "directional", analytical = analytical, numerical = numerical,
      relative.error = relative.error, tolerance = ctl$gradient.check.tol,
      passed = passed, predictor = if (analytical.predictor) "analytical" else "central"
    )
    if (!passed && analytical.predictor) {
      analytical.predictor <- FALSE
      fallback.value <- tryCatch(sum(gradientFunction(p0) * direction),
                                 error = function(e) NA_real_)
      fallback.error <- abs(fallback.value - numerical) / max(1, abs(numerical))
      fallback.passed <- is.finite(fallback.error) &&
        fallback.error <= ctl$gradient.check.tol
      gradient.check$fallback <- list(
        predictor = "central", analytical = fallback.value,
        numerical = numerical, relative.error = fallback.error,
        tolerance = ctl$gradient.check.tol, passed = fallback.passed
      )
      passed <- fallback.passed
      gradient.check$passed <- passed
      gradient.check$predictor.used <- if (passed) "central" else "none"
    } else {
      gradient.check$predictor.used <- if (passed) gradient.check$predictor else "none"
    }
    if (!passed) {
      mod_msg_warn("LMS-cat gradient validation failed; reverting to numerical optimizer derivatives.")
      gradientFunction <- NULL
      use.gradient <- FALSE
    }
  }
  if (!use.gradient || !analytical.predictor) adaptiveGradientFunction <- NULL

  constructAdaptiveRules <- function(p, update) {
    beginStage("adaptation", sprintf("recomputing adaptive quadrature (update %d)",
                                     update))
    theta.adapt <- theta0
    theta.adapt[active] <- p[seq_along(active)]
    delta.adapt <- p[length(active) + seq_along(threshold.active)]
    theta.adapt[threshold.active] <- delta.adapt
    filled.adapt <- fillModel(model, theta.adapt, fillPhi = TRUE, method = "lms")
    ans <- vector("list", ng)
    patterns.done <- 0L
    patterns.total <- sum(vapply(compressed, function(x) NROW(x$data), integer(1L)))
    for (g in seq_len(ng)) {
      ans[[g]] <- lapply(seq_len(NROW(compressed[[g]]$data)), function(i) {
        row <- compressed[[g]]$data[i, , drop = FALSE]
        log.kernel <- function(nodes) {
          vapply(seq_len(NROW(nodes)), function(q) {
            one <- list(nodes = matrix(nodes[q, ], nrow = 1L), weights = 1,
                        log.adjust = 0, method = "gh", signed = FALSE,
                        dimension = NCOL(nodes), nodes.total = 1L)
            lmsCatLogLikGroup(filled.adapt$models[[g]], row, ordered,
                              threshold.info, delta.adapt, one,
                              reductions[[g]], g, 1)
          }, numeric(1L))
        }
        previous <- if (is.null(adaptive.rules)) NULL else
          adaptive.rules[[g]][[i]]$adaptation$mode
        adapted <- lmsCatAdaptiveRule(rule, log.kernel, start = previous)
        patterns.done <<- patterns.done + 1L
        if (verbose && (patterns.done == 1L || patterns.done == patterns.total ||
                        patterns.done %% max(1L, ceiling(patterns.total / 20)) == 0L)) {
          clearConsoleLine()
          printf("lms-cat [adaptation] patterns=%d/%d elapsed=%.1fs",
                 patterns.done, patterns.total,
                 proc.time()[["elapsed"]] - progress$started)
        }
        adapted
      })
    }
    if (verbose) printf("\n")
    ans
  }

  algorithm <- toupper(ctl$algorithm %||% "NLMINB")
  mod_stopif(!algorithm %in% c("NLMINB", "EMA", "EM", "QN"),
             paste0("`integration.control$algorithm` must be `NLMINB`, ",
                    "`EMA`, `EM`, or `QN`."))
  adaptive.frequency <- as.integer(ctl$adaptive.frequency %||% 3L)
  mod_stopif(!is.finite(adaptive.frequency) || adaptive.frequency < 1L,
             "`adaptive.frequency` must be a positive integer.")
  if (rule$signed && algorithm %in% c("EMA", "EM")) {
    mod_msg_warn("EM/EMA requires nonnegative quadrature weights; using analytical QN for the signed sparse-grid rule.")
    algorithm <- "QN"
  }
  detailsAt <- function(p, score = TRUE) {
    if (score && !is.null(ema.cache$par) && identical(c(p), c(ema.cache$par)))
      return(ema.cache$state)
    if (score) progress$score.evaluations <- progress$score.evaluations + 1L
    else progress$light.evaluations <- progress$light.evaluations + 1L
    theta <- theta0; theta[active] <- p[seq_along(active)]
    delta <- p[length(active) + seq_along(threshold.active)]
    theta[threshold.active] <- delta
    filled <- tryCatch(fillModel(model, theta, fillPhi = TRUE, method = "lms"),
                       error = function(e) NULL)
    if (is.null(filled)) return(NULL)
    scores <- lapply(seq_len(ng), function(g) lmsCatLogLikGroup(
      filled$models[[g]], compressed[[g]]$data, ordered, threshold.info,
      delta, rule, reductions[[g]], g, compressed[[g]]$weights,
      compiled = TRUE, score = score, details = TRUE,
      posterior.details = score
    ))
    if (any(!vapply(scores, function(x) is.finite(x$logLik), logical(1L))))
      return(NULL)
    state <- list(logLik = sum(vapply(scores, `[[`, numeric(1L), "logLik")),
                  scores = scores, theta = theta, delta = delta, filled = filled)
    if (score) {
      ema.cache$par <- c(p)
      ema.cache$state <- state
    }
    state
  }
  expectedComplete <- function(candidate, reference) {
    sum(vapply(seq_len(ng), function(g) {
      w <- compressed[[g]]$weights %||% rep(1, NROW(reference$scores[[g]]$posterior))
      lmsCatFrozenQCpp(reference$scores[[g]]$posterior,
                       candidate$scores[[g]]$log.kernel, w)
    }, numeric(1L)))
  }
  emMap <- function(p) {
    reference <- detailsAt(p)
    if (is.null(reference)) return(list(par = p, state = NULL, accepted = FALSE))
    score <- -gradientFunction(p)
    if (any(!is.finite(score)) || !any(abs(score) > 0))
      return(list(par = p, state = reference, accepted = FALSE))
    direction <- score / max(abs(score)) * ctl$ema.max.step
    q0 <- expectedComplete(reference, reference)
    alpha <- 1
    while (alpha >= ctl$ema.min.step) {
      trial <- pmin(upper, pmax(lower, p + alpha * direction))
      candidate <- detailsAt(trial, score = FALSE)
      if (!is.null(candidate) && candidate$logLik >= reference$logLik - 1e-10 &&
          expectedComplete(candidate, reference) >= q0 - 1e-10)
        return(list(par = trial, state = candidate, accepted = TRUE))
      alpha <- alpha / 2
    }
    list(par = p, state = reference, accepted = FALSE)
  }

  control <- list(iter.max = as.integer(max.iter %||% 500L),
                  eval.max = max(200L, 2L * as.integer(max.iter %||% 500L)),
                  abs.tol = convergence.abs %||% 1e-5,
                  trace = 0)
  beginStage("optimization", "optimizing with the base integration rule")
  adaptation.updates <- 0L
  if (algorithm == "NLMINB") {
    p.current <- p0
    old.logLik <- -objective(p.current)
    converged.nlminb <- FALSE
    iterations.nlminb <- 0L
    last.step <- NULL
    for (iteration in seq_len(as.integer(max.iter %||% 500L))) {
      iterations.nlminb <- iteration
      adapted <- isTRUE(ctl$adapt) && identical(rule$method, "gh") &&
        iteration %% adaptive.frequency == 0L
      if (adapted) {
        adaptation.updates <- adaptation.updates + 1L
        adaptive.rules <- constructAdaptiveRules(p.current, adaptation.updates)
        # A quadrature refresh changes the numerical approximation.  Establish
        # a fresh baseline before assessing the improvement from this M step.
        old.logLik <- -objective(p.current)
        beginStage("optimization", "continuing with refreshed adaptive quadrature")
      }
      step.control <- control
      step.control$iter.max <- 1L
      step.control$eval.max <- max(100L, 2L * length(p.current) + 20L)
      step <- stats::nlminb(
        p.current, objective,
        gradient = if (is.null(adaptive.rules)) gradientFunction else
          adaptiveGradientFunction,
        lower = lower, upper = upper, control = step.control
      )
      new.logLik <- -step$objective
      deltaLL <- new.logLik - old.logLik
      relDeltaLL <- deltaLL / max(1, abs(old.logLik))
      updateStatusLog(iteration, if (adapted) "QN-A" else "QN", new.logLik,
                      deltaLL, relDeltaLL, verbose)
      p.current <- step$par
      last.step <- step
      converged.nlminb <- is.finite(new.logLik) &&
        (abs(deltaLL) < (convergence.abs %||% 1e-5) ||
           abs(relDeltaLL) < ctl$rel.tol)
      if (isTRUE(ctl$adapt) && identical(rule$method, "gh") &&
          adaptation.updates == 0L) converged.nlminb <- FALSE
      old.logLik <- new.logLik
      if (converged.nlminb) break
    }
    if (verbose) cat("\n")
    est <- list(par = p.current, objective = -old.logLik,
                iterations = iterations.nlminb,
                convergence = if (converged.nlminb) 0L else 1L,
                message = if (converged.nlminb) "converged" else
                  (last.step$message %||% "iteration limit reached"))
  } else if (algorithm == "QN" || is.null(gradientFunction)) {
    est <- stats::nlminb(p0, objective, gradient = gradientFunction,
                         lower = lower, upper = upper, control = control)
  } else {
    current <- detailsAt(p0)
    p.current <- p0
    old.logLik <- -Inf
    converged.em <- FALSE
    iterations.em <- 0L
    for (iteration in seq_len(as.integer(max.iter %||% 500L))) {
      iterations.em <- iteration
      first <- emMap(p.current)
      proposal <- first
      mode <- "EM"
      if (algorithm == "EMA" && first$accepted) {
        second <- emMap(first$par)
        proposal <- second
        if (second$accepted) {
          r <- first$par - p.current
          v <- (second$par - first$par) - r
          denom <- sum(v^2)
          if (is.finite(denom) && denom > 0) {
            alpha <- -sqrt(sum(r^2) / denom)
            alpha <- max(-ctl$ema.max.extrapolation, min(-1, alpha))
            square <- pmin(upper, pmax(lower,
              p.current - 2 * alpha * r + alpha^2 * v))
            accelerated <- emMap(square)
            if (!is.null(accelerated$state) && accelerated$state$logLik >=
                second$state$logLik - 1e-10) {
              proposal <- accelerated
              mode <- "QN"
            }
          }
        }
      }
      if (is.null(proposal$state)) break
      new.logLik <- proposal$state$logLik
      deltaLL <- new.logLik - old.logLik
      relDeltaLL <- if (is.finite(old.logLik)) deltaLL / abs(old.logLik) else Inf
      updateStatusLog(iteration, mode, new.logLik, deltaLL, relDeltaLL, verbose)
      p.current <- proposal$par
      current <- proposal$state
      converged.em <- is.finite(old.logLik) &&
        (abs(deltaLL) < (convergence.abs %||% 1e-5) ||
           abs(relDeltaLL) < ctl$rel.tol)
      old.logLik <- new.logLik
      if (converged.em || !proposal$accepted) break
    }
    if (verbose) cat("\n")
    est <- list(par = p.current, objective = -current$logLik,
                iterations = iterations.em,
                convergence = if (converged.em) 0L else 1L,
                message = if (converged.em) "converged" else "stopped")
  }
  qmc.diagnostics <- NULL
  if (identical(rule$method, "qmc")) {
    repeat {
      qmc.diagnostics <- qmcDiagnostics(est$par, rule)
      enough <- !is.finite(qmc.diagnostics$standard.error) ||
        isTRUE(qmc.diagnostics$reliable) ||
        rule$points >= as.integer(ctl$qmc.max.points) || !isTRUE(ctl$qmc.adapt)
      if (enough) break
      next.points <- min(2L * rule$points, as.integer(ctl$qmc.max.points))
      if (next.points <= rule$points) break
      beginStage("QMC refinement",
                 sprintf("increasing randomized Sobol points from %d to %d (logLik SE %.4g)",
                         rule$points, next.points, qmc.diagnostics$standard.error))
      metadata <- rule[c("requested", "selection.reason")]
      rule <- lmsCatQmcRule(next.points, dims[1L], ctl$qmc.replicates, ctl$seed)
      rule[names(metadata)] <- metadata
      est <- stats::nlminb(est$par, objective, gradient = gradientFunction,
                           lower = lower, upper = upper,
                           control = control)
    }
  }
  theta <- theta0; theta[active] <- est$par[seq_along(active)]
  delta <- est$par[length(active) + seq_along(threshold.active)]
  theta[threshold.active] <- delta
  converged <- est$convergence == 0L && is.finite(est$objective)

  out <- finalizeModelEstimatesDA(
    model = model, theta = theta, method = "lms", data = data,
    logLik = -est$objective, iterations = est$iterations,
    converged = converged, optimizer = "nlminb", calc.se = FALSE,
    includeStartModel = TRUE, startModel = model
  )
  # Continuous-response nuisance rows do not belong to the ordinal model.
  pt <- as.data.frame(out$parTable)
  nuisance.rows <- (pt$lhs %in% ordered & pt$op == "~1") |
    (pt$op == "~~" & (pt$lhs %in% ordered | pt$rhs %in% ordered))
  pt <- pt[!nuisance.rows, , drop = FALSE]

  V <- NULL
  if (calc.se) {
    beginStage("standard errors", "calculating the observed numerical Hessian")
    hessian.gradient <- if (!use.gradient) NULL else if (is.null(adaptive.rules))
      gradientFunction else adaptiveGradientFunction
    H <- tryCatch(stats::optimHess(est$par, objective, gr = hessian.gradient),
                  error = function(e) NULL)
    V <- if (!is.null(H)) tryCatch(solve((H + t(H)) / 2), error = function(e) NULL) else NULL
  }
  se.theta <- rep(NA_real_, length(theta))
  if (!is.null(V)) se.theta[active] <- sqrt(pmax(diag(V)[seq_along(active)], 0))
  if (!is.null(V)) se.theta[threshold.active] <- sqrt(pmax(
    diag(V)[length(active) + seq_along(threshold.active)], 0
  ))
  free.labels <- names(out$coefs.free)
  if (length(free.labels) == length(theta)) {
    idx <- match(getParTableLabels(pt), free.labels)
    use <- !is.na(idx)
    pt$std.error[use] <- se.theta[idx[use]]
  }

  threshold.rows <- NULL; threshold.values <- numeric()
  for (spec in threshold.info$specs) {
    tau <- lmsCatDeltaToThresholds(delta, spec)
    threshold.values <- c(threshold.values, stats::setNames(tau, spec$labels))
    threshold.rows <- rbind(threshold.rows, data.frame(
      lhs = spec$variable, op = "|", rhs = paste0("t", seq_along(tau)),
      label = "", group = spec$group, est = tau, std.error = NA_real_
    ))
  }
  if (!is.null(V) && length(delta)) {
    vd <- V[length(active) + seq_along(delta), length(active) + seq_along(delta), drop = FALSE]
    J <- lmsCatThresholdJacobian(delta, threshold.info)
    st <- sqrt(pmax(diag(J %*% vd %*% t(J)), 0))
    threshold.rows$std.error <- st
  }
  missing.cols <- setdiff(names(pt), names(threshold.rows))
  threshold.rows[missing.cols] <- NA
  threshold.rows <- threshold.rows[names(pt)]
  pt <- addZStatsParTable(rbind(pt, threshold.rows))

  out$method.requested <- requested.method
  out$method <- "lms-jv"
  out$link <- "probit"
  out$theta.base <- theta
  out$theta <- theta
  out$coefs.all <- c(out$coefs.all, threshold.values)
  out$parTable <- modsemParTable(sortParTableDA(pt, model = out$model))
  out$logLik <- -est$objective
  out$type.se <- if (is.null(V)) "none" else "observed-hessian"
  out$vcov.optim <- V
  vcov.delta <- if (!is.null(V) && length(delta))
    V[length(active) + seq_along(delta),
      length(active) + seq_along(delta), drop = FALSE] else NULL
  out$model <- lmsCatAttachThresholdMatrices(out$model, threshold.info, delta,
                                             vcov.delta = vcov.delta)
  integration.method <- if (is.null(adaptive.rules)) rule$method else "aghq"
  out$integration <- list(method = integration.method, dimension = dims[1L],
                          nodes.total = rule$nodes.total, control = ctl,
                          adaptation.updates = adaptation.updates,
                          adaptive.frequency = adaptive.frequency,
                          requested = integration,
                          selection.reason = rule$selection.reason,
                          explicit = reductions[[1L]]$explicit,
                          analytic = reductions[[1L]]$analytic,
                          source.dimension = reductions[[1L]]$full.dimension)
  out$integration$response.patterns <- lapply(compressed, function(x) {
    list(compressed = x$compressed, original = x$n.original %||% NROW(x$data),
         patterns = x$n.patterns %||% NROW(x$data))
  })
  diagnostics <- list()
  if (identical(rule$method, "qmc")) {
    qmc.diagnostics <- qmcDiagnostics(est$par, rule)
    diagnostics$qmc <- qmc.diagnostics
  }
  if (identical(rule$method, "sparse") && rule$level > 1L) {
    lower.rule <- lmsCatSparseRule(rule$level - 1L, dims[1L], ctl$max.nodes)
    lower.logLik <- baseLogLik(est$par, lower.rule)
    diagnostics$sparse <- list(
      family = rule$family, level = rule$level,
      comparison.level = rule$level - 1L,
      logLik = out$logLik, comparison.logLik = lower.logLik,
      absolute.change = abs(out$logLik - lower.logLik)
    )
  }
  if (identical(rule$method, "gh") && ctl$order > 3L) {
    comparison.order <- max(3L, as.integer(ctl$order) - 2L)
    comparison.rule <- lmsCatTensorRule(comparison.order, dims[1L])
    current.base.logLik <- baseLogLik(est$par, rule)
    comparison.logLik <- baseLogLik(est$par, comparison.rule)
    diagnostics$gauss.hermite <- list(
      order = as.integer(ctl$order), comparison.order = comparison.order,
      logLik = current.base.logLik, comparison.logLik = comparison.logLik,
      absolute.change = abs(current.base.logLik - comparison.logLik)
    )
  }
  out$integration$diagnostics <- diagnostics
  out$info.quad <- list(dim = dims[1L], nodes.dim = ctl$order,
                        nodes.total = rule$nodes.total, method = integration.method)
  out$args$ordered <- ordered
  out$args$integration <- integration
  out$args$integration.control <- ctl
  out$optimization.progress <- list(
    evaluations = progress$evaluations,
    score.evaluations = progress$score.evaluations,
    lightweight.evaluations = progress$light.evaluations,
    adaptation.updates = adaptation.updates,
    best.logLik = progress$best,
    elapsed = proc.time()[["elapsed"]] - progress$started
  )
  out$optimization.gradient <- list(
    used = use.gradient, method = if (use.gradient)
      if (all.categorical) paste(
        "exact ordinal score with",
        if (analytical.predictor) "analytical" else "central",
        "predictor sensitivities"
      ) else paste(
        "exact mixed-response score with",
        if (analytical.predictor) "analytical" else "central",
        "mean and covariance sensitivities"
      ) else "optimizer numerical",
    check = gradient.check,
    hessian = if (calc.se && use.gradient && (!is.null(gradientFunction) ||
                                                !is.null(adaptiveGradientFunction)))
      "finite difference of score" else if (calc.se) "finite difference of objective" else "not calculated"
  )
  out$optimization.algorithm <- algorithm
  if (verbose) {
    clearConsoleLine()
    printf("lms-jv: %s after %d iterations (logLik=%.6f, elapsed=%.1fs)\n",
           if (converged) "converged" else "stopped without convergence",
           est$iterations, out$logLik, out$optimization.progress$elapsed)
  }
  class(out) <- c("modsem_da", "modsem")
  out
}
