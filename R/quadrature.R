# Numerical integration rules for the LMS approach.
#
# A *spec* describes which rule to build; it never holds nodes. Rules are
# materialised on demand by `buildQuadRule()`, so a k-dimensional product grid
# only ever exists when it is about to be evaluated. This matters because a
# product rule costs `nodes ^ k` rows, and building one eagerly at model
# construction (as earlier versions did) wastes that memory even when the rule
# is immediately replaced by an adaptive one.
#
# Two axes are orthogonal:
#
#   integration  "gh"   Gauss-Hermite product rule
#                "rect" equispaced rectangular rule, as in Mplus INTEGRATION=STANDARD
#                "mc"   fixed Monte-Carlo draws
#
#   adaptive     "none"  one fixed rule, shared by every observation
#                "quasi" one rule per group, adapted to the aggregate posterior
#                "full"  one rule per observation, from its mode and Hessian
#
# `adaptive = "full"` is not implemented here: it builds the corresponding
# non-adaptive rule and hands it to the backend, which transforms it per
# observation. That keeps the two axes independent.

INTEGRATION_TYPES <- c("gh", "rect", "mc")
ADAPTIVE_TYPES    <- c("none", "quasi", "full")

INTEGRATION_LABELS <- c(gh = "Gauss-Hermite", rect = "Rectangular",
                        mc = "Monte-Carlo")
ADAPTIVE_LABELS <- c(none = "None (fixed rule)",
                     quasi = "Quasi (one rule per group)",
                     full = "Full (one rule per observation)")

# Default nodes keyed on (integration, adaptive). For `gh` and `rect` this is
# the number of nodes per latent dimension; for `mc` it is the total number of
# draws, whose cost does not grow with dimension.
#
# The shared-rule defaults match Mplus: 15 points per dimension for both
# product rules, 500 draws in total for Monte-Carlo. A per-observation rule is
# centred on each observation's own posterior and needs far fewer points.
QUAD_NODE_DEFAULTS <- list(
  gh   = c(none = 15,  quasi = 15,  full = 5),
  rect = c(none = 15,  quasi = 15,  full = 5),
  mc   = c(none = 500, quasi = 500, full = 250)
)

# Which (backend, integration, adaptive) combinations are supported. The legacy
# LMS backend has no per-observation machinery and its E-step is built around a
# product rule, so it offers the two common rules only.
QUAD_AVAILABLE <- list(
  legacy = list(gh = c("none", "quasi"), rect = c("none", "quasi"), mc = character()),
  graph  = list(gh = ADAPTIVE_TYPES, rect = ADAPTIVE_TYPES, mc = c("none", "full"))
)


quadSpec <- function(integration = "gh",
                     adaptive = "quasi",
                     nodes = NULL,
                     k = 0L,
                     quad.range = Inf,
                     rect.range = 5,
                     adaptive.frequency = 3L,
                     adaptive.quad.tol = 1e-12,
                     seed = 1234L,
                     max.nodes = 1e6) {
  integration <- match.arg(integration, INTEGRATION_TYPES)
  adaptive    <- match.arg(adaptive, ADAPTIVE_TYPES)

  if (is.null(nodes)) nodes <- QUAD_NODE_DEFAULTS[[integration]][[adaptive]]

  if (quad.range < 0) {
    mod_msg_warn("`quad.range` should be positive, using `-quad.range` instead!\n")
    quad.range <- -quad.range
  }

  structure(list(
    integration = integration,
    adaptive    = adaptive,
    m           = as.numeric(nodes),
    k           = as.integer(k),
    quad.range  = quad.range,
    rect.range  = rect.range,
    a           = -quad.range,
    b           = quad.range,
    adaptive.frequency = adaptive.frequency,
    adaptive.quad.tol  = adaptive.quad.tol,
    seed        = seed,
    max.nodes   = max.nodes
  ), class = "modsem_quad_spec")
}


# A model with no dimensions to integrate (a purely linear LMS model) has a
# single degenerate node, so there is nothing for any adaptation to do.
quadIsDegenerate <- function(spec)
  !isTRUE(spec$k >= 1L) || !isTRUE(spec$m >= 1)


isAdaptiveQuad <- function(spec)
  !identical(spec$adaptive, "none") && !quadIsDegenerate(spec)


setAdaptiveQuad <- function(spec, adaptive) {
  spec$adaptive <- match.arg(adaptive, ADAPTIVE_TYPES)
  spec
}


# `nodes` counts draws for Monte-Carlo and nodes per dimension otherwise.
quadTotalNodes <- function(spec) {
  if (identical(spec$integration, "mc")) as.numeric(spec$m)
  else as.numeric(spec$m) ^ as.numeric(spec$k)
}


checkQuadBudget <- function(spec, total = quadTotalNodes(spec), what = NULL) {
  if (is.finite(total) && total <= spec$max.nodes) return(invisible(TRUE))
  label <- if (is.null(what)) "requested rule" else what
  mod_msg_stop(sprintf(paste0(
    "The %s needs %.0f integration points in %d latent dimensions, which ",
    "exceeds the limit of %.0f.\n",
    "Reduce `nodes`, use `adaptive = \"quasi\"` to prune the rule, ",
    "`adaptive = \"full\"` for a small per-observation rule, or ",
    "`integration = \"mc\"` whose cost does not grow with dimension."),
    label, total, spec$k, spec$max.nodes))
}


checkQuadAvailable <- function(spec, backend = "legacy") {
  available <- QUAD_AVAILABLE[[backend]][[spec$integration]]
  if (spec$adaptive %in% available) return(invisible(TRUE))

  alternatives <- if (!length(available))
    sprintf("`integration = \"%s\"` is not available for the %s LMS backend",
            spec$integration, backend)
  else sprintf("the %s LMS backend supports `adaptive = %s` with `integration = \"%s\"`",
               backend, paste0("\"", available, "\"", collapse = " or "),
               spec$integration)

  hint <- if (identical(backend, "legacy") &&
              spec$adaptive %in% QUAD_AVAILABLE$graph[[spec$integration]])
    "\nUse `lms.backend = \"graph\"` for this combination." else ""

  mod_msg_stop(sprintf(
    "`integration = \"%s\"` with `adaptive = \"%s\"` is not supported: %s.%s",
    spec$integration, spec$adaptive, alternatives, hint))
}


# One-dimensional rules ------------------------------------------------------
#
# Every product rule is assembled from per-dimension node/weight vectors, which
# may differ in length once an adaptive rule has pruned them. Weights are
# normalised per dimension so their tensor product integrates a constant to one.

gaussHermite1d <- function(m) {
  if (m <= 1L) return(list(n = 0, w = 1))
  single <- fastGHQuad::gaussHermiteData(m)
  list(n = sqrt(2) * single$x, w = single$w * pi ^ (-1 / 2))
}


# Trapezoid weights on an equispaced grid: the two endpoints each cover half an
# interval. Mplus calls this "rectangular (trapezoid)" integration. Normalising
# by the sum makes the rule integrate a constant exactly.
#
# On an unbounded, Gaussian-decaying integrand this rule is not the O(h^2)
# scheme its name suggests: Poisson summation gives an aliasing error of
# 2*exp(-2*pi^2/h^2), so it converges exponentially in 1/h. Accuracy is then
# governed by whichever is larger, that term or the truncation error
# 2*(1 - pnorm(range)) -- which is why `rect.range` matters as much as `nodes`.
rectangularWeights <- function(nodes) {
  weights <- stats::dnorm(nodes)
  ends <- c(1L, length(weights))
  if (length(weights) > 1L) weights[ends] <- 0.5 * weights[ends]
  weights / sum(weights)
}


rectangular1d <- function(m, range = 5) {
  if (m <= 1L) return(list(n = 0, w = 1))
  nodes <- seq(-range, range, length.out = m)
  list(n = nodes, w = rectangularWeights(nodes))
}


quadNodes1d <- function(spec, m = spec$m) {
  switch(spec$integration,
         gh   = gaussHermite1d(m),
         rect = rectangular1d(m, range = spec$rect.range),
         mod_msg_stop(sprintf("`%s` has no one-dimensional rule.", spec$integration)))
}


# `expand.grid()` varies the first column fastest, so the matching weight
# product is a kronecker over the reversed list.
tensorProductRule <- function(nodes, weights) {
  n <- as.matrix(do.call(expand.grid, nodes))
  dimnames(n) <- NULL
  # `kronecker()` returns a 1-d array, which does not recycle against a matrix.
  list(n = n, w = as.numeric(Reduce(kronecker, rev(weights))))
}


withSeed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    state <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", state, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}


# Rule builders --------------------------------------------------------------
#
# Each builder returns `list(n, w)`; `buildQuadRule()` adds the bookkeeping.
# The registry below is the seam for changing an adaptation strategy: swapping
# the dense-candidate approach back in for `gh$quasi`, say, is a one-line change
# here and touches nothing else.

quadRuleFixed <- function(spec, ...) {
  if (spec$k < 1L) return(list(n = matrix(0), w = 1))
  checkQuadBudget(spec, what = sprintf("fixed %s rule", spec$integration))
  single <- quadNodes1d(spec)
  tensorProductRule(replicate(spec$k, single$n, simplify = FALSE),
                    replicate(spec$k, single$w, simplify = FALSE))
}


quadRuleMonteCarlo <- function(spec, ...) {
  samples <- as.integer(spec$m)
  mod_stopif(samples < 2L, "Monte-Carlo integration requires at least two draws.")
  checkQuadBudget(spec, total = samples, what = "Monte-Carlo rule")
  if (spec$k < 1L) return(list(n = matrix(0), w = 1))
  nodes <- withSeed(spec$seed,
                    matrix(stats::rnorm(samples * spec$k), samples, spec$k))
  list(n = nodes, w = rep(1 / samples, samples))
}


# Gauss-Hermite, quasi-adaptive: prune one dimension at a time on a 1-D grid
# and only then form the product of the surviving nodes. A dense k-dimensional
# candidate grid is never built, which is what lets this scale past three or
# four latent dimensions.
quadRuleGaussHermiteQuasi <- function(spec, density, previous = NULL, ...) {
  if (spec$k < 1L) return(list(n = matrix(0), w = 1))
  m.ceil <- previous$m.ceil
  if (is.null(m.ceil) || !all(is.finite(m.ceil)))
    m.ceil <- if (spec$k > 1L) spec$m else
      round(estMForNodesInRange(spec$m, a = -5, b = 5))

  rule <- adaptiveGaussQuadrature(
    fun = density, a = spec$a, b = spec$b, m = spec$m, k = spec$k,
    m.ceil = m.ceil, tol = spec$adaptive.quad.tol, ...
  )
  checkQuadBudget(spec, total = NROW(rule$n), what = "quasi-adaptive rule")
  rule
}


# Rectangular, quasi-adaptive: a rectangular rule's one real lever is its
# range, so rather than dropping nodes we shrink the interval per dimension to
# where the aggregate posterior actually has mass and respace the same number
# of points inside it. Pruning would waste the lever; see `?modsem_da`.
quadRuleRectangularQuasi <- function(spec, density, previous = NULL, ...) {
  if (spec$k < 1L) return(list(n = matrix(0), w = 1))
  checkQuadBudget(spec, what = "quasi-adaptive rectangular rule")

  probe.m <- max(4L * as.integer(spec$m), 64L)
  nodes <- weights <- vector("list", spec$k)

  for (j in seq_len(spec$k)) {
    probe <- rectangular1d(probe.m, range = spec$rect.range)
    grid <- matrix(0, probe.m, spec$k)
    grid[, j] <- probe$n

    trimmed <- tryCatch({
      pruned <- pruneQuadratureNodes(
        quadw = probe$w, quadn = grid, quadf = density(grid, ...),
        a = spec$a, b = spec$b, tol = spec$adaptive.quad.tol
      )
      range(pruned$quadn[, j])
    }, error = function(e) c(-spec$rect.range, spec$rect.range))

    lower <- min(trimmed)
    upper <- max(trimmed)
    if (!is.finite(lower) || !is.finite(upper) || upper - lower <= 0) {
      lower <- -spec$rect.range
      upper <- spec$rect.range
    }

    axis <- if (spec$m <= 1L) 0 else seq(lower, upper, length.out = spec$m)
    nodes[[j]] <- axis
    weights[[j]] <- if (length(axis) > 1L) rectangularWeights(axis) else 1
  }

  out <- tensorProductRule(nodes, weights)
  out$range <- vapply(nodes, range, numeric(2L))
  out
}


QUAD_BUILDERS <- list(
  gh   = list(none = quadRuleFixed,      quasi = quadRuleGaussHermiteQuasi),
  rect = list(none = quadRuleFixed,      quasi = quadRuleRectangularQuasi),
  mc   = list(none = quadRuleMonteCarlo, quasi = NULL)
)


# Build the rule a spec describes.
#
# `density(nodes, ...)` must return an N by Q matrix of densities and is only
# needed for `adaptive = "quasi"`. `adaptive = "full"` returns the underlying
# common rule; transforming it per observation is the backend's job, since only
# the backend can find each observation's mode and Hessian.
buildQuadRule <- function(spec, density = NULL, previous = NULL, ...) {
  adaptive <- if (identical(spec$adaptive, "full") ||
                  quadIsDegenerate(spec)) "none" else spec$adaptive
  builder <- QUAD_BUILDERS[[spec$integration]][[adaptive]]
  mod_stopif(is.null(builder), sprintf(
    "`integration = \"%s\"` has no `adaptive = \"%s\"` rule.",
    spec$integration, spec$adaptive))
  mod_stopif(identical(adaptive, "quasi") && is.null(density),
             "A quasi-adaptive rule needs a density function.")

  rule <- builder(spec, density = density, previous = previous, ...)

  rule$Q <- NROW(rule$n)
  rule$k <- spec$k
  rule$m <- spec$m
  rule$integration <- spec$integration
  rule$adaptive <- spec$adaptive
  rule$adaptive.frequency <- spec$adaptive.frequency
  rule$quad.range <- spec$quad.range
  rule$spec <- spec
  # A rule shared by every observation permits the sufficient-statistic M-step;
  # a per-observation rule does not. `full` rules are marked by the backend
  # once it has packed them.
  rule$common <- TRUE
  rule$packed <- FALSE
  rule
}


# Legacy quasi-adaptive Gauss-Hermite ----------------------------------------
#
# Prune a denser one-dimensional rule until the number of surviving nodes lands
# near the requested order, then take the tensor product of the pruned axes.

adaptiveGaussQuadrature <- function(fun,
                                    a = -7,
                                    b = 7,
                                    m = 32,
                                    m.ceil = m + m / 2,
                                    k = 1,
                                    iter.max = 10,
                                    node.max = 500,
                                    tol = 1e-12,
                                    mdiff.tol = 2,
                                    secondary.pruning = TRUE,
                                    ...) {
  if (k == 0 || m == 0)
    return(list(n = matrix(0), w = 1, f = NA, m = 1, k = 1))

  mod_stopif(tol >= 1 || tol < 0,
         "`adaptive.quad.tol` must be in the boundary `[0, 1)`")

  if (is.null(m.ceil) || all(is.na(m.ceil)))
    m.ceil <- m + m / 2

  if (k <= 1) {
    out <- adaptiveGaussQuadratureK(
      fun = fun, a = a, b = b, m = m, m.ceil = m.ceil,
      k = 1, K = k, iter = 1, iter.max = iter.max,
      node.max = node.max, tol = tol, mdiff.tol = mdiff.tol, ...
    )

    return(out)
  }

  a      <- rep(a, length.out = k)
  b      <- rep(b, length.out = k)
  m.ceil <- rep(m.ceil, length.out = k)

  NODES     <- vector("list", k)
  WEIGHTS   <- vector("list", k)
  new.ceils <- integer(k)

  for (i in seq_len(k)) {
    QUAD <- adaptiveGaussQuadratureK(
      fun = fun, a = a[i], b = b[i], m = m, m.ceil = m.ceil[i],
      k = i, K = k, iter = 1, iter.max = iter.max,
      node.max = node.max, tol = tol, mdiff.tol = mdiff.tol, ...
    )

    NODES[[i]]     <- QUAD$n[, i]
    WEIGHTS[[i]]   <- QUAD$w
    new.ceils[[i]] <- QUAD$m.ceil
  }

  product <- tensorProductRule(NODES, WEIGHTS)
  quadn <- product$n
  quadw <- product$w
  quadf <- fun(quadn, ...)
  weighted <- sweep(quadf, 2, quadw, "*")
  integral <- quadratureLogLik(weighted)
  error <- 0
  error.abs <- 0
  error.rel <- 0

  if (secondary.pruning) {
    pruned <- pruneQuadratureNodes(
      quadw = quadw, quadn = quadn, quadf = quadf,
      a = a, b = b, tol = tol
    )

    quadw     <- pruned$quadw
    quadn     <- pruned$quadn
    quadf     <- pruned$quadf
    integral  <- pruned$I.cur
    error     <- pruned$I.err
    error.abs <- pruned$I.err.abs
    error.rel <- pruned$I.err.rel
  }

  list(
    n = quadn,
    w = quadw,
    F = quadf,
    m.ceil = new.ceils,
    error = error,
    error.abs = error.abs,
    error.rel = error.rel,
    integral = integral
  )
}


adaptiveGaussQuadratureK <- function(fun,
                                     a = -7,
                                     b = 7,
                                     m = 32,
                                     m.ceil = m + m / 2,
                                     k = 1,
                                     K = 1,
                                     iter = 1,
                                     iter.max = 10,
                                     node.max = 500,
                                     tol = 1e-12,
                                     mdiff.tol = 2,
                                     ...) {
  if (k == 0 || m == 0)
    return(list(n = matrix(0), w = 1, f = NA, m = 1, k = 1))

  if (is.null(m.ceil) || is.na(m.ceil) || m.ceil <= 0)
    m.ceil <- round(estMForNodesInRange(m, a = -5, b = 5))

  single <- gaussHermite1d(m.ceil)

  if (K > 1) { # more computationally efficient/scalable
    quadn <- matrix(0, nrow = m.ceil, ncol = K)
    quadn[,k] <- single$n
  } else quadn <- matrix(single$n, ncol = 1L)

  quadf <- fun(quadn, ...)
  quadw <- single$w

  pruned <- pruneQuadratureNodes(
    quadw = quadw, quadn = quadn, quadf = quadf,
    a = a, b = b, tol = tol
  )

  quadw     <- pruned$quadw
  quadn     <- pruned$quadn
  quadf     <- pruned$quadf
  I.err     <- pruned$I.err
  I.full    <- pruned$I.full
  I.cur     <- pruned$I.cur
  I.err.abs <- pruned$I.err.abs
  I.err.rel <- pruned$I.err.rel

  lower <- min(quadn)
  upper <- max(quadn)

  diff.m <- NROW(quadn) - m

  # Calculate next ceiling
  if (abs(diff.m) > 2 * mdiff.tol) {
    # If the difference is large, use a rough estimate
    new.ceil <- round(estMForNodesInRange(k = m, a = lower, b = upper))

    # else do a step change
  } else new.ceil <- m.ceil - diff.m

  OKNextCeil <- new.ceil <= node.max && new.ceil != m.ceil
  converged <- abs(diff.m) <= mdiff.tol
  OKIter <- iter < iter.max

  if (!converged && OKIter && OKNextCeil) {

    return(adaptiveGaussQuadratureK(
      fun = fun,
      a = a,
      b = b,
      k = k,
      K = K,
      m = m,
      m.ceil = new.ceil,
      iter = iter + 1,
      tol = tol,
      iter.max = iter.max,
      node.max = node.max,
      mdiff.tol = mdiff.tol,
      ...
    ))
  }

  mod_warnif_immediate(
    iter >= iter.max,
    paste0("Max iterations reached fitting quasi-adaptive quadrature...\n",
           sprintf("Iter: %d, total: %d, target: %d, kept: %d, discarded: %d",
                   iter, m.ceil, m, NROW(quadn), m.ceil - NROW(quadn))),
    .newline = TRUE
  )

  list(
    n = quadn,
    w = quadw,
    F = quadf,
    k = k,
    m = nrow(quadn) ^ (1 / k),

    m.ceil   = m.ceil,
    iter     = iter,
    error    = I.err,
    error.abs = I.err.abs,
    error.rel = I.err.rel,
    integral = I.cur
  )
}


estGHNodesInRange <- function(m, a, b, scale = TRUE) {
  mod_stopif(!is.numeric(m) || length(m) != 1 || m <= 0,
         "'m' must be a single positive number")
  mod_stopif(!is.numeric(a) || length(a) != 1 || !is.numeric(b) || length(b) != 1,
         "'a' and 'b' must be numeric scalars")

  if (scale) {
    a <- a / sqrt(2)
    b <- b / sqrt(2)
  }
  # Ensure interval is ordered
  lo <- min(a, b)
  hi <- max(a, b)

  # Scaling factor from the semicircle law support
  scale <- sqrt(2 * m)

  # Normalize interval endpoints to t = x / sqrt(2m)
  t_lo <- lo / scale
  t_hi <- hi / scale

  # Clamp to the theoretical support of the density, [-1, 1]
  t_lo <- pmax(pmin(t_lo,  1), -1)
  t_hi <- pmax(pmin(t_hi,  1), -1)

  # Define the cumulative function F(t) = t*sqrt(1 - t^2) + arcsin(t)
  F <- \(t) t * sqrt(pmax(0, 1 - t^2)) + asin(t)

  (m / pi) * (F(t_hi) - F(t_lo))
}


estMForNodesInRange <- function(k, a, b,
                                lower = 1,
                                upper = NULL,
                                tol = 1e-6,
                                maxiter = 100,
                                scale = TRUE) {
  mod_stopif(!is.numeric(k) || length(k) != 1 || k <= 0,
         "'k' must be a single positive number")
  mod_stopif(!is.numeric(a) || length(a) != 1 || !is.numeric(b) || length(b) != 1,
         "'a' and 'b' must be numeric scalars")

  if (scale) {
    a <- a / sqrt(2)
    b <- b / sqrt(2)
  }

  # Ensure interval
  lo <- min(a, b)
  hi <- max(a, b)

  # Define f(m) = estimate_gh_nodes(m) - k
  f <- function(m) estGHNodesInRange(m, lo, hi, scale = FALSE) - k

  # Determine upper bound if not provided by doubling
  if (is.null(upper)) {
    upper <- lower * 2
    iter <- 0
    while (f(upper) < 0 && iter < maxiter) {
      upper <- upper * 2
      iter <- iter + 1
    }
    mod_stopif(f(upper) < 0, "Unable to bracket root: increase maxiter or provide a larger 'upper'.")
  }

  # Use uniroot for inversion
  root <- stats::uniroot(f, lower = lower, upper = upper, tol = tol)

  root$root
}


quadratureLogLik <- function(weighted) {
  rs <- rowSums(weighted)

  mod_warnif(any(rs < 0, na.rm = TRUE),
             "Found negative quadrature node contributions, this is likely a bug!")

  rs[!is.finite(rs)] <- 0
  rs.safe <- pmax(rs, .Machine$double.xmin)
  sum(log(rs.safe))
}


quadratureRelativeError <- function(error, reference) {
  if (!is.finite(reference) || reference == 0)
    return(NA_real_)

  abs(error) / abs(reference)
}


pruneQuadratureNodes <- function(quadw, quadn, quadf, a, b, tol) {
  mod_stopif(!is.numeric(quadw), "`quadw` must be numeric.")

  weight_vec <- as.numeric(quadw)
  n.nodes <- NCOL(quadf)
  mod_stopif(length(weight_vec) != n.nodes,
         "`quadw` must have the same length as the number of quadrature nodes.")

  n.input <- NROW(quadn)

  # A node whose density is not finite carries no usable information about
  # where to place the rule, so treat it as empty rather than letting the NaN
  # propagate into every downstream comparison. Ordinal kernels underflow
  # readily -- they are a product of one categorical probability per indicator
  # -- and the E-step itself, which works in log space, is unaffected.
  quadf[!is.finite(quadf)] <- 0

  # precompute weighted information to drop empty nodes early
  weighted <- sweep(quadf, 2, weight_vec, "*")
  I.full <- quadratureLogLik(weighted)

  zeroInfoNodes <- colSums(weighted) <= .Machine$double.xmin

  lower_vec <- rep(a, length.out = NCOL(quadn))
  upper_vec <- rep(b, length.out = NCOL(quadn))
  lower_mat <- matrix(lower_vec, nrow = NROW(quadn), ncol = NCOL(quadn), byrow = TRUE)
  upper_mat <- matrix(upper_vec, nrow = NROW(quadn), ncol = NCOL(quadn), byrow = TRUE)
  nodesOutside <- rowSums((quadn < lower_mat) | (quadn > upper_mat)) > 0

  isValidNode <- !(zeroInfoNodes | nodesOutside)

  quadn <- quadn[isValidNode, , drop = FALSE]
  quadf <- quadf[, isValidNode, drop = FALSE]
  weight_vec <- weight_vec[isValidNode]
  weighted <- weighted[, isValidNode, drop = FALSE]

  # Calculate per node contributions
  rs <- rowSums(weighted)

  # guard against log(0) / division by 0
  rs.safe <- pmax(rs, .Machine$double.xmin)

  # B[j,i] = log1p( - A[j,i] / rs[j] )
  #
  # An observation whose density underflows to zero at every node leaves
  # `rs.safe` at the smallest representable double, so the ratio blows up and
  # `log1p()` returns NaN. Such a row carries no information about which nodes
  # to keep, so its contribution is zero rather than a propagating NaN. This
  # happens readily for ordinal indicators, where the kernel is a product of
  # one categorical probability per indicator.
  ratio <- sweep(weighted, 1, rs.safe, "/")
  ratio[!is.finite(ratio)] <- 0
  ratio <- pmin(pmax(ratio, 0), 1 - .Machine$double.eps)
  B <- log1p(-ratio)
  I.base <- quadratureLogLik(weighted)
  I.subvec <- I.base + colSums(B)

  contributions <- abs(abs(I.subvec) - abs(I.base))
  contributions[!is.finite(contributions)] <- 0

  # identify nodes trivially safe to remove
  contrib.rank <- order(abs(contributions))
  cumulative   <- cumsum(contributions[contrib.rank])
  is.removable <- abs(cumulative) < tol * abs(I.base)

  # reverse ordering of is.removable
  is.removable <- is.removable[order(contrib.rank)]
  removable    <- which(is.removable)

  mod_warnif(sum(contributions[removable]) > tol * abs(I.base),
         paste0("Something went wrong when pruning nodes:\n",
         "More information than expected was lost!\n"), .newline = TRUE)

  mod_stopif(length(removable) >= NROW(quadn), "Cannot remove all nodes!")

  if (length(removable) > 0) {
    quadn <- quadn[-removable, , drop = FALSE]
    quadf <- quadf[, -removable, drop = FALSE]
    weighted <- weighted[, -removable, drop = FALSE]
    weight_vec <- weight_vec[-removable]
  }

  I.cur <- quadratureLogLik(weighted)
  I.err <- I.full - I.cur
  I.err.abs <- abs(I.err)
  I.err.rel <- quadratureRelativeError(I.err, I.full)

  n.output  <- NROW(quadn)
  n.removed <- n.input - n.output

  list(
    quadw   = weight_vec,
    quadn   = quadn,
    quadf   = quadf,
    I.full  = I.full,
    I.cur   = I.cur,
    I.err   = I.err,
    I.err.abs = I.err.abs,
    I.err.rel = I.err.rel,
    n.in    = n.input,
    n.out   = n.output,
    n.rm    = n.removed
  )
}
