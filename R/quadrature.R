# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000). Also stores information
# regarding adaptive quadrature.
quadrature <- function(m, k,
                       adaptive = FALSE,
                       quad.range = Inf,
                       adaptive.frequency = 3,
                       ...) {
  if (quad.range < 0) {
    warning2("`quad.range` should be positive, using `-quad.range` instead!\n")
    quad.range <- -quad.range
  }

  a <- -quad.range
  b <- quad.range

  if (k == 0 || m == 0) {
    return(list(
      n = matrix(0),
      w = 1,
      k = 0,
      m = m,
      a = a,
      b = b,
      adaptive = FALSE,
      quad.range = quad.range,
      adaptive.frequency = adaptive.frequency
    ))
  }

  singleDimGauss <- fastGHQuad::gaussHermiteData(m)
  nodes <- singleDimGauss$x
  weights <- singleDimGauss$w

  nodes_list   <- replicate(k, nodes, simplify = FALSE)
  weight_list  <- replicate(k, weights, simplify = FALSE)
  nodes   <- do.call(expand.grid, nodes_list) |> as.matrix()
  weights <- Reduce(kronecker, rev(weight_list))

  nodes <- sqrt(2) * nodes
  weights <- weights * pi ^ (-k/2)

  list(
    n = nodes,
    w = weights,
    k = k,
    m = m,
    a = a,
    b = b,
    quad.range = quad.range,
    adaptive = adaptive,
    adaptive.frequency = adaptive.frequency
  )
}


adaptiveGaussQuadrature <- function(fun,
                                    a = -7,
                                    b = 7,
                                    m = 32,
                                    m.ceil = m + m / 2,
                                    k = 1,
                                    iter.max = 10,
                                    node.max = 1000,
                                    tol = 1e-12,
                                    mdiff.tol = 2,
                                    secondary.pruning = TRUE,
                                    ...) {
  if (k == 0 || m == 0)
    return(list(n = matrix(0), w = 1, f = NA, m = 1, k = 1))

  stopif(tol >= 1 || tol < 0,
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

  quadn <- do.call(expand.grid, NODES) |> as.matrix()
  quadw <- Reduce(kronecker, rev(WEIGHTS))
  quadf <- fun(quadn, ...)

  if (secondary.pruning) {
    pruned <- pruneQuadratureNodes(
      quadw = quadw, quadn = quadn, quadf = quadf,
      a = a, b = b, tol = tol
    )

    quadw  <- pruned$quadw
    quadn  <- pruned$quadn
    quadf  <- pruned$quadf
  }

  list(n = quadn,
       w = quadw,
       F = quadf,
       m.ceil = new.ceils)
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
                                     node.max = 1000,
                                     tol = 1e-12,
                                     mdiff.tol = 2,
                                     ...) {
  if (k == 0 || m == 0)
    return(list(n = matrix(0), w = 1, f = NA, m = 1, k = 1))

  if (is.null(m.ceil) || is.na(m.ceil) || m.ceil <= 0)
    m.ceil <- round(estMForNodesInRange(m, a = -5, b = 5))

  quad  <- quadrature(m = m.ceil, k = 1)

  if (K > 1) { # more computationally efficient/scalable
    quadn <- matrix(0, nrow = m.ceil, ncol = K)
    quadn[,k] <- quad$n
  } else quadn <- quad$n

  quadf <- fun(quadn, ...)
  quadw <- quad$w

  pruned <- pruneQuadratureNodes(
    quadw = quadw, quadn = quadn, quadf = quadf,
    a = a, b = b, tol = tol
  )

  quadw  <- pruned$quadw
  quadn  <- pruned$quadn
  quadf  <- pruned$quadf
  I.err  <- pruned$I.err
  I.full <- pruned$I.full
  I.cur  <- pruned$I.cur

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

  warnif(iter >= iter.max,
         "Max iterations reached fitting quasi-adaptive quadrature...\n",
         sprintf("Iter %d, total: %d, target: %d, kept: %d, discarded: %d",
                 iter, m.ceil, m, NROW(quadn), m.ceil - NROW(quadn)),
         .newline = TRUE)

  list(n = quadn,
       w = quadw,
       F = quadf,
       k = k,
       m = nrow(quadn) ^ (1 / k),

       m.ceil   = m.ceil,
       iter     = iter,
       error    = I.err,
       integral = I.cur)
}


estGHNodesInRange <- function(m, a, b, scale = TRUE) {
  stopif(!is.numeric(m) || length(m) != 1 || m <= 0,
         "'m' must be a single positive number")
  stopif(!is.numeric(a) || length(a) != 1 || !is.numeric(b) || length(b) != 1,
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
  stopif(!is.numeric(k) || length(k) != 1 || k <= 0,
         "'k' must be a single positive number")
  stopif(!is.numeric(a) || length(a) != 1 || !is.numeric(b) || length(b) != 1,
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
    stopif(f(upper) < 0, "Unable to bracket root: increase maxiter or provide a larger 'upper'.")
  }

  # Use uniroot for inversion
  root <- stats::uniroot(f, lower = lower, upper = upper, tol = tol)

  root$root
}


pruneQuadratureNodes <- function(quadw, quadn, quadf, a, b, tol) {
  stopif(!is.numeric(quadw), "`quadw` must be numeric.")

  weight_vec <- as.numeric(quadw)
  n.nodes <- NCOL(quadf)
  stopif(length(weight_vec) != n.nodes,
         "`quadw` must have the same length as the number of quadrature nodes.")

  n.input <- NROW(quadn)

  # precompute weighted information to drop empty nodes early
  weighted <- sweep(quadf, 2, weight_vec, "*")
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
  warnif(any(rs < 0), "Found negative quadrature node contributions, this is likely a bug!")
  rs.safe <- pmax(rs, .Machine$double.xmin)

  I.full <- sum(log(rs.safe))

  # B[j,i] = log1p( - A[j,i] / rs[j] )
  B <- log1p(-sweep(weighted, 1, rs.safe, "/"))
  I.subvec <- I.full + colSums(B)

  contributions <- abs(abs(I.subvec) - abs(I.full))

  # identify nodes trivially safe to remove
  contrib.rank <- order(abs(contributions))
  cumulative   <- cumsum(contributions[contrib.rank])
  is.removable <- abs(cumulative) < tol * abs(I.full)

  # reverse ordering of is.removable
  is.removable <- is.removable[order(contrib.rank)]
  removable    <- which(is.removable)

  warnif(sum(contributions[removable]) > tol * abs(I.full),
         "Something went wrong when pruning nodes:\n",
         "More information than expected was lost!\n", .newline = TRUE)

  stopif(length(removable) >= NROW(quadn), "Cannot remove all nodes!")

  I.cur <- I.full
  I.err <- 0

  if (length(removable) > 0) {
    quadn <- quadn[-removable, , drop = FALSE]
    quadf <- quadf[, -removable, drop = FALSE]
    weighted <- weighted[, -removable, drop = FALSE]
    weight_vec <- weight_vec[-removable]

    I.cur  <- I.full - sum(contributions[removable])
    I.err  <- I.full - I.cur
  }

  n.output  <- NROW(quadn)
  n.removed <- n.input - n.output

  list(
    quadw   = weight_vec,
    quadn   = quadn,
    quadf   = quadf,
    I.full  = I.full,
    I.cur   = I.cur,
    I.err   = I.err,
    n.in    = n.input,
    n.out   = n.output,
    n.rm    = n.removed
  )
}
