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

  nodes <- lapply(seq_len(k), function(k) nodes) |>
    expand.grid() |> as.matrix()
  weights <- lapply(seq_len(k), function(k) weights) |>
    expand.grid() |> apply(MARGIN = 1, prod)

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
                                    collapse = \(x) x, # function to collapse the results, if 'relevant'
                                    a = -7, 
                                    b = 7, 
                                    m = 32,
                                    m.ceil = m + m / 2,
                                    k = 1, 
                                    iter = 1,
                                    iter.max = 40,
                                    node.max = 1000,
                                    tol = 1e-12,
                                    mdiff.tol = 2,
                                    ...) {
  if (k == 0 || m == 0) return(list(n = matrix(0), w = 1, f = NA, m = 1, k = 1))

  m.target <- m ^ k

  quad <- quadrature(m = m.ceil, k = k)
  quadn <- quad$n
  quadf <- fun(quadn, ...)
  quadw <- matrix(quad$w, ncol = length(quad$w), nrow = NROW(quadf), byrow = TRUE)

  integral.full <- collapse(quadf * quadw)

  zeroInfoNodes <- apply(quadw * quadf, MARGIN=2, FUN = \(x) sum(x) <= .Machine$double.xmin * 2)
  nodesOutside  <- apply(quadn, MARGIN = 1, FUN = \(x) any(x < a | x > b))
  isValidNode   <- !(zeroInfoNodes | nodesOutside)

  quadn <- quadn[isValidNode, , drop = FALSE]
  quadf <- quadf[, isValidNode, drop = FALSE]
  quadw <- quadw[, isValidNode, drop = FALSE]
  
  # vector of full-integral for thresholding
  I.full <- integral.full

  # This loop is in practive quite superfluous... I've commented it out for now
  # let's see how it works for a while

  # repeat {
    # current integral and error
    I.cur <- collapse(quadf * quadw)
    err   <- abs(I.cur - I.full)
    m.cur <- nrow(quadn)
    m.start  <- m.cur

    # if we're within tolerance or cannot drop below m.start, stop
  #  if (err > tol * abs(I.full) || m.cur <= 2) break

    # compute per-node contributions in one pass:
    # c[j] = contribution of node j to the integral
    contributions <- numeric(length=NCOL(quadf))
    for (i in seq_len(NCOL(quadf))) {
      I.sub <- collapse(quadf[, -i, drop=FALSE] * quadw[, -i, drop=FALSE])
      contributions[[i]] <- abs(abs(I.sub) - abs(I.full))
    }

    # identify nodes trivially safe to remove
    removable <- which(abs(contributions) < tol * abs(I.full))

    if (length(removable) > 0) {
      # update everything by dropping them all at once
      quadn <- quadn[-removable, , drop = FALSE]
      quadf <- quadf[, -removable, drop = FALSE]
      quadw <- quadw[, -removable, drop = FALSE]
      # subtract their total from the current integral
      I.cur  <- I.cur - sum(contributions[removable])
  #    next
    }

  #  # This never does anything, might make it do something later
  #  # I.e., this criterion is just as strict as `removable`
  #  # so removable will remove any valid candidates...

  #  # if no bulk‐removable nodes, try peeling off the single smallest
  #  # rank.idx <- order(abs(contributions))
  #  # for (j in rank.idx) {
  #  #   if (m.cur <= 2) break
  #  #   # would removing j be acceptable?
  #  #   if (abs(contributions[j]) < err) {
  #  #     quadn <- quadn[-j, , drop = FALSE]
  #  #     quadf <- quadf[, -j, drop = FALSE]
  #  #     quadw <- quadw[, -j, drop = FALSE]
  #  #     I.cur  <- I.cur - contributions[j]
  #  #     m.cur  <- m.cur - 1
  #  #     # update error and restart removal loop
  #  #     err <- abs(I.cur - I.full)
  #  #     break
  #  #   } else break
  #  # }

  #  # if no single removal was made, exit
  #  if (err >= tol * abs(I.full) || m.cur == m.start) break
  #}

  lower  <- min(quadn)
  upper <- max(quadn)
  diff.m <- round(nrow(quadn) ^ (1 / k)) - m

  # Calculate next ceiling
  if (abs(diff.m) > 2 * mdiff.tol) {
    # If the difference is large, use a rough estimate
    new.ceil <- round(estMForNodesInRange(k = m, a = lower, b = upper))

    # else do a step change
  } else new.ceil <- m.ceil - diff.m

  OKMaxNodes <- m.ceil > node.max && new.ceil >= m.ceil
  OKDiff <- abs(diff.m) <= mdiff.tol
  OKIter <- iter < iter.max

  if (!OKDiff && OKIter && !OKMaxNodes) {

    return(adaptiveGaussQuadrature(
      fun = fun, 
      collapse = collapse, 
      a = a, 
      b = b, 
      k = k, 
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
         "Maximum number of iterations reached when calculating adaptive quadrature!\n",
         "Try lowering the number of nodes, or increase `quad.range` argument!")

  # only need to calculate this before returning the final version
  integral.reduced <- collapse(quadf * quadw)
  error <- integral.reduced - integral.full

  list(n = quadn, 
       w = quadw[1, , drop = TRUE], 
       f = quadf, 
       k = k, 
       m = nrow(quadn) ^ (1 / k), 

       m.ceil   = m.ceil,
       iter     = iter,
       error    = error,
       integral = integral.reduced)
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
