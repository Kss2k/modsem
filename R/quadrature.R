# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000). Also stores information
# regarding adaptive quadrature.
quadrature <- function(m, k, 
                       cut = Inf, 
                       adaptive = FALSE, 
                       m.start = 4, 
                       ...) {
  a <- -cut
  b <- cut

  if (k == 0 || m == 0) {
    return(list(
      n = matrix(0), 
      w = 1, 
      k = 0, 
      m = m,
      a = a, 
      b = b, 
      cut = cut,
      m.start = m.start, 
      adaptive = FALSE
    ))
  }
  

  if (is.finite(cut)) {
    singleDimGauss <- finiteGaussQuadrature(m = m, k = k, a = a, b = b)
    nodes <- singleDimGauss$nodes
    weights <- singleDimGauss$weights

  } else {
    singleDimGauss <- fastGHQuad::gaussHermiteData(m)
    nodes <- singleDimGauss$x
    weights <- singleDimGauss$w

    nodes <- lapply(seq_len(k), function(k) nodes) |>
      expand.grid() |> as.matrix()
    weights <- lapply(seq_len(k), function(k) weights) |>
      expand.grid() |> apply(MARGIN = 1, prod)
    
    nodes <- sqrt(2) * nodes
    weights <- weights * pi ^ (-k/2)
  }

  list(
    n = nodes,
    w = weights,
    k = k, 
    m = m,
    a = a, 
    b = b, 
    cut = cut,
    adaptive = adaptive
  )
}


finiteGaussQuadrature <- function(a, b, m = 10, k = 1) {
  if (k == 0 || m == 0) {
    return(list(nodes = matrix(0), weights = 1, intersects = matrix(0), m = m, k = k))
  }

  if (length(a) == 1L) a <- rep_len(a, k)
  if (length(b) == 1L) b <- rep_len(b, k)
  stopifnot(length(a) == k, length(b) == k)

  build_1d <- function(a1, b1, m1) {
    intervals <- (b1 - a1) / m1
    intersects <- (seq_len(m1 + 1L) - 1L) * intervals + a1
    nodes1d <- weights1d <- numeric(m1)
    for (j in seq_len(m1)) {
      x1 <- intersects[[j]]
      x2 <- intersects[[j + 1L]]
      nodes1d[[j]]   <- (x1 + x2) / 2
      weights1d[[j]] <- stats::pnorm(x2) - stats::pnorm(x1)
    }
    list(nodes = nodes1d, weights = weights1d, intersects = intersects)
  }

  per_dim <- Map(build_1d, a, b, m)

  grid   <- expand.grid(lapply(per_dim, `[[`, "nodes"))
  w_list <- expand.grid(lapply(per_dim, `[[`, "weights"))
  weights <- apply(w_list, 1L, prod)
  nodes   <- as.matrix(grid)

  list(nodes = nodes,
       weights = weights,
       intersects = lapply(per_dim, `[[`, "intersects"),
       m = m, k = k)
}


collapseQuadratureSectors <- function(sectors, collapse = \(x) x) {
  values <- NULL

  for (quad in sectors) {
    f <- quad$f
    w <- matrix(quad$weights, ncol = length(quad$weights), nrow = NROW(f),
                byrow = TRUE)

    values <- cbind(values, w * f)
  }

  collapse(values)
}


adaptiveGaussQuadrature <- function(fun, 
                                    collapse = \(x) x, # function to collapse the results, if 'relevant'
                                    a = -7, 
                                    b = 7, 
                                    m = 32,
                                    m.ceil = 150,
                                    k = 1, 
                                    tol = 1e-12, ...) {
  stopif(k > 1, "Adaptive quadrature not implemented for multiple dimensions!")

  if (k == 0 || m == 0) {
    return(list(integral = 0, n = matrix(0), w = 1, f = NA, m = 0, k = k))
  }

  m.target <- m ^ k

  quad <- quadrature(m = m.ceil, k = k)
  quadn <- quad$n
  quadf <- fun(quadn, ...)
  quadw <- matrix(quad$w, ncol = length(quad$w), nrow = NROW(quadf), byrow = TRUE)

  integral        <- collapse(quadf * quadw)

  zeroInfoNodes <- apply(quadw * quadf, MARGIN=2, FUN = \(x) sum(x) <= .Machine$double.xmin * 2)

  quadn <- quadn[!zeroInfoNodes, , drop = FALSE]
  quadf <- quadf[ , !zeroInfoNodes, drop = FALSE]
  quadw <- quadw[ , !zeroInfoNodes, drop = FALSE]
  
  while (TRUE) {
    integral.current <- collapse(quadf * quadw)
    error <- abs(integral.current - integral)
    m.current <- nrow(quadn)

    if (error < tol * integral && m.current > m.start) break

    min.contrib <- NA 
    min.index <- NA

    for (i in seq_len(m.current)) {
      quadsub_n <- quadn[-i, , drop = FALSE]
      quadsub_f <- quadf[-i, , drop = FALSE]
      quadsub_w <- quadw[-i, , drop = FALSE]

      contrib <- abs(integral.current - collapse(quadsub_f * quadsub_w))

      isValid    <- contrib < tol * integral
      isSmallest <- is.na(min.contrib) || contrib < min.contrib

      if (isValid && isSmallest) {
        min.contrib <- min(contrib)
        min.index <- i
      }
    }

    if (is.na(min.index)) break

    quadn <- quadn[-min.index, , drop = FALSE]
    quadf <- quadf[, -min.index, drop = FALSE]
    quadw <- quadw[, -min.index, drop = FALSE]
  }


  diff.m <- round(m.target - nrow(quadn)) ^ (1 / k)

  if (abs(diff.m) > 5) {
    new.ceil <- m.ceil + diff.m
    return(adaptiveGaussQuadrature(
      fun = fun, 
      collapse = collapse, 
      a = a, 
      b = b, 
      m = m, 
      m.ceil = new.ceil, 
      k = k, 
      tol = tol, 
      ...
    ))
  }

  list(integral = integral.current, 
       error = error,
       n = quadn, 
       w = quadw[1, , drop = TRUE], 
       f = quadf, 
       m = m.current ^ (1 / k), 
       m.ceil = m.ceil,
       k = k)
}
