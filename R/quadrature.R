# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000)
quadrature <- function(m, k, cut = Inf) {
  if (k == 0 || m == 0) return(list(n = matrix(0), w = 1, k = 0, m = m))
  singleDimGauss <- fastGHQuad::gaussHermiteData(m)

  nodes <- singleDimGauss$x
  weights <- singleDimGauss$w

  select <- abs(nodes) < cut
  nodes <- nodes[select]
  weights  <- weights[select]
  m <- length(weights)

  nodes <- lapply(seq_len(k), function(k) nodes) |>
    expand.grid() |> as.matrix()
  weights <- lapply(seq_len(k), function(k) weights) |>
    expand.grid() |> apply(MARGIN = 1, prod)
  list(n = nodes * sqrt(2), w = weights * pi ^ (-k/2), k = k, m = m)
}


finiteGaussQuadrature <- function(a, b, m = 10) {
  intervals <- (b - a) / (m)
  intersects <- (seq_len(m + 1) -1) * intervals + a

  intersects
  nodes <- numeric(m)
  weights <- numeric(m)

  for (i in seq_len(m)) {
    x1 <- intersects[[i]]
    x2 <- intersects[[i + 1]]

    nodes[[i]] <- (x1 + x2) / 2
    weights[[i]] <- stats::pnorm(x2) - stats::pnorm(x1) # Normal distribution CDF
  }

  list(nodes = nodes, 
       weights = weights,
       intersects = intersects)
}



adaptiveGaussQuadrature <- function(fun, a = -7, b = 7, m.start = 4, m.max = 32, tol = 1e-6, ...) {
  if (m.start >= m.max) {
    quad <- finiteGaussQuadrature(a, b, m.max)
    integral <- sum(quad$weights * fun(quad$nodes, ...))

    return(list(integral = integral, n = quad$nodes, w = quad$weights))
  }

  quad.l <- finiteGaussQuadrature(a, b, m.start)
  quad.h <- finiteGaussQuadrature(a, b, m.start * 2)

  integral.l <- sum(quad.l$weights * fun(quad.l$nodes, ...))
  integral.h <- sum(quad.h$weights * fun(quad.h$nodes, ...))

  E <- abs(integral.h - integral.l) / abs(integral.h)
  
  # cat(sprintf("a: %f, b: %f, m.start: %d, m.max: %.f, E: %.6f\n", a, b, m.start, m.max, E))
  
  if (E < tol || m.start >= m.max) {
    return(list(integral = integral.l, n = quad.l$nodes, w = quad.l$weights))
  }

  c <- (a + b) / 2

  Ea <- dnorm(a) * fun(a, ...)
  Eb <- dnorm(b) * fun(b, ...)

  if (Ea < Eb) {
    a_left <- a
    b_left <- c
    a_right <- c
    b_right <- b
  } else {
    a_left <- c
    b_left <- b
    a_right <- a
    b_right <- c
  }

  left <- adaptiveGaussQuadrature(fun, a_left, b_left, m.start + 1, round(m.max / 2), tol, ...)
  
  left_m <- length(left$n)
  right <- adaptiveGaussQuadrature(fun, a_right, b_right, m.start + 1, m.max - left_m, tol, ...)

  if (abs(sum(left$w) + sum(right$w) - pnorm(b) + pnorm(a)) > 1e-12)  {
    browser()
    stop("Weights sum to incorrect value after recursion. Check the function or bounds.")
  }

  list(integral = left$integral + right$integral,
       n = c(left$n, right$n),
       w = c(left$w, right$w))
}


# q <- adaptiveGaussQuadrature(fun=\(x) dnorm(4*x) + 1, a = -7, b = 7, m.max=64)
# plot(q$nodes, q$weights)
# q$weights |> sum()
# diff(sort(q$nodes))
# 
# o <- order(q$nodes)
# nodes <- q$nodes[o]
# weights <- q$weights[o]
# plot(nodes, cumsum(weights))
