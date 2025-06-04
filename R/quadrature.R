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


finiteGaussQuadrature <- function(a, b, m = 10, k = 1) {
  intervals <- (b - a) / (m)
  intersects <- (seq_len(m + 1) -1) * intervals + a

  nodes <- numeric(m)
  weights <- numeric(m)

  for (i in seq_len(m)) {
    x1 <- intersects[[i]]
    x2 <- intersects[[i + 1]]

    nodes[[i]] <- (x1 + x2) / 2
    weights[[i]] <- stats::pnorm(x2) - stats::pnorm(x1) # Normal distribution CDF
  }

  nodes <- lapply(seq_len(k), function(k) nodes) |>
    expand.grid() |> as.matrix()
  weights <- lapply(seq_len(k), function(k) weights) |>
    expand.grid() |> apply(MARGIN = 1, prod)

  list(nodes = nodes, 
       weights = weights,
       intersects = intersects,
       m = m, k = k)
}



adaptiveGaussQuadrature <- function(fun, a = -7, b = 7, m.start = 4, m.max = 32, tol = 1e-6,
                                    k = 1, total.integral = NULL, ...) {
  stopif(k > 1, "Adaptive quadrature for k > 1 is not supported yet. Use fixed instead!")

  if (m.start >= m.max) {
    quad <- finiteGaussQuadrature(a=a, b=b, m=m.max, k=k)
    f <- fun(quad$nodes, ...)
    integral <-  sum(quad$weights * f)

    return(
      list(
         integral = integral, 
         f = f,
         n = quad$nodes, 
         w = quad$weights, 
         m = quad$m,
         k = k
      )
    )
  }

  quad.l <- finiteGaussQuadrature(a=a, b=b, m=m.start, k=k)
  quad.h <- finiteGaussQuadrature(a=a, b=b, m=m.start * 2, k=k)

  f.l <- fun(quad.l$nodes, ...)
  f.h <- fun(quad.h$nodes, ...)
  integral.l <- sum(quad.l$weights * f.l)
  integral.h <- sum(quad.h$weights * f.h)

  if (is.null(total.integral)) total.integral <- integral.h

  E <- abs(integral.h - integral.l) / abs(total.integral)
 
  if (E < tol || m.start >= m.max) {
    return(
      list(
         integral = integral.l, 
         f = f.l, 
         n = quad.l$nodes, 
         w = quad.l$weights, 
         m = quad.l$m,
         k = k
      )
    )
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

  left <- adaptiveGaussQuadrature(fun=fun, a=a_left, b=b_left, m.start=m.start,
                                  m.max=round(m.max / 2), tol=tol, 
                                  total.integral=total.integral, k = k, ...)

  right <- adaptiveGaussQuadrature(fun=fun, a=a_right, b=b_right, m.start=m.start, 
                                   m.max=m.max - left$m, tol=tol, 
                                   total.integral=total.integral, k = k, ...)

  list(
    integral = left$integral + right$integral,
    f = c(left$f, right$f),
    m = left$m + right$m,
    n = rbind(left$n, right$n),
    w = c(left$w, right$w),
    k = k
  )
}


q <- adaptiveGaussQuadrature(fun=\(x) dnorm(4*x) + 1, a = -7, b = 7, m.max=64)
q <- adaptiveGaussQuadrature(fun=\(x) mvtnorm::dmvnorm(x, log=TRUE), a = -7, b = 7, m.max=64, k = 2)
q <- adaptiveGaussQuadrature(fun=\(x) 1, a = -7, b = 7, m.max=64, k = 2)
# plot(q$nodes, q$weights)
# q$weights |> sum()
# diff(sort(q$nodes))
# 
# o <- order(q$nodes)
# nodes <- q$nodes[o]
# weights <- q$weights[o]
# plot(nodes, cumsum(weights))
