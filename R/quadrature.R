# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000). Also stores information
# regarding adaptive quadrature.
quadrature <- function(m, k, 
                       cut = Inf, 
                       adaptive = FALSE, 
                       a = -7, 
                       b = 7, 
                       m.start = 4, 
                       adaptive.quad.tol = 1e-4,
                       ...) {
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

  singleDimGauss <- fastGHQuad::gaussHermiteData(m)

  nodes <- singleDimGauss$x
  weights <- singleDimGauss$w

  select <- abs(nodes) < cut
  nodes <- nodes[select]
  weights  <- weights[select]
  m <- if (!adaptive) length(weights) else m

  nodes <- lapply(seq_len(k), function(k) nodes) |>
    expand.grid() |> as.matrix()
  weights <- lapply(seq_len(k), function(k) weights) |>
    expand.grid() |> apply(MARGIN = 1, prod)

  if (is.finite(cut)) {
    a <- -cut
    b <- cut
  }

  list(
    n = nodes * sqrt(2), 
    w = weights * pi ^ (-k/2), 
    k = k, 
    m = m,
    a = a, 
    b = b, 
    cut = cut,
    m.start = m.start, 
    adaptive = adaptive,
    adaptive.quad.tol = adaptive.quad.tol
  )
}
