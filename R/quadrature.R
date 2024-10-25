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
