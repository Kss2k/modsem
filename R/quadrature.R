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
                                    iter.max = 20,
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

  diff.m <- round(nrow(quadn) ^ (1 / k)) - m

  if (diff.m != 0 && iter < iter.max) {
    new.ceil <- m.ceil - diff.m
    return(adaptiveGaussQuadrature(
      fun = fun, 
      collapse = collapse, 
      a = a, 
      b = b, 
      k = k, 
      m = m, 
      m.ceil = new.ceil, 
      iter = iter + 1,
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
