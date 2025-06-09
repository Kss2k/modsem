# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000). Also stores information
# regarding adaptive quadrature.
quadrature <- function(m, k, 
                       cut = Inf, 
                       adaptive = FALSE, 
                       a = -7, 
                       b = 7, 
                       m.start = 4, 
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
  

  if (is.finite(cut)) {
    a <- -cut
    b <- cut
    singleDimGauss <- finiteGaussQuadrature(m = m, k = k, a = -7, b = 7)
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


adaptiveGaussQuadrature <- function(fun, a = -7, b = 7, m = 32,
                                    k = 1, ...) {
  stopif(k > 1, "Adaptive quadrature not implemented for multiple dimensions!")

  if (k == 0 || m == 0) {
    return(list(integral = 0, n = matrix(0), w = 1, f = NA, m = 0, k = k))
  }

  nsectors <- round(m / 2) # in the base case all sectors get two nodes each
  step <- (b - a) / 12 
  lower <- a + (seq_len(nsectors) - 1) * step
  upper <- a + seq_len(nsectors) * step
  
  sectorsLow  <- vector("list", nsectors)
  sectorsHigh <- vector("list", nsectors)

  for (i in seq_len(nsectors)) {
    low <- finiteGaussQuadrature(m = 1, a = lower[i], b = upper[i], k = k)
    high <- finiteGaussQuadrature(m = 5, a = lower[i], b = upper[i], k = k)
    
    low$d <- fun(low$nodes, ...)
    low$f <- low$weights * low$d
    
    high$d <- fun(high$nodes, ...)
    high$f <- high$weights * high$d

    sectorsLow[[i]]  <- low
    sectorsHigh[[i]] <- high
  }


  available <- m - nsectors
  sector <- 1
  while (available) {
    maxdiff       <- NA
    maxdiffSector <- NA

    for (i in seq_len(nsectors)) {
      high <- sectorsHigh[[i]]
      low  <- sectorsLow[[i]]

      diff <- abs(sum(high$f) - sum(low$f))

      if (is.na(maxdiff) || diff > maxdiff) {
        maxdiff       <- diff
        maxdiffSector <- i
      }
    }

    stopif(is.na(maxdiff), "Unexpected NA in maxdiff!")

    newM <- sectorsLow[[maxdiffSector]]$m + 1
    newQ <- finiteGaussQuadrature(m = newM, k = k,
                                  a = lower[maxdiffSector], 
                                  b = upper[maxdiffSector])

    newQ$d <- fun(newQ$nodes, ...)
    newQ$f <- newQ$weights * newQ$d

    sectorsLow[[maxdiffSector]] <- newQ
    available <- available - 1
  }

  n <- NULL
  w <- NULL
  f <- NULL
  d <- NULL
  
  for (quad in sectorsLow) {
    n <- rbind(n, quad$nodes)
    w <- c(w, quad$w)
    f <- c(f, quad$f)
    d <- c(d, quad$d)
  }

  list(n = n, w = w, f = f, d = d, k = k, m = m)
}
