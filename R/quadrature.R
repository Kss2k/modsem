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
      weights1d[[j]] <- pnorm(x2) - pnorm(x1)
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


adaptiveGaussQuadrature <- function(fun, a = -7, b = 7, m.start = 4, m.max = 32,
                                    tol = 1e-6, k = 1, total.integral = NULL, ...) {
  stopif(k > 1, "Adaptive quadrature not implemented for multiple dimensions!")

  if (k == 0 || m.max == 0) {
    return(list(integral = 0, n = matrix(0), w = 1, f = NA, m = 0, k = k))
  }

  if (m.start >= m.max) {
    quad <- finiteGaussQuadrature(a = a, b = b, m = m.max, k = k)
    f <- fun(quad$nodes, ...)
    integral <- sum(quad$weights * f)
    return(list(integral = integral, f = f, n = quad$nodes, w = quad$weights, m = quad$m, k = k))
  }

  quad.l <- finiteGaussQuadrature(a = a, b = b, m = m.start, k = k)
  quad.h <- finiteGaussQuadrature(a = a, b = b, m = m.start * 2, k = k)
  f.l <- fun(quad.l$nodes, ...)
  f.h <- fun(quad.h$nodes, ...)
  integral.l <- sum(quad.l$weights * f.l)
  integral.h <- sum(quad.h$weights * f.h)

  if (is.null(total.integral)) total.integral <- integral.h

  E <- abs(integral.h - integral.l) / abs(total.integral)

  if (E < tol || m.start >= m.max) {
    return(list(integral = integral.l, f = f.l, n = quad.l$nodes, w = quad.l$weights, m = quad.l$m, k = k))
  }

  c_mid <- (a + b) / 2
  Ea <- dnorm(a) * fun(matrix(a), ...)
  Eb <- dnorm(b) * fun(matrix(b), ...)

  if (Ea < Eb) {
    a_left <- a; b_left <- c_mid; a_right <- c_mid; b_right <- b
  } else {
    a_left <- c_mid; b_left <- b; a_right <- a; b_right <- c_mid
  }

  left <- adaptiveGaussQuadrature(fun, a_left, b_left, m.start, round(m.max / 2), tol, 1, total.integral, ...)
  right <- adaptiveGaussQuadrature(fun, a_right, b_right, m.start, m.max - left$m, tol, 1, total.integral, ...)

  list(integral = left$integral + right$integral, f = c(left$f, right$f), m = left$m + right$m,
       n = rbind(left$n, right$n), w = c(left$w, right$w), k = k)
}


# This function is a work in progress...
# laplaceAdaptiveQuadrature <- function(fun, init, L = 6, m = 12, tol = 1e-6, max.expand = 3, ...) {
#   stopifnot(is.numeric(init))
#   k <- length(init)
#
#   # --- 1. locate mode ------------------------------------------------------
#   nlogf <- function(theta) {
#     val <- fun(matrix(theta, nrow = 1L), ...)
#     if (val <= 0 || !is.finite(val)) return(Inf)
#     -log(val)
#   }
#   opt <- optim(init, nlogf, hessian = TRUE, control = list(fnscale = 1))
#   if (opt$convergence != 0)
#     warning("optim did not converge; result may be inaccurate")
#   mu <- opt$par
#
#   # --- 2. Hessian & covariance --------------------------------------------
#   H <- opt$hessian
#   if (is.null(H) || any(is.na(H)) || any(!is.finite(H))) {
#     H <- numDeriv::hessian(function(x) -log(fun(matrix(x, nrow = 1L), ...)), mu)
#   }
#   Sigma <- tryCatch(solve(H), error = function(e) {
#     warning("Hessian not invertible; using identity covariance")
#     diag(1, k)
#   })
#   sds <- sqrt(diag(Sigma))
#
#   # helper to compute integral on [a,b] for given L ------------------------
#   integrate_rect <- function(Lmult) {
#     a <- mu - Lmult * sds
#     b <- mu + Lmult * sds
#     quad <- finiteGaussQuadrature(a, b, m = m, k = k)
#     list(val = sum(quad$weights * fun(quad$nodes, ...)), quad = quad, a = a, b = b)
#   }
#
#   base <- integrate_rect(L)      # core region
#   total <- base$val
#   Lcurr <- L
#
#   # --- 3. expand region adaptively ----------------------------------------
#   for (step in seq_len(max.expand)) {
#     Lnext <- Lcurr * 1.5           # geometric expansion
#     outer <- integrate_rect(Lnext)
#     incr  <- outer$val - total
#     if (abs(incr) / abs(outer$val) < tol) {
#       total <- outer$val; base <- outer; break
#     }
#     total <- outer$val; Lcurr <- Lnext; base <- outer
#   }
#
#   with(base, {
#     list(integral = val, mu = mu, Sigma = Sigma, L = Lcurr, k = k,
#          a = a, b = b, n = quad$nodes, w = quad$weights)
#   })
# }