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
    # reverse order, if sorted in reverse
    ab <- c(a1, b1)
    a1 <- min(ab)
    b1 <- max(ab)

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


adaptiveGaussQuadrature <- function(fun, a = -7, b = 7, 
                                    m = 32, rec.max = 4, 
                                    m.start = 4, rel.tol = 1e-6,
                                    k = 1, type = "symmetric", 
                                    total.integral = NULL, 
                                    ...) {
  switch(type,
      symmetric = adaptiveGaussQuadratureSymmetric(
        fun = fun, a = a, b = b, m = m, rel.tol = rel.tol, 
        k = k, rec.max = rec.max, ...
      ),
      asymmetric = adaptiveGaussQuadratureASymmetric(
        fun = fun, a = a, b = b, m.max = m, m.start = m.start,
        tol = rel.tol, k = k, total.integral = total.integral, ...
      )
  )
}


adaptiveGaussQuadratureSymmetric <- function(fun, a = -7, b = 7,
                                             m = 32, rec.max = 4,
                                             rel.tol = 1e-10, k = 1,
                                             corners = TRUE, 
                                             slabs = TRUE, ...) {
  errors   <- NULL
  m.target <- m
  OK       <- FALSE
  iter     <- 0
  ai       <- if (length(a) == 1) rep(a, k) else a
  bi       <- if (length(b) == 1) rep(b, k) else b

  final.nodes   <- NULL
  final.weights <- NULL
  final.f       <- NULL

  while (!OK && (iter <- iter+1) <= rec.max && m.target >= 1) {
    quad    <- finiteGaussQuadrature(m = m.target, k = k, a = ai, b = bi)
    nodes   <- quad$nodes
    weights <- quad$weights
    mk      <- nrow(nodes)
    f       <- fun(nodes, ...)
    density <- weights * abs(f) # take the absolute density
    integral <- sum(density) # + sum(errors)

    o <- order(density) # sort by relative contribution to the integral

    density <- density[o]
    nodes   <- nodes[o, , drop=FALSE]
    weights <- weights[o]
    f       <- f[o]

    cumintegral <- cumsum(density)
    ok <- cumintegral > rel.tol * integral # should maybe divide my rec.max...
                                           # currently addressed by adding errors to integral
    OK <- all(ok)

    if (OK || iter >= rec.max) {
      final.nodes   <- rbind(final.nodes, nodes)
      final.weights <- c(final.weights, weights)
      final.f       <- c(final.f, f)

      break

    } else if (slabs) {
      errors <- c(errors, sum(density[!ok]))

      nodes   <- nodes[ok, , drop=FALSE]
      weights <- weights[ok]
      f       <- f[ok]

      aj <- apply(nodes, MARGIN=2, FUN=min)
      bj <- apply(nodes, MARGIN=2, FUN=max)

      slabs <- buildEdgeSlabs(ai = ai, aj = aj, bi = bi, bj = bj, corners = corners)
      A <- slabs$A
      B <- slabs$B

      nsectors <- NROW(A) # k * 2
      expected <- if (corners && k > 1) 2 * k + k ^ 2  else 2 * k

      warnif(nsectors != expected, "Unexpected number of quadrature quadrants to fill in the tail!\n")
      
      tail.nodes   <- matrix(NA, nrow=nsectors, ncol=k)
      tail.weights <- numeric(nsectors)
      tail.f       <- rep(NA, nsectors) # not worth the extra computation to fill in...

      for (i in seq_len(nrow(A))) { # for k > 1, we don't fill in the corner blocks, as they
                           # have lower densities..., thus we only fill in two extra 
                           # tail nodes per dimension (k)
        ak <- A[i, ]  
        bk <- B[i, ]

        tailquad <- finiteGaussQuadrature(m = 1, k = k, a = ak, b = bk)

        tail.nodes[i, ] <- tailquad$n
        tail.weights[[i]] <- tailquad$w
      }

      final.nodes   <- rbind(final.nodes, tail.nodes)
      final.weights <- c(final.weights, tail.weights)
      final.f       <- c(final.f, tail.f)

      m.target <- m.target - if (corners && k > 1) 2 + k else 2
      ai <- aj
      bi <- bj
    }
  }

  approx.error <- sum(errors)

  # Order before sending out
  for (d in seq_len(k)) {
    o <- order(final.nodes[, k])

    final.nodes   <- final.nodes[o, , drop = FALSE]
    final.weights <- final.weights[o]
    final.f       <- final.f[o]
  }
  
  list(
    n = final.nodes, 
    w = final.weights, 
    f = final.f, 
    a = a,
    b = b,
    k = 2, 
    m = NROW(final.nodes) ^ (1 / k), # may not be an integer...
    integral = sum(final.weights * final.f, na.rm=TRUE),
    error.abs = approx.error,
    error.rel = approx.error / integral
  )
}


buildEdgeSlabs <- function(ai, aj, bi, bj, corners = FALSE) {
  stopifnot(is.numeric(ai),  is.numeric(aj),  is.numeric(bi),  is.numeric(bj),
            length(ai) == length(aj),
            length(ai) == length(bi),
            length(ai) == length(bj),
            identical(dim(ai), NULL))

  k <- length(ai)
  if (k == 1L) corners <- FALSE

  n_slabs   <- 2L * k
  n_corners <- if (corners) 2L^k else 0L
  n_total   <- n_slabs + n_corners

  A <- matrix(NA_real_, n_total, k)
  B <- matrix(NA_real_, n_total, k)

  row_counter <- 0L
  add_row <- function(Ar, Br) {
    row_counter <<- row_counter + 1L          # advance the global counter
    A[row_counter, ] <<- Ar
    B[row_counter, ] <<- Br
  }

  for (d in seq_len(k)) {

    ## negative-end slab on axis d
    A_neg <- bj ; B_neg <- aj
    A_neg[d] <- ai[d]
    add_row(A_neg, B_neg)

    ## positive-end slab on axis d
    A_pos <- bj ; B_pos <- aj
    B_pos[d] <- bi[d]
    add_row(A_pos, B_pos)
  }

  if (corners) {
    sign_tbl <- as.matrix(expand.grid(replicate(k, c(-1, 1), simplify = FALSE)))
    for (s in seq_len(nrow(sign_tbl))) {

      signs <- sign_tbl[s, ]
      Ar    <- bj        # start from inner–positive, then overwrite extremes
      Br    <- aj

      neg_axes        <- which(signs == -1)
      pos_axes        <- which(signs ==  1)
      Ar[neg_axes]    <- ai[neg_axes]
      Br[pos_axes]    <- bi[pos_axes]

      add_row(Ar, Br)
    }
  }

  list(A = A, B = B)
}



adaptiveGaussQuadratureASymmetric <- function(fun, a = -7, b = 7, m.start = 4, m.max = 32,
                                              tol = 1e-6, k = 1, total.integral = NULL, ...) {
  stopif(k > 1, "Non-symmetric adaptive quadrature not implemented for multiple dimensions!\n",
         'Try `adaptive.quad=FALSE` or `adaptive.quad.type="symmetric" instead!')

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
  Ea <- stats::dnorm(a) * fun(matrix(a), ...)
  Eb <- stats::dnorm(b) * fun(matrix(b), ...)

  if (Ea < Eb) {
    a_left <- a; b_left <- c_mid; a_right <- c_mid; b_right <- b
  } else {
    a_left <- c_mid; b_left <- b; a_right <- a; b_right <- c_mid
  }

  left <- adaptiveGaussQuadratureASymmetric(
    fun = fun, a = a_left, b = b_left, m.start =  m.start, m.max = round(m.max / 2), 
    tol = tol, k = k, total.integral = total.integral, ...
  )
  right <- adaptiveGaussQuadratureASymmetric(
    fun = fun, a = a_right, b = b_right, m.start = m.start, m.max = m.max - left$m, 
    tol = tol, k = k, total.integral = total.integral, ...
  )

  list(integral = left$integral + right$integral, f = c(left$f, right$f), m = left$m + right$m,
       n = rbind(left$n, right$n), w = c(left$w, right$w), k = k)
}


# Not finished yet...
# adaptive_quadrature <- function(
#   log_integrand,          # function(z, data, …) -> scalar log f(z)
#   z, w,                   # prototype nodes (m×d) and weights (length m)
#   init     = NULL,        # optional optimiser starting point (length d)
#   grad     = numDeriv::grad,
#   hess     = numDeriv::hessian,
#   control  = list(),
#   ...
# ) {
#   ## ---- coerce shapes -------------------------------------------------
#   if (is.null(dim(z)))                # 1-D rule supplied as a vector
#     z <- matrix(z, ncol = 1)
# 
#   m <- nrow(z)
#   k <- ncol(z)
# 
#   stopifnot(length(w) == m)
# 
#   ## ---- 1. mode & Hessian --------------------------------------------
#   if (is.null(init))
#     init <- rep(0, k)
# 
#   opt <- optim(
#     par      = init,
#     fn       = function(theta) -log_integrand(theta, ...),
#     gr       = function(theta) -grad(log_integrand, theta, ...),
#     method   = "BFGS",
#     hessian  = TRUE,
#     control  = modifyList(list(reltol = 1e-12), control)
#   )
#   if (opt$convergence != 0)
#     warning("optim did not converge; results may be unreliable")
# 
#   mu <- opt$par
#   # H  <- -opt$hessian               # ∇² log f  (negative-definite at the mode)
#   # Hessian is already negative...
#   H  <- opt$hessian               # ∇² log f  (negative-definite at the mode)
# 
#   ## ---- 2. Laplace covariance & transform ----------------------------
#   Sigma <- solve(H)                # (-H)⁻¹
#   L     <- t(chol(Sigma))          # lower-triangular s.t. LLᵀ = Σ
#   detL  <- abs(det(L))             # |det L|  (positive by construction)
# 
#   ## ---- 3. Affine map each node --------------------------------------
#   z_new <- sweep(z %*% t(L), 2, mu, FUN = "+")   # m×d
#   w_new <- w * detL                              # length m
# 
#   list(n = z_new * sqrt(2),
#        w = w_new,
#        mode    = mu,
#        hessian = H,
#        k = k, m = m)
# }
# 
# 
# log_integrand_lms <- function(z, data, model, add_prior = TRUE, ...) {
# 
#   ## 1. Model-specific mean & covariance at this z ---------------------
#   mu    <- muLmsCpp(   model = model, z = z)
#   Sigma <- sigmaLmsCpp(model = model, z = z)
# 
#   ## 2. Contribution of the observed data ------------------------------
#   ##    (dmvnorm() returns a vector of log-densities, one per row)
#   ll_data <- mvtnorm::dmvnorm(data,
#                               mean  = mu,
#                               sigma = Sigma,
#                               log   = TRUE)
# 
#   ## 3. Optional prior on z (here: standard multivariate normal) -------
#   ll_prior <- if (add_prior) -0.5 * crossprod(z) else 0
# 
#   ## 4. Total log-integrand (scalar) -----------------------------------
#   sum(ll_data) + ll_prior          # this is what adaptive_quadrature() needs
# }
