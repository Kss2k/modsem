# -----------------------------------------------------------------------------
#  Adaptive LMS utilities and EM helpers (R)                                   
#  ---------------------------------------------------------------------------
#  This source file augments the original fixed‑grid LMS implementation with  
#  mean–variance adaptive Gauss–Hermite quadrature (AGHQ) à la Liu & Pierce   
#  (1994) – the same strategy used by Mplus when EM‑estimating latent‑variable
#  models.  All new code lives in clearly‑separated helper functions so the    
#  original C++ back‑end (muLmsCpp, sigmaLmsCpp, …) can be reused unchanged.  
# -----------------------------------------------------------------------------

#' Null‑coalescing operator (internal)
`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# 1.  Standard tensor‑product Gauss–Hermite rule ------------------------------
# -----------------------------------------------------------------------------
#  Returns a list with
#    $S … M×k matrix of standard GH nodes  (rows = nodes, cols = dimensions)
#    $w … length‑M vector of probability weights  (sum = 1)
#  where M = m^k and m is the number of nodes per dimension.
# -----------------------------------------------------------------------------

gh_rule <- function(k, m = 5L) {
  if (k < 1L)        stop("'k' must be ≥ 1")
  if (m < 1L)        stop("'m' must be ≥ 1")
  r1 <- statmod::gauss.quad.prob(m, dist = "normal")  # N(0,1) – weights sum to 1
  grid <- as.matrix(do.call(expand.grid, replicate(k, r1$nodes, simplify = FALSE)))
  wMat <- as.matrix(do.call(expand.grid, replicate(k, r1$weights, simplify = FALSE)))
  list(S = grid, w = apply(wMat, 1L, prod))            # length‑M
}

# -----------------------------------------------------------------------------
# 2.  Case‑specific AGHQ transformation ---------------------------------------
# -----------------------------------------------------------------------------
#  aghq_adapt(S, w, log_f)
#  • S, w  … standard nodes & weights from gh_rule()
#  • log_f … function(z) returning log f(y_i, z | θ)
#  Returns a list with adaptive nodes Z (same dim as S) and weights w_new.
# -----------------------------------------------------------------------------

aghq_adapt <- function(S, w, log_f,
                       tol = 1e-8, maxit = 50L, attempt_PD = TRUE) {
  k   <- ncol(S)
  # 2a.  Posterior mode via optimisation of –log f
  opt <- optim(rep(0, k), fn = function(z) -log_f(z), method = "BFGS",
               hessian = TRUE, control = list(reltol = tol, maxit = maxit))
  zhat <- opt$par
  H    <- opt$hessian

  # 2b.  Cholesky of (H^{-1} / 2)
  Sigma <- tryCatch(solve(H) / 2, error = function(e) NULL)
  if (is.null(Sigma) && attempt_PD) {
    # Numerical quirks: pull towards nearest PD matrix
    Sigma <- as.matrix(Matrix::nearPD(solve(H) / 2)$mat)
  }
  L  <- chol(Sigma, pivot = TRUE)

  # 2c.  Transform nodes and weights
  Z  <- sweep(S %*% t(L), 2L, zhat, `+`)        # M×k adaptive nodes
  w_new <- as.numeric(w * det(L))               # length‑M  (sum != 1)
  list(Z = Z, w = w_new)
}


estepLms <- function(model, theta, data, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  adaptive  <- isTRUE(model$quad$adaptive)
  k         <- model$quad$k
  m         <- model$quad$m %||% 5L   # sensible default

  if (adaptive) {
    stdGH <- gh_rule(k, m)            # reused for all cases
    S     <- stdGH$S
    w_std <- stdGH$w
  } else {
    V     <- as.matrix(model$quad$n)
    w_std <- model$quad$w
  }

  N <- nrow(data)
  P_list <- vector("list", N)        # posterior weights per case
  V_list <- vector("list", N)        # adaptive nodes per case
  w_list <- vector("list", N)        # adaptive weights per case
  density <- numeric(N)               # f(y_i)

  # Per‑case loop (can be parallelised)
  cat("\nCalculating Grid!...\n")
  for (i in seq_len(N)) {
    if (adaptive) {
      # log f(y_i, z | θ)
      log_fzi <- function(z) {
        mu    <- muLmsCpp   (model = modFilled, z = z)
        sigma <- sigmaLmsCpp(model = modFilled, z = z)
        mvtnorm::dmvnorm(data[i, ], mean = mu, sigma = sigma, log = TRUE) +
          sum(dnorm(z, log = TRUE))                         # N(0, I) prior
      }
      ad   <- aghq_adapt(S, w_std, log_fzi)
      Z_i  <- ad$Z
      w_i  <- ad$w
    } else {
      Z_i <- V
      w_i <- w_std
    }

    m_i <- nrow(Z_i)
    p_i <- numeric(m_i)
    for (j in seq_len(m_i)) {
      mu_ij    <- muLmsCpp   (model = modFilled, z = Z_i[j, ])
      sigma_ij <- sigmaLmsCpp(model = modFilled, z = Z_i[j, ])
      p_i[j]   <- mvtnorm::dmvnorm(data[i, ], mean = mu_ij, sigma = sigma_ij) *
                  w_i[j]
    }
    density[i] <- sum(p_i)
    P_list[[i]] <- p_i / density[i]   # normalised posterior probs
    V_list[[i]] <- Z_i
    w_list[[i]] <- w_i
  }

  observedLogLik <- sum(log(density))

  # -------------------------------------------------------------------------
  #  Weighted moments  (needed by M‑step)  ----------------------------------
  # -------------------------------------------------------------------------
  wMeans <- wCovs <- tGamma <- vector("list", length = N)
  for (i in seq_len(N)) {
    p    <- P_list[[i]]                # length m_i
    Z_i  <- V_list[[i]]
    tGamma[[i]] <- sum(p)
    wMeans[[i]] <- colSums(Z_i * p) / tGamma[[i]]
    Xc   <- sweep(Z_i, 2L, wMeans[[i]], `-`)
    wCovs[[i]] <- t(Xc) %*% (Xc * p) / tGamma[[i]]
  }

  cacheXptr <- QuadGrid2XPtr(list(P = P_list, V = V_list, w = w_list))

  list(P = P_list, mean = wMeans, cov = wCovs, tgamma = tGamma,
       V = V_list, w = w_list, obsLL = observedLogLik,
       adaptive = adaptive, xptrQuadCache = cacheXptr)
}


logLikLms <- function(theta, model, data, P, sign = -1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  sign * completeLogLikLmsCpp(modelR=modFilled, data=data, xpCache=P$xptrQuadCache)
} 


obsLogLikLms <- function(theta, model, data, P, sign = 1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")
  sign * obsLogLikLmsCpp(modelR=modFilled, data=data, quadR=P)
}
