# Script for testing the accuracy of different quadrature integration
# specifications
devtools::load_all()

matrix3 <- \(...) matrix(c(...), ncol = 3, byrow = TRUE)

Psi <- matrix3(
  1.2, 0.7, 0.2,
  0.7, 1.8, 0.4,
  0.2, 0.4, 0.9
)

A <- chol(Psi)

Lambda <- matrix3(
  0.9, 0.0, 0.0,
  1.2, 0.0, 0.0,
  0.7, 0.0, 0.0,
  0.0, 1.4, 0.0,
  0.0, 2.1, 0.0,
  0.0, 0.9, 0.0,
  0.0, 0.0, 0.4,
  0.0, 0.0, 0.9,
  0.0, 0.0, 0.7
)

Theta <- diag(c(
  0.4, 0.2, 0.4,
  0.2, 0.9, 1.1,
  0.2, 0.4, 0.5
))


scale <- 0.5
ThetaScaled <- Theta * scale


Sigma <- Lambda %*% Psi %*% t(Lambda) + ThetaScaled

n <- 2000

set.seed(293482)
X <- mvtnorm::rmvnorm(n, mean = rep(0, 9), sigma = Sigma)

density <- function(z) {
  if (!is.matrix(z))
    z <- matrix(z, nrow = 1L)

  D <- matrix(0, nrow = nrow(X), ncol = nrow(z))
  Z <- (z %*% A) %*% t(Lambda)

  for (i in seq_len(nrow(Z)))
    D[,i] <- mvtnorm::dmvnorm(X, mean = Z[i, ], sigma = Theta, log = FALSE)

  D 
}


sumIntegral <- function(densities, w = NULL) {
  if (!is.matrix(densities))
    densities <- matrix(densities, ncol = 1L)

  if (is.null(w))
    w <- rep(1, ncol(densities))

  sum(log(densities %*% w))
}


logIntegral <- function(z, w) {
  sumIntegral(density(z), w)
}


# Examples of the integral using fixed Gauss-Hermite quadrature
gh05  <- quadRuleFixed(quadSpec(nodes = 5, k = 3, integration = "gh"))
gh10  <- quadRuleFixed(quadSpec(nodes = 10, k = 3, integration = "gh"))
gh15  <- quadRuleFixed(quadSpec(nodes = 15, k = 3, integration = "gh"))
gh20  <- quadRuleFixed(quadSpec(nodes = 20, k = 3, integration = "gh", max.nodes = Inf))
gh100 <- quadRuleFixed(quadSpec(nodes = 100, k = 3, integration = "gh", max.nodes = Inf))
gh200 <- quadRuleFixed(quadSpec(nodes = 200, k = 3, integration = "gh", max.nodes = Inf))

# really slow but it should serve as an initial "true" value
keep <- apply(gh200$n, MARGIN = 1L, FUN = \(x) all(abs(x) <= 5.5)) 
logIntegral(gh200$n[keep,], gh200$w[keep])
true <- sum(mvtnorm::dmvnorm(X, sigma = Sigma, log = TRUE))
#> [1] -24622.6

benchmarkQuadRules <- function(...) {
  # we should probably add some code for benchmarking how long it takes to
  # "adapt" the different adaptive quadratures as well
 
  quadRules <- list(...)
  
  n <- length(quadRules)
  out <- data.frame(
    rule  = names(quadRules),
    ll    = numeric(n),
    err   = numeric(n),
    sign  = integer(n),
    nodes = integer(n),
    time  = numeric(n) # seconds
  )

  for (i in seq_along(quadRules)) {
    rule <- quadRules[[i]]

    time <- system.time({
      ll <- logIntegral(rule$n, rule$w)
    })

    out$ll[i]    <- ll
    out$err[i]   <- abs(true - ll)
    out$sign[i]  <- sign(true - ll)
    out$nodes[i] <- length(rule$w)
    out$time[i]  <- time[[1L]] # user time
  }

  out
}


# The adaptive rules need the density to adapt to. `buildQuadRule()` passes
# `...` straight through to it, and expects an N-by-Q matrix, which is exactly
# what `density()` above returns.
adaptTo <- function(nodes, ...) density(nodes)

rule <- function(integration, adaptive, nodes, k = 3, ...) {
  spec <- quadSpec(integration = integration, adaptive = adaptive,
                   nodes = nodes, k = k, max.nodes = Inf, ...)
  buildQuadRule(spec, density = adaptTo)
}


# Common rules only: `adaptive = "full"` builds one rule per observation and so
# cannot be evaluated as a single node set here. Extending the harness to cover
# it means integrating per row of X rather than once for the whole sample.
benchmarkQuadRules(
  `GH   fixed  m=5`     = rule("gh",   "none",  5),
  `GH   fixed  m=10`    = rule("gh",   "none",  10),
  `GH   fixed  m=15`    = rule("gh",   "none",  15),
  `GH   fixed  m=20`    = rule("gh",   "none",  20),
  `GH   fixed  m=24`    = rule("gh",   "none",  24),

  `GH   quasi  m=10`    = rule("gh",   "quasi", 10),
  `GH   quasi  m=15`    = rule("gh",   "quasi", 15),
  `GH   quasi  m=24`    = rule("gh",   "quasi", 24),

  `Rect fixed  m=15 r5` = rule("rect", "none",  15, rect.range = 5),
  `Rect fixed  m=15 r3` = rule("rect", "none",  15, rect.range = 3),
  `Rect fixed  m=24 r5` = rule("rect", "none",  24, rect.range = 5),

  `Rect quasi  m=15`    = rule("rect", "quasi", 15),
  `Rect quasi  m=24`    = rule("rect", "quasi", 24),

  `MC   fixed  n=500`   = rule("mc",   "none",  500),
  `MC   fixed  n=2000`  = rule("mc",   "none",  2000),
  `MC   fixed  n=8000`  = rule("mc",   "none",  8000)
)
