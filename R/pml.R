# Pairwise composite likelihood (PML) for ordinal LMS models.
#
# An estimator in its own right, built on the LMS model: it reads the same DA
# matrices and the same interaction structure, but replaces the full-information
# likelihood with a sum over indicator PAIRS.
#
# WHY. Full information integrates over every latent variable, so the graph
# backend costs O(N * m^k * J) with k the number of latent variables. PML's
# integral runs only over the latents needed to linearise the interaction --
# usually one, at most two -- and for pairs of ordinal indicators the data
# reduces to contingency tables, so N leaves the per-evaluation cost entirely.
#
# THE STRUCTURE THAT MAKES IT WORK. `sortXis()` already chooses a set of
# "nonlinear" xis covering every interaction term (a minimum vertex cover of the
# interaction graph). Conditional on those, EVERY product term is linear in the
# remaining latents -- which is exactly what `kronZ.t() * Oxx` does in the legacy
# LMS kernel. So given a conditioning node the whole model is linear-Gaussian,
# the underlying variables x* are jointly normal, and each pair probability is a
# closed-form bivariate-normal rectangle.
#
#   P(x_j = c, x_k = d) = sum_q w_q [ Phi2(a_c, b_d) - Phi2(a_{c-1}, b_d)
#                                     - Phi2(a_c, b_{d-1}) + Phi2(a_{c-1}, b_{d-1}) ]
#
# PROBIT ONLY. Under a logit link the formulation is conditional independence
# given the latents, so a pair probability would still require integrating over
# them FOR EVERY PAIR -- strictly worse than full information. There is also no
# bivariate logistic consistent with a linear factor structure. Residual
# variances are therefore fixed at 1, not pi^2/3.
#
# WHAT IS VALIDATED. Only k = 0 (no interaction), against lavaan's
# `estimator = "PML"`. The conditioning path below is written for the general
# case but is NOT yet verified; see test_pml.R.


# Quadrature over the conditioning latents. k = 0 is the linear model: a single
# node of weight 1, so the node loop below runs once and the rule disappears.
pmlQuadRule <- function(k, m = 30L) {
  if (k < 1L) return(list(n = matrix(0, 1L, 0L), w = 1))
  single <- fastGHQuad::gaussHermiteData(m)
  one <- list(n = sqrt(2) * single$x, w = single$w * pi^(-1 / 2))
  nodes <- as.matrix(expand.grid(rep(list(one$n), k)))
  dimnames(nodes) <- NULL
  list(n = nodes, w = as.numeric(Reduce(kronecker, rep(list(one$w), k))))
}


# Mean and covariance of the full latent vector (xi, eta), conditional on the
# first `k` xis being held at `znode`.
#
# With those held fixed the structural equation is linear, so this is ordinary
# Gaussian conditioning followed by the usual reduced-form solve. At k = 0 the
# conditioning is vacuous and this returns the unconditional moments.
pmlLatentMoments <- function(M, numXis, numEtas, znode = numeric(0)) {
  k <- length(znode)
  Phi <- M$A %*% t(M$A)
  Binv <- solve(diag(numEtas) - M$gammaEta)
  muXi <- as.vector(M$beta0)
  Gx <- M$gammaXi
  alpha <- as.vector(M$alpha)
  covZetaXi <- if (is.null(M$covZetaXi) || !length(M$covZetaXi))
    matrix(0, numEtas, numXis) else M$covZetaXi

  if (k > 0L) {
    ci <- seq_len(k); fi <- setdiff(seq_len(numXis), ci)

    # The rule's nodes are standard normal INNOVATIONS, matching LMS's
    # xi = beta0 + A u with A lower triangular. Conditioning on u_1..k is
    # therefore conditioning on xi_1..k, and this is the xi value meant.
    znode <- muXi[ci] + as.vector(M$A[ci, ci, drop = FALSE] %*% znode)

    # The product terms become linear in the free xis once the conditioning
    # ones are fixed: this is `kronZ' Oxx` from the legacy kernel.
    kronZ <- Reduce(kronecker, rep(list(znode), 1L))
    Gx <- Gx + matrix(kronZ, nrow = 1L) %*% M$omegaXiXi[seq_len(k), , drop = FALSE]

    Scc <- Phi[ci, ci, drop = FALSE]
    Sfc <- Phi[fi, ci, drop = FALSE]
    solved <- Sfc %*% solve(Scc)
    muFree <- muXi[fi] + as.vector(solved %*% (znode - muXi[ci]))
    covFree <- Phi[fi, fi, drop = FALSE] - solved %*% t(Sfc)

    # Re-embed at full width, with the conditioned block degenerate.
    muXi <- replace(muXi, fi, muFree); muXi[ci] <- znode
    Phi <- matrix(0, numXis, numXis)
    Phi[fi, fi] <- covFree
  }

  muEta <- as.vector(Binv %*% (alpha + Gx %*% muXi))
  covXiEta <- Phi %*% t(Gx) %*% t(Binv) + t(covZetaXi) %*% t(Binv)
  covEta <- Binv %*% (Gx %*% Phi %*% t(Gx) + Gx %*% t(covZetaXi) +
                        covZetaXi %*% t(Gx) + M$psi) %*% t(Binv)

  list(mean = c(muXi, muEta),
       cov = rbind(cbind(Phi, covXiEta), cbind(t(covXiEta), covEta)))
}


# Underlying-variable moments. Every indicator routes through lambdaX/tauX,
# whose columns span the xis followed by the etas.
pmlImpliedMoments <- function(M, latent) {
  mean <- as.vector(M$tauX) + as.vector(M$lambdaX %*% latent$mean)
  cov <- M$lambdaX %*% latent$cov %*% t(M$lambdaX) + M$thetaDelta
  list(mean = mean, sd = sqrt(diag(cov)), cov = cov)
}


# Bivariate normal CDF tolerant of infinite arguments, which the outer rows and
# columns of every rectangle need.
pmlBivariate <- function(a, b, rho) {
  # `pbivnorm` is a stopgap: it is an R-level call at ~240ns, and the intended
  # implementation is a Drezner-Wesolowsky routine in C++ where rho is constant
  # across a pair's corners and the setup amortises. Suggested, not imported,
  # until that lands.
  mod_stopif(!requireNamespace("pbivnorm", quietly = TRUE),
             "PML requires the `pbivnorm` package.")
  out <- numeric(length(a))
  lowA <- is.infinite(a) & a < 0
  lowB <- is.infinite(b) & b < 0
  highA <- is.infinite(a) & a > 0
  highB <- is.infinite(b) & b > 0
  out[highA & highB] <- 1
  use <- highA & !highB & !lowB
  out[use] <- stats::pnorm(b[use])
  use <- highB & !highA & !lowA
  out[use] <- stats::pnorm(a[use])
  finite <- !lowA & !lowB & !highA & !highB
  if (any(finite))
    out[finite] <- pbivnorm::pbivnorm(a[finite], b[finite], rho = rho)
  out[lowA | lowB] <- 0
  out
}


# Cell probabilities for one pair, at one node. Every cell reads from the same
# grid of corners, so the corners are computed once and differenced.
pmlPairProbabilities <- function(implied, thresholds, j, k) {
  a <- c(-Inf, (thresholds[[j]] - implied$mean[j]) / implied$sd[j], Inf)
  b <- c(-Inf, (thresholds[[k]] - implied$mean[k]) / implied$sd[k], Inf)
  rho <- implied$cov[j, k] / (implied$sd[j] * implied$sd[k])
  grid <- expand.grid(a = a, b = b)
  corner <- matrix(pmlBivariate(grid$a, grid$b, rho), length(a), length(b))
  corner[-1L, -1L] - corner[-nrow(corner), -1L] -
    corner[-1L, -ncol(corner)] + corner[-nrow(corner), -ncol(corner)]
}


# Bivariate contingency tables, formed once. This is what removes N from the
# per-evaluation cost: the objective only ever sees counts.
pmlPairTables <- function(data, ordered, categories) {
  codes <- vapply(data[ordered], as.integer, integer(NROW(data)))
  pairs <- t(utils::combn(length(ordered), 2L))
  tables <- lapply(seq_len(NROW(pairs)), function(p) {
    j <- pairs[p, 1L]; k <- pairs[p, 2L]
    table(factor(codes[, j], levels = seq_len(categories[[j]])),
          factor(codes[, k], levels = seq_len(categories[[k]])))
  })
  list(pairs = pairs, tables = tables)
}


# Which etas carry an interaction, directly or through a structural path.
#
# omegaXiXi and omegaEtaXi are stacked per-eta blocks of `numXis` rows each, the
# layout `kron(Ieta, muXi)' Omega` reads. An eta is DIRECTLY affected if its own
# block is nonzero, and affectedness then propagates forwards along gammaEta:
# anything regressed on an affected eta is itself affected.
#
# The criterion is "downstream of an interaction", NOT "exogenous". In the TPB
# model `INT ~ ATT + SN + PBC` is linear, so INT is clean even though it is
# endogenous; and if `BEH ~ INT` with an interaction in INT's equation, BEH is
# affected even though it carries no omega of its own.
pmlAffectedEtas <- function(M, numXis, numEtas) {
  if (numEtas < 1L) return(logical(0L))
  nonzeroBlock <- function(Omega, i) {
    if (is.null(Omega) || !length(Omega)) return(FALSE)
    rows <- (i - 1L) * numXis + seq_len(numXis)
    any(Omega[rows, , drop = FALSE] != 0)
  }
  affected <- vapply(seq_len(numEtas), function(i)
    nonzeroBlock(M$omegaXiXi, i) || nonzeroBlock(M$omegaEtaXi, i), logical(1L))

  repeat {
    # unname(): gammaEta carries dimnames, and a named/unnamed mismatch would
    # make `identical()` never fire.
    grown <- affected |
      unname(rowSums(abs(M$gammaEta[, affected, drop = FALSE])) > 0)
    if (identical(grown, affected)) break
    affected <- grown
  }
  affected
}


# Indicators that load on no affected eta. Their underlying variables are
# unconditionally normal, so a pair of them needs no quadrature.
pmlCleanIndicators <- function(M, numXis, numEtas) {
  affected <- pmlAffectedEtas(M, numXis, numEtas)
  if (!any(affected)) return(rep(TRUE, NROW(M$lambdaX)))
  unname(rowSums(abs(M$lambdaX[, numXis + which(affected), drop = FALSE])) == 0)
}


# The pairwise composite log-likelihood.
#
# Conditional on a node the model is linear-Gaussian, so EVERY pair is a
# bivariate-normal rectangle there. But a pair whose indicators are untouched by
# any interaction is bivariate normal UNCONDITIONALLY as well -- mixing those
# conditional normals back over the node distribution returns exactly the
# unconditional pair. Running such a pair through the quadrature is therefore not
# only wasted work but strictly worse ACCURACY: it approximates a quantity the
# closed form already delivers exactly. Those pairs are computed once, outside
# the loop.
#
# `sign = -1` gives a minimisation objective, matching the DA convention.
pmlObjective <- function(M, numXis, numEtas, thresholds, tables, rule,
                         sign = -1) {
  pairs <- tables$pairs
  clean <- pmlCleanIndicators(M, numXis, numEtas)
  hoisted <- which(clean[pairs[, 1L]] & clean[pairs[, 2L]])
  integrated <- setdiff(seq_len(NROW(pairs)), hoisted)
  probability <- vector("list", NROW(pairs))

  if (length(hoisted)) {
    implied <- pmlImpliedMoments(M, pmlLatentMoments(M, numXis, numEtas))
    for (p in hoisted)
      probability[[p]] <- pmlPairProbabilities(implied, thresholds,
                                               pairs[p, 1L], pairs[p, 2L])
  }

  if (length(integrated)) for (q in seq_len(NROW(rule$n))) {
    latent <- pmlLatentMoments(M, numXis, numEtas, as.vector(rule$n[q, ]))
    implied <- pmlImpliedMoments(M, latent)
    for (p in integrated) {
      cell <- rule$w[[q]] *
        pmlPairProbabilities(implied, thresholds, pairs[p, 1L], pairs[p, 2L])
      probability[[p]] <- if (q == 1L) cell else probability[[p]] + cell
    }
  }

  total <- 0
  for (p in seq_len(NROW(pairs)))
    total <- total + sum(tables$tables[[p]] * log(pmax(probability[[p]], 1e-12)))
  sign * total
}


# Probit identification for ordered indicators, and the threshold matrix.
#
# The scale of an underlying variable is not identified alongside free
# thresholds, so it has to be fixed. PML uses the THETA parameterisation, which
# is also what lavaan's `parameterization = "theta"` does:
#   * indicator intercept fixed at 0,
#   * residual variance fixed at 1 (probit; there is no logit option -- see the
#     header),
#   * every threshold free, parameterised through `thresholdDelta` so ordering
#     holds by construction.
#
# Starting values come from the observed cumulative proportions via
# `orderedThresholdSpec()`.
pmlPrepareOrdered <- function(model, data, ordered) {
  if (!length(ordered)) return(model)
  info <- orderedThresholdSpec(data, ordered, link = "probit")
  starts <- numeric()

  for (g in seq_along(model$models)) {
    submodel <- model$models[[g]]
    indicators <- rownames(submodel$matrices$lambdaX)
    index <- match(intersect(ordered, indicators), indicators)
    if (!length(index)) next

    submodel$matrices$tauX[index, 1L] <- 0
    submodel$matrices$thetaDelta[index, ] <- 0
    submodel$matrices$thetaDelta[, index] <- 0
    diag(submodel$matrices$thetaDelta)[index] <- 1
    submodel$labelMatrices$tauX[index, 1L] <- ""
    submodel$labelMatrices$thetaDelta[index, ] <- ""
    submodel$labelMatrices$thetaDelta[, index] <- ""

    specs <- Filter(function(x) identical(x$group, g), info$specs)
    width <- max(c(0L, vapply(specs, function(x) x$K - 1L, integer(1L))))
    delta <- matrix(NaN, length(indicators), width,
                    dimnames = list(indicators, paste0("t", seq_len(width))))
    for (spec in specs)
      delta[spec$variable, seq_len(spec$K - 1L)] <- info$delta[spec$index]

    submodel$matrices$thresholdDelta <- delta
    submodel$matrices$thresholds <- thresholdDeltaToThresholdMatrix(delta)
    submodel$labelMatrices$thresholdDelta <- matrix(
      "", NROW(delta), NCOL(delta), dimnames = dimnames(delta))
    submodel$labelMatrices$thresholds <- submodel$labelMatrices$thresholdDelta

    free <- is.finite(delta)
    values <- delta[free]
    names(values) <- getParamNamesMatrix(delta, "thresholdDelta")[free]
    if (g > 1L) names(values) <- sprintf("%s.g%d", names(values), g)
    starts <- c(starts, values)
    submodel$matrices$thresholdDelta[free] <- NA_real_
    submodel$info$ordered <- intersect(ordered, indicators)
    model$models[[g]] <- submodel
  }

  old <- model$theta
  params <- createTheta(model, parTable.in = model$parTable)
  shared <- intersect(names(old), names(params$theta))
  params$theta[shared] <- old[shared]
  params$theta[names(starts)] <- starts
  model$params[names(params)] <- params
  model$theta <- params$theta
  model$params$bounds <- getParamBounds(model)
  model$info$ordered <- ordered
  model
}


# Thresholds as a per-indicator list, dropping the NaN padding that ragged rows
# carry. Indicators with no thresholds (continuous) get `numeric(0)`.
pmlThresholdList <- function(M) {
  if (is.null(M$thresholds) || !NCOL(M$thresholds))
    return(rep(list(numeric(0)), NROW(M$lambdaX)))
  lapply(seq_len(NROW(M$thresholds)), function(i) {
    row <- M$thresholds[i, ]
    as.numeric(row[is.finite(row)])
  })
}
