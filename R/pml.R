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


# Underlying-variable moments, restricted to the ordered indicators.
#
# The moments themselves come from `LMSModel::muSigma()` in src/pml.cpp -- the
# same kernel the full-information likelihood uses -- so there is no second
# definition of the model's conditional distribution to drift out of step.
# `rows` maps the ordered indicators onto rows of lambdaX, which puts everything
# downstream (pairs, thresholds, the clean/integrated split) in one index space.
pmlSelectMoments <- function(moments, rows) {
  cov <- moments$cov[rows, rows, drop = FALSE]
  list(mean = moments$mean[rows], cov = cov, sd = sqrt(diag(cov)))
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
#
# This reads STRUCTURE, not values, so it must be given the UNFILLED matrices:
# a free entry is NA there and counts as nonzero. Reading a filled model would
# call everything clean whenever omega starts at zero, which is exactly where
# the optimiser starts.
pmlAffectedEtas <- function(M, numXis, numEtas) {
  if (numEtas < 1L) return(logical(0L))
  nonzero <- function(x) is.na(x) | x != 0
  nonzeroBlock <- function(Omega, i) {
    if (is.null(Omega) || !length(Omega)) return(FALSE)
    rows <- (i - 1L) * numXis + seq_len(numXis)
    any(nonzero(Omega[rows, , drop = FALSE]))
  }
  affected <- vapply(seq_len(numEtas), function(i)
    nonzeroBlock(M$omegaXiXi, i) || nonzeroBlock(M$omegaEtaXi, i), logical(1L))

  repeat {
    # unname(): gammaEta carries dimnames, and a named/unnamed mismatch would
    # make `identical()` never fire.
    grown <- affected | unname(rowSums(
      nonzero(M$gammaEta[, affected, drop = FALSE])) > 0)
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
  loadings <- M$lambdaX[, numXis + which(affected), drop = FALSE]
  unname(rowSums(is.na(loadings) | loadings != 0) == 0)
}


# Split the pairs into those the node loop has to see and those it does not.
# Computed once, from the unfilled model, and held fixed for the whole fit.
pmlPartition <- function(M, numXis, numEtas, pairs, rows) {
  clean <- pmlCleanIndicators(M, numXis, numEtas)[rows]
  hoisted <- which(clean[pairs[, 1L]] & clean[pairs[, 2L]])
  list(hoisted = hoisted,
       integrated = setdiff(seq_len(NROW(pairs)), hoisted))
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
pmlObjective <- function(submodel, thresholds, tables, rule, partition, rows,
                         sign = -1) {
  pairs <- tables$pairs
  probability <- vector("list", NROW(pairs))

  if (length(partition$hoisted)) {
    implied <- pmlSelectMoments(pmlUnconditionalMomentsCpp(submodel), rows)
    for (p in partition$hoisted)
      probability[[p]] <- pmlPairProbabilities(implied, thresholds,
                                               pairs[p, 1L], pairs[p, 2L])
  }

  if (length(partition$integrated)) {
    moments <- pmlMomentsCpp(submodel, rule$n)
    for (q in seq_along(moments)) {
      implied <- pmlSelectMoments(moments[[q]], rows)
      for (p in partition$integrated) {
        cell <- rule$w[[q]] *
          pmlPairProbabilities(implied, thresholds, pairs[p, 1L], pairs[p, 2L])
        probability[[p]] <- if (q == 1L) cell else probability[[p]] + cell
      }
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
