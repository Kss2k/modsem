# Pairwise composite likelihood estimation.
#
# Same model as LMS, different objective: instead of integrating the full
# likelihood over every latent, sum log-probabilities over indicator PAIRS. See
# R/pml.R for why that is cheap and what it costs.
#
# There is no E/M loop -- the composite objective is minimised directly -- so
# this is much closer to `estQml` in shape than to `emLms`.

estPml <- function(model,
                   data = NULL,
                   ordered = NULL,
                   verbose = FALSE,
                   calc.se = TRUE,
                   max.iter = 500L,
                   optimizer = "nlminb",
                   epsilon = 1e-6,
                   FIM = "observed",
                   OFIM.hessian = FALSE,
                   EFIM.S = 3e4,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   R.max = 1e6,
                   cr1s = NULL,
                   nodes = 30L,
                   ...) {
  submodel <- model$models[[1L]]
  mod_stopif(model$info$n.groups > 1L,
             "PML does not support multiple groups yet.")
  mod_stopif(!length(submodel$info$ordered %||% character()),
             "PML currently requires ordered indicators.")

  ordered <- submodel$info$ordered
  tables <- pmlPairTables(data, ordered)

  # The conditioning set is whatever `sortXis()` marked nonlinear: one variable
  # per interaction term, chosen to cover them all. k = 0 is the linear model,
  # where the rule collapses to a single node and no integration happens.
  k <- submodel$quad$k %||% 0L
  rule <- pmlQuadRule(k, m = nodes)

  # `rows` maps each ordered indicator onto its row of lambdaX. Everything the
  # objective touches -- pairs, thresholds, the clean/integrated split -- is
  # indexed against `ordered`, and this is the only place that translates.
  rows <- match(ordered, rownames(submodel$matrices$lambdaX))
  mod_stopif(anyNA(rows), "Ordered indicators not found in the model matrices.")

  # From the UNFILLED matrices, so a free omega counts as nonzero even though it
  # starts at zero, and the split cannot shift underneath the optimiser.
  partition <- pmlPartition(submodel$matrices, submodel$info$numXis,
                            submodel$info$numEtas, tables$pairs, rows)

  # Everything the objective and gradient need, resolved once. The kernel takes
  # 0-based indices; everything above is 1-based.
  plan <- list(
    rule = rule, tables = tables, rows = rows,
    rows0 = rows - 1L, pairs0 = tables$pairs - 1L,
    hoisted0 = partition$hoisted - 1L,
    integrated0 = partition$integrated - 1L,
    locations = pmlParamLocations(model),
    thresholds = pmlThresholdLocations(model))

  # nlminb reconstructs a gradient from ~p+1 objective evaluations when it is
  # not given one, and each of those is a full pass over every pair, cell and
  # node. Supplying the analytic gradient collapses that to a single pass, so
  # both are computed together and cached on theta.
  cache <- new.env(parent = emptyenv())
  cache$theta <- NULL
  cache$last <- NULL

  evaluate <- function(theta) {
    if (!is.null(cache$theta) && identical(cache$theta, theta)) return(cache$value)
    value <- tryCatch(pmlObjectiveGradient(model, theta, plan, sign = -1),
                      error = function(e) NULL)
    if (is.null(value) || !is.finite(value$objective) ||
        anyNA(value$gradient) || any(!is.finite(value$gradient))) {
      # An infeasible point needs a gradient that POINTS BACK. Returning zero
      # tells nlminb the region is flat and it may settle there; reusing the
      # last good gradient keeps a descent direction pointing at feasible
      # ground, and the large objective makes the line search retreat.
      value <- list(objective = .Machine$double.xmax^(1 / 4),
                    gradient = cache$last %||% rep(0, length(theta)))
    } else {
      cache$last <- value$gradient
    }
    cache$theta <- theta
    cache$value <- value
    value
  }

  bounds <- pmlBounds(model$params$bounds, plan$locations)
  start <- model$theta
  fit <- suppressWarnings(stats::nlminb(
    start = start,
    objective = function(theta) evaluate(theta)$objective,
    gradient = function(theta) evaluate(theta)$gradient,
    lower = bounds$lower, upper = bounds$upper,
    control = list(iter.max = max.iter, eval.max = 2L * max.iter,
                   trace = if (isTRUE(verbose)) 1L else 0L)
  ))

  # Standard errors from the Hessian of a COMPOSITE likelihood are wrong: they
  # need the Godambe sandwich H^-1 J H^-1, and J is the covariance of dependent
  # pair scores. Until that exists, SEs are refused rather than reported
  # incorrectly. `logLik` is NA for the same reason -- the composite objective
  # is not a likelihood, so AIC/BIC built on it would not mean anything.
  mod_warnif(isTRUE(calc.se),
    paste0("PML standard errors require a Godambe sandwich, which is not ",
           "implemented; returning estimates without standard errors."))

  finalizeModelEstimatesDA(
    model = model, theta = fit$par, method = "pml", data = data,
    logLik = NA_real_, iterations = fit$iterations,
    converged = fit$convergence == 0L, optimizer = optimizer,
    calc.se = FALSE, FIM = FIM, OFIM.hessian = OFIM.hessian,
    EFIM.S = EFIM.S, EFIM.parametric = EFIM.parametric,
    robust.se = robust.se, epsilon = epsilon, cr1s = cr1s, R.max = R.max,
    verbose = verbose
  )
}
