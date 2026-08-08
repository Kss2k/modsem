lmsGraphLaplaceEvaluation <- function(model, theta,
                                      order = c("laplace", "laplace2"),
                                      previous = NULL,
                                      link = c("logit", "probit"),
                                      tolerance = 1e-8,
                                      correction.floor = 1e-8) {
  order <- match.arg(order)
  link <- match.arg(link)
  filled <- fillModel(model = model, theta = theta, method = "lms")
  groups <- vector("list", model$info$n.groups)
  total <- 0
  individual <- numeric()

  for (g in seq_len(model$info$n.groups)) {
    submodel <- filled$models[[g]]
    submodel$quad <- model$models[[g]]$quad
    submodel$quad$m <- 1L
    submodel$quad$adaptive <- TRUE
    old <- if (!is.null(previous)) previous[[g]] else NULL
    workspace <- old$workspace %||% lmsGraphWorkspace(submodel)
    rule <- suppressWarnings(lmsGraphAdaptiveRule(
      submodel, previous = old, link = link, tolerance = tolerance,
      workspace = workspace
    ))
    M <- submodel$matrices
    if (is.null(M$thresholds))
      M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
    ordered <- match(submodel$info$ordered %||% character(),
                     rownames(M$lambdaX)) - 1L
    ordered <- ordered[!is.na(ordered)]
    sampling.weights <- submodel$data$weights
    if (is.null(sampling.weights)) sampling.weights <- rep(1, submodel$data$n)
    if (order == "laplace2") {
      approximation <- lmsGraphLaplace2Cpp(
        M, rule$modes, submodel$info$numXis, submodel$info$numEtas,
        submodel$data$data.split, submodel$data$colidx0, ordered,
        identical(link, "logit"), 1e-6, correction.floor,
        ThreadEnv$n.threads %||% 1L
      )
      contribution <- approximation$second
    } else {
      aggregate <- lmsGraphPstepWorkspaceCpp(
        M, rule$n, submodel$info$numXis, submodel$info$numEtas,
        workspace, rule$w, sampling.weights, identical(link, "logit"),
        ThreadEnv$n.threads %||% 1L
      )
      contribution <- aggregate$logDensity
      approximation <- list(
        first = contribution, second = contribution,
        correction = rep(0, length(contribution)),
        adjusted = rule$curvature.adjusted
      )
    }
    value <- sum(sampling.weights * contribution)
    individual <- c(individual, sampling.weights * contribution)
    total <- total + value
    rule$workspace <- workspace
    groups[[g]] <- list(
      quad = rule, approximation = approximation, logLik = value,
      workspace = workspace
    )
  }
  list(logLik = total, individual = individual, groups = groups,
       quad = lapply(groups, `[[`, "quad"))
}


lmsGraphLaplaceInformation <- function(theta, model, P,
                                       epsilon = 1e-5) {
  integration <- P$integration %||% "laplace2"
  link <- P$P_GROUPS[[1L]]$link %||% "logit"
  base <- lmsGraphLaplaceEvaluation(
    model, theta, integration, P$quad, link
  )
  step <- epsilon * pmax(1, abs(theta))
  scores <- matrix(0, length(base$individual), length(theta))
  for (j in seq_along(theta)) {
    plus <- theta
    plus[j] <- plus[j] + step[j]
    shifted <- lmsGraphLaplaceEvaluation(
      model, plus, integration, base$quad, link
    )
    scores[, j] <- (shifted$individual - base$individual) / step[j]
  }
  information <- crossprod(scores)
  dimnames(information) <- list(names(theta), names(theta))
  information
}


estLmsGraphLaplace <- function(model,
                               order = c("laplace", "laplace2"),
                               link = c("logit", "probit"),
                               optimizer = c("nlminb", "L-BFGS-B"),
                               max.iter = 500L,
                               convergence.rel = 1e-8,
                               calc.se = FALSE,
                               verbose = FALSE,
                               adaptive.quad.tol = 1e-8,
                               ...) {
  order <- match.arg(order)
  link <- match.arg(link)
  optimizer <- match.arg(optimizer)
  bounds <- model$params$bounds %||% getParamBounds(model)
  initial <- lmsGraphLaplaceEvaluation(
    model, model$theta, order, NULL, link, adaptive.quad.tol
  )
  previous <- initial$quad
  evaluations <- 0L
  objective <- function(theta) {
    value <- tryCatch(
      lmsGraphLaplaceEvaluation(
        model, theta, order, previous, link, adaptive.quad.tol
      ), error = identity
    )
    evaluations <<- evaluations + 1L
    if (inherits(value, "error") || !is.finite(value$logLik))
      return(.Machine$double.xmax / 1e100)
    previous <<- value$quad
    -value$logLik
  }
  start.time <- proc.time()[["elapsed"]]
  if (optimizer == "nlminb") {
    fit <- stats::nlminb(
      start = model$theta, objective = objective,
      lower = bounds$lower, upper = bounds$upper,
      control = list(iter.max = as.integer(max.iter),
                     eval.max = max(1000L, 10L * as.integer(max.iter)),
                     rel.tol = convergence.rel)
    )
    theta <- fit$par
    iterations <- fit$iterations
    converged <- fit$convergence == 0L
    message <- fit$message
  } else {
    fit <- stats::optim(
      par = model$theta, fn = objective, method = "L-BFGS-B",
      lower = bounds$lower, upper = bounds$upper,
      control = list(maxit = as.integer(max.iter), factr = max(
        1, convergence.rel / .Machine$double.eps
      ))
    )
    theta <- fit$par
    iterations <- fit$counts[["function"]]
    converged <- fit$convergence == 0L
    message <- fit$message %||% ""
  }
  final <- lmsGraphLaplaceEvaluation(
    model, theta, order, previous, link, adaptive.quad.tol
  )
  information.P <- if (isTRUE(calc.se)) pstepLmsGraph(
    model, theta, lastQuad = final$quad, recalcQuad = FALSE,
    link = link, adaptive.quad.tol = adaptive.quad.tol,
    integration = order
  ) else NULL
  out <- finalizeModelEstimatesDA(
    model = model, theta = theta, method = "lms",
    data = lapply(model$models, `[[`, "data"), logLik = final$logLik,
    iterations = iterations, converged = converged,
    optimizer = paste0("DIRECT-", optimizer, "-", toupper(order)),
    calc.se = calc.se, FIM = "observed", OFIM.hessian = TRUE,
    EFIM.S = 0, EFIM.parametric = FALSE, robust.se = FALSE,
    epsilon = 1e-6, cr1s = FALSE, R.max = 0, verbose = verbose,
    P = information.P, includeStartModel = TRUE, startModel = model
  )
  out$integration <- order
  out$optimizer.message <- message
  out$function.evaluations <- evaluations
  out$elapsed.time <- proc.time()[["elapsed"]] - start.time
  out$laplace <- list(
    groups = lapply(final$groups, function(group) list(
      modes = group$quad$modes,
      approximation = group$approximation,
      logLik = group$logLik
    ))
  )
  out
}
