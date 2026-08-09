lmsGraphLaplaceEvaluation <- function(model, theta,
                                      order = c("laplace", "laplace2"),
                                      previous = NULL,
                                      link = c("logit", "probit"),
                                      tolerance = 1e-8,
                                      correction.floor = 1e-8,
                                      fixed.modes = NULL) {
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
    if (is.null(fixed.modes)) {
      rule <- suppressWarnings(lmsGraphAdaptiveRule(
        submodel, previous = old, link = link, tolerance = tolerance,
        workspace = workspace
      ))
    } else {
      rule <- old
      rule$modes <- fixed.modes[[g]]
    }
    M <- submodel$matrices
    if (is.null(M$thresholds))
      M$thresholds <- matrix(NaN, NROW(M$lambdaX), 0L)
    ordered <- match(submodel$info$ordered %||% character(),
                     rownames(M$lambdaX)) - 1L
    ordered <- ordered[!is.na(ordered)]
    sampling.weights <- submodel$data$weights
    if (is.null(sampling.weights)) sampling.weights <- rep(1, submodel$data$n)
    if (order == "laplace2" || !is.null(fixed.modes)) {
      approximation <- lmsGraphLaplace2Cpp(
        M, rule$modes, submodel$info$numXis, submodel$info$numEtas,
        submodel$data$data.split, submodel$data$colidx0, ordered,
        identical(link, "logit"), 1e-6, correction.floor,
        ThreadEnv$n.threads %||% 1L
      )
      contribution <- if (order == "laplace2")
        approximation$second else approximation$first
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


lmsGraphLaplaceImplicitGradient <- function(theta, model,
                                             order = c("laplace", "laplace2"),
                                             previous = NULL,
                                             link = c("logit", "probit"),
                                             parameter.epsilon = 1e-5,
                                             mode.epsilon = 1e-4,
                                             tolerance = 1e-8,
                                             base.evaluation = NULL) {
  order <- match.arg(order)
  link <- match.arg(link)
  base <- base.evaluation %||% lmsGraphLaplaceEvaluation(
    model, theta, order, previous, link, tolerance
  )
  modes <- lapply(base$groups, function(group) group$quad$modes)
  mode.derivative <- vector("list", length(base$groups))
  mode.inverse <- vector("list", length(base$groups))

  for (g in seq_along(base$groups)) {
    current <- modes[[g]]
    d <- NCOL(current)
    derivative <- matrix(0, NROW(current), d)
    for (k in seq_len(d)) {
      step <- mode.epsilon * pmax(1, abs(current[, k]))
      plus.modes <- minus.modes <- modes
      plus.modes[[g]][, k] <- current[, k] + step
      minus.modes[[g]][, k] <- current[, k] - step
      plus <- lmsGraphLaplaceEvaluation(
        model, theta, order, base$quad, link, tolerance,
        fixed.modes = plus.modes
      )
      minus <- lmsGraphLaplaceEvaluation(
        model, theta, order, base$quad, link, tolerance,
        fixed.modes = minus.modes
      )
      use <- seq_len(NROW(current)) + sum(vapply(
        base$groups[seq_len(g - 1L)], function(x) NROW(x$quad$modes), integer(1L)
      ))
      derivative[, k] <- (plus$individual[use] - minus$individual[use]) /
        (2 * step)
    }
    mode.derivative[[g]] <- derivative
    inverse <- vector("list", NROW(current))
    for (i in seq_len(NROW(current))) {
      H <- base$groups[[g]]$approximation$modeHessian[, , i]
      H <- .5 * (H + t(H))
      decomposition <- eigen(H, symmetric = TRUE)
      values <- pmax(decomposition$values, 1e-6)
      inverse[[i]] <- decomposition$vectors %*%
        diag(1 / values, length(values)) %*% t(decomposition$vectors)
    }
    mode.inverse[[g]] <- inverse
  }

  parameter.step <- parameter.epsilon * pmax(1, abs(theta))
  parameter.variants <- vector("list", 2L * length(theta))
  for (j in seq_along(theta)) {
    plus.theta <- minus.theta <- theta
    plus.theta[j] <- plus.theta[j] + parameter.step[j]
    minus.theta[j] <- minus.theta[j] - parameter.step[j]
    parameter.variants[[2L * j - 1L]] <- fillModel(
      model = model, theta = plus.theta, method = "lms"
    )
    parameter.variants[[2L * j]] <- fillModel(
      model = model, theta = minus.theta, method = "lms"
    )
  }
  parameter.batch <- vector("list", length(base$groups))
  for (g in seq_along(base$groups)) {
    submodel <- parameter.variants[[1L]]$models[[g]]
    matrices <- lapply(parameter.variants, function(variant) {
      current <- variant$models[[g]]$matrices
      if (is.null(current$thresholds))
        current$thresholds <- matrix(NaN, NROW(current$lambdaX), 0L)
      current
    })
    ordered <- match(submodel$info$ordered %||% character(),
                     rownames(matrices[[1L]]$lambdaX)) - 1L
    ordered <- ordered[!is.na(ordered)]
    parameter.batch[[g]] <- lmsGraphLaplaceFixedBatchCpp(
      matrices, modes[[g]], submodel$info$numXis, submodel$info$numEtas,
      submodel$data$data.split, submodel$data$colidx0, ordered,
      identical(order, "laplace2"), identical(link, "logit"),
      1e-6, 1e-8, ThreadEnv$n.threads %||% 1L
    )
  }
  gradient <- numeric(length(theta))
  for (j in seq_along(theta)) {
    implicit <- 0
    partial <- 0
    for (g in seq_along(base$groups)) {
      batch <- parameter.batch[[g]]
      weights <- parameter.variants[[1L]]$models[[g]]$data$weights
      if (is.null(weights)) weights <- rep(1, NROW(modes[[g]]))
      plus.contribution <- batch$contribution[, 2L * j - 1L]
      minus.contribution <- batch$contribution[, 2L * j]
      partial <- partial + sum(weights *
        (plus.contribution - minus.contribution) / (2 * parameter.step[j]))
      plus.gradient <- batch$modeGradient[, , 2L * j - 1L, drop = FALSE]
      minus.gradient <- batch$modeGradient[, , 2L * j, drop = FALSE]
      dim(plus.gradient) <- dim(minus.gradient) <-
        c(NROW(modes[[g]]), NCOL(modes[[g]]))
      mixed <- (plus.gradient - minus.gradient) / (2 * parameter.step[j])
      n <- NROW(mixed)
      for (i in seq_len(n)) {
        mode.change <- -drop(mode.inverse[[g]][[i]] %*% mixed[i, ])
        implicit <- implicit + sum(mode.derivative[[g]][i, ] * mode.change)
      }
    }
    gradient[j] <- partial + implicit
  }
  names(gradient) <- names(theta)
  list(logLik = base$logLik, gradient = gradient, evaluation = base)
}


lmsGraphLaplaceImplicitHessian <- function(theta, model,
                                            order = c("laplace", "laplace2"),
                                            previous = NULL,
                                            link = c("logit", "probit"),
                                            epsilon = 1e-4,
                                            tolerance = 1e-8) {
  order <- match.arg(order)
  link <- match.arg(link)
  lmsGraphHessianFromGradient(theta, function(value) {
    lmsGraphLaplaceImplicitGradient(
      value, model, order, previous, link, tolerance = tolerance
    )$gradient
  }, epsilon)
}


lmsGraphBoundedLbfgs <- function(start, evaluate, lower, upper,
                                 max.iter = 500L, convergence.rel = 1e-8,
                                 memory = 10L, max.line.search = 20L,
                                 max.step = 1,
                                 verbose = FALSE) {
  project <- function(x) pmin(upper, pmax(lower, x))
  projected.gradient <- function(x, gradient) {
    out <- gradient
    tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(x))
    out[x <= lower + tolerance & gradient > 0] <- 0
    out[x >= upper - tolerance & gradient < 0] <- 0
    out
  }
  direction.lbfgs <- function(gradient, s.history, y.history) {
    count <- length(s.history)
    if (!count) return(-gradient)
    q <- gradient
    alpha <- numeric(count)
    rho <- numeric(count)
    for (i in rev(seq_len(count))) {
      rho[i] <- 1 / sum(y.history[[i]] * s.history[[i]])
      alpha[i] <- rho[i] * sum(s.history[[i]] * q)
      q <- q - alpha[i] * y.history[[i]]
    }
    last <- count
    scale <- sum(s.history[[last]] * y.history[[last]]) /
      sum(y.history[[last]]^2)
    r <- pmax(.Machine$double.eps, scale) * q
    for (i in seq_len(count)) {
      beta <- rho[i] * sum(y.history[[i]] * r)
      r <- r + s.history[[i]] * (alpha[i] - beta)
    }
    -r
  }

  theta <- project(start)
  state <- evaluate(theta, NULL)
  if (!isTRUE(state$valid))
    stop("The initial Laplace QN evaluation is not finite.", call. = FALSE)
  s.history <- y.history <- list()
  converged <- FALSE
  message <- "iteration limit reached"
  iterations <- 0L
  evaluations <- 1L

  for (iteration in seq_len(as.integer(max.iter))) {
    iterations <- iteration
    pg <- projected.gradient(theta, state$gradient)
    gradient.norm <- max(abs(pg))
    if (isTRUE(verbose))
      cat(sprintf("QN %4d  objective %.10g  |projected gradient| %.4g\n",
                  iteration - 1L, state$objective, gradient.norm))
    if (gradient.norm <= convergence.rel * pmax(1, abs(state$objective))) {
      converged <- TRUE
      message <- "projected gradient convergence"
      break
    }
    direction <- direction.lbfgs(pg, s.history, y.history)
    direction[(theta <= lower & direction < 0) |
              (theta >= upper & direction > 0)] <- 0
    slope <- sum(state$gradient * direction)
    if (!is.finite(slope) || slope >= 0) {
      s.history <- y.history <- list()
      direction <- -pg
      slope <- -sum(pg^2)
    }
    if (is.finite(max.step) && max(abs(direction)) > max.step)
      direction <- direction * (max.step / max(abs(direction)))
    slope <- sum(state$gradient * direction)
    positive <- direction > 0 & is.finite(upper)
    negative <- direction < 0 & is.finite(lower)
    maximum.step <- Inf
    if (any(positive)) maximum.step <- min(maximum.step,
      min((upper[positive] - theta[positive]) / direction[positive]))
    if (any(negative)) maximum.step <- min(maximum.step,
      min((lower[negative] - theta[negative]) / direction[negative]))
    step <- min(1, maximum.step)
    accepted <- FALSE
    trial <- NULL
    for (line.iteration in seq_len(as.integer(max.line.search))) {
      trial.theta <- project(theta + step * direction)
      if (identical(trial.theta, theta)) break
      trial <- evaluate(trial.theta, state$evaluation$quad)
      evaluations <- evaluations + 1L
      if (isTRUE(trial$valid) && trial$objective <=
          state$objective + 1e-4 * step * slope) {
        accepted <- TRUE
        break
      }
      step <- step * .5
    }
    if (!accepted) {
      message <- "line search failed"
      break
    }
    s <- trial.theta - theta
    y <- trial$gradient - state$gradient
    curvature <- sum(s * y)
    if (is.finite(curvature) && curvature >
        1e-10 * sqrt(sum(s^2) * sum(y^2))) {
      s.history <- c(s.history, list(s))
      y.history <- c(y.history, list(y))
      if (length(s.history) > memory) {
        s.history <- s.history[-1L]
        y.history <- y.history[-1L]
      }
    }
    relative.change <- abs(state$objective - trial$objective) /
      pmax(1, abs(state$objective))
    theta <- trial.theta
    state <- trial
    trial.gradient.norm <- max(abs(projected.gradient(theta, state$gradient)))
    if (relative.change <= convergence.rel && trial.gradient.norm <=
        sqrt(convergence.rel) * pmax(1, abs(state$objective))) {
      converged <- TRUE
      message <- "relative objective convergence"
      break
    }
  }
  if (isTRUE(verbose)) {
    pg <- projected.gradient(theta, state$gradient)
    cat(sprintf("QN final objective %.10g  |projected gradient| %.4g\n",
                state$objective, max(abs(pg))))
  }
  list(par = theta, objective = state$objective, state = state,
       iterations = iterations, evaluations = evaluations,
       convergence = if (converged) 0L else 1L, message = message)
}


estLmsGraphLaplace <- function(model,
                               order = c("laplace", "laplace2"),
                               link = c("logit", "probit"),
                               optimizer = c("custom", "nlminb", "L-BFGS-B"),
                               max.iter = 500L,
                               max.step = 1,
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
  cache <- new.env(parent = emptyenv())
  cache$theta <- NULL
  cache$evaluation <- NULL
  evaluate <- function(theta) {
    if (!is.null(cache$theta) && identical(theta, cache$theta))
      return(cache$evaluation)
    value <- tryCatch(
      lmsGraphLaplaceEvaluation(
        model, theta, order, previous, link, adaptive.quad.tol
      ), error = identity
    )
    evaluations <<- evaluations + 1L
    if (!inherits(value, "error") && is.finite(value$logLik)) {
      previous <<- value$quad
      cache$theta <- theta
      cache$evaluation <- value
    }
    value
  }
  objective <- function(theta) {
    value <- evaluate(theta)
    if (inherits(value, "error") || !is.finite(value$logLik))
      return(.Machine$double.xmax / 1e100)
    -value$logLik
  }
  gradient <- function(theta) {
    value <- evaluate(theta)
    if (inherits(value, "error") || !is.finite(value$logLik))
      return(rep(0, length(theta)))
    result <- tryCatch(lmsGraphLaplaceImplicitGradient(
      theta, model, order, value$quad, link,
      tolerance = adaptive.quad.tol, base.evaluation = value
    ), error = identity)
    if (inherits(result, "error") || any(!is.finite(result$gradient)))
      return(lmsGraphGradient(theta, objective, 1e-5))
    -result$gradient
  }
  start.time <- proc.time()[["elapsed"]]
  if (optimizer == "custom") {
    initial.available <- TRUE
    evaluate.fg <- function(theta, accepted.quad) {
      if (initial.available && is.null(accepted.quad) &&
          identical(theta, model$theta)) {
        value <- initial
        initial.available <<- FALSE
      } else {
        value <- tryCatch(lmsGraphLaplaceEvaluation(
          model, theta, order, accepted.quad, link, adaptive.quad.tol
        ), error = identity)
        evaluations <<- evaluations + 1L
      }
      if (inherits(value, "error") || !is.finite(value$logLik))
        return(list(valid = FALSE, objective = Inf,
                    gradient = rep(NA_real_, length(theta))))
      derivative <- tryCatch(lmsGraphLaplaceImplicitGradient(
        theta, model, order, value$quad, link,
        tolerance = adaptive.quad.tol, base.evaluation = value
      ), error = identity)
      if (inherits(derivative, "error") ||
          any(!is.finite(derivative$gradient)))
        return(list(valid = FALSE, objective = Inf,
                    gradient = rep(NA_real_, length(theta))))
      list(valid = TRUE, objective = -value$logLik,
           gradient = -derivative$gradient, evaluation = value)
    }
    fit <- lmsGraphBoundedLbfgs(
      model$theta, evaluate.fg, bounds$lower, bounds$upper,
      max.iter, convergence.rel, max.step = max.step, verbose = verbose
    )
    theta <- fit$par
    iterations <- fit$iterations
    converged <- fit$convergence == 0L
    message <- fit$message
    previous <- fit$state$evaluation$quad
  } else if (optimizer == "nlminb") {
    fit <- stats::nlminb(
      start = model$theta, objective = objective,
      gradient = gradient,
      lower = bounds$lower, upper = bounds$upper,
      control = list(iter.max = as.integer(max.iter),
                     eval.max = max(1000L, 10L * as.integer(max.iter)),
                     rel.tol = convergence.rel,
                     trace = if (isTRUE(verbose)) 1L else 0L)
    )
    theta <- fit$par
    iterations <- fit$iterations
    converged <- fit$convergence == 0L
    message <- fit$message
  } else {
    fit <- stats::optim(
      par = model$theta, fn = objective, gr = gradient, method = "L-BFGS-B",
      lower = bounds$lower, upper = bounds$upper,
      control = list(maxit = as.integer(max.iter), factr = max(
        1, convergence.rel / .Machine$double.eps
      ), trace = if (isTRUE(verbose)) 1L else 0L, REPORT = 1L)
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
    optimizer = paste0("QN-", optimizer, "-", toupper(order)),
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
