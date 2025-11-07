emGsem <- function(model,
                   algorithm = c("EMA", "EM"),
                   em.control = list(),
                   verbose = FALSE,
                   convergence.abs = 1e-4,
                   convergence.rel = 1e-10,
                   max.iter = 500,
                   max.step = 1,
                   control = list(),
                   calc.se = TRUE,
                   FIM = "observed",
                   OFIM.hessian = FALSE,
                   EFIM.S = 3e4,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   epsilon = 1e-6,
                   optimizer = "nlminb",
                   R.max = 1e6,
                   adaptive.quad = FALSE,
                   quad.range = -Inf,
                   adaptive.quad.tol = 1e-12,
                   nodes = 24,
                   cr1s = TRUE,
                   # new knobs for FS/Louis
                   fs.matrix = c("Iobs","Icom"),        # use Iobs (Louis) or fall back to Icom
                   fs.fd.scheme = c("forward","central"),
                   fs.fd.epsilon = 1e-6,
                   fs.jitter.mult = sqrt(.Machine$double.eps),
                   ...) {

  theta.lower  <- model$params$bounds$lower
  theta.upper  <- model$params$bounds$upper
  bounds.all   <- c(theta.lower, theta.upper)

  tryCatch({
    logLikNew <- -Inf
    logLikOld <- -Inf
    thetaNew  <- model$theta

    iterations <- 0L
    run <- TRUE

    testSimpleGradient <- !model$params$gradientStruct$hasCovModel
    lastQuad <- NULL
    adaptiveQuad <- isTRUE(model$models[[1L]]$quad$adaptive)
    adaptiveFreq <- model$models[[1L]]$quad$adaptive.frequency
    if (is.null(adaptiveFreq) || adaptiveFreq < 1L)
      adaptiveFreq <- 1L

    while (run) {
      iterations <- iterations + 1L
      logLikOld  <- logLikNew
      thetaOld   <- thetaNew
      recalcQuad <- adaptiveQuad && iterations %% adaptiveFreq == 0L

      # E-step at thetaOld
      P <- estepGsem(model = model, theta = thetaOld, lastQuad = lastQuad,
                     recalcQuad = recalcQuad, adaptive.quad.tol = adaptive.quad.tol)
      lastQuad <- P$QUAD

      logLikNew  <- P$obsLL
      deltaLL    <- logLikNew - logLikOld
      relDeltaLL <- if (is.finite(logLikOld)) deltaLL / abs(logLikOld) else Inf

      updateStatusLog(iterations, "EM", logLikNew, deltaLL, relDeltaLL, verbose)

      converged <- (abs(deltaLL) < convergence.abs) ||
                   (abs(relDeltaLL) < convergence.rel)

      if (iterations >= max.iter || converged) break

      if (deltaLL < -1e-8) {
        if (verbose) cat("\n")
        warning2(sprintf("Loglikelihood decreased by %.2g", deltaLL))
      }

      mstep <- mstepGsem(model = model, P = P, theta = thetaOld,
                         max.step = max.step, control = control, ...)
      if (!any(is.na(mstep$par))) thetaNew <- mstep$par
    }

    if (verbose) cat("\n")
    warnif(iterations >= max.iter, "Maximum iterations reached!\n",
           "Consider tweaking these parameters:\n",
           formatParameters(convergence.abs, convergence.rel, algorithm,
                            max.step, max.iter, nodes, adaptive.quad,
                            adaptive.quad.tol, quad.range))

  })
  return(mstep)
}

resolve_gsem_threads <- function() {
  ncores <- ThreadEnv$n.threads
  if (is.null(ncores) || !is.finite(ncores) || ncores < 1L)
    1L
  else
    as.integer(ncores)
}

gsem_extract_nodes_matrix <- function(quad, k = NULL) {
  if (is.null(k)) {
    if (!is.null(quad$k)) k <- quad$k
    else if (!is.null(quad$nodes_matrix)) k <- ncol(quad$nodes_matrix)
    else k <- 0L
  }

  if (!is.null(quad$nodes_matrix))
    return(quad$nodes_matrix)

  nodes <- quad$n
  if (is.null(nodes))
    return(matrix(numeric(), nrow = 0L, ncol = k))

  if (is.list(nodes)) {
    nNodes <- length(nodes)
    if (!k)
      return(matrix(numeric(), nrow = nNodes, ncol = 0L))

    out <- vapply(nodes, function(mat) {
      cols <- ncol(mat)
      start <- max(1L, cols - k + 1L)
      mat[1L, start:cols]
    }, numeric(k))
    nodes_mat <- t(out)
  } else {
    nNodes <- if (is.null(dim(nodes))) 1L else nrow(nodes)
    if (!k)
      return(matrix(numeric(), nrow = nNodes, ncol = 0L))

    cols <- ncol(nodes)
    start <- max(1L, cols - k + 1L)
    nodes_mat <- matrix(nodes[1L, start:cols], nrow = 1L)
  }

  storage.mode(nodes_mat) <- "double"
  nodes_mat
}

gsem_extract_weight_vector <- function(quad) {
  if (!is.null(quad$weights_vec))
    return(quad$weights_vec)

  w <- quad$w
  if (is.null(w))
    return(numeric())

  if (is.matrix(w))
    as.numeric(w[1L, , drop = TRUE])
  else
    as.numeric(w)
}

gsem_nodes_matrix_to_list <- function(nodes_matrix, nObs) {
  if (!NROW(nodes_matrix))
    return(vector("list", 0L))

  lapply(seq_len(nrow(nodes_matrix)), function(i) {
    matrix(nodes_matrix[i, , drop = FALSE], nrow = nObs,
           ncol = ncol(nodes_matrix), byrow = TRUE)
  })
}

gsem_expand_nodes <- function(submodel, nodes_matrix) {
  nObs <- submodel$info$N
  node_list <- gsem_nodes_matrix_to_list(nodes_matrix, nObs)

  omega <- submodel$matrices$OMEGA
  if (!NROW(omega) || !length(node_list))
    return(node_list)

  xis <- submodel$info$xis
  etas <- submodel$info$etas
  latent_names <- c(xis, etas)
  if (!length(latent_names))
    return(node_list)

  psi <- submodel$matrices$psi
  psi_block <- psi[latent_names, latent_names, drop = FALSE]
  alpha <- submodel$matrices$alpha
  alpha_latent <- if (nrow(alpha)) alpha[latent_names, 1L, drop = TRUE] else rep(0, length(latent_names))

  expand_quadrature_cpp(
    Z_list      = node_list,
    cholPsi     = chol(psi_block),
    gamma       = submodel$matrices$gamma,
    omega       = omega,
    alphaLatent = alpha_latent,
    xiNames     = xis,
    etaNames    = etas,
    latentNames = latent_names
  )
}

gsem_build_quad <- function(submodel, quadTemplate, nodes_matrix, weights_vec) {
  if (is.null(nodes_matrix)) {
    k <- quadTemplate$k
    if (is.null(k)) k <- 0L
    nodes_matrix <- matrix(numeric(), nrow = 0L, ncol = k)
  }

  nodes_matrix <- as.matrix(nodes_matrix)
  nObs <- submodel$info$N
  nNodes <- nrow(nodes_matrix)

  if (is.null(weights_vec) || !length(weights_vec)) {
    weights_vec <- if (nNodes) rep(1 / nNodes, nNodes) else numeric()
  }

  if (length(weights_vec) && length(weights_vec) != nNodes)
    stop2("Number of weights does not match number of quadrature nodes.")

  weights_matrix <- if (length(weights_vec)) {
    matrix(weights_vec, nrow = nObs, ncol = length(weights_vec), byrow = TRUE)
  } else matrix(0, nrow = nObs, ncol = 0L)

  expanded <- gsem_expand_nodes(submodel, nodes_matrix)

  quadTemplate$n <- expanded
  quadTemplate$w <- weights_matrix
  quadTemplate$nodes_matrix <- nodes_matrix
  quadTemplate$weights_vec <- weights_vec

  quadTemplate
}

gsem_density_from_nodes <- function(submodel, nodes_matrix, ncores) {
  quadTmp <- submodel$quad
  expanded <- gsem_expand_nodes(submodel, nodes_matrix)
  nObs <- submodel$info$N
  quadTmp$n <- expanded
  quadTmp$w <- matrix(1, nrow = nObs, ncol = length(expanded), byrow = TRUE)

  tmp <- submodel
  tmp$quad <- quadTmp

  P_Step_GSEM_Group(tmp, normalized = FALSE, ncores = ncores)
}

adaptive_quadrature_gsem <- function(submodel, quadTemplate, lastQuad, tol, ncores) {
  cat("\nFITTING ADPATIVE QUADRATURE...")
  k <- quadTemplate$k
  if (is.null(k) || !length(k)) {
    base_cols <- if (!is.null(quadTemplate$nodes_matrix)) ncol(quadTemplate$nodes_matrix) else 0L
    k <- base_cols
  }
  if (!length(k) || !k)
    return(list(
      nodes_matrix = quadTemplate$nodes_matrix,
      weights_vec = quadTemplate$weights_vec,
      m.ceil = quadTemplate$m.ceil
    ))

  m <- quadTemplate$m
  a <- quadTemplate$a
  b <- quadTemplate$b

  m.ceil <- NULL
  if (!is.null(lastQuad$m.ceil))
    m.ceil <- lastQuad$m.ceil
  else if (!is.null(quadTemplate$m.ceil))
    m.ceil <- quadTemplate$m.ceil
  else if (k > 1)
    m.ceil <- rep(m, length.out = k)
  else
    m.ceil <- round(estMForNodesInRange(m, a = -5, b = 5))

  dens_fun <- function(node_matrix) {
    gsem_density_from_nodes(submodel, node_matrix, ncores = ncores)
  }

  quad_adapt <- adaptiveGaussQuadrature(
    fun = dens_fun,
    collapse = \(x) sum(log(rowSums(x))),
    a = rep(a, length.out = k),
    b = rep(b, length.out = k),
    m = m,
    m.ceil = m.ceil,
    k = k,
    tol = tol
  )

  cat(" DONE!\n")
  list(
    nodes_matrix = quad_adapt$n,
    weights_vec = quad_adapt$w,
    m.ceil = quad_adapt$m.ceil
  )
}

prepare_gsem_quadrature <- function(submodel, lastQuad = NULL, recalcQuad = FALSE,
                                    adaptive.quad.tol = 1e-10, ncores = 1L) {
  quad <- submodel$quad
  k <- quad$k

  nodes_matrix <- NULL
  weights_vec <- NULL

  if (!is.null(lastQuad)) {
    nodes_matrix <- if (!is.null(lastQuad$nodes_matrix)) lastQuad$nodes_matrix else gsem_extract_nodes_matrix(lastQuad, k)
    weights_vec <- if (!is.null(lastQuad$weights_vec)) lastQuad$weights_vec else gsem_extract_weight_vector(lastQuad)
    if (!is.null(lastQuad$m.ceil))
      quad$m.ceil <- lastQuad$m.ceil
  }

  if (isTRUE(quad$adaptive)) {
    if (is.null(nodes_matrix) || recalcQuad) {
      adapt <- adaptive_quadrature_gsem(submodel, quad, lastQuad, adaptive.quad.tol, ncores)
      nodes_matrix <- adapt$nodes_matrix
      weights_vec <- adapt$weights_vec
      quad$m.ceil <- adapt$m.ceil
    }
  } else if (is.null(nodes_matrix)) {
    nodes_matrix <- gsem_extract_nodes_matrix(quad, k)
    weights_vec <- gsem_extract_weight_vector(quad)
  }

  if (is.null(nodes_matrix))
    nodes_matrix <- gsem_extract_nodes_matrix(quad, k)
  if (is.null(weights_vec))
    weights_vec <- gsem_extract_weight_vector(quad)

  gsem_build_quad(submodel, quad, nodes_matrix, weights_vec)
}


estepGsem <- function(model, theta, lastQuad = NULL, recalcQuad = FALSE, adaptive.quad.tol = 1e-10) {
  modFilled <- fillModelGsem(model = model, theta = theta)

  ncores <- resolve_gsem_threads()
  QUAD <- vector("list", length = length(modFilled$models))
  for (g in seq_along(modFilled$models)) {
    submodel <- modFilled$models[[g]]
    lastQuad.g <- if (!is.null(lastQuad) && length(lastQuad) >= g) lastQuad[[g]] else NULL
    quadPrepared <- prepare_gsem_quadrature(
      submodel = submodel,
      lastQuad = lastQuad.g,
      recalcQuad = recalcQuad,
      adaptive.quad.tol = adaptive.quad.tol,
      ncores = ncores
    )
    modFilled$models[[g]]$quad <- quadPrepared
    QUAD[[g]] <- quadPrepared
  }

  P <- P_Step_GSEM(modFilled, normalized = FALSE, ncores = ncores)

  density <- rowSums(P)
  P       <- P / density
  obsLL   <- sum(log(density))

  list(P = P, QUAD = QUAD, obsLL = obsLL)
}


Q_Gsem <- function(theta, model, P_Step) {
  modFilled <- fillModelGsem(model = model, theta = theta)
  # This should probably be refactored to C++ code
  for (g in seq_along(modFilled$models)) { # This should probably be refactored to C++ code
    quad.g <- P_Step$QUAD[[g]]
    modFilled$models[[g]]$quad <- quad.g
  }

  Q_GSEM(modFilled, P = P_Step$P, ncores = resolve_gsem_threads())
}


fdgrad <- function(x, .f, ...,  eps = 1e-5) {
  g <- numeric(length(x))

  f0 <- .f(x, ...)
  for (i in seq_along(x)) {
    xi <- x
    xi[i] <- x[i] + eps

    fi <- .f(xi, ...)
    g[i] <- (fi - f0) / eps
  }

  g
}


# fastOneStepNLMINB <- function(start, objective, gradient = fdgrad, lower = NULL, upper = NULL, ...) {
#   direction <- -fdgrad(x = start, .f = objective, ...)
# 
#   f0 <- objective(start, ...)
# 
#   success <- FALSE
# 
#   if (!is.null(lower) || !is.null(upper)) {
#     if (is.null(lower)) lower <- -Inf
#     if (is.null(upper)) upper <- Inf
# 
#     truncatePars <- function(x) {
#       x[x < lower] <- lower[x < lower]
#       x[x > upper] <- upper[x > upper]
#       x
#     }
# 
#   } else truncatePars <- \(x) x
# 
#   alpha <- 1
#   while (alpha > 1e-5) {
#     test <- truncatePars(start + alpha * direction)
#     f1 <- objective(test, ...)
# 
#     if (!is.nan(f1) && !is.na(f1) && f1 < f0) {
#       success <- TRUE
#       break
#     }
# 
#     alpha <- alpha / 2L
#   }
# 
#   list(
#     objective = if (success) f1   else f0,
#     par       = if (success) test else start,
#     success   = success
#   )
# }


gradQ_Gsem <- function(theta, model, P_Step) {
  modFilled <- fillModelGsem(model = model, theta = theta)
  for (g in seq_along(modFilled$models)) { # This should probably be refactored to C++ code
    quad.g <- P_Step$QUAD[[g]]
    modFilled$models[[g]]$quad <- quad.g
  }

  Grad_Q_GSEM(modelR = modFilled, P = P_Step$P, ncores = resolve_gsem_threads())[names(theta)]
}


mstepGsem <- function(theta, model, P_Step, max.step = 1L, control = list(),
                      lower = NULL, upper = NULL) {
  control$iter.max <- max.step

  gradient  <- \(theta) -gradQ_Gsem(theta = theta, model = model, P_Step = P_Step)
  objective <- \(theta) -Q_Gsem(theta = theta, model = model, P_Step = P_Step)
  # gradient2  <- \(theta) setNames(fdgrad(x = theta, .f = objective, eps = 1e-7), names(theta))

  # g1 <- gradient(theta)
  # g2 <- gradient2(theta)

  # cat("\nAnalytical Gradient:\n")
  # print(g1)
  # cat("Numerical Gradient:\n")
  # print(g2)
  # cat("Difference:\n")
  # print(round(g2 - g1, 4))
  # stop("Stopping after grad check")
  # timeExpr(
  # fit2 <- fastOneStepNLMINB(start = theta, objective = objective,
  #                           lower = lower, upper = upper)
  # )

  fit <- nlminb(start = theta, objective = objective, control = control,
                lower = lower, upper = upper, gradient = gradient)

  fit
}
