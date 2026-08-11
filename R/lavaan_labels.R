combineLavLabels <- function(lavLabelsCov, lavLabelsMain, currentLabels, g = 1L) {
  lavLabels <- c(lavLabelsCov, lavLabelsMain)

  if (g > 1L)
    lavLabels <- stats::setNames(sprintf("%s.g%d", lavLabels, g),
                                 nm = names(lavLabels))

  finalLabels <- currentLabels
  finalLabels[finalLabels %in% names(lavLabels)] <-
    lavLabels[names(lavLabels) %in% finalLabels]
  finalLabels
}


createLavLabels <- function(matrices, subset, etas, parTable.in = NULL) {
  lambdaX      <- createLabelsMatrix(matrices$lambdaX, op = "=~")
  lambdaY      <- createLabelsMatrix(matrices$lambdaY, op = "=~")
  thetaDelta   <- createLabelsMatrix(matrices$thetaDelta, op = "~~")
  thetaEpsilon <- createLabelsMatrix(matrices$thetaEpsilon, op = "~~")
  W            <- createLabelsMatrix(matrices$W, op = "<~")
  T            <- createLabelsMatrix(matrices$T, op = "~~")
  phi          <- createLabelsMatrix(matrices$phi, op = "~~")
  A            <- createLabelsMatrix(matrices$A, op = "~~")
  covZetaXi    <- createLabelsMatrix(matrices$covZetaXi, op = "~~", first = "cols")
  psi          <- createLabelsMatrix(matrices$psi, op = "~~")
  tauX         <- createLabelsMatrix(matrices$tauX, op = "~", first = "rows")
  tauY         <- createLabelsMatrix(matrices$tauY, op = "~", first = "rows")
  alpha        <- createLabelsMatrix(matrices$alpha, op = "~", first = "rows")
  beta0        <- createLabelsMatrix(matrices$beta0, op = "~", first = "rows")
  gammaXi      <- createLabelsMatrix(matrices$gammaXi, op = "~", first = "rows")
  gammaEta     <- createLabelsMatrix(matrices$gammaEta, op = "~", first = "rows")
  omegaXiXi    <- createLabelsOmega(matrices$omegaXiXi, parTable.in = parTable.in)
  omegaEtaXi   <- createLabelsOmega(matrices$omegaEtaXi, parTable.in = parTable.in)
  thresholdDelta <- createLabelsThreshold(matrices$thresholdDelta)
  omega        <- createLabelsProductOmega(matrices$omega,
                                           parTable.in = parTable.in)

  labels <- c(
    lambdaX = lambdaX,
    lambdaY = lambdaY,
    tauX = tauX,
    tauY = tauY,
    thetaDelta = thetaDelta,
    thetaEpsilon = thetaEpsilon,
    W = W,
    T = T,
    phi = phi,
    A = A,
    covZetaXi = covZetaXi,
    psi = psi,
    alpha = alpha,
    beta0 = beta0,
    gammaXi = gammaXi,
    gammaEta = gammaEta,
    omegaXiXi = omegaXiXi,
    omegaEtaXi = omegaEtaXi,
    thresholdDelta = thresholdDelta,
    omega = omega
  )

  labels[subset]
}


createLavLabelsCov <- function(matrices, subset) {
  if (is.null(matrices)) return(NULL)

  phi      <- createLabelsMatrix(matrices$phi, op = "~~")
  psi      <- createLabelsMatrix(matrices$psi, op = "~~")
  gammaXi  <- createLabelsMatrix(matrices$gammaXi, op = "~", first = "rows")
  gammaEta <- createLabelsMatrix(matrices$gammaEta, op = "~", first = "rows")

  labels <- c("phi" = phi,
              "psi" = psi,
              "gammaXi" = gammaXi,
              "gammaEta" = gammaEta)

  labels[subset]
}


createLabelsMatrix <- function(X, op = "~", first = "cols") {
  labels <- character(0L)
  rows   <- rownames(X)
  cols   <- colnames(X)

  # this is ugly, but... we have to read by cols first
  getLabel <- switch(first, cols = function(col, row) paste0(col, op, row),
                     rows = function(col, row) paste0(row, op, col))

  for (i in seq_len(NCOL(X))) for (j in seq_len(NROW(X))) {
    labels <- c(labels, getLabel(col = cols[[i]], row = rows[[j]]))
  }

  labels
}


createLabelsOmega <- function(X, parTable.in = NULL) {
  C <- \(x, y) sprintf("%s:%s", x, y)
  getIntTerm <- function(lhs, rhs) {
    xz <- C(lhs, rhs)
    zx <- C(rhs, lhs)

    if (is.null(parTable.in))         xz
    else if (xz %in% parTable.in$rhs) xz
    else                              zx
  }

  rows   <- rownames(X)
  cols   <- colnames(X)
  labels <- character(0L)

  for (i in seq_len(ncol(X))) for (j in seq_len(nrow(X))) {
    eta_i <- stringr::str_split_i(rows[[j]], pattern = "~", i = 1) # y~x -> y
    lhs_i <- stringr::str_split_i(rows[[j]], pattern = "~", i = 2) # y~x -> x
    rhs_i <- cols[[i]]
    xz    <- getIntTerm(lhs_i, rhs_i)

    labels <- c(labels, sprintf("%s~%s", eta_i, xz))
  }

  labels
}


getLavCoefs <- function(model, theta, method) {
  fullTheta <- getTransformationsTheta(model, theta, method)
  fullNames  <- names(fullTheta)
  thetaNames <- names(theta)

  if (!is.null(fullNames) && !is.null(thetaNames)) {
    isFree <- fullNames %in% thetaNames
  } else if (!is.null(fullNames) && is.null(thetaNames)) {
    isFree <- logical(length(fullTheta))
  } else {
    isFree <- seq_along(fullTheta) %in% seq_along(theta)
  }

  lavLabels <- model$params$lavLabels

  if (!is.null(lavLabels) && length(lavLabels) == length(fullTheta)) {
    names(fullTheta) <- lavLabels
  } else if (!is.null(fullNames)) {
    names(fullTheta) <- fullNames
  }

  list(all = fullTheta, free = fullTheta[isFree])
}


# Labels for the graph backend's product coefficients. Far simpler than the
# Kronecker case: the product name is already the right-hand side. Column-major
# to match `as.vector()`, as every other label constructor is.
createLabelsProductOmega <- function(X, parTable.in = NULL) {
  if (is.null(X) || !length(X)) return(character(0L))
  rows <- rownames(X)
  cols <- colnames(X)
  labels <- character(0L)

  key <- function(product)
    paste(sort(lmsGraphProductFactors(product)), collapse = ":")

  for (i in seq_len(ncol(X))) for (j in seq_len(nrow(X))) {
    product <- cols[[i]]
    # Report the factor ordering the user wrote rather than the canonical one,
    # so `Z:X` in the syntax does not come back as `X:Z`.
    if (!is.null(parTable.in)) {
      candidates <- parTable.in$rhs[parTable.in$op == "~" &
                                      parTable.in$lhs == rows[[j]]]
      same <- candidates[vapply(candidates, function(r)
        identical(key(r), key(product)), logical(1L))]
      if (length(same)) product <- same[[1L]]
    }
    labels <- c(labels, sprintf("%s~%s", rows[[j]], product))
  }
  labels
}


# Thresholds are stored as an indicators-by-categories matrix with NaN where a
# category does not exist. Absent entries still need a slot so the label vector
# stays aligned with `allModelValues`.
createLabelsThreshold <- function(X) {
  if (is.null(X) || !length(X)) return(character(0L))
  rows <- rownames(X)
  labels <- character(0L)
  for (i in seq_len(ncol(X))) for (j in seq_len(nrow(X)))
    labels <- c(labels, sprintf("%s|t%d", rows[[j]], i))
  labels
}
