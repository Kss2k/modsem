#' @describeIn modsem_predict
#' Computes (optionally standardised) factor scores via the
#'   regression method using the baseline model unless \code{H0 = FALSE}.
#'
#' @param object \code{\link{modsem_da}} object
#' @param standardized Logical. If \code{TRUE}, return standardized factor scores.
#' @param newdata Compute factor scores based on a different dataset, than the one used in the model estimation.
#' @export
modsem_predict.modsem_da <- function(object,
                                     newdata = NULL,
                                     method = c("EBM", "ML"),
                                     type = c("lv", "ov", "all"),
                                     standardized = FALSE,
                                     ...) {
  modsemPredictDA(
    object = object,
    newdata = newdata,
    method = method
  )$Eta
}


modsemPredictDA <- function(object, newdata = NULL, method = c("EBM", "ML")) {
  model    <- object$model
  submodels <- model$models

  if (is.null(newdata)) {
    DATA <- object$data
  } else {
    DATA <- prepNewdataDA(object, newdata)
  }

  ETA <- vector("list", length = length(submodels))
  YY  <- vector("list", length = length(submodels))

  for (g in seq_along(submodels)) {

    ETA[[g]] <- modsemPredictEtaDA_Group(
      submodel = submodels[[g]],
      data     = DATA[[g]],
      method   = method
    )

    YY[[g]] <- YFromEta(
      Eta      = ETA[[g]],
      matrices = submodels[[g]]$matrices
    )
  }

  list(
    Eta = modsemMatrix(do.call(rbind, ETA), is.public = TRUE),
    Y   = modsemMatrix(do.call(rbind, YY), is.public = TRUE)
  )
}


modsemPredictEtaDA_Group <- function(submodel, data, method = c("EBM", "ML")) {
  method <- match.arg(toupper(method), c("EBM", "ML"))

  # TODO:
  #  We can handle composites (theoretically) by converting them to observed
  #  variables, and modifying data and matrices$lambda/matrices$theta accordingly.
  #  This can be done in modsemPredictPreProcessMatrices().
  stopif(isTRUE(submodel$info$hasComposites),
    "modsem_predict() does not work with composite variables (yet)!"
  )

  # data
  data.split <- data$data.split
  patterns   <- data$patterns
  rowidx     <- data$rowidx
  n.total    <- data$n

  # matrices
  matrices <- modsemPredictPreProcessMatrices(submodel$matrices, data = data)
  xptr     <- modelMatrixCacheCpp(matrices)

  .f <- switch(method,
    ML  = logLikFromZetaMLCpp,
    EBM = logLikFromZetaEBMCpp,
    stop2("Unrecognized method: ", method, "!")
  )

  k   <- NCOL(matrices$phi) + NCOL(matrices$psi)
  dim <- c(colnames(matrices$phi), colnames(matrices$psi))

  Eta <- matrix(NA_real_, nrow = n.total, ncol = k, dimnames = list(NULL, dim))

  for (p in seq_along(data.split)) {
    data.p    <- data.split[[p]]
    pattern.p <- patterns[p, , drop = TRUE]
    idx.p     <- which(pattern.p) - 1L

    n     <- NROW(data.p)
    start <- rep(0, k)
    Eta.p <- matrix(NA_real_, nrow = n, ncol = k)

    for (i in seq_len(n)) {
      y <- data.p[i, , drop = TRUE]

      opt_i <- nlminb(
        start     = start,
        objective = \(zeta) -.f(zeta = zeta, y = y, xptr = xptr, idx = idx.p),
        gradient  = NULL
      )

      Eta.p[i, ] <- impliedEtaFromZetaCpp(zeta = opt_i$par, xptr = xptr)
    }

    Eta[rowidx[[p]], ] <- Eta.p
  }

  Eta
}


YFromEta <- function(Eta, matrices) {
  LambdaX <- matrices$lambdaX
  LambdaY <- matrices$lambdaY
  tauX    <- matrices$tauX
  tauY    <- matrices$tauY

  Lambda <- diagPartitionedMat(LambdaX, LambdaY)
  tau    <- c(tauX, tauY)

  sweep(Eta %*% t(Lambda), MARGIN = 2, STATS = tau, FUN = "+")
}


prepNewdataDA <- function(object, newdata) {
  newdata    <- as.data.frame(newdata)
  newdata    <- addStructOVColumnsDA(newdata, object)  # add TEMP_OV cols
  submodels  <- object$model$models
  group.var  <- object$args$group
  n.groups   <- length(submodels)

  if (!is.null(group.var)) {
    stopif(!group.var %in% colnames(newdata),
           sprintf("Grouping variable '%s' not found in `newdata`.", group.var))

    group.levels <- object$model$info$group.levels
    stopif(is.null(group.levels),
           "Model has a grouping variable but group levels could not be determined.")

    group.values <- as.character(newdata[[group.var]])

    DATA <- lapply(seq_len(n.groups), function(g) {
      cols <- colnames(submodels[[g]]$data$data.full)

      missing.cols <- setdiff(cols, colnames(newdata))
      stopif(length(missing.cols),
             "Missing columns in `newdata`:\n  ", paste(missing.cols, collapse = ", "))

      rows.g <- group.values == group.levels[[g]]
      mat.g  <- as.matrix(newdata[rows.g, cols, drop = FALSE])
      storage.mode(mat.g) <- "double"

      patternizeMissingDataFIML(mat.g)
    })

  } else {
    cols <- colnames(submodels[[1L]]$data$data.full)

    missing.cols <- setdiff(cols, colnames(newdata))
    stopif(length(missing.cols),
           "Missing columns in `newdata`:\n  ", paste(missing.cols, collapse = ", "))

    mat <- as.matrix(newdata[, cols, drop = FALSE])
    storage.mode(mat) <- "double"

    DATA <- list(patternizeMissingDataFIML(mat))
  }

  DATA
}


addStructOVColumnsDA <- function(newdata, object) {
  group.info <- object$model$info$group.info
  structovs  <- group.info$structovs
  ovIntTerms <- group.info$ovIntTerms

  if (!length(structovs)) return(newdata)

  # OV interaction terms: compute product columns (e.g. x1:x2 -> x1___INTERACTION___x2)
  for (ovInt in ovIntTerms) {
    vars     <- stringr::str_split_1(ovInt, ":")
    ovIntNew <- stringr::str_replace_all(ovInt, ":", OP_OV_INT)

    missing <- setdiff(vars, colnames(newdata))
    stopif(length(missing),
           "Missing columns in `newdata` required for interaction term '", ovInt, "':\n  ",
           paste(missing, collapse = ", "))

    newdata[[ovIntNew]] <- apply(newdata[, vars, drop = FALSE], MARGIN = 1L, FUN = prod)
  }

  # Plain structural OVs: check they are present, then add TEMP_OV copy
  ovIntNew_names  <- stringr::str_replace_all(ovIntTerms, ":", OP_OV_INT)
  plainStructovs <- setdiff(structovs, ovIntNew_names)

  missing <- setdiff(plainStructovs, colnames(newdata))
  stopif(length(missing),
    "Missing columns in `newdata`:\n  ", paste(missing, collapse = ", ")
  )

  for (ov in structovs) {
    newdata[[paste0(TEMP_OV_PREFIX, ov)]] <- newdata[[ov]]
  }

  newdata
}


modsemPredictPreProcessMatrices <- function(matrices, data, pct.fill = 0.01) {
  # There are two goals here:
  #   1. Handle observed variables with zeros in diag(Theta)
  #   2. Handle composite variables
  # For point 2 we convert composite variables to observed variables using the weights
  # Point 2 is not implemented (yet).
  # For point 2 we will also need to modify and return data

  W <- matrices$W
  is.composite <- apply(W, MARGIN = 2L, FUN = \(x) any(x != 0))

  if (any(is.composite)) {
    stop2("modsem_predict() does not work with composite variables (yet)!")
  }

  fillTheta <- function(Theta) {
    if (!NROW(Theta) || !NCOL(Theta))
      return(Theta)

    nm <- colnames(Theta)
    bad.idx <- which(diag(Theta) <= 0 | is.na(diag(Theta)))

    for (idx in bad.idx) {
      col <- nm[idx]
      var <- var(data$data.full[,col], na.rm = TRUE)

      Theta[idx, idx] <- pct.fill * var
    }

    Theta
  }

  matrices$thetaDelta <- fillTheta(matrices$thetaDelta)
  matrices$thetaEpsilon <- fillTheta(matrices$thetaEpsilon)

  matrices
}
