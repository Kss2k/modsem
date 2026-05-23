#' @describeIn modsem_predict
#' Computes factor scores or model-implied observed-variable scores for a
#' \code{\link{modsem_da}} model via MAP optimisation.
#' @export
modsem_predict.modsem_da <- function(object,
                                     newdata = NULL,
                                     method = c("EBM", "ML", "Bartlett", "Regression"),
                                     type = c("lv", "ov", "all"),
                                     standardized = FALSE,
                                     ...) {
  method <- match.arg(
    arg        = toupper(method)[[1L]],
    choices    = c("EBM", "ML", "BARTLETT", "REGRESSION"),
    several.ok = FALSE
  )

  if      (method == "BARTLETT")   method <- "EBM"
  else if (method == "REGRESSION") method <- "ML"

  predicted <- modsemPredictDA(
    object  = object,
    newdata = newdata,
    method  = method
  )

  Eta <- predicted$Eta
  Y   <- predicted$Y

  type <- match.arg(
    arg        = tolower(type)[[1L]],
    choices    = c("lv", "ov", "all"),
    several.ok = FALSE
  )

  out <- switch(type,
    lv  = Eta,
    ov  = Y,
    all = cbind(Eta, Y[,setdiff(colnames(Y), colnames(Eta)), drop = FALSE])
  )

  if (standardized)
    out <- modsemMatrix(scale(out))

  out
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
    if (is.null(DATA[[g]])) next

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

  # collect
  Eta <- do.call(rbind, ETA)
  Y   <- do.call(rbind, YY)

  # Restore original row order when newdata was split by group.
  # Each non-NULL DATA[[g]] carries an "orig.idx" attribute recording its rows'
  # positions in newdata; unlist() collects them in group-processing order so
  # order() gives the permutation back to input order.
  orig.idx <- unlist(lapply(DATA, attr, "orig.idx"))

  if (length(orig.idx)) {
    ord <- order(orig.idx)
    Eta <- Eta[ord, , drop = FALSE]
    Y   <- Y[ord, , drop = FALSE]
  }

  list(
    Eta = modsemMatrix(Eta, is.public = TRUE),
    Y   = modsemMatrix(Y, is.public = TRUE)
  )
}


modsemPredictEtaDA_Group <- function(submodel, data, method = c("EBM", "ML")) {
  method <- match.arg(toupper(method), c("EBM", "ML"))

  # Convert composite LVs to observed scores upfront: replaces their multi-
  # indicator blocks with a single score column and rebuilds the patternized
  # data to match. This is analogous to how TEMP_OV wrapping handles structural
  # OVs -- the composite score becomes a perfect single-indicator of the LV.
  processed  <- convertCompositesToObs(submodel$matrices, data)
  data       <- processed$data

  # data
  data.split <- data$data.split
  patterns   <- data$patterns
  rowidx     <- data$rowidx
  n.total    <- data$n

  # matrices: fill any zero/NA theta diagonals (TEMP_OV, single-indicator composites)
  matrices <- modsemPredictPreProcessMatrices(processed$matrices, data = data)
  xptr     <- modelMatrixCacheCpp(matrices)

  .f <- switch(method,
    ML  = logLikFromZetaMLCpp,
    EBM = logLikFromZetaEBMCpp,
    stop2("Unrecognized method: ", method, "!")
  )

  .g <- switch(method,
    ML  = gradLogLikFromZetaMLCpp,
    EBM = gradLogLikFromZetaEBMCpp,
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

      tryCatch({
        opt_i <- nlminb(
          start     = start,
          objective = \(zeta) -.f(zeta = zeta, y = y, xptr = xptr, idx = idx.p),
          gradient  = \(zeta) -.g(zeta = zeta, y = y, xptr = xptr, idx = idx.p)
        )

        Eta.p[i, ] <- impliedEtaFromZetaCpp(zeta = opt_i$par, xptr = xptr)

      }, error = function(e) {
        warning2(sprintf("nlminb() failed for pattern %i, row %i!", p, i))
        NULL
      })
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
      if (!any(rows.g)) return(NULL)

      mat.g  <- as.matrix(newdata[rows.g, cols, drop = FALSE])
      storage.mode(mat.g) <- "double"

      data.g <- patternizeMissingDataFIML(mat.g)
      attr(data.g, "orig.idx") <- which(rows.g)
      data.g
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
  ovIntNewNms    <- stringr::str_replace_all(ovIntTerms, ":", OP_OV_INT)
  plainStructOVs <- setdiff(structovs, ovIntNewNms)

  missing <- setdiff(plainStructOVs, colnames(newdata))
  stopif(length(missing),
    "Missing columns in `newdata`:\n  ", paste(missing, collapse = ", ")
  )

  for (ov in structovs) {
    newdata[[paste0(TEMP_OV_PREFIX, ov)]] <- newdata[[ov]]
  }

  newdata
}


modsemPredictPreProcessMatrices <- function(matrices, data, pct.fill = 0.01) {
  # Fill zero/NA diagonal entries in Theta with pct.fill * empirical variance.
  # Handles TEMP_OV indicators (theta = 0 by construction) and composite score
  # columns added by convertCompositesToObs (also initialised to 0).
  fillTheta <- function(Theta) {
    if (!NROW(Theta) || !NCOL(Theta))
      return(Theta)

    nm      <- colnames(Theta)
    bad.idx <- which(diag(Theta) <= 0 | is.na(diag(Theta)))

    for (idx in bad.idx) {
      col <- nm[idx]
      Theta[idx, idx] <- pct.fill * var(data$data.full[, col], na.rm = TRUE)
    }

    Theta
  }

  matrices$thetaDelta   <- fillTheta(matrices$thetaDelta)
  matrices$thetaEpsilon <- fillTheta(matrices$thetaEpsilon)

  matrices
}


convertCompositesToObs <- function(matrices, data) {
  # For each composite LV, replace its multi-indicator block in the data with a
  # single composite-score column (score = W' * indicators). The modified lambda
  # gets a loading of 1 from the LV to its score column; theta is set to 0
  # (fillTheta adds a small regularising value afterward). This is the same
  # pattern as TEMP_OV wrapping for structural OVs.
  W.mat <- matrices$W

  has.composites <- NCOL(W.mat) > 0L &&
                    any(apply(W.mat, MARGIN = 2L, FUN = \(x) any(x != 0)))

  if (!has.composites)
    return(list(matrices = matrices, data = data))

  compLVs  <- colnames(W.mat)[apply(W.mat, 2L, \(x) any(x != 0))]
  compInds <- rownames(W.mat)[apply(W.mat, 1L, \(x) any(x != 0))]

  xiInds  <- rownames(matrices$lambdaX) %||% character(0)
  etaInds <- rownames(matrices$lambdaY) %||% character(0)
  xis     <- colnames(matrices$lambdaX) %||% character(0)
  etas    <- colnames(matrices$lambdaY) %||% character(0)

  xiCompLVs  <- intersect(xis,     compLVs)
  etaCompLVs <- intersect(etas,    compLVs)
  xiRegInds  <- setdiff(xiInds,    compInds)
  etaRegInds <- setdiff(etaInds,   compInds)

  xiScoreNms  <- if (length(xiCompLVs))  paste0(".PRED_COMP__", xiCompLVs)  else character(0)
  etaScoreNms <- if (length(etaCompLVs)) paste0(".PRED_COMP__", etaCompLVs) else character(0)

  # Compute composite scores. W.mat rows and dataFull columns are both indexed
  # by indicator names, so named indexing gives correct alignment.
  # NA propagates naturally when any component indicator is missing.
  dataFull  <- data$data.full
  W.comp    <- W.mat[, compLVs, drop = FALSE]
  shared    <- intersect(rownames(W.comp), colnames(dataFull))
  compScore <- dataFull[, shared, drop = FALSE] %*% W.comp[shared, , drop = FALSE]
  colnames(compScore) <- c(xiScoreNms, etaScoreNms)

  # New data matrix: xi regular inds, xi scores, eta regular inds, eta scores.
  # This ordering matches the new lambdaX/lambdaY row ordering below.
  newDataFull <- cbind(
    dataFull[, xiRegInds,    drop = FALSE],
    compScore[, xiScoreNms,  drop = FALSE],
    dataFull[, etaRegInds,   drop = FALSE],
    compScore[, etaScoreNms, drop = FALSE]
  )

  # lambdaX: keep regular-indicator rows; add a row per xi composite with loading 1
  lX.reg  <- matrices$lambdaX[xiRegInds, , drop = FALSE]
  lX.comp <- matrix(0, nrow = length(xiCompLVs), ncol = length(xis),
                    dimnames = list(xiScoreNms, xis))
  for (i in seq_along(xiCompLVs))
    lX.comp[xiScoreNms[i], xiCompLVs[i]] <- 1

  # lambdaY: same for eta composites
  lY.reg  <- matrices$lambdaY[etaRegInds, , drop = FALSE]
  lY.comp <- matrix(0, nrow = length(etaCompLVs), ncol = length(etas),
                    dimnames = list(etaScoreNms, etas))
  for (i in seq_along(etaCompLVs))
    lY.comp[etaScoreNms[i], etaCompLVs[i]] <- 1

  # thetaDelta / thetaEpsilon: drop composite-indicator rows/cols; add score
  # rows/cols initialised to 0 (fillTheta regularises them)
  newThetaDelta   <- expandTheta(matrices$thetaDelta,   keep = xiRegInds,  add = xiScoreNms)
  newThetaEpsilon <- expandTheta(matrices$thetaEpsilon, keep = etaRegInds, add = etaScoreNms)

  # tauX / tauY: tauX is stored as a 1-column matrix; extract as named vectors.
  # tau for a composite score = W' * tau_indicators (weighted sum of intercepts).
  tauX.full <- setNames(as.vector(matrices$tauX), rownames(matrices$tauX) %||% character(0))
  tauY.full <- setNames(as.vector(matrices$tauY), rownames(matrices$tauY) %||% character(0))

  tau.xiComp  <- as.vector(t(W.mat[xiInds,  xiCompLVs,  drop = FALSE]) %*% tauX.full[xiInds])
  tau.etaComp <- as.vector(t(W.mat[etaInds, etaCompLVs, drop = FALSE]) %*% tauY.full[etaInds])

  newTauX <- c(tauX.full[xiRegInds],  setNames(tau.xiComp,  xiScoreNms))
  newTauY <- c(tauY.full[etaRegInds], setNames(tau.etaComp, etaScoreNms))

  matrices$lambdaX      <- rbind(lX.reg, lX.comp)
  matrices$lambdaY      <- rbind(lY.reg, lY.comp)
  matrices$thetaDelta   <- newThetaDelta
  matrices$thetaEpsilon <- newThetaEpsilon
  matrices$tauX         <- newTauX
  matrices$tauY         <- newTauY

  list(matrices = matrices, data = patternizeMissingDataFIML(newDataFull))
}


expandTheta <- function(Theta, keep, add) {
  thetaKept <- Theta[keep, keep, drop = FALSE]
  n.add     <- length(add)
  if (!n.add) return(thetaKept)

  n.keep <- nrow(thetaKept)
  nm     <- c(keep, add)
  out    <- matrix(0, nrow = n.keep + n.add, ncol = n.keep + n.add,
                   dimnames = list(nm, nm))
  out[keep, keep] <- thetaKept
  out
}
