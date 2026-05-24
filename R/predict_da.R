#' @describeIn modsem_predict
#' Computes factor scores or model-implied observed-variable scores for a
#' \code{\link{modsem_da}} model via MAP optimisation.
#' @export
modsem_predict.modsem_da <- function(object,
                                     newdata = NULL,
                                     method = c("EBM", "ML", "Bartlett", "Regression"),
                                     type = c("lv", "ov", "all"),
                                     standardized = FALSE,
                                     se = FALSE,
                                     drop.list.single.group = TRUE,
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
    method  = method,
    se      = se,
    drop.list.single.group = drop.list.single.group
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
    all = predicted$All
  )

  if (standardized) {
    predAttrs <- attributes(out)[c("se", "acov")]

    scaled <- scale(out)
    sc     <- attr(scaled, "scaled:scale")
    out    <- modsemMatrix(scaled)

    if (!is.null(predAttrs$se))
      attr(out, "se") <- sweep(predAttrs$se, MARGIN = 2L,
                               STATS = sc[colnames(predAttrs$se)], FUN = "/")
    if (!is.null(predAttrs$acov))
      attr(out, "acov") <- predictDA_RescaleAcov(predAttrs$acov, sc)
  }

  out
}


modsemPredictDA <- function(object, newdata = NULL, method = c("EBM", "ML"),
                            se = FALSE, drop.list.single.group = TRUE) {
  model    <- object$model
  submodels <- model$models

  if (is.null(newdata)) {
    DATA <- object$data
  } else {
    DATA <- prepNewdataDA(object, newdata)
  }

  ETA <- vector("list", length = length(submodels))
  YY  <- vector("list", length = length(submodels))
  ALL <- vector("list", length = length(submodels))

  for (g in seq_along(submodels)) {
    if (is.null(DATA[[g]])) next

    ETA[[g]] <- modsemPredictEtaDA_Group(
      submodel = submodels[[g]],
      data     = DATA[[g]],
      method   = method,
      se       = se
    )

    YY[[g]] <- YFromEta(
      Eta      = ETA[[g]],
      matrices = submodels[[g]]$matrices
    )

    ALL[[g]] <- cbind(
      ETA[[g]],
      YY[[g]][, setdiff(colnames(YY[[g]]), colnames(ETA[[g]])), drop = FALSE]
    )

    if (se) {
      matrices.g <- submodels[[g]]$matrices
      Lambda.g   <- predictDALambda(matrices.g)
      A.eta.g    <- attr(ETA[[g]], "acov")

      ySe <- predictDAEmptySE(
        n = nrow(YY[[g]]),
        vars = colnames(YY[[g]]),
        npatterns = length(A.eta.g)
      )
      yAcov <- vector("list", length = length(A.eta.g))

      allSe <- predictDAEmptySE(
        n = nrow(ALL[[g]]),
        vars = colnames(ALL[[g]]),
        npatterns = length(A.eta.g)
      )
      allAcov <- vector("list", length = length(A.eta.g))

      yKeep <- match(colnames(ALL[[g]])[-seq_len(ncol(ETA[[g]]))], colnames(YY[[g]]))
      Aall  <- rbind(diag(ncol(ETA[[g]])), Lambda.g[yKeep, , drop = FALSE])
      rownames(Aall) <- colnames(ALL[[g]])

      pattern.g <- attr(ETA[[g]], "pattern")

      for (p in seq_along(A.eta.g)) {
        if (is.null(A.eta.g[[p]])) next

        rows.p <- pattern.g == p
        Vy <- Lambda.g %*% A.eta.g[[p]] %*% t(Lambda.g)
        dimnames(Vy) <- list(colnames(YY[[g]]), colnames(YY[[g]]))

        Va <- Aall %*% A.eta.g[[p]] %*% t(Aall)
        dimnames(Va) <- list(colnames(ALL[[g]]), colnames(ALL[[g]]))

        yAcov[[p]]   <- Vy
        ySe[predictDASeRows(rows.p, length(A.eta.g)), ] <- matrix(
          predictDASqrtDiag(Vy),
          nrow = if (length(A.eta.g) == 1L) 1L else sum(rows.p),
          ncol = ncol(ySe),
          byrow = TRUE
        )

        allAcov[[p]] <- Va
        allSe[predictDASeRows(rows.p, length(A.eta.g)), ] <- matrix(
          predictDASqrtDiag(Va),
          nrow = if (length(A.eta.g) == 1L) 1L else sum(rows.p),
          ncol = ncol(allSe),
          byrow = TRUE
        )
      }

      names(yAcov) <- names(A.eta.g)
      names(allAcov) <- names(A.eta.g)

      attr(YY[[g]], "se")   <- ySe
      attr(YY[[g]], "acov") <- yAcov
      attr(attr(YY[[g]], "acov"), "patterns") <- attr(A.eta.g, "patterns")
      attr(YY[[g]], "pattern") <- pattern.g

      attr(ALL[[g]], "se")   <- allSe
      attr(ALL[[g]], "acov") <- allAcov
      attr(attr(ALL[[g]], "acov"), "patterns") <- attr(A.eta.g, "patterns")
      attr(ALL[[g]], "pattern") <- pattern.g
    }
  }

  # collect
  Eta <- do.call(rbind, ETA)
  Y   <- do.call(rbind, YY)
  All <- do.call(rbind, ALL)

  if (se) {
    SE.eta <- predictDACollectSEAttr(ETA, "se")
    SE.y   <- predictDACollectSEAttr(YY,  "se")
    SE.all <- predictDACollectSEAttr(ALL, "se")

    ACOV.eta <- predictDACollectGroupAttr(ETA, "acov")
    ACOV.y   <- predictDACollectGroupAttr(YY,  "acov")
    ACOV.all <- predictDACollectGroupAttr(ALL, "acov")

    group.names <- model$info$group.levels %||% paste0("group", seq_along(ACOV.eta))
    group.names <- group.names[seq_along(ACOV.eta)]

    names(ACOV.eta) <- names(ACOV.y) <- names(ACOV.all) <- group.names

  }

  # Restore original row order when newdata was split by group.
  # Each non-NULL DATA[[g]] carries an "orig.idx" attribute recording its rows'
  # positions in newdata; unlist() collects them in group-processing order so
  # order() gives the permutation back to input order.
  orig.idx <- unlist(lapply(DATA, attr, "orig.idx"))

  if (length(orig.idx)) {
    ord <- order(orig.idx)
    Eta <- Eta[ord, , drop = FALSE]
    Y   <- Y[ord, , drop = FALSE]
    All <- All[ord, , drop = FALSE]

    if (se) {
      if (nrow(SE.eta) == nrow(Eta)) SE.eta <- SE.eta[ord, , drop = FALSE]
      if (nrow(SE.y) == nrow(Y))     SE.y   <- SE.y[ord, , drop = FALSE]
      if (nrow(SE.all) == nrow(All)) SE.all <- SE.all[ord, , drop = FALSE]
    }
  }

  if (se) {
    Eta <- modsemMatrix(Eta, is.public = TRUE)
    Y   <- modsemMatrix(Y, is.public = TRUE)
    All <- modsemMatrix(All, is.public = TRUE)

    attr(Eta, "se")           <- modsemMatrix(SE.eta, is.public = TRUE)
    attr(Eta, "acov")         <- predictDAFormatGroupVcovList(
      ACOV.eta, drop.list.single.group = drop.list.single.group
    )

    attr(Y, "se")   <- modsemMatrix(SE.y, is.public = TRUE)
    attr(Y, "acov") <- predictDAFormatGroupVcovList(
      ACOV.y, drop.list.single.group = drop.list.single.group
    )

    attr(All, "se")   <- modsemMatrix(SE.all, is.public = TRUE)
    attr(All, "acov") <- predictDAFormatGroupVcovList(
      ACOV.all, drop.list.single.group = drop.list.single.group
    )

  } else {
    Eta <- modsemMatrix(Eta, is.public = TRUE)
    Y   <- modsemMatrix(Y, is.public = TRUE)
    All <- modsemMatrix(All, is.public = TRUE)

  }

  list(Eta = Eta, Y = Y, All = All)
}


modsemPredictEtaDA_Group <- function(submodel, data, method = c("EBM", "ML"), se = FALSE) {
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
    mod_msg_stop(paste0("Unrecognized method: ", method, "!"))
  )

  .g <- switch(method,
    ML  = gradLogLikFromZetaMLCpp,
    EBM = gradLogLikFromZetaEBMCpp,
    mod_msg_stop(paste0("Unrecognized method: ", method, "!"))
  )

  .h <- switch(method,
    ML  = hessLogLikFromZetaMLCpp,
    EBM = hessLogLikFromZetaEBMCpp,
    mod_msg_stop(paste0("Unrecognized method: ", method, "!"))
  )

  k   <- NCOL(matrices$phi) + NCOL(matrices$psi)
  dim <- c(colnames(matrices$phi), colnames(matrices$psi))

  Eta <- matrix(NA_real_, nrow = n.total, ncol = k, dimnames = list(NULL, dim))

  if (se) {
    SE <- predictDAEmptySE(
      n = n.total,
      vars = colnames(Eta),
      npatterns = length(data.split)
    )
    ACOV <- vector("list", length(data.split))
    PATTERN <- integer(n.total)
  }

  for (p in seq_along(data.split)) {
    data.p    <- data.split[[p]]
    pattern.p <- patterns[p, , drop = TRUE]
    idx.p     <- which(pattern.p) - 1L

    n     <- NROW(data.p)
    start <- rep(0, k)
    Eta.p <- matrix(NA_real_, nrow = n, ncol = k)

    if (se) {
      VCOV.p  <- vector("list", n)
      SE.p <- predictDAEmptySE(
        n = n,
        vars = colnames(Eta),
        npatterns = length(data.split)
      )
    }

    for (i in seq_len(n)) {
      y <- data.p[i, , drop = TRUE]

      tryCatch({
        opt_i <- stats::nlminb(
          start     = start,
          objective = \(zeta) -.f(zeta = zeta, y = y, xptr = xptr, idx = idx.p),
          gradient  = \(zeta) -.g(zeta = zeta, y = y, xptr = xptr, idx = idx.p)
        )

        Eta.p[i, ] <- impliedEtaFromZetaCpp(zeta = opt_i$par, xptr = xptr)

        if (se) {
          H <- .h(
            zeta = opt_i$par,
            y    = y,
            xptr = xptr,
            idx  = idx.p
          )

          info.zeta <- -H$Hessian
          dimnames(info.zeta) <- list(dim, dim)

          vcov.zeta <- solveFIM(info.zeta)
          dimnames(vcov.zeta) <- list(dim, dim)

          J.eta <- jacobianEtaFromZeta(
            zeta = opt_i$par,
            xptr = xptr
          )
          dimnames(J.eta) <- list(dim, dim)

          vcov.eta <- J.eta %*% vcov.zeta %*% t(J.eta)
          dimnames(vcov.eta) <- list(dim, dim)

          VCOV.p[[i]]  <- vcov.eta
        }

      }, error = function(e) {
        mod_msg_warn(sprintf("nlminb() failed for pattern %i, row %i!", p, i))
        NULL
      })
    }

    Eta[rowidx[[p]], ] <- Eta.p

    if (se) {
      ACOV[[p]]  <- predictDAMeanMatrix(VCOV.p)

      SE.p[TRUE]  <- rep(predictDASqrtDiag(ACOV[[p]]), each = nrow(SE.p))

      SE[predictDASeRows(rowidx[[p]], length(data.split)), ] <- SE.p
      PATTERN[rowidx[[p]]]   <- p
    }
  }

  if (se) {
    pattern.names <- paste0("pattern", seq_along(ACOV))
    names(ACOV) <- pattern.names
    attr(ACOV, "patterns") <- predictDAPatternTable(patterns, pattern.names)

    attr(Eta, "se")               <- SE
    attr(Eta, "acov")             <- ACOV
    attr(Eta, "pattern")          <- PATTERN
  }

  Eta
}


jacobianEtaFromZeta <- function(zeta, xptr, relStep = 1e-4, minAbsPar = 0.0) {
  p <- length(zeta)
  J <- matrix(0, nrow = p, ncol = p)

  incr <- pmax(abs(zeta), minAbsPar) * relStep
  incr[incr == 0] <- relStep

  for (j in seq_len(p)) {
    zp <- zeta
    zm <- zeta
    zp[j] <- zp[j] + incr[j]
    zm[j] <- zm[j] - incr[j]

    J[, j] <- (
      impliedEtaFromZetaCpp(zeta = zp, xptr = xptr) -
      impliedEtaFromZetaCpp(zeta = zm, xptr = xptr)
    ) / (2 * incr[j])
  }

  J
}


predictDALambda <- function(matrices) {
  diagPartitionedMat(matrices$lambdaX, matrices$lambdaY)
}


predictDACollectMatrixAttr <- function(x, which) {
  mats <- lapply(x, attr, which = which)
  mats <- mats[lengths(mats) > 0L]
  do.call(rbind, mats)
}


predictDACollectSEAttr <- function(x, which) {
  mats <- lapply(x, attr, which = which)
  mats <- mats[lengths(mats) > 0L]

  if (all(vapply(mats, nrow, integer(1)) == 1L)) {
    out <- do.call(rbind, mats)
    rownames(out) <- names(mats)
    return(out)
  }

  do.call(rbind, mats)
}


predictDACollectGroupAttr <- function(x, which) {
  out <- lapply(x, attr, which = which)
  out[lengths(out) > 0L]
}


predictDAFormatGroupVcovList <- function(x, drop.list.single.group = TRUE) {
  out <- lapply(x, \(g) {
    if (is.null(g)) return(NULL)

    if (is.matrix(g))
      return(modsemMatrix(g, symmetric = TRUE, is.public = TRUE))

    pattern.table <- attr(g, "patterns")

    if (length(g) == 1L)
      return(modsemMatrix(g[[1L]], symmetric = TRUE, is.public = TRUE))

    out.g <- lapply(g, \(v) {
      if (is.null(v)) return(NULL)
      modsemMatrix(v, symmetric = TRUE, is.public = TRUE)
    })

    attr(out.g, "patterns") <- pattern.table
    out.g
  })

  if (drop.list.single.group && length(out) == 1L) out[[1L]] else out
}


predictDAMeanMatrix <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(NULL)

  out <- Reduce(`+`, x) / length(x)
  dimnames(out) <- dimnames(x[[1L]])
  out
}


predictDAEmptySE <- function(n, vars, npatterns) {
  matrix(
    NA_real_,
    nrow = if (npatterns == 1L) 1L else n,
    ncol = length(vars),
    dimnames = list(NULL, vars)
  )
}


predictDASeRows <- function(rows, npatterns) {
  if (npatterns == 1L) TRUE else rows
}


predictDAPatternTable <- function(patterns, pattern.names) {
  patterns <- as.matrix(patterns)
  rownames(patterns) <- pattern.names

  colnames(patterns) <- stringr::str_replace_all(
    string = stringr::str_remove_all(colnames(patterns), pattern = TEMP_OV_PREFIX),
    pattern = OP_REPLACEMENTS_INV
  )

  patterns
}


predictDASqrtDiag <- function(x) {
  if (is.null(x)) return(rep(NA_real_, 0L))

  v <- diag(x)
  v[v < 0] <- NA_real_
  sqrt(v)
}


predictDA_RescaleAcov <- function(acov, sc) {
  if (is.matrix(acov)) {
    d <- 1 / sc[rownames(acov)]
    return(sweep(sweep(acov, 1L, d, "*"), 2L, d, "*"))
  }

  patterns.attr <- attr(acov, "patterns")
  out <- lapply(acov, \(v) if (is.null(v)) NULL else predictDA_RescaleAcov(v, sc))
  attr(out, "patterns") <- patterns.attr
  out
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
    mod_stopif(!group.var %in% colnames(newdata),
           sprintf("Grouping variable '%s' not found in `newdata`.", group.var))

    group.levels <- object$model$info$group.levels
    mod_stopif(is.null(group.levels),
           "Model has a grouping variable but group levels could not be determined.")

    group.values <- as.character(newdata[[group.var]])

    DATA <- lapply(seq_len(n.groups), function(g) {
      cols <- colnames(submodels[[g]]$data$data.full)

      missing.cols <- setdiff(cols, colnames(newdata))
      mod_stopif(length(missing.cols),
             paste0("Missing columns in `newdata`:\n  ", paste(missing.cols, collapse = ", ")))

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
    mod_stopif(length(missing.cols),
           paste0("Missing columns in `newdata`:\n  ", paste(missing.cols, collapse = ", ")))

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
    mod_stopif(length(missing),
           paste0("Missing columns in `newdata` required for interaction term '", ovInt, "':\n  ",
           paste(missing, collapse = ", ")))

    newdata[[ovIntNew]] <- apply(newdata[, vars, drop = FALSE], MARGIN = 1L, FUN = prod)
  }

  # Plain structural OVs: check they are present, then add TEMP_OV copy
  ovIntNewNms    <- stringr::str_replace_all(ovIntTerms, ":", OP_OV_INT)
  plainStructOVs <- setdiff(structovs, ovIntNewNms)

  missing <- setdiff(plainStructOVs, colnames(newdata))
  mod_stopif(length(missing),
    paste0("Missing columns in `newdata`:\n  ", paste(missing, collapse = ", "))
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
      Theta[idx, idx] <- pct.fill * stats::var(data$data.full[, col], na.rm = TRUE)
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
  newThetaDelta   <- expandTheta(
    matrices$thetaDelta, keep = xiRegInds,  add = xiScoreNms
  )
  newThetaEpsilon <- expandTheta(
    matrices$thetaEpsilon, keep = etaRegInds, add = etaScoreNms
  )

  # tauX / tauY: tauX is stored as a 1-column matrix; extract as named vectors.
  # tau for a composite score = W' * tau_indicators (weighted sum of intercepts).
  tauX.full <- stats::setNames(
    as.vector(matrices$tauX), rownames(matrices$tauX) %||% character(0)
  )
  tauY.full <- stats::setNames(
    as.vector(matrices$tauY), rownames(matrices$tauY) %||% character(0)
  )

  tau.xiComp  <- as.vector(
    t(W.mat[xiInds,  xiCompLVs,  drop = FALSE]) %*% tauX.full[xiInds]
  )
  tau.etaComp <- as.vector(
    t(W.mat[etaInds, etaCompLVs, drop = FALSE]) %*% tauY.full[etaInds]
  )

  newTauX <- c(tauX.full[xiRegInds],  stats::setNames(tau.xiComp,  xiScoreNms))
  newTauY <- c(tauY.full[etaRegInds], stats::setNames(tau.etaComp, etaScoreNms))

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
