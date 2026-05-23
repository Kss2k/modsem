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

  list(Eta = do.call(rbind, ETA), Y = do.call(rbind, YY))
}


modsemPredictEtaDA_Group <- function(submodel, data, method = c("EBM", "ML")) {
  method <- match.arg(toupper(method), c("EBM", "ML"))

  matrices   <- submodel$matrices
  data.split <- data$data.split
  patterns   <- data$patterns
  rowidx     <- data$rowidx
  n.total    <- data$n
  xptr       <- modelMatrixCacheCpp(matrices)

  .f <- switch(method,
    ML  = logLikFromZetaMLCpp,
    EBM = logLikFromZetaEBMCpp,
    stop2("Unrecognized method: ", method, "!")
  )

  k   <- NCOL(matrices$phi) + NCOL(matrices$psi)
  dim <- c(colnames(matrices$phi), colnames(matrices$psi))

  Eta.out <- matrix(NA_real_, nrow = n.total, ncol = k, dimnames = list(NULL, dim))

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

    Eta.out[rowidx[[p]], ] <- Eta.p
  }

  Eta.out
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
