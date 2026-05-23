modsemPredictDA <- function(object, newdata = NULL, method = c("EBM", "ML")) {
  model <- object$model
  submodels <- model$models

  if (is.null(newdata)) {
    data <- object$data
  } else {
    stop2("newdata argument is not implemented (yet)!")
  }

  ETA <- vector("list", length = length(submodels))
  YY  <- vector("list", length = length(submodels))

  for (g in seq_along(submodels)) {

    ETA[[g]] <- modsemPredictEtaDA_Group(
      submodel = submodels[[g]],
      data     = data[[g]],
      method   = method
    )

    YY[[g]] <- YFromEta(
      Eta       = ETA[[g]],
      matrices = submodels[[g]]$matrices
    )
  }

  list(Eta = do.call(rbind, ETA), Y = do.call(rbind, YY))
}


modsemPredictY_DA_Group <- function(object, newdata = NULL, method = c("EBM", "ML")) {
  Eta <- modsemPredictEtaDA(
    object  = object,
    newdata = newdata,
    method  = method
  )

  matrices <-
  LambdaX <- matrices$lambdaX
  LambdaY <- matrices$lambdaY

  Lambda <- diagPartitionedMat(LambdaX, LambdaY)
  Eta %*% t(Lambda)
}


modsemPredictEtaDA_Group <- function(submodel, data, method = c("EBM", "ML")) {
  method <- match.arg(toupper(method), c("EBM", "ML"))

  matrices <- submodel$matrices
  data.split <- data$data.split
  patterns <- data$patterns
  xptr <- modelMatrixCacheCpp(matrices)

  .f <- switch(method,
    ML  = logLikFromZetaMLCpp,
    EBM = logLikFromZetaEBMCpp,
    stop2("Unrecognized method: ", method, "!")
  )

  ETA <- vector("list", length = length(patterns))

  for (p in seq_along(data.split)) {
    data.p <- data.split[[p]]
    pattern.p <- patterns[p,,drop=TRUE]
    idx.p <- which(pattern.p) - 1

    n <- NROW(data.p)
    k <- NCOL(matrices$phi) + NCOL(matrices$psi)
    dim <- c(colnames(matrices$phi), colnames(matrices$psi))

    start <- rep(0, k)
    Eta <- matrix(
      NA_real_, nrow = n, ncol = k,
      dimnames = list(NULL, dim)
    )

    for (i in seq_len(n)) {
      y <- data.p[i,, drop=TRUE]

      objective_i <- \(zeta)
        -.f(zeta = zeta, y = y, xptr = xptr, idx = idx.p)

      opt_i <- nlminb(
        start = start,
        objective = objective_i,
        gradient = NULL # for now
      )

      Eta[i,] <- impliedEtaFromZetaCpp(
        zeta = opt_i$par, xptr = xptr
      )
    }

    ETA[[p]] <- Eta
  }

  do.call(rbind, ETA)
}


YFromEta <- function(Eta, matrices) {
  LambdaX <- matrices$lambdaX
  LambdaY <- matrices$lambdaY

  Lambda <- diagPartitionedMat(LambdaX, LambdaY)
  Eta %*% t(Lambda)
}
