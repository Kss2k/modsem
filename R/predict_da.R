modsemPredictEtaDA <- function(object, newdata = NULL, method = c("EBM", "ML")) {
  model <- object$model
  submodels <- model$models

  if (is.null(newdata)) {
    data <- object$data
  } else {
    stop2("newdata argument is not implemented (yet)!")
  }

  Eta <- vector("list", length = length(submodels))
  for (g in seq_along(submodels)) {

    Eta[[g]] <- modsemPredictEtaDA_Group(
      submodel = submodels[[g]],
      data     = data[[g]],
      method   = method
    )

  }

  do.call(rbind, Eta)
}

modsemPredictEtaDA_Group <- function(submodel, data, method = c("EBM", "ML")) {
  method <- match.arg(toupper(method), c("EBM", "ML"))

  matrices <- submodel$matrices
  data.split <- data$data.split
  patterns <- data$patterns

  for (p in seq_along(data.split)) {
    data.p <- data.split[[p]]
    pattern.p <- patterns[p,,drop=TRUE]
  
    n <- NROW(data.p)
    k <- NCOL(data.p)

    Zeta0 <- matrix(0, nrow = n, ncol = k)
    
    Eta <- predictedLatentScoresCpp(Zeta0, matrices = matrices)
    Y <- predictedObservedScores(Eta, matrices = matrices)
    Theta <- data.p - Theta
  }
}


predictedObservedScores <- function(Eta, matrices) {
  LambdaX <- matrices$lambdaX
  LambdaY <- matrices$lambdaY

  Lambda <- diagPartitionedMat(LambdaX, LambdaY)
  Eta %*% t(Lambda)
}
ovLogLikelihoods <- function(Zeta) NULL
lvLogLikelihoods <- function(Zeta) NULL
