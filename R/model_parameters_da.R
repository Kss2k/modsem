# Global variables 
namesParMatrices <- c("lambdaX", "lambdaY", "gammaXi", "gammaEta", 
                      "thetaDelta", "thetaEpsilon", "phi", "A",
                      "psi", "tauX", "tauY", "alpha", "beta0", "omegaEtaXi", 
                      "omegaXiXi")
namesParMatricesCov <- c("gammaXi", "gammaEta", "A", "psi", "phi")


createTheta <- function(model, start = NULL) {
  listThetaCov <- createThetaCovModel(model$covModel)
  thetaCov     <- listThetaCov$theta
  lavLabelsCov <- listThetaCov$labels
  thetaLabel   <- createThetaLabel(model$labelMatrices, 
                                   model$covModel$labelMatrices,
                                   model$constrExprs)
  totalThetaLabel <- calcThetaLabel(thetaLabel, model$constrExprs)

  M  <- model$matrices
  lambdaX      <- as.vector(M$lambdaX)
  lambdaY      <- as.vector(M$lambdaY)
  thetaDelta   <- as.vector(M$thetaDelta)
  thetaEpsilon <- as.vector(M$thetaEpsilon)
  phi          <- as.vector(M$phi)
  A            <- as.vector(M$A)
  psi          <- as.vector(M$psi)
  tauX         <- as.vector(M$tauX)
  tauY         <- as.vector(M$tauY)
  alpha        <- as.vector(M$alpha)
  beta0        <- as.vector(M$beta0)
  gammaXi      <- as.vector(M$gammaXi)
  gammaEta     <- as.vector(M$gammaEta)
  omegaXiXi    <- as.vector(M$omegaXiXi)
  omegaEtaXi   <- as.vector(M$omegaEtaXi)

  allModelValues <- c("lambdaX" = lambdaX,
                      "lambdaY" = lambdaY,
                      "tauX" = tauX,
                      "tauY" = tauY,
                      "thetaDelta" = thetaDelta,
                      "thetaEpsilon" = thetaEpsilon,
                      "phi" = phi,
                      "A" = A,
                      "psi" = psi,
                      "alpha" = alpha,
                      "beta0" = beta0,
                      "gammaXi" = gammaXi,
                      "gammaEta" = gammaEta,
                      "omegaXiXi" = omegaXiXi,
                      "omegaEtaXi" = omegaEtaXi)
  
  lavLabelsMain <- createLavLabels(M, subset=is.na(allModelValues))
  lavLabels <- c(lavLabelsCov, lavLabelsMain) 

  thetaMain <- allModelValues[is.na(allModelValues)]
  if (is.null(start)) {
    thetaMain <- vapply(thetaMain, FUN.VALUE = vector("numeric", 1L),
                        FUN = function(x) stats::runif(1))
  }
  
  theta <- c(thetaLabel, thetaCov, thetaMain)

  list(theta = theta, lenThetaMain = length(thetaMain),
       lenThetaLabel = length(thetaLabel),
       totalLenThetaLabel = length(totalThetaLabel),
       lenThetaCov = length(thetaCov), lavLabels = lavLabels)
}


createThetaCovModel <- function(covModel, start = NULL) {
  M <- covModel$matrices

  phi      <- as.vector(M$phi)
  A        <- as.vector(M$A)
  psi      <- as.vector(M$psi)
  alpha    <- as.vector(M$alpha)
  gammaXi  <- as.vector(M$gammaXi)
  gammaEta <- as.vector(M$gammaEta)
  thetaCov <- c("phi" = phi,
                "A" = A,
                "psi" = psi,
                "gammaXi" = gammaXi,
                "gammaEta" = gammaEta)
  
  lavLabelsCov <- createLavLabelsCov(M, subset = is.na(thetaCov))
  thetaCov <- thetaCov[is.na(thetaCov)]

  if (is.null(start)) {
   thetaCov <- vapply(thetaCov, FUN.VALUE = vector("numeric", 1L),
                   FUN = function(x) stats::runif(1))
  }
  
  list(theta = thetaCov, labels = lavLabelsCov)
}


fillModel <- function(model, theta, fillPhi = FALSE, method = "lms") {
  if (is.null(names(theta))) names(theta) <- names(model$theta)

  # labeled parameters
  thetaLabel <- NULL
  if (model$totalLenThetaLabel > 0) {
    if (model$lenThetaLabel > 0) {
      thetaLabel <- theta[seq_len(model$lenThetaLabel)] 
      theta <- theta[-seq_len(model$lenThetaLabel)]
    } 
    thetaLabel <- calcThetaLabel(thetaLabel, model$constrExprs)
  }

  # cov model
  thetaCov <- NULL
  thetaMain <- theta
  if (model$lenThetaCov > 0) {
    thetaCov <- theta[seq_len(model$lenThetaCov)]
    thetaMain <- theta[-seq_len(model$lenThetaCov)]
  }

  model$covModel <- fillCovModel(model$covModel, thetaCov, thetaLabel,
                                 fillPhi = fillPhi, method = method)
  model$matrices <- fillMainModel(model, thetaMain, thetaLabel,
                                  fillPhi = fillPhi, method = method)
  model
}


fillMainModel <- function(model, theta, thetaLabel, fillPhi = FALSE, 
                          method = "lms") {
  xis      <- model$info$xis
  numXis   <- model$info$numXis
  numEtas  <- model$info$numEtas
  M        <- model$matrices
  covModel <- model$covModel

  lMatrices <- model$labelMatrices[namesParMatrices]
  pMatrices <- M[namesParMatrices]
  M[namesParMatrices] <- fillMatricesLabels(pMatrices, lMatrices, thetaLabel)

  if (!is.null(model$covModel$matrices)) {
    M$phi <- M$A <- expectedCovModel(covModel, method = method, sortedXis = xis)
  } else if (method == "lms") {
    M$A <- fillNA_Matrix(M$A, theta = theta, pattern = "^A[0-9]*$")
  } else if (method == "qml") {
    M$phi <- fillSymmetric(M$phi, fetch(theta, "^phi"))
  }

  M$lambdaX      <- fillNA_Matrix(M$lambdaX, theta = theta, pattern = "^lambdaX")
  M$lambdaY      <- fillNA_Matrix(M$lambdaY, theta = theta, pattern = "^lambdaY")
  M$thetaDelta   <- fillSymmetric(M$thetaDelta, fetch(theta, "^thetaDelta"))
  M$thetaEpsilon <- fillSymmetric(M$thetaEpsilon, fetch(theta, "thetaEpsilon"))
  M$psi          <- fillSymmetric(M$psi, fetch(theta, "^psi"))
  M$tauX         <- fillNA_Matrix(M$tauX, theta = theta, pattern = "^tauX")
  M$tauY         <- fillNA_Matrix(M$tauY, theta = theta, pattern = "^tauY")
  M$alpha        <- fillNA_Matrix(M$alpha, theta = theta, pattern = "^alpha")
  M$beta0        <- fillNA_Matrix(M$beta0, theta = theta, pattern = "^beta0")
  M$gammaEta     <- fillNA_Matrix(M$gammaEta, theta = theta, pattern = "^gammaEta")
  M$gammaXi      <- fillNA_Matrix(M$gammaXi, theta = theta, pattern = "^gammaXi")
  M$omegaXiXi    <- fillNA_Matrix(M$omegaXiXi, theta = theta, pattern = "^omegaXiXi")
  M$omegaEtaXi   <- fillNA_Matrix(M$omegaEtaXi, theta = theta, pattern = "^omegaEtaXi")
  
  if (fillPhi) M$phi <- M$A %*% t(M$A)
  M
}


fillCovModel <- function(covModel, theta, thetaLabel, fillPhi = FALSE, 
                         method = "lms") {
  if (is.null(names(theta))) names(theta) <- names(covModel$theta)
  if (is.null(covModel$matrices)) return(NULL)
  M <- covModel$matrices 

  lMatrices <- covModel$labelMatrices[namesParMatricesCov]
  pMatrices <- M[namesParMatricesCov]
  M[namesParMatricesCov] <- fillMatricesLabels(pMatrices, lMatrices, thetaLabel)

  M$psi      <- fillSymmetric(M$psi, fetch(theta, "^psi"))
  M$gammaEta <- fillNA_Matrix(M$gammaEta, theta = theta, pattern = "^gammaEta")
  M$gammaXi  <- fillNA_Matrix(M$gammaXi, theta = theta, pattern = "^gammaXi")
  
  if (method == "lms") {
    M$A <- fillSymmetric(M$A, fetch(theta, "^A[0-9]+"))
  } else if (method == "qml") {
    M$phi <- fillSymmetric(M$phi, fetch(theta, "^phi"))
  }

  if (fillPhi) M$phi <- M$A %*% t(M$A)
  
  covModel$matrices <- M 
  covModel
}


fillNA_Matrix <- function(X, theta, pattern) {
  X[is.na(X)] <- fetch(theta, pattern)
  X
}


fillSymmetric <- function(mat, values) {
  mat[is.na(mat)] <- values
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}


# Set bounds for parameters to (0, Inf)
getParamBounds <- function(model, lowest = 1e-6) {
  namePattern <- paste0("lambdaX[0-9]*$|lambdaY[0-9]*$|",
                        "thetaDelta[0-9]*$|thetaEpsilon[0-9]*$|",
                        "phi[0-9]*$|psi[0-9]*$")
  lower <- rep(-Inf, model$freeParams)
  upper <- rep(Inf, model$freeParams)
  names(lower) <- names(upper) <- names(model$theta)
  lower[grepl(namePattern, names(lower))] <- lowest
  list(lower = lower, upper = upper)
}


checkStartingParams <- function(start, model) {
  if (length(start) != length(model$theta)) {
    stop2("The length of the starting parameters does not match the number of parameters in the model")
  }
  if (is.null(names(start))) {
    names(start) <- names(model$theta)
  }
  if (!all(names(start) %in% names(model$theta))) {
    stop2("The names of the starting parameters do not match the names of the parameters in the model")
  }

  NULL
}
