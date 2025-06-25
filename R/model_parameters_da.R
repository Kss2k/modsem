# Global variables
namesParMatrices <- c("lambdaX", "lambdaY", "gammaXi", "gammaEta",
                      "thetaDelta", "thetaEpsilon", "phi", "A",
                      "psi", "tauX", "tauY", "alpha", "beta0", "omegaEtaXi",
                      "omegaXiXi")
namesParMatricesCov <- c("gammaXi", "gammaEta", "A", "psi", "phi")


createTheta <- function(model, start = NULL) {
  etas <- model$info$etas

  listThetaCov <- createThetaCovModel(model$covModel)
  thetaCov     <- listThetaCov$theta
  lavLabelsCov <- listThetaCov$lavLabels
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

  lavLabelsMain <- createLavLabels(M, subset = is.na(allModelValues),
                                   etas = etas)

  thetaMain <- allModelValues[is.na(allModelValues)]
  thetaMain <- fillThetaIfStartNULL(start = start, theta = thetaMain)
  theta     <- c(thetaLabel, thetaCov, thetaMain)

  allLabels <- names(c(totalThetaLabel, thetaCov, thetaMain))
  lavLabels <- combineLavLabels(lavLabelsMain = lavLabelsMain,
                                lavLabelsCov = lavLabelsCov,
                                currentLabels = allLabels)
  
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
  thetaCov <- fillThetaIfStartNULL(start = start, theta = thetaCov)

  list(theta = thetaCov, lavLabels = lavLabelsCov)
}


fillThetaIfStartNULL <- function(start, theta) {
  if (!is.null(start)) return(theta)
  vapply(theta, FUN = function(x) stats::runif(1),
         FUN.VALUE = vector("numeric", 1L))
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
    thetaLabel <- suppressWarnings(calcThetaLabel(thetaLabel, model$constrExprs))
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
  if (is.null(covModel$matrices)) return(covModel)
  M <- covModel$matrices

  lMatrices <- covModel$labelMatrices[namesParMatricesCov]
  pMatrices <- M[namesParMatricesCov]
  M[namesParMatricesCov] <- fillMatricesLabels(pMatrices, lMatrices, thetaLabel)

  M$psi      <- fillSymmetric(M$psi, fetch(theta, "^psi"))
  M$gammaEta <- fillNA_Matrix(M$gammaEta, theta = theta, pattern = "^gammaEta")
  M$gammaXi  <- fillNA_Matrix(M$gammaXi, theta = theta, pattern = "^gammaXi")

  if (method == "lms") {
    M$A <- fillNA_Matrix(M$A, theta = theta, pattern = "^A[0-9]+")
  } else if (method == "qml") {
    M$phi <- fillSymmetric(M$phi, fetch(theta, "^phi"))
  }

  if (fillPhi) M$phi <- M$A %*% t(M$A)

  covModel$matrices <- M
  covModel
}


fillNA_Matrix <- function(X, theta, pattern) {
  X[is.na(X) & !is.nan(X)] <- fetch(theta, pattern)
  X
}


fillSymmetric <- function(mat, values) {
  mat[is.na(mat) & !is.nan(mat)] <- values
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}


# Set bounds for parameters to (0, Inf)
getParamBounds <- function(model, lowest = 0, varParams=NULL) {
  lower <- rep(-Inf, model$freeParams)
  upper <- rep(Inf, model$freeParams)
  names(lower) <- names(upper) <- names(model$theta)
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


calcPhiTheta <- function(theta, model, method) {
  if (method != "lms") return(theta)
  filledModel <- fillModel(theta = theta, model = model, method = method,
                           fillPhi = TRUE)

  if (!is.null(model$covModel$matrices)) {
    matEst <- filledModel$covModel$matrices
    matLab <- model$covModel$labelMatrices
    matNA  <- model$covModel$matrices
  } else {
    matEst <- filledModel$matrices
    matNA  <- model$matrices
    matLab <- model$labelMatrices
  }

  vals   <- as.vector(matEst$phi[is.na(matNA$A) & !is.nan(matNA$A)])
  labels <- as.vector(matLab$A)

  if (any(labels != "")) {
    allVals <- as.vector(matEst$phi)
    labVals <- allVals[labels != ""]
    labels  <- labels[labels != ""]
    theta[labels] <- labVals
  }

  theta[grepl("^A[0-9]+$", names(theta))] <- vals
  theta
}


LMS_BLOCKS = list(
  lambdaX      = 0,
  lambdaY      = 1,
  tauX         = 2,
  tauY         = 3,
  thetaDelta   = 4,
  thetaEpsilon = 5,
  A            = 6,
  psi          = 7,
  alpha        = 8, 
  beta0        = 9,
  gammaXi      = 10,
  gammaEta     = 11, 
  omegaXiXi    = 12,
  omegaEtaXi   = 13,
  phi          = NA
)


getParamNamesMatrix <- function(mat, matname) {
  if (is.character(mat)) {
    c(mat)

  } else {
    M <- list()
    M[[matname]] <- mat
    names(unlist(M))
  }
}


getParamLocationsMatrices <- function(matrices, isFree = is.na) {

  locations <- data.frame(param = NULL, block = NULL, row = NULL, col = NULL)
  for (blockname in names(matrices)) {
    X <- matrices[[blockname]]
    n <- nrow(X) 
    m <- ncol(X)

    if (!any(isFree(X))) next

    params <- getParamNamesMatrix(mat = X, matname = blockname)
    block  <- LMS_BLOCKS[[blockname]]
    rowidx <- matrix(seq_len(n) - 1, nrow = n, ncol = m, byrow = FALSE)
    colidx <- matrix(seq_len(m) - 1, nrow = n, ncol = m, byrow = TRUE)

    params <- params[isFree(X)]
    rowidx <- rowidx[isFree(X)]
    colidx <- colidx[isFree(X)]

    locationsBlock <- data.frame(
      param = params,
      block = block,
      row   = rowidx,
      col   = colidx
    )

    locations <- rbind(locations, locationsBlock)
  }

  locations
}


getGradientStruct <- function(model, theta) {
  tryCatch(
    getGradientStructSimple(model = model, theta = theta),
    error = function(e) {
      warning2("Failed to compute gradient structure: ", e$message)
      
      list(
        locations   = NULL, 
        Jacobian    = NULL,
        nlinDerivs  = NULL,
        evalTheta   = NULL,
        hasCovModel = TRUE, # may not be true, but we should behave as if it is
        isNonLinear = TRUE  # may not be true, but we should behave as if it is
      )
    }
  )
}


getGradientStructSimple <- function(model, theta) {
  hasCovModel <- !is.null(model$covModel$matrices)

  if (hasCovModel) {
    out <- list(
      locations   = NULL, 
      Jacobian    = NULL,
      nlinDerivs  = NULL,
      evalTheta   = NULL,
      hasCovModel = TRUE, 
      isNonLinear = TRUE  # may not be true, but we should behave as if it is
    )

    return(out)
  }

  parTable <- model$parTable

  isConstraint <- parTable$op %in% CONSTRAINT_OPS
  constraints  <- parTable[isConstraint, ]
  restParTable <- parTable[!isConstraint, ]
  constraints  <- constraints[constraints$lhs %in% restParTable$mod, ]

  derivatives <- list()
  for (i in seq_len(NROW(constraints))) {
    constrVar <- constraints[i, "lhs"]
    constrEq  <- constraints[i, "rhs"]

    derivatives[[constrVar]] <- derivateConstraint(constrEq)
  }

  isLinear <- vapply(derivatives, FUN.VALUE = logical(1L), FUN = is.atomic)

  linDerivs  <- derivatives[isLinear]
  nlinDerivs <- derivatives[!isLinear]
  evalTheta  <- \(theta) c(theta, suppressWarnings(calcThetaLabel(theta, constraints))) # This could be made a bit better

  locations <- rbind(
    getParamLocationsMatrices(model$matrices, isFree=is.na),
    getParamLocationsMatrices(model$labelMatrices, isFree=\(x) x != "")
  )

  k <- nrow(locations)
  m <- length(theta)

  locations  <- locations[sample(k), ]
  param.full <- locations$param
  param.part <- names(theta)

  ordering <- structure(seq_along(theta), names = param.part)
  ordering <- ordering[param.full]
  
  locations  <- locations[order(ordering), ]
  param.full <- locations$param

  Jacobian <- matrix(0, nrow = m, ncol = k,  
                     dimnames = list(param.part, param.full))

  for (par in param.full) {
    match.full <- param.full == par
    match.part <- param.part == par

    Jacobian[match.part, match.full] <- 1
  }

  for (dep in names(linDerivs)) {
    deriv <- linDerivs[[dep]]

    for (indep in names(deriv)) {
      match.full <- param.full == dep
      match.part <- param.part == indep
      Jacobian[match.part, match.full] <- deriv[[indep]]
    }
  }

  list(
    locations   = locations, 
    Jacobian    = Jacobian,
    nlinDerivs  = nlinDerivs,
    evalTheta   = evalTheta,
    hasCovModel = hasCovModel,
    isNonLinear = length(nlinDerivs) > 1
  )
}


derivateConstraint <- function(constr) {
  f <- stats::formula(paste0("~", constr))
  eq <- Deriv::Deriv(f)

  if (is.null(names(eq))) names(eq) <- all.vars(f)

  eq
}
