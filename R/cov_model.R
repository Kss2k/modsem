# functions for computing the constrained covariance matrix, 
# based on causal relationships. This can make the lms method more flexible, 
# as you can split the model into a non-linear, and linear part. allowing 
# you to use (normally distributed) endogenous variables as non-normal
# as of now the mean-structure is excluded


# global variables 
paramMatricesCov <- c("gammaXi", "gammaEta", "A", "psi", "phi")


# functions
covModel <- function(syntax = NULL, method = "lms", parTable = NULL) {
  if (!is.null(syntax)) parTable <- modsemify(syntax)
  if (is.null(parTable)) {
    return(list(matrices = NULL, freeParams = 0, 
                theta = NULL, syntax = NULL, parTable = NULL))
  }

  etas <- getSortedEtas(parTable, isLV = FALSE, checkAny = TRUE)
  numEtas <- length(etas)
  xis <- getXis(parTable, checkAny = TRUE, isLV = FALSE)
  numXis <- length(xis)

  # Gamma
  listGammaXi <- constructGamma(etas, xis, parTable = parTable)
  gammaXi <- listGammaXi$numeric 
  labelGammaXi <- listGammaXi$label 

  listGammaEta <- constructGamma(etas, etas, parTable = parTable)
  gammaEta <- listGammaEta$numeric 
  labelGammaEta <- listGammaEta$label

  # covariance matrices
  listPsi <- constructPsi(etas, parTable = parTable)
  psi <- listPsi$numeric
  labelPsi <- listPsi$label

  listPhi <- constructPhi(xis, method = method, parTable = parTable)
  phi <- listPhi$numeric
  labelPhi <- listPhi$label

  listA <- constructA(xis, method = method, parTable = parTable)
  A <- listA$numeric
  labelA <- listA$label

  matrices <- list(
    gammaXi = gammaXi,
    gammaEta = gammaEta,
    A = A,
    psi = psi,
    phi = phi) 

  labelMatrices <- list(
    gammaXi = labelGammaXi,
    gammaEta = labelGammaEta,
    A = labelA,
    psi = labelPsi,
    phi = labelPhi)


  model <- list(info =
                list(etas = etas,
                     numEtas = numEtas,
                     xis = xis,
                     numXis = numXis),
                matrices = matrices,
                labelMatrices = labelMatrices,
                syntax = syntax,
                parTable = parTable)
   
  listTheta <- createThetaCovModel(model)
  model$freeParams <- listTheta$freeParams
  model$lenLabelTheta <- listTheta$lenLabelTheta
  model$theta <- listTheta$thetaCov

  model
}


fillCovModel <- function(covModel, theta, fillPhi = FALSE, method = "lms") {
  if (is.null(names(theta))) names(theta) <- names(covModel$theta)
  if (is.null(covModel$matrices)) return(NULL)
  matrices <- covModel$matrices 
  
  if (covModel$lenLabelTheta == 0) {
    thetaLabelCov <- NULL
    thetaCov <- theta
  } else {
    thetaCov <- theta[-seq_len(covModel$lenLabelTheta)]
    thetaLabelCov <- theta[seq_len(covModel$lenLabelTheta)]
    matrices <- fillMatricesLabels(matrices[paramMatricesCov], 
                                   covModel$labelMatrices[paramMatricesCov], 
                                   thetaLabelCov)
  }

  matrices$psi[is.na(matrices$psi)] <-
    fetch(thetaCov, "^psi[0-9]*$")
  matrices$gammaEta[is.na(matrices$gammaEta)] <-
    fetch(thetaCov, "gammaEta[0-9]*$")
  matrices$gammaXi[is.na(matrices$gammaXi)] <-
    fetch(thetaCov, "gammaXi[0-9]*$")
  if (method == "lms") {
    matrices$A[is.na(matrices$A)] <-
      fetch(thetaCov, "^A[0-9]*$")
  } else if (method == "qml") {
    matrices$phi <- fillSymmetric(matrices$phi, fetch(thetaCov, "^phi[0-9]*$"))
  }

  if (fillPhi) matrices$phi <- matrices$A %*% t(matrices$A)
  covModel$matrices <- matrices 
  covModel
}


createThetaCovModel <- function(covModel, start = NULL) {
  matrices <- covModel$matrices
  labelThetaCov <- createLabelTheta(covModel$labelMatrices)

  set.seed(123)
  phi <- as.vector(matrices$phi)
  A <- as.vector(matrices$A)
  psi <- as.vector(matrices$psi)
  alpha <- as.vector(matrices$alpha)
  gammaXi <- as.vector(matrices$gammaXi)
  gammaEta <- as.vector(matrices$gammaEta)
  thetaCov <- c("phi" = phi,
             "A" = A,
             "psi" = psi,
             "gammaXi" = gammaXi,
             "gammaEta" = gammaEta)
  thetaCov <- thetaCov[is.na(thetaCov)]
  if (is.null(start)) {
   thetaCov <- vapply(thetaCov, FUN.VALUE = vector("numeric", 1L),
                   FUN = function(x) stats::runif(1))
  }

  list(thetaCov = c(labelThetaCov, thetaCov),
       freeParams = length(thetaCov) + length(labelThetaCov),
       lenLabelTheta = length(labelThetaCov))
}


countFreeCovModel <- function(matrices) {
  vapply(matrices, FUN.VALUE = integer(1L), 
         FUN = function(x) sum(is.na(x))) |> sum()
}


expectedCovModel <- function(model, method = "lms", sortedXis) {
  gammaXi <- model$matrices$gammaXi
  gammaEta <- model$matrices$gammaEta

  if (method == "lms") {
    A <- model$matrices$A
    phi <- A %*% t(A)
  } else if (method == "qml") {
    phi <- model$matrices$phi
  }
  psi <- model$matrices$psi
  
  Binv <- solve(diag(nrow(gammaEta)) - gammaEta)
  covEtaEta <- Binv %*% (gammaXi %*% phi %*% t(gammaXi) + psi) %*% t(Binv)
  covEtaXi <- Binv %*% gammaXi %*% phi
  sigma <- rbind(cbind(covEtaEta, covEtaXi),
                 cbind(t(covEtaXi), phi))
  sigma <- sigma[sortedXis, sortedXis]

  if (method == "lms") {
    sigma <- tryCatch(t(chol(sigma)), 
                      error = function(e) {
                        sigma[TRUE] <- NaN
                        sigma
                      })
  }
  sigma
}


covModelToParTable <- function(model, method = "lms") {
  matricesEst <- model$covModel$matrices
  matricesSE <- model$covModelSE$matrices
  matricesNA <- model$covModelNA$matrices
  if (is.null(matricesEst) || is.null(matricesNA) ||
      is.null(matricesSE)) return(NULL)
  etas <- model$info$etas
  numXis <- model$info$numXis
  parTable <- NULL

  # coefficients Structural Model 
  newRows <- matrixToParTable(matricesNA$gammaXi,
                              matricesEst$gammaXi,
                              matricesSE$gammaXi, 
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)
  
  newRows <- matrixToParTable(matricesNA$gammaEta,
                              matricesEst$gammaEta,
                              matricesSE$gammaEta, 
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)
  
  if (method == "lms") {
    phiNA <- matricesNA$A
    phiEst <- matricesEst$phi
    phiSE <- matricesSE$A
  } else if (method == "qml") {
    phiNA <- matricesNA$phi
    phiEst <- matricesEst$phi
    phiSE <- matricesSE$phi
  }

  newRows <- matrixToParTable(phiNA,
                              phiEst,
                              phiSE,
                              op = "~~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$psi,
                              matricesEst$psi,
                              matricesSE$psi, 
                              op = "~~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

  # return
  parTable
}
