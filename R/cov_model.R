# functions for computing the constrained covariance matrix, 
# based on causal relationships. This can make the lms method more flexible, 
# as you can split the model into a non-linear, and linear part. allowing 
# you to use (normally distributed) endogenous variables as non-normal
# as of now the mean-structure is excluded
covModel <- function(syntax = NULL, method = "lms", parTable = NULL) {
  if (is.null(parTable) && !is.null(syntax)) parTable <- modsemify(syntax)
  if (is.null(parTable)) {
    return(list(matrices = NULL, freeParams = 0, info = NULL,
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

  model
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
  matricesLabel <- model$covModel$labelMatrices
  
  if (is.null(matricesEst) || is.null(matricesNA)) return(NULL)
  if (is.null(matricesSE)) matricesSE <- matricesNA

  etas <- model$info$etas
  numXis <- model$info$numXis
  parTable <- NULL

  # coefficients Structural Model 
  newRows <- matrixToParTable(matricesNA$gammaXi,
                              matricesEst$gammaXi,
                              matricesSE$gammaXi, 
                              matricesLabel$gammaXi,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)
  
  newRows <- matrixToParTable(matricesNA$gammaEta,
                              matricesEst$gammaEta,
                              matricesSE$gammaEta, 
                              matricesLabel$gammaEta,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)
  
  if (method == "lms") {
    phiNA <- matricesNA$A
    phiEst <- matricesEst$phi
    phiSE <- matricesSE$A
    phiLabel <- matricesLabel$A
  } else if (method == "qml") {
    phiNA <- matricesNA$phi
    phiEst <- matricesEst$phi
    phiSE <- matricesSE$phi
    phiLabel <- matricesLabel$phi
  }

  newRows <- matrixToParTable(phiNA,
                              phiEst,
                              phiSE,
                              phiLabel,
                              op = "~~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$psi,
                              matricesEst$psi,
                              matricesSE$psi, 
                              matricesLabel$psi,
                              op = "~~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

  parTable <- lapplyDf(parTable, FUN = function(x) replace(x, x == -999, NA))
  # return
  parTable
}
