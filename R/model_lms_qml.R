# Functions for Specifiying lms and qml model. modsem(method = c("lms", "qml"))
# Last updated: 06.06.2024


specifyModelLmsQml <- function(syntax, data = NULL, method = "lms", m = 16,
                               cov_syntax = NULL, double = FALSE) {
  parTable <- modsem::modsemify(syntax)
  # endogenous variables (etas)
  etas <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  numEtas <- length(etas)
  
  indsEtas <- getIndsLVs(parTable, etas)
  numIndsEtas <- vapply(indsEtas, FUN.VALUE = vector("integer", 1L),
                        FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  
  # exogenouts variables (xis) and interaction terms 
  intTerms <- getIntTerms(parTable)
  varsInts <- getVarsInts(intTerms)
  allVarsInInts <- unique(unlist(varsInts))
  xis <- getXis(parTable, checkAny = TRUE)
  numXis <- length(xis)

  omegaAndSortedXis <- sortXisConstructOmega(xis, varsInts, etas, intTerms,
                                             method = method, double = double)
  xis <- omegaAndSortedXis$sortedXis # get sorted xis according to interaction terms

  indsXis <- getIndsLVs(parTable, xis)
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis <- unlist(indsXis)
  numAllIndsXis <- length(allIndsXis)

  # measurement model x
  lambdaX <- constructLambda(xis, indsXis)
  tauX <- constructTau(xis, indsXis)
  thetaDelta <- constructTheta(xis, indsXis)

  # measurement model y
  lambdaY <- constructLambda(etas, indsEtas)
  tauY <- constructTau(etas, indsEtas)
  thetaEpsilon <- constructTheta(etas, indsEtas)

  # structural model 
  Ieta <- diag(numEtas) # used for (B^-1 = (Ieta - gammaEta)^-1)
  gammaXi <- constructGamma(etas, xis, parTable = parTable)
  gammaEta <- constructGamma(etas, etas, parTable = parTable)

  # covariance matrices
  psi <- constructPsi(etas)
  phi <- constructPhi(xis, method = method, cov_syntax = cov_syntax)
  A <- constructA(xis, method = method, cov_syntax = cov_syntax)

  # mean etas
  alpha <- constructAlpha(etas)

  # quadratic terms 
  omegaEtaXi <- omegaAndSortedXis$omegaEtaXi
  omegaXiXi <- omegaAndSortedXis$omegaXiXi

  # matrices for scaling variables in qml
  selectScalingY <- selectScalingY(lambdaY, method = method)
  selectBetaRows <- selectBetaRows(lambdaY, method = method)
  emptyR <- constructR(etas, indsEtas, lambdaY, method = method)
  scalingInds <- getScalingInds(allIndsEtas, R = emptyR, method = method)
  selectThetaEpsilon <- selectThetaEpsilon(lambdaY, thetaEpsilon, 
                                           scalingInds, method = method)
  subThetaEpsilon <- constructSubThetaEpsilon(indsEtas, thetaEpsilon, 
                                              scalingInds, method = method)

  # list of matrices
  matrices <- list(
    lambdaX = lambdaX,
    lambdaY = lambdaY,
    gammaXi = gammaXi,
    gammaEta = gammaEta,
    thetaDelta = thetaDelta,
    thetaEpsilon = thetaEpsilon,
    phi = phi,
    A = A,
    Ieta = Ieta,
    psi = psi,
    tauX = tauX,
    tauY = tauY,
    alpha = alpha,
    omegaEtaXi = omegaEtaXi,
    omegaXiXi = omegaXiXi,
    selectScalingY = selectScalingY,
    selectThetaEpsilon = selectThetaEpsilon,
    selectBetaRows = selectBetaRows,
    emptyR = emptyR,
    subThetaEpsilon = subThetaEpsilon)

  k <- omegaAndSortedXis$k
  quad <- quadrature(m, k)

  model <- list(info =
                list(etas = etas,
                     numEtas = numEtas,
                     indsEtas = indsEtas,
                     allIndsEtas = allIndsEtas,
                     xis = xis,
                     varsInts = varsInts,
                     numXis = numXis,
                     indsXis = indsXis,
                     allIndsXis = allIndsXis,
                     scalingInds = scalingInds,
                     kOmegaEta = getK_NA(omegaEtaXi)),
                quad = quad,
                matrices = matrices,
                syntax = syntax,
                parTable = parTable,
                covModel = covModel(cov_syntax, method = method))

  model$theta <- createParamVector(model)
  model$info$bounds <- getParamBounds(model)
  
  model$data <- cleanAndSortData(data, allIndsXis, allIndsEtas)
  model$info$N <- NROW(model$data)

  model
}


createParamVector <- function(model, start = NULL) {
  thetaCovModel <- model$covModel$theta
  set.seed(123)
  matrices <- model$matrices
  lambdaX <- as.vector(matrices$lambdaX)
  lambdaY <- as.vector(matrices$lambdaY)
  thetaDelta <- as.vector(matrices$thetaDelta)
  thetaEpsilon <- as.vector(matrices$thetaEpsilon)
  phi <- as.vector(matrices$phi)
  A <- as.vector(matrices$A)
  psi <- as.vector(matrices$psi)
  tauX <- as.vector(matrices$tauX)
  tauY <- as.vector(matrices$tauY)
  alpha <- as.vector(matrices$alpha)
  gammaXi <- as.vector(matrices$gammaXi)
  gammaEta <- as.vector(matrices$gammaEta)
  omgeaXiXi <- as.vector(matrices$omegaXiXi)
  omegaEtaXi <- as.vector(matrices$omegaEtaXi)
  thetaMainModel <- c("lambdaX" = lambdaX,
                      "lambdaY" = lambdaY,
                      "tauX" = tauX,
                      "tauY" = tauY,
                      "thetaDelta" = thetaDelta,
                      "thetaEpsilon" = thetaEpsilon,
                      "phi" = phi,
                      "A" = A,
                      "psi" = psi,
                      "alpha" = alpha,
                      "gammaXi" = gammaXi,
                      "gammaEta" = gammaEta,
                      "omegaXiXi" = omgeaXiXi,
                      "omegaEtaXi" = omegaEtaXi)
  thetaMainModel <- thetaMainModel[is.na(thetaMainModel)]
  if (is.null(start)) {
   thetaMainModel <- vapply(thetaMainModel, FUN.VALUE = vector("numeric", 1L),
                            FUN = function(x) stats::runif(1))
  }
  c(thetaCovModel, thetaMainModel)
}


fillModel <- function(model, theta, fillPhi = FALSE, method = "lms") {
  if (is.null(names(theta))) names(theta) <- names(model$theta)
  # cov model
  if (model$covModel$freeParams == 0) {
    thetaCovModel <- NULL
    thetaMainModel <- theta
  } else {
    thetaCovModel <- theta[seq_len(model$covModel$freeParams)]
    thetaMainModel <- theta[-seq_len(model$covModel$freeParams)]
  }

  model$covModel <- fillCovModel(model$covModel, thetaCovModel, 
                                 fillPhi = fillPhi, method = method)
  model$matrices <- fillMainModel(model, thetaMainModel, fillPhi = fillPhi,
                                  method = method)
  model
}


fillMainModel <- function(model, thetaMainModel, fillPhi = FALSE, 
                          method = "lms") {
  xis <- model$info$xis
  numXis <- model$info$numXis
  numEtas <- model$info$numEtas
  matrices <- model$matrices
  matrices$lambdaX[is.na(matrices$lambdaX)] <-
    fetch(thetaMainModel, "lambdaX[0-9]*$")
  matrices$lambdaY[is.na(matrices$lambdaY)] <-
    fetch(thetaMainModel, "lambdaY[0-9]*$")
  matrices$thetaDelta[is.na(matrices$thetaDelta)] <-
    fetch(thetaMainModel, "thetaDelta[0-9]*")
  matrices$thetaEpsilon[is.na(matrices$thetaEpsilon)] <-
    fetch(thetaMainModel, "thetaEpsilon[0-9]*")

  if (method == "lms") {
    if (!is.null(model$covModel$matrices)) {
      matrices$A <- expectedCovModel(model$covModel, method = "lms", 
                                     sortedXis = xis)
    } else {
      matrices$A[is.na(matrices$A)] <- fetch(thetaMainModel, "^A[0-9]*$")
    }
  } else if (method == "qml"){
    if (!is.null(model$covModel$matrices)) {
      matrices$phi <- expectedCovModel(model$covModel, method = "qml",
                                       sortedXis = xis)
    } else {
      matrices$phi <- fillSymmetric(matrices$phi, fetch(thetaMainModel, 
                                                        "^phi[0-9]*$"))
    }
  }

  matrices$psi[is.na(matrices$psi)] <-
    fetch(thetaMainModel, "^psi[0-9]*$")
  matrices$tauX[is.na(matrices$tauX)] <-
    fetch(thetaMainModel, "tauX[0-9]*$")
  matrices$tauY[is.na(matrices$tauY)] <-
    fetch(thetaMainModel, "tauY[0-9]*$")
  matrices$alpha[is.na(matrices$alpha)] <-
    fetch(thetaMainModel, "alpha[0-9]*$")
  matrices$gammaEta[is.na(matrices$gammaEta)] <-
    fetch(thetaMainModel, "gammaEta[0-9]*$")
  matrices$gammaXi[is.na(matrices$gammaXi)] <-
    fetch(thetaMainModel, "gammaXi[0-9]*$")
  matrices$omegaXiXi[is.na(matrices$omegaXiXi)] <-
    fetch(thetaMainModel, "omegaXiXi[0-9]*$")
  matrices$omegaEtaXi[is.na(matrices$omegaEtaXi)] <-
    fetch(thetaMainModel, "omegaEtaXi[0-9]*$")
  if (fillPhi) matrices$phi <- matrices$A %*% t(matrices$A)
  matrices
}


# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000)
quadrature <- function(m, k) {
  if (k == 0) return(list(n = matrix(0), w = 1, k = 0, m = m))
  singleDimGauss <- gaussquad::hermite.h.quadrature.rules(m)[[m]]
  nodes <- lapply(seq_len(k), function(k) singleDimGauss$x) |>
    expand.grid() |> as.matrix()
  weights <- lapply(seq_len(k), function(k) singleDimGauss$w) |>
    expand.grid() |> apply(MARGIN = 1, prod)
  list(n = nodes * sqrt(2), w = weights * pi ^ (-k/2), k = k, m = m)
}


# Set bounds for parameters to (0, Inf)
getParamBounds <- function(model) {
  namePattern <- paste0("lambdaX[0-9]*$|lambdaY[0-9]*$|",
                        "thetaDelta[0-9]*$|thetaEpsilon[0-9]*$|",
                        "phi[0-9]*$|psi[0-9]*$")
  lower <- rep(-Inf, countFreeParams(model))
  upper <- rep(Inf, countFreeParams(model))
  names(lower) <- names(upper) <- names(model$theta)
  lower[grepl(namePattern, names(lower))] <- 0
  list(lower = lower, upper = upper)
}


matrixToParTable <- function(matrixNA, matrixEst, matrixSE, 
                             op = "=~", rowsLhs = TRUE) {
  parTable <- NULL
  if (!rowsLhs) {
    matrixNA <- t(matrixNA)
    matrixEst <- t(matrixEst)
    matrixSE <- t(matrixSE)
  }
  for (lhs in rownames(matrixNA)) {
    for (rhs in colnames(matrixNA)) {
      if (is.na(matrixNA[lhs, rhs])) {
        newRow <- data.frame(lhs = lhs, op = op, 
                             rhs = rhs,
                             est = matrixEst[lhs, rhs],
                             se = matrixSE[lhs, rhs])
        parTable <- rbind(parTable, newRow)
      }
    }
  }
  parTable 
}


omegaToParTable <- function(omegaNA, omegaEst, omegaSE, etas) {
  numEtas <- length(etas) 
  subNrow <- nrow(omegaNA) %/% numEtas
  subSeqRows <- seq_len(subNrow) 
  lastRow <- 0 
  lastCol <- 0
  rowNames <- rownames(omegaNA)
  colNames <- colnames(omegaNA)
  parTable <- NULL
  for (eta_i in seq_len(numEtas)) {
    for (i in subSeqRows + (eta_i - 1) * subNrow) {
      for (j in seq_len(ncol(omegaNA))) {
        if (!is.na(omegaNA[i, j])) next
        intTerm <- paste0(rowNames[[i]], ":", colNames[[j]])
        newRow <- data.frame(lhs = etas[[eta_i]], 
                             op = "~", rhs = intTerm,
                             est = omegaEst[i, j],
                             se = omegaSE[i, j]) 
        parTable <- rbind(parTable, newRow)
      }
    }
  }
  parTable
}


finalModelToParTable <- function(finalModel, method = "lms") {
  matricesEst <- finalModel$matrices
  matricesSE <- finalModel$matricesSE
  matricesNA <- finalModel$matricesNA
  etas <- finalModel$info$etas
  numXis <- finalModel$info$numXis
  parTable <- NULL

  # Coefficients Measurement Model 
  matricesNA$lambdaX[matricesEst$lambdaX == 1] <- NA
  matricesSE$lambdaX[matricesEst$lambdaX == 1] <- NA
  newRows <- matrixToParTable(matricesNA$lambdaX,
                              matricesEst$lambdaX,
                              matricesSE$lambdaX, 
                              op = "=~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

  matricesNA$lambdaY[matricesEst$lambdaY == 1] <- NA
  matricesSE$lambdaY[matricesEst$lambdaY == 1] <- NA
  newRows <- matrixToParTable(matricesNA$lambdaY,
                              matricesEst$lambdaY,
                              matricesSE$lambdaY, 
                              op = "=~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

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
  
  # interaction effects
  newRows <- omegaToParTable(matricesNA$omegaXiXi,
                             matricesEst$omegaXiXi, 
                             matricesSE$omegaXiXi,
                             etas = etas)
  parTable <- rbind(parTable, newRows)

  newRows <- omegaToParTable(matricesNA$omegaEtaXi,
                             matricesEst$omegaEtaXi, 
                             matricesSE$omegaEtaXi,
                             etas = etas)
  parTable <- rbind(parTable, newRows)

  # Means Measurement Model 
  tauXNA <- matricesNA$tauX
  tauXEst <- matricesEst$tauX
  tauXSE <- matricesSE$tauX
  colnames(tauXNA) <- colnames(tauXEst) <- colnames(tauXSE) <- "1"
  newRows <- matrixToParTable(tauXNA,
                              tauXEst,
                              tauXSE,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  tauYNA <- matricesNA$tauY
  tauYEst <- matricesEst$tauY
  tauYSE <- matricesSE$tauY
  colnames(tauYNA) <- colnames(tauYEst) <- colnames(tauYSE) <- "1"
  newRows <- matrixToParTable(tauYNA,
                              tauYEst,
                              tauYSE,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  # Residual (co) variances Measurement Model 
  newRows <- matrixToParTable(matricesNA$thetaDelta,
                              matricesEst$thetaDelta,
                              matricesSE$thetaDelta, 
                              op = "~~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$thetaEpsilon,
                              matricesEst$thetaEpsilon,
                              matricesSE$thetaEpsilon, 
                              op = "~~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  # (Co) variances Structural Model 
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
