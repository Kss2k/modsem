# Functions for Specifiying lms and qml model. modsem(method = c("lms", "qml"))
# Last updated: 06.06.2024


# Global variables 
paramMatrices <- c("lambdaX", "lambdaY", "gammaXi", "gammaEta", 
                   "thetaDelta", "thetaEpsilon", "phi", "A",
                   "psi", "tauX", "tauY", "alpha", "omegaEtaXi", 
                   "omegaXiXi")


# Functions
specifyModelLmsQml <- function(syntax = NULL, data = NULL, method = "lms", m = 16,
                               cov_syntax = NULL, double = FALSE, 
                               parTable = NULL, parTableCovModel = NULL) {
  if (!is.null(syntax)) parTable <- modsem::modsemify(syntax)
  if (is.null(parTable)) stop("No parTable found")

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
  listLambdaX <- constructLambda(xis, indsXis, parTable = parTable)
  lambdaX <- listLambdaX$numeric 
  labelLambdaX <- listLambdaX$label 

  listTauX <- constructTau(xis, indsXis, parTable = parTable)
  tauX <- listTauX$numeric 
  labelTauX <- listTauX$label

  listThetaDelta <- constructTheta(xis, indsXis, parTable = parTable)
  thetaDelta <- listThetaDelta$numeric
  labelThetaDelta <- listThetaDelta$label

  # measurement model y
  listLambdaY <- constructLambda(etas, indsEtas, parTable = parTable)
  lambdaY <- listLambdaY$numeric 
  labelLambdaY <- listLambdaY$label 

  listTauY <- constructTau(etas, indsEtas, parTable = parTable)
  tauY <- listTauY$numeric 
  labelTauY <- listTauY$label 
  
  listThetaEpsilon <- constructTheta(etas, indsEtas, parTable = parTable)
  thetaEpsilon <- listThetaEpsilon$numeric 
  labelThetaEpsilon <- listThetaEpsilon$label

  # structural model 
  Ieta <- diag(numEtas) # used for (B^-1 = (Ieta - gammaEta)^-1)
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

  listPhi <- constructPhi(xis, method = method, cov_syntax = cov_syntax,
                          parTable = parTable)
  phi <- listPhi$numeric
  labelPhi <- listPhi$label

  listA <- constructA(xis, method = method, cov_syntax = cov_syntax,
                  parTable = parTable)
  A <- listA$numeric
  labelA <- listA$label

  # mean etas
  listAlpha <- constructAlpha(etas, parTable = parTable)
  alpha <- listAlpha$numeric
  labelAlpha <- listAlpha$label

  # quadratic terms 
  listOmegaEtaXi <- omegaAndSortedXis$omegaEtaXi
  omegaEtaXi <- listOmegaEtaXi$numeric
  labelOmegaEtaXi <- listOmegaEtaXi$label

  listOmegaXiXi <- omegaAndSortedXis$omegaXiXi
  omegaXiXi <- listOmegaXiXi$numeric
  labelOmegaXiXi <- listOmegaXiXi$label

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

  labelMatrices <- list(
    lambdaX = labelLambdaX,
    lambdaY = labelLambdaY,
    gammaXi = labelGammaXi,
    gammaEta = labelGammaEta,
    thetaDelta = labelThetaDelta,
    thetaEpsilon = labelThetaEpsilon,
    phi = labelPhi,
    A = labelA,
    psi = labelPsi,
    tauX = labelTauX,
    tauY = labelTauY,
    alpha = labelAlpha,
    omegaEtaXi = labelOmegaEtaXi,
    omegaXiXi = labelOmegaXiXi)

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
                labelMatrices = labelMatrices,
                syntax = syntax,
                cov_syntax = cov_syntax,
                parTable = parTable,
                covModel = covModel(cov_syntax, method = method,
                                    parTable = parTableCovModel))

  listTheta <- createTheta(model)
  model$freeParams = length(listTheta$theta)
  model$theta <- listTheta$theta
  model$lenLabelMain <- listTheta$lenLabelMain

  model$info$bounds <- getParamBounds(model)
  
  model$data <- cleanAndSortData(data, allIndsXis, allIndsEtas)
  model$info$N <- NROW(model$data)

  model
}


createTheta <- function(model, start = NULL) {
  thetaCov <- model$covModel$theta
  labelThetaMain <- createLabelTheta(model$labelMatrices)

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
  omegaXiXi <- as.vector(matrices$omegaXiXi)
  omegaEtaXi <- as.vector(matrices$omegaEtaXi)
  thetaMain <- c("lambdaX" = lambdaX,
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
                      "omegaXiXi" = omegaXiXi,
                      "omegaEtaXi" = omegaEtaXi)
  thetaMain <- thetaMain[is.na(thetaMain)]
  if (is.null(start)) {
   thetaMain <- vapply(thetaMain, FUN.VALUE = vector("numeric", 1L),
                            FUN = function(x) stats::runif(1))
  }

  thetaMain <- c(labelThetaMain, thetaMain)
  theta <- c(thetaCov, thetaMain)
  list(theta = theta,
       lenLabelMain = length(labelThetaMain))
}


fillModel <- function(model, theta, fillPhi = FALSE, method = "lms") {
  if (is.null(names(theta))) names(theta) <- names(model$theta)
  # cov model
  if (model$covModel$freeParams == 0) {
    thetaCov <- NULL
    thetaMain <- theta
  } else {
    thetaCov <- theta[seq_len(model$covModel$freeParams)]
    thetaMain <- theta[-seq_len(model$covModel$freeParams)]
  }

  model$covModel <- fillCovModel(model$covModel, thetaCov, 
                                 fillPhi = fillPhi, method = method)
  model$matrices <- fillMainModel(model, thetaMain, fillPhi = fillPhi,
                                  method = method)
  model
}


fillMainModel <- function(model, theta, fillPhi = FALSE, 
                          method = "lms") {
  xis <- model$info$xis
  numXis <- model$info$numXis
  numEtas <- model$info$numEtas
  matrices <- model$matrices
  if (model$lenLabelMain == 0) {
    thetaLabelMain <- NULL
    thetaMain <- theta
  } else {
    thetaMain <- theta[-seq_len(model$lenLabelMain)]
    thetaLabelMain <- theta[seq_len(model$lenLabelMain)]
    matrices[paramMatrices] <- 
      fillMatricesLabels(matrices[paramMatrices], 
                         model$labelMatrices[paramMatrices], 
                         thetaLabelMain)
  }

  matrices$lambdaX[is.na(matrices$lambdaX)] <-
    fetch(thetaMain, "lambdaX[0-9]*$")
  matrices$lambdaY[is.na(matrices$lambdaY)] <-
    fetch(thetaMain, "lambdaY[0-9]*$")
  matrices$thetaDelta <- 
    fillSymmetric(matrices$thetaDelta, 
                  fetch(thetaMain, "thetaDelta[0-9]*$"))
  matrices$thetaEpsilon <- 
    fillSymmetric(matrices$thetaEpsilon, 
                  fetch(thetaMain, "thetaEpsilon[0-9]*$"))

  if (method == "lms") {
    if (!is.null(model$covModel$matrices)) {
      matrices$A <- expectedCovModel(model$covModel, method = "lms", 
                                     sortedXis = xis)
    } else {
      matrices$A[is.na(matrices$A)] <- fetch(thetaMain, "^A[0-9]*$")
    }
  } else if (method == "qml"){
    if (!is.null(model$covModel$matrices)) {
      matrices$phi <- expectedCovModel(model$covModel, method = "qml",
                                       sortedXis = xis)
    } else {
      matrices$phi <- fillSymmetric(matrices$phi, 
                                    fetch(thetaMain, "^phi[0-9]*$"))
    }
  }

  matrices$psi <- fillSymmetric(matrices$psi, 
                                fetch(thetaMain, "^psi[0-9]*$"))

  matrices$tauX[is.na(matrices$tauX)] <-
    fetch(thetaMain, "tauX[0-9]*$")
  matrices$tauY[is.na(matrices$tauY)] <-
    fetch(thetaMain, "tauY[0-9]*$")
  matrices$alpha[is.na(matrices$alpha)] <-
    fetch(thetaMain, "alpha[0-9]*$")
  matrices$gammaEta[is.na(matrices$gammaEta)] <-
    fetch(thetaMain, "gammaEta[0-9]*$")
  matrices$gammaXi[is.na(matrices$gammaXi)] <-
    fetch(thetaMain, "gammaXi[0-9]*$")
  matrices$omegaXiXi[is.na(matrices$omegaXiXi)] <-
    fetch(thetaMain, "omegaXiXi[0-9]*$")
  matrices$omegaEtaXi[is.na(matrices$omegaEtaXi)] <-
    fetch(thetaMain, "omegaEtaXi[0-9]*$")
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
  lower <- rep(-Inf, model$freeParams)
  upper <- rep(Inf, model$freeParams)
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
                             std.error = matrixSE[lhs, rhs])
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
                             std.error = omegaSE[i, j]) 
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

  # Intercepts
  newRows <- matrixToParTable(matricesNA$tauX,
                              matricesEst$tauX,
                              matricesSE$tauX,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$tauY,
                              matricesEst$tauY,
                              matricesSE$tauY,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$alpha,
                              matricesEst$alpha,
                              matricesSE$alpha,
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

  # remove -999
  parTable <- lapplyDf(parTable, FUN = function(x) replace(x, x == -999, NA))
  # return
  parTable
}
