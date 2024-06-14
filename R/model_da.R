# Functions for Specifiying lms and qml model. modsem(method = c("lms", "qml"))
# Last updated: 09.06.2024


# Global variables 
paramMatrices <- c("lambdaX", "lambdaY", "gammaXi", "gammaEta", 
                   "thetaDelta", "thetaEpsilon", "phi", "A",
                   "psi", "tauX", "tauY", "alpha", "omegaEtaXi", 
                   "omegaXiXi")


# Functions
specifyModelDA <- function(syntax = NULL, 
                               data = NULL, 
                               method = "lms", 
                               m = 16,
                               cov.syntax = NULL, 
                               double = FALSE, 
                               parTable = NULL, 
                               parTableCovModel = NULL,
                               auto.constraints = TRUE, 
                               create.theta = TRUE,
                               mean.observed = TRUE,
                               standardize.inp = FALSE, 
                               standardize.out = FALSE
                               ) {
  if (!is.null(syntax)) parTable <- modsemify(syntax)
  if (is.null(parTable)) stop2("No parTable found")

  # endogenous variables (etas)model
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
  listLambdaX <- constructLambda(xis, indsXis, parTable = parTable,
                                 auto.constraints = auto.constraints)
  lambdaX <- listLambdaX$numeric 
  labelLambdaX <- listLambdaX$label 

  listTauX <- constructTau(xis, indsXis, parTable = parTable,
                           mean.observed = mean.observed)
  tauX <- listTauX$numeric 
  labelTauX <- listTauX$label

  listThetaDelta <- constructTheta(xis, indsXis, parTable = parTable,
                                   auto.constraints = auto.constraints)
  thetaDelta <- listThetaDelta$numeric
  thetaLabelDelta <- listThetaDelta$label

  # measurement model y
  listLambdaY <- constructLambda(etas, indsEtas, parTable = parTable,
                                 auto.constraints = auto.constraints)
  lambdaY <- listLambdaY$numeric 
  labelLambdaY <- listLambdaY$label 

  listTauY <- constructTau(etas, indsEtas, parTable = parTable,
                           mean.observed = mean.observed)
  tauY <- listTauY$numeric 
  labelTauY <- listTauY$label 
  
  listThetaEpsilon <- constructTheta(etas, indsEtas, parTable = parTable,
                                     auto.constraints = auto.constraints)
  thetaEpsilon <- listThetaEpsilon$numeric 
  thetaLabelEpsilon <- listThetaEpsilon$label

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

  listPhi <- constructPhi(xis, method = method, cov.syntax = cov.syntax,
                          parTable = parTable)
  phi <- listPhi$numeric
  labelPhi <- listPhi$label

  listA <- constructA(xis, method = method, cov.syntax = cov.syntax,
                  parTable = parTable)
  A <- listA$numeric
  labelA <- listA$label

  # mean etas
  listAlpha <- constructAlpha(etas, parTable = parTable, 
                              auto.constraints = auto.constraints,
                              mean.observed = mean.observed)
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
  
  covModel <- covModel(cov.syntax, method = method, parTable = parTableCovModel)

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
    thetaDelta = thetaLabelDelta,
    thetaEpsilon = thetaLabelEpsilon,
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
                cov.syntax = cov.syntax,
                parTable = parTable,
                covModel = covModel)


  model$constrExprs <- getConstrExprs(parTable, model$covModel$parTable)
 
  if (create.theta) {
    listTheta <- create.theta(model)
    model$freeParams <- length(listTheta$theta)
    model$lenThetaMain <- listTheta$lenThetaMain
    model$lenThetaCov <- listTheta$lenThetaCov
    model$lenThetaLabel <- listTheta$lenThetaLabel
    model$totalLenThetaLabel <- listTheta$totalLenThetaLabel
    model$theta <- listTheta$theta
    model$info$bounds <- getParamBounds(model)
  }

  model$data <- cleanAndSortData(data, allIndsXis, allIndsEtas)
  model$info$N <- NROW(model$data)

  checkModel(model = model, covModel = covModel) 

  model
}


create.theta <- function(model, start = NULL) {
  set.seed(123)
  thetaCov <- create.thetaCovModel(model$covModel)
  thetaLabel <- create.thetaLabel(model$labelMatrices, 
                                 model$covModel$labelMatrices,
                                 model$constrExprs)
  totalThetaLabel <- calcThetaLabel(thetaLabel, model$constrExprs)

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

  theta <- c(thetaLabel, thetaCov, thetaMain)
  list(theta = theta,
       lenThetaMain = length(thetaMain),
       lenThetaLabel = length(thetaLabel),
       totalLenThetaLabel = length(totalThetaLabel),
       lenThetaCov = length(thetaCov))
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
  xis <- model$info$xis
  numXis <- model$info$numXis
  numEtas <- model$info$numEtas
  matrices <- model$matrices

  matrices[paramMatrices] <- 
    fillMatricesLabels(matrices[paramMatrices], 
                       model$labelMatrices[paramMatrices], 
                       thetaLabel)

  matrices$lambdaX[is.na(matrices$lambdaX)] <-
    fetch(theta, "lambdaX[0-9]*$")
  matrices$lambdaY[is.na(matrices$lambdaY)] <-
    fetch(theta, "lambdaY[0-9]*$")
  matrices$thetaDelta <- 
    fillSymmetric(matrices$thetaDelta, 
                  fetch(theta, "thetaDelta[0-9]*$"))
  matrices$thetaEpsilon <- 
    fillSymmetric(matrices$thetaEpsilon, 
                  fetch(theta, "thetaEpsilon[0-9]*$"))

  if (method == "lms") {
    if (!is.null(model$covModel$matrices)) {
      matrices$A <- expectedCovModel(model$covModel, method = "lms", 
                                     sortedXis = xis)
    } else {
      matrices$A[is.na(matrices$A)] <- fetch(theta, "^A[0-9]*$")
    }
  } else if (method == "qml"){
    if (!is.null(model$covModel$matrices)) {
      matrices$phi <- expectedCovModel(model$covModel, method = "qml",
                                       sortedXis = xis)
    } else {
      matrices$phi <- fillSymmetric(matrices$phi, 
                                    fetch(theta, "^phi[0-9]*$"))
    }
  }

  matrices$psi <- fillSymmetric(matrices$psi, 
                                fetch(theta, "^psi[0-9]*$"))

  matrices$tauX[is.na(matrices$tauX)] <-
    fetch(theta, "tauX[0-9]*$")
  matrices$tauY[is.na(matrices$tauY)] <-
    fetch(theta, "tauY[0-9]*$")
  matrices$alpha[is.na(matrices$alpha)] <-
    fetch(theta, "alpha[0-9]*$")
  matrices$gammaEta[is.na(matrices$gammaEta)] <-
    fetch(theta, "gammaEta[0-9]*$")
  matrices$gammaXi[is.na(matrices$gammaXi)] <-
    fetch(theta, "gammaXi[0-9]*$")
  matrices$omegaXiXi[is.na(matrices$omegaXiXi)] <-
    fetch(theta, "omegaXiXi[0-9]*$")
  matrices$omegaEtaXi[is.na(matrices$omegaEtaXi)] <-
    fetch(theta, "omegaEtaXi[0-9]*$")
  if (fillPhi) matrices$phi <- matrices$A %*% t(matrices$A)
  matrices
}


# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000)
quadrature <- function(m, k) {
  if (k == 0 || m == 0) return(list(n = matrix(0), w = 1, k = 0, m = m))
  singleDimGauss <- gaussquad::hermite.h.quadrature.rules(m)[[m]]
  nodes <- lapply(seq_len(k), function(k) singleDimGauss$x) |>
    expand.grid() |> as.matrix()
  weights <- lapply(seq_len(k), function(k) singleDimGauss$w) |>
    expand.grid() |> apply(MARGIN = 1, prod)
  list(n = nodes * sqrt(2), w = weights * pi ^ (-k/2), k = k, m = m)
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


matrixToParTable <- function(matrixNA, matrixEst, matrixSE, matrixLabel,
                             op = "=~", rowsLhs = TRUE) {
  parTable <- NULL
  if (!rowsLhs) {
    matrixNA <- t(matrixNA)
    matrixEst <- t(matrixEst)
    matrixSE <- t(matrixSE)
    matrixLabel <- t(matrixLabel)
  }
  for (lhs in rownames(matrixNA)) {
    for (rhs in colnames(matrixNA)) {
      if (!is.na(matrixNA[lhs, rhs]) && matrixLabel[lhs, rhs] == "") next
      newRow <- data.frame(lhs = lhs, 
                           op = op, 
                           rhs = rhs,
                           label = matrixLabel[lhs, rhs],
                           est = matrixEst[lhs, rhs],
                           std.error = matrixSE[lhs, rhs])
      parTable <- rbind(parTable, newRow)
    }
  }
  parTable 
}


omegaToParTable <- function(omegaNA, omegaEst, omegaSE, omegaLabel, etas) {
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
        if (!is.na(omegaNA[i, j]) && omegaLabel[i, j] == "") next
        intTerm <- paste0(rowNames[[i]], ":", colNames[[j]])
        newRow <- data.frame(lhs = etas[[eta_i]], 
                             op = "~", 
                             rhs = intTerm,
                             label = omegaLabel[i, j],
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
  matricesLabel <- finalModel$labelMatrices
  etas <- finalModel$info$etas
  numXis <- finalModel$info$numXis
  parTable <- NULL

  # Coefficients Measurement Model 
  newRows <- matrixToParTable(matricesNA$lambdaX,
                              matricesEst$lambdaX,
                              matricesSE$lambdaX, 
                              matricesLabel$lambdaX,
                              op = "=~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$lambdaY,
                              matricesEst$lambdaY,
                              matricesSE$lambdaY, 
                              matricesLabel$lambdaY,
                              op = "=~",
                              rowsLhs = FALSE)
  parTable <- rbind(parTable, newRows)

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
  
  # interaction effects
  newRows <- omegaToParTable(matricesNA$omegaXiXi,
                             matricesEst$omegaXiXi, 
                             matricesSE$omegaXiXi,
                             matricesLabel$omegaXiXi,
                             etas = etas)
  parTable <- rbind(parTable, newRows)

  newRows <- omegaToParTable(matricesNA$omegaEtaXi,
                             matricesEst$omegaEtaXi, 
                             matricesSE$omegaEtaXi,
                             matricesLabel$omegaEtaXi,
                             etas = etas)
  parTable <- rbind(parTable, newRows)

  # Intercepts
  newRows <- matrixToParTable(matricesNA$tauX,
                              matricesEst$tauX,
                              matricesSE$tauX,
                              matricesLabel$tauX,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$tauY,
                              matricesEst$tauY,
                              matricesSE$tauY,
                              matricesLabel$tauY,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$alpha,
                              matricesEst$alpha,
                              matricesSE$alpha,
                              matricesLabel$alpha,
                              op = "~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  # Residual (co) variances Measurement Model 
  newRows <- matrixToParTable(matricesNA$thetaDelta,
                              matricesEst$thetaDelta,
                              matricesSE$thetaDelta, 
                              matricesLabel$thetaDelta,
                              op = "~~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$thetaEpsilon,
                              matricesEst$thetaEpsilon,
                              matricesSE$thetaEpsilon, 
                              matricesLabel$thetaEpsilon,
                              op = "~~",
                              rowsLhs = TRUE)
  parTable <- rbind(parTable, newRows)

  # (Co) variances Structural Model 
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
  parTable
}
