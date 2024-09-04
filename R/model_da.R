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
                           createTheta = TRUE,
                           mean.observed = TRUE,
                           standardize.inp = FALSE, 
                           standardize.out = FALSE,
                           checkModel = TRUE,
                           quad.range = Inf) {
  if (!is.null(syntax)) parTable <- modsemify(syntax)
  if (is.null(parTable)) stop2("No parTable found")

  # additions to lavaan-syntax for optimizer
  lavOptimizerSyntaxAdditions <- ""

  # endogenous variables (etas)model
  etas <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  numEtas <- length(etas)
  
  indsEtas <- getIndsLVs(parTable, etas)
  numIndsEtas <- vapply(indsEtas, FUN.VALUE = vector("integer", 1L),
                        FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  
  # exogenouts variables (xis) and interaction terms 
  intTerms <- getIntTermRows(parTable)
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
  lavOptimizerSyntaxAdditions <- paste0(lavOptimizerSyntaxAdditions, 
                                        listTauX$syntaxAdditions)

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
  lavOptimizerSyntaxAdditions <- paste0(lavOptimizerSyntaxAdditions, 
                                        listTauY$syntaxAdditions)
  
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

  # mean xis 
  listBeta0 <- constructAlpha(xis, parTable = parTable, 
                              auto.constraints = auto.constraints,
                              mean.observed = mean.observed)
  beta0 <- listBeta0$numeric
  labelBeta0 <- listBeta0$label

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
  fullR <- constructFullR(etas, indsEtas, lambdaY, method = method)

  latentEtas <- getLatentEtasQml(indsEtas, method = method)
  colsU <- getColsU(etas, indsEtas, lambdaY, method = method)

  fullL2 <- constructFullL2(colsU, etas = etas, method = method)
  selectSubL2 <- getSelectSubL2(fullL2, colsU = colsU, latentEtas = latentEtas, 
                                method = method)
  fullSigma2ThetaEpsilon <- constructFullSigma2ThetaEpsilon(psi, method = method)
  selectSubSigma2ThetaEpsilon <- 
    getSelectSubSigma2ThetaEpsilon(fullSigma2ThetaEpsilon, latentEtas = latentEtas, 
                                   method = method)
  fullU <- constructFullU(fullL2 = fullL2, N = NROW(data), etas = etas, method = method)

  scalingInds <- getScalingInds(indsEtas, R = emptyR, latentEtas = latentEtas, method = method)
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
    beta0 = beta0,
    omegaEtaXi = omegaEtaXi,
    omegaXiXi = omegaXiXi,
    selectScalingY = selectScalingY,
    selectThetaEpsilon = selectThetaEpsilon,
    selectBetaRows = selectBetaRows,
    emptyR = emptyR,
    fullR = fullR,
    fullSigma2ThetaEpsilon = fullSigma2ThetaEpsilon, 
    selectSubSigma2ThetaEpsilon = selectSubSigma2ThetaEpsilon,
    fullL2 = fullL2,
    selectSubL2 = selectSubL2,
    fullU = fullU,
    colsU = colsU,
    colsR = colnames(emptyR),
    rowsR = rownames(emptyR),
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
    beta0 = labelBeta0,
    omegaEtaXi = labelOmegaEtaXi,
    omegaXiXi = labelOmegaXiXi)

  k <- omegaAndSortedXis$k
  quad <- quadrature(m, k, cut = quad.range)

  model <- list(info =
                list(etas = etas,
                     latentEtas = latentEtas,
                     numEtas = numEtas,
                     indsEtas = indsEtas,
                     allIndsEtas = allIndsEtas,
                     xis = xis,
                     varsInts = varsInts,
                     numXis = numXis,
                     indsXis = indsXis,
                     allIndsXis = allIndsXis,
                     scalingInds = scalingInds,
                     kOmegaEta = getK_NA(omegaEtaXi),
                     lavOptimizerSyntaxAdditions = lavOptimizerSyntaxAdditions),
                quad = quad,
                matrices = matrices,
                labelMatrices = labelMatrices,
                syntax = syntax,
                cov.syntax = cov.syntax,
                parTable = parTable,
                covModel = covModel)


  model$constrExprs <- getConstrExprs(parTable, model$covModel$parTable)
 
  if (createTheta) {
    listTheta         <- createTheta(model)
    model             <- c(model, listTheta)
    model$freeParams  <- length(listTheta$theta)
    model$info$bounds <- getParamBounds(model)
  }

  model$data <- cleanAndSortData(data, allIndsXis, allIndsEtas)
  model$info$N <- NROW(model$data)

  if (checkModel) checkModel(model = model, covModel = covModel, method = method)

  model
}


matrixToParTable <- function(matrixNA, matrixEst, matrixSE, matrixLabel,
                             op = "=~", rowsLhs = TRUE) {
  if (!rowsLhs) {
    matrixNA <- t(matrixNA)
    matrixEst <- t(matrixEst)
    matrixSE <- t(matrixSE)
    matrixLabel <- t(matrixLabel)
  }
  
  parTable <- NULL
  for (lhs in rownames(matrixEst)) {
    for (rhs in colnames(matrixEst)) {
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
  rowNames <- rownames(omegaEst)
  colNames <- colnames(omegaEst)

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


mainModelToParTable <- function(finalModel, method = "lms") {
  matricesEst <- finalModel$matrices
  matricesSE <- finalModel$matricesSE
  matricesNA <- finalModel$matricesNA
  matricesLabel <- finalModel$labelMatrices
  
  if (is.null(matricesSE)) matricesSE <- matricesNA
  
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
  
  newRows <- matrixToParTable(matricesNA$beta0,
                              matricesEst$beta0,
                              matricesSE$beta0,
                              matricesLabel$beta0,
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



modelToParTable <- function(model, method = "lms") {
  rbind(mainModelToParTable(model, method = method),
        covModelToParTable(model, method = method))
}
