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
                           quad.range = Inf,
                           adaptive.quad = FALSE,
                           adaptive.frequency = 3,
                           impute.na = FALSE,
                           orthogonal.y = FALSE) {
  if (!is.null(syntax)) parTable <- modsemify(syntax)
  stopif(is.null(parTable), "No parTable found")
  
  checkParTableDA(parTable)
  # additions to lavaan-syntax for optimizer
  lavOptimizerSyntaxAdditions <- ""

  # endogenous variables (etas)model
  etas    <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  numEtas <- length(etas)

  indsEtas    <- getIndsLVs(parTable, etas)
  numIndsEtas <- vapply(indsEtas, FUN.VALUE = vector("integer", 1L),
                        FUN = length)
  allIndsEtas    <- unique(unlist(indsEtas))
  numAllIndsEtas <- length(allIndsEtas)

  # exogenouts variables (xis) and interaction terms
  intTerms      <- getIntTermRows(parTable)
  varsInts      <- getVarsInts(intTerms)
  allVarsInInts <- unique(unlist(varsInts))
  xis           <- getXis(parTable, checkAny = TRUE)
  numXis        <- length(xis)
  
  omegaAndSortedXis <- sortXisConstructOmega(xis, varsInts, etas, intTerms,
                                             method = method, double = double)
  xis <- omegaAndSortedXis$sortedXis # get sorted xis according to interaction terms
  nonLinearXis <- omegaAndSortedXis$nonLinearXis

  indsXis    <- getIndsLVs(parTable, xis)
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis    <- unique(unlist(indsXis))
  numAllIndsXis <- length(allIndsXis)

  # clean data
  data <- cleanAndSortData(data, allIndsXis, allIndsEtas, impute.na = impute.na)

  # measurement model x
  listLambdaX <- constructLambda(xis, indsXis, parTable = parTable,
                                 auto.constraints = auto.constraints)
  lambdaX      <- listLambdaX$numeric
  labelLambdaX <- listLambdaX$label

  listTauX <- constructTau(xis, indsXis, parTable = parTable,
                           mean.observed = mean.observed)
  tauX      <- listTauX$numeric
  labelTauX <- listTauX$label
  lavOptimizerSyntaxAdditions <- paste0(lavOptimizerSyntaxAdditions,
                                        listTauX$syntaxAdditions)

  listThetaDelta <- constructTheta(xis, indsXis, parTable = parTable,
                                   auto.constraints = auto.constraints)
  thetaDelta      <- listThetaDelta$numeric
  thetaLabelDelta <- listThetaDelta$label

  # measurement model y
  listLambdaY <- constructLambda(etas, indsEtas, parTable = parTable,
                                 auto.constraints = auto.constraints)
  lambdaY      <- listLambdaY$numeric
  labelLambdaY <- listLambdaY$label

  listTauY <- constructTau(etas, indsEtas, parTable = parTable,
                           mean.observed = mean.observed)
  tauY      <- listTauY$numeric
  labelTauY <- listTauY$label
  lavOptimizerSyntaxAdditions <- paste0(lavOptimizerSyntaxAdditions,
                                        listTauY$syntaxAdditions)

  listThetaEpsilon <- constructTheta(etas, indsEtas, parTable = parTable,
                                     auto.constraints = auto.constraints)
  thetaEpsilon      <- listThetaEpsilon$numeric
  thetaLabelEpsilon <- listThetaEpsilon$label

  # structural model
  Ieta         <- diag(numEtas) # used for (B^-1 = (Ieta - gammaEta)^-1)
  listGammaXi  <- constructGamma(etas, xis, parTable = parTable)
  gammaXi      <- listGammaXi$numeric
  labelGammaXi <- listGammaXi$label

  listGammaEta  <- constructGamma(etas, etas, parTable = parTable)
  gammaEta      <- listGammaEta$numeric
  labelGammaEta <- listGammaEta$label

  # covariance matrices
  listPsi  <- constructPsi(etas, parTable = parTable, orthogonal.y = orthogonal.y)
  psi      <- listPsi$numeric
  labelPsi <- listPsi$label

  listPhi <- constructPhi(xis, method = method, cov.syntax = cov.syntax,
                          parTable = parTable)
  phi      <- listPhi$numeric
  labelPhi <- listPhi$label

  listA <- constructA(xis, method = method, cov.syntax = cov.syntax,
                  parTable = parTable)
  A      <- listA$numeric
  labelA <- listA$label

  # mean etas
  listAlpha <- constructAlpha(etas, parTable = parTable,
                              mean.observed = mean.observed) 
  alpha      <- listAlpha$numeric                            
  labelAlpha <- listAlpha$label

  # mean xis
  listBeta0 <- constructAlpha(xis, parTable = parTable,
                              mean.observed = mean.observed)
  beta0      <- listBeta0$numeric                   
  labelBeta0 <- listBeta0$label

  # quadratic terms
  listOmegaEtaXi  <- omegaAndSortedXis$omegaEtaXi
  omegaEtaXi      <- listOmegaEtaXi$numeric
  labelOmegaEtaXi <- listOmegaEtaXi$label

  listOmegaXiXi  <- omegaAndSortedXis$omegaXiXi
  omegaXiXi      <- listOmegaXiXi$numeric
  labelOmegaXiXi <- listOmegaXiXi$label

  # matrices for scaling variables in qml
  selectScalingY <- selectScalingY(lambdaY, method = method)
  selectBetaRows <- selectBetaRows(lambdaY, method = method)
  emptyR <- constructR(etas, indsEtas, lambdaY, method = method)
  fullR  <- constructFullR(etas, indsEtas, lambdaY, method = method)

  latentEtas <- getLatentEtasQml(indsEtas, method = method)
  colsU      <- getColsU(etas, indsEtas, lambdaY, method = method)

  fullL2      <- constructFullL2(colsU, etas = etas, method = method)
  selectSubL2 <- getSelectSubL2(fullL2, colsU = colsU, latentEtas = latentEtas,
                                method = method)
  fullSigma2ThetaEpsilon <- constructFullSigma2ThetaEpsilon(psi, method = method)
  selectSubSigma2ThetaEpsilon <-
    getSelectSubSigma2ThetaEpsilon(fullSigma2ThetaEpsilon, latentEtas = latentEtas,
                                   method = method)
  fullU <- constructFullU(fullL2 = fullL2, N = NROW(data), etas = etas, method = method)

  scalingInds <- getScalingInds(indsEtas, R = emptyR, latentEtas = latentEtas,
                                method = method)
  selectThetaEpsilon <- selectThetaEpsilon(lambdaY, thetaEpsilon,
                                           scalingInds, method = method)
  subThetaEpsilon <- constructSubThetaEpsilon(indsEtas, thetaEpsilon,
                                              scalingInds, method = method)

  covModel <- covModel(cov.syntax, method = method, parTable = parTableCovModel,
                       xis.main = xis, parTable.main = parTable)

  # list of matrices
  matrices <- list(
    lambdaX      = lambdaX,
    lambdaY      = lambdaY,
    gammaXi      = gammaXi,
    gammaEta     = gammaEta,
    thetaDelta   = thetaDelta,
    thetaEpsilon = thetaEpsilon,
    phi          = phi,
    A            = A,
    Ieta         = Ieta,
    psi          = psi,
    tauX         = tauX,
    tauY         = tauY,
    alpha        = alpha,
    beta0        = beta0,
    omegaEtaXi   = omegaEtaXi,
    omegaXiXi    = omegaXiXi,

    selectScalingY     = selectScalingY,
    selectThetaEpsilon = selectThetaEpsilon,
    selectBetaRows     = selectBetaRows,

    emptyR = emptyR,
    fullR  = fullR,

    fullSigma2ThetaEpsilon      = fullSigma2ThetaEpsilon,
    selectSubSigma2ThetaEpsilon = selectSubSigma2ThetaEpsilon,

    fullL2      = fullL2,
    selectSubL2 = selectSubL2,

    fullU = fullU,
    colsU = colsU,
    colsR = colnames(emptyR),
    rowsR = rownames(emptyR),

    subThetaEpsilon = subThetaEpsilon)

  labelMatrices <- list(
    lambdaX      = labelLambdaX,
    lambdaY      = labelLambdaY,
    gammaXi      = labelGammaXi,
    gammaEta     = labelGammaEta,
    thetaDelta   = thetaLabelDelta,
    thetaEpsilon = thetaLabelEpsilon,

    phi   = labelPhi,
    A     = labelA,
    psi   = labelPsi,
    tauX  = labelTauX,
    tauY  = labelTauY,
    alpha = labelAlpha,
    beta0 = labelBeta0,

    omegaEtaXi = labelOmegaEtaXi,
    omegaXiXi  = labelOmegaXiXi)

  k <- omegaAndSortedXis$k
  quad <- quadrature(m, k, quad.range = quad.range, adaptive = adaptive.quad, 
                     adaptive.frequency = adaptive.frequency)

  model <- list(
    info = list(
      N             = NROW(data),
      ncol          = NCOL(data),
      xis           = xis,
      etas          = etas,
      numXis        = numXis,
      numEtas       = numEtas,
      indsXis       = indsXis,
      indsEtas      = indsEtas,
      allIndsXis    = allIndsXis,
      allIndsEtas   = allIndsEtas,
      varsInts      = varsInts,
      latentEtas    = latentEtas,
      scalingInds   = scalingInds,
      kOmegaEta     = getK_NA(omegaEtaXi, labelOmegaEtaXi),
      nonLinearXis  = nonLinearXis,
      mean.observed = mean.observed,
      lavOptimizerSyntaxAdditions = lavOptimizerSyntaxAdditions
    ),

    data          = data,
    quad          = quad,
    matrices      = matrices,
    labelMatrices = labelMatrices,
    syntax        = syntax,
    cov.syntax    = cov.syntax,
    parTable      = parTable,
    covModel      = covModel
  )

  model$constrExprs <- getConstrExprs(parTable, model$covModel$parTable)
  if (createTheta) {
    listTheta         <- createTheta(model)
    model             <- c(model, listTheta)
    model$freeParams  <- length(listTheta$theta)
    model$info$bounds <- getParamBounds(model, varParams=listTheta$diagFreeParams)
    model$gradientStruct  <- getGradientStruct(model, theta = model$theta)
  }

  if (checkModel) 
    checkModel(model = model, covModel = covModel, method = method)

  model
}


matrixToParTable <- function(matrixNA, matrixEst, matrixSE, matrixLabel,
                             op = "=~", rowsLhs = TRUE, symmetric = FALSE) {
  if (symmetric) {
    matrixNA[upper.tri(matrixNA)]       <- 0
    matrixLabel[upper.tri(matrixLabel)] <- ""
  }

  if (!rowsLhs) {
    matrixNA    <- t(matrixNA)
    matrixEst   <- t(matrixEst)
    matrixSE    <- t(matrixSE)
    matrixLabel <- t(matrixLabel)
  }

  parTable <- NULL
  for (lhs in rownames(matrixEst)) {
    for (rhs in colnames(matrixEst)) {
      if (!is.na(matrixNA[lhs, rhs]) && matrixLabel[lhs, rhs] == "") next
      newRow <- data.frame(lhs = lhs, op = op, rhs = rhs,
                           label = matrixLabel[lhs, rhs],
                           est = matrixEst[lhs, rhs],
                           std.error = matrixSE[lhs, rhs])
      parTable <- rbind(parTable, newRow)
    }
  }
  parTable
}


interceptsToParTable <- function(matrixNA, matrixEst, matrixSE, matrixLabel) {
  parTable <- matrixToParTable(matrixNA, matrixEst, matrixSE, matrixLabel,
                               op = "~1", rowsLhs = TRUE)
  if (!is.null(parTable)) parTable$rhs <- ""
  parTable
}


omegaToParTable <- function(omegaNA, omegaEst, omegaSE, omegaLabel) {
  rows <- rownames(omegaEst)
  cols <- colnames(omegaEst)

  parTable <- NULL
  for (row in rows) for (col in cols) {
    if (!is.na(omegaNA[row, col]) && omegaLabel[row, col] == "") next
    eta     <- getEtaRowLabelOmega(row)
    x       <- getXiRowLabelOmega(row)
    intTerm <- paste0(x, ":", col)

    newRow <- data.frame(lhs = eta, op = "~", rhs = intTerm,
                         label = omegaLabel[row, col], est = omegaEst[row, col],
                         std.error = omegaSE[row, col])
    parTable <- rbind(parTable, newRow)
  }
  parTable
}


mainModelToParTable <- function(finalModel, method = "lms") {
  matricesEst   <- finalModel$matrices
  matricesSE    <- finalModel$matricesSE
  matricesNA    <- finalModel$matricesNA
  matricesLabel <- finalModel$labelMatrices

  if (is.null(matricesSE)) matricesSE <- matricesNA

  etas     <- finalModel$info$etas
  numXis   <- finalModel$info$numXis
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
                             matricesLabel$omegaXiXi)
  parTable <- rbind(parTable, newRows)

  newRows <- omegaToParTable(matricesNA$omegaEtaXi,
                             matricesEst$omegaEtaXi,
                             matricesSE$omegaEtaXi,
                             matricesLabel$omegaEtaXi)
  parTable <- rbind(parTable, newRows)

  # Intercepts
  newRows <- interceptsToParTable(matricesNA$tauX,
                                  matricesEst$tauX,
                                  matricesSE$tauX,
                                  matricesLabel$tauX)
  parTable <- rbind(parTable, newRows)

  newRows <- interceptsToParTable(matricesNA$tauY,
                                  matricesEst$tauY,
                                  matricesSE$tauY,
                                  matricesLabel$tauY)
  parTable <- rbind(parTable, newRows)

  newRows <- interceptsToParTable(matricesNA$alpha,
                                  matricesEst$alpha,
                                  matricesSE$alpha,
                                  matricesLabel$alpha)
  parTable <- rbind(parTable, newRows)

  newRows <- interceptsToParTable(matricesNA$beta0,
                                  matricesEst$beta0,
                                  matricesSE$beta0,
                                  matricesLabel$beta0)
  parTable <- rbind(parTable, newRows)

  # Residual (co) variances Measurement Model
  newRows <- matrixToParTable(matricesNA$thetaDelta,
                              matricesEst$thetaDelta,
                              matricesSE$thetaDelta,
                              matricesLabel$thetaDelta,
                              op = "~~", rowsLhs = TRUE,
                              symmetric = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$thetaEpsilon,
                              matricesEst$thetaEpsilon,
                              matricesSE$thetaEpsilon,
                              matricesLabel$thetaEpsilon,
                              op = "~~", rowsLhs = TRUE,
                              symmetric = TRUE)
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
                              rowsLhs = FALSE,
                              symmetric = TRUE)
  parTable <- rbind(parTable, newRows)

  newRows <- matrixToParTable(matricesNA$psi,
                              matricesEst$psi,
                              matricesSE$psi,
                              matricesLabel$psi,
                              op = "~~",
                              rowsLhs = FALSE,
                              symmetric = TRUE)
  parTable <- rbind(parTable, newRows)

  parTable <- lapplyDf(parTable, FUN = function(x) replace(x, x == -999, NA))
  parTable
}


customParamsToParTable <- function(model, coefs, se) {
  parTable <- model$parTable
  custom   <- parTable[parTable$op == ":=", ]

  if (!NROW(custom$lhs)) return(NULL)
  parTable <- NULL
  for (i in seq_len(NROW(custom))) {
    lhs <- custom[i, "lhs"]
    rhs <- custom[i, "rhs"]

    newRow <- data.frame(lhs = lhs, op = ":=", rhs = rhs,
                         label = lhs, est = coefs[[lhs]],
                         std.error = se[[lhs]])
    parTable <- rbind(parTable, newRow)
  }
  parTable
}


modelToParTable <- function(model, coefs = NULL, se = NULL, method = "lms") {
  parTable <- rbind(mainModelToParTable(model, method = method),
                    covModelToParTable(model, method = method))

  if (!is.null(coefs) && !is.null(se) && !is.null(names(se))) {
    parTable <- rbind(parTable, customParamsToParTable(model, coefs, se))

    # this is ugly but should work...
    # due to how values are read from the matrices, std.errors are overwritten
    # by the custom parameter-values (e.g., 'X=~a*x1; a==1.2' results in a std.error of 1.2, when it should be 0)
    isLabelled <- parTable$label != ""
    labels     <- parTable[isLabelled, "label"]
    parTable[isLabelled, "std.error"] <- se[labels]
    # if the std.error of a labelled parameter is 0, it is invariant, and should be NA
    # NB: there is a very small chance that a std.error of 0 is caused by a rounding error
    parTable[isLabelled & parTable$std.error == 0, "std.error"] <- NA
  }

  parTable
}
