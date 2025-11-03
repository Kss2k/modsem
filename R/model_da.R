# Functions
specifModelDA_Group <- function(syntax = NULL,
                                data = NULL,
                                method = "lms",
                                m = 16,
                                cov.syntax = NULL,
                                double = FALSE,
                                parTable = NULL,
                                parTableCovModel = NULL,
                                auto.fix.first = TRUE,
                                auto.fix.single = TRUE,
                                createTheta = TRUE,
                                mean.observed = TRUE,
                                standardize.inp = FALSE,
                                standardize.out = FALSE,
                                checkModel = TRUE,
                                quad.range = Inf,
                                adaptive.quad = FALSE,
                                adaptive.frequency = 3,
                                missing = "complete",
                                orthogonal.x = FALSE,
                                orthogonal.y = FALSE,
                                cluster = NULL,
                                sampling.weights = NULL) {
  if (is.null(parTable) && !is.null(syntax)) parTable <- modsemify(syntax)
  stopif(is.null(parTable), "No parTable found")

  checkParTableDA(parTable, method = method)
  # additions to lavaan-syntax for optimizer
  lavOptimizerSyntaxAdditions <- ""

  # General Information
  higherOrderLVs <- getHigherOrderLVs(parTable)
  indsHigherOrderLVs <- getIndsLVs(parTable, lVs = higherOrderLVs, isOV = FALSE)
  ovs <- getOVs(parTable)

  # endogenous variables (etas)
  etas    <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  numEtas <- length(etas)

  indsEtas    <- getIndsLVs(parTable, lVs = etas, isOV = TRUE, ovs = ovs)
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

  indsXis    <- getIndsLVs(parTable, lVs = xis, isOV = TRUE, ovs = ovs)
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis    <- unique(unlist(indsXis))
  numAllIndsXis <- length(allIndsXis)

  # clean data
  data.cleaned <- prepDataModsemDA(
    data             = data,
    allIndsXis       = allIndsXis,
    allIndsEtas      = allIndsEtas,
    missing          = missing,
    cluster          = cluster,
    sampling.weights = sampling.weights
  )

  # measurement model x
  if (method == "qml") {
    listLambdaX <- constructLambda(xis, indsXis, parTable = parTable,
                                   auto.fix.first = auto.fix.first)
  } else {
    listLambdaX <- constructLambda(c(xis, etas), c(indsXis, indsEtas),
                                   parTable = parTable,
                                   auto.fix.first = auto.fix.first)
  }

  lambdaX      <- listLambdaX$numeric
  labelLambdaX <- listLambdaX$label

  if (method == "qml") {
    listTauX <- constructTau(xis, indsXis, parTable = parTable,
                             mean.observed = mean.observed)
  } else {
    listTauX <- constructTau(c(xis, etas), c(indsXis, indsEtas),
                             parTable = parTable,
                             mean.observed = mean.observed)
  }

  tauX      <- listTauX$numeric
  labelTauX <- listTauX$label
  lavOptimizerSyntaxAdditions <- paste0(lavOptimizerSyntaxAdditions,
                                        listTauX$syntaxAdditions)

  checkResCovX_Y(parTable = parTable, allIndsXis = allIndsXis,
                 allIndsEtas = allIndsEtas, method = method)

  if (method == "qml") {
    listThetaDelta <- constructTheta(xis, indsXis, parTable = parTable,
                                     auto.fix.single = auto.fix.single)
  } else {
    listThetaDelta <- constructTheta(c(xis, etas), c(indsXis, indsEtas),
                                     parTable = parTable,
                                     auto.fix.single = auto.fix.single)
  }

  thetaDelta      <- listThetaDelta$numeric
  thetaLabelDelta <- listThetaDelta$label

  # measurement model y
  if (method == "qml") {
    listLambdaY <- constructLambda(etas, indsEtas, parTable = parTable,
                                   auto.fix.first = auto.fix.first)
  } else {
    listLambdaY <- constructLambda(NULL, NULL, parTable = parTable,
                                   auto.fix.first = auto.fix.first)
  }

  lambdaY      <- listLambdaY$numeric
  labelLambdaY <- listLambdaY$label

  if (method == "qml") {
    listTauY <- constructTau(etas, indsEtas, parTable = parTable,
                             mean.observed = mean.observed)
  } else {
    listTauY <- constructTau(NULL, NULL, parTable = parTable,
                             mean.observed = mean.observed)
  }

  tauY      <- listTauY$numeric
  labelTauY <- listTauY$label
  lavOptimizerSyntaxAdditions <- paste0(lavOptimizerSyntaxAdditions,
                                        listTauY$syntaxAdditions)

  if (method == "qml") {
    listThetaEpsilon <- constructTheta(etas, indsEtas, parTable = parTable,
                                       auto.fix.single = auto.fix.single)
  } else {
    listThetaEpsilon <- constructTheta(NULL, NULL, parTable = parTable,
                                       auto.fix.single = auto.fix.single)
  }

  thetaEpsilon      <- listThetaEpsilon$numeric
  thetaLabelEpsilon <- listThetaEpsilon$label

  # structural model
  Ieta         <- diag(numEtas) # used for (B^-1 = (Ieta - gammaEta)^-1)
  listGammaXi  <- constructGamma(etas, xis, parTable = parTable,
                                 auto.fix.first = auto.fix.first)
  gammaXi      <- listGammaXi$numeric
  labelGammaXi <- listGammaXi$label

  listGammaEta  <- constructGamma(etas, etas, parTable = parTable,
                                  auto.fix.first = auto.fix.first)
  gammaEta      <- listGammaEta$numeric
  labelGammaEta <- listGammaEta$label

  # covariance matrices
  listPsi  <- constructPsi(etas, parTable = parTable, orthogonal.y = orthogonal.y)
  psi      <- listPsi$numeric
  labelPsi <- listPsi$label

  listPhi <- constructPhi(xis, method = method, cov.syntax = cov.syntax,
                          parTableCovModel = parTableCovModel,
                          parTable = parTable, orthogonal.x = orthogonal.x)
  phi      <- listPhi$numeric
  labelPhi <- listPhi$label

  listA <- constructA(xis, method = method, cov.syntax = cov.syntax,
                      parTableCovModel = parTableCovModel,
                      parTable = parTable, orthogonal.x = orthogonal.x)
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
  fullU <- constructFullU(fullL2 = fullL2, N = data.cleaned$n, etas = etas, method = method)

  scalingInds <- getScalingInds(indsEtas, R = emptyR, latentEtas = latentEtas,
                                method = method)

  # 1 = truly latent variables, 2 = latent variables with single indicators
  selectThetaEpsilon1 <- selectThetaEpsilon1(indsEtas, thetaEpsilon,
                                             scalingInds, method = method)
  selectThetaEpsilon2 <- selectThetaEpsilon2(indsEtas, thetaEpsilon,
                                             scalingInds, method = method)
  subThetaEpsilon1 <- constructSubThetaEpsilon1(indsEtas, thetaEpsilon,
                                                scalingInds, method = method)
  subThetaEpsilon2 <- constructSubThetaEpsilon2(indsEtas, thetaEpsilon,
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

    selectScalingY      = selectScalingY,
    selectThetaEpsilon1 = selectThetaEpsilon1,
    selectThetaEpsilon2 = selectThetaEpsilon2,
    selectBetaRows      = selectBetaRows,

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

    subThetaEpsilon1 = subThetaEpsilon1,
    subThetaEpsilon2 = subThetaEpsilon2)

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
      N             = data.cleaned$n,
      ncol          = data.cleaned$k,
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

      has.interaction    = NROW(intTerms) > 0L,
      higherOrderLVs     = higherOrderLVs,
      indsHigherOrderLVs = indsHigherOrderLVs,

      lavOptimizerSyntaxAdditions = lavOptimizerSyntaxAdditions
    ),

    data          = data.cleaned,
    data.raw      = data,
    quad          = quad,
    matrices      = matrices,
    labelMatrices = labelMatrices,
    syntax        = syntax,
    cov.syntax    = cov.syntax,
    parTable      = parTable,
    covModel      = covModel,
    lavaan.fit    = NULL
  )

  if (checkModel)
    preCheckModel(model = model, covModel = covModel, method = method,
                  missing = missing)

  model
}


specifyModelDA <- function(..., group.info, createTheta = TRUE) {
  args <- list(...)

  n.groups <- group.info$n.groups
  stopif(n.groups < 1L, "Invalid grouping structure supplied.")

  parTable    <- group.info$parTable
  parTableCov <- group.info$parTableCov
  group.col   <- parTable$group

  stopif(is.null(group.col) || max(group.col) != n.groups,
         "Number of group-specific parameter tables does not match number of groups.")

  submodels <- vector("list", length = n.groups)

  for (g in seq_len(n.groups)) {
    args.g <- args

    if (!is.null(group.info$data))
      args.g$data <- group.info$data[group.info$indices[[g]], , drop = FALSE]
    else
      args.g <- addNamedNullField(args.g, field = "data")

    args.g$parTable <- parTable[parTable$group == g, , drop = FALSE]

    if (!is.null(parTableCov))
      args.g$parTableCovModel <- parTableCov[parTableCov$group == g, , drop = FALSE]
    else
      args.g <- addNamedNullField(args.g, field = "parTableCovModel")

    submodel.g <- do.call(specifModelDA_Group, args.g)
    submodel.g$info$group <- group.info$levels[[g]]
    submodel.g$info$n.groups <- 1L

    submodels[[g]] <- submodel.g
  }

  model <- list(
    models   = submodels,
    syntax   = args$syntax,
    data.raw = group.info$data.raw,
    parTable = parTable,
    info     = list(
      n.groups      = n.groups,
      group.levels  = group.info$levels,
      group.info    = group.info,

      # Constants across groups
      xis           = submodels[[1L]]$info$xis,
      etas          = submodels[[1L]]$info$etas,
      numXis        = submodels[[1L]]$info$numXis,
      numEtas       = submodels[[1L]]$info$numEtas,
      indsXis       = submodels[[1L]]$info$indsXis,
      indsEtas      = submodels[[1L]]$info$indsEtas,
      allIndsXis    = submodels[[1L]]$info$allIndsXis,
      allIndsEtas   = submodels[[1L]]$info$allIndsEtas,
      varsInts      = submodels[[1L]]$info$varsInts,
      latentEtas    = submodels[[1L]]$info$latentEtas,
      scalingInds   = submodels[[1L]]$info$scalingInds,
      kOmegaEta     = submodels[[1L]]$info$kOmegaEta,
      nonLinearXis  = submodels[[1L]]$info$nonLinearXis,
      mean.observed = submodels[[1L]]$info$mean.observed,

      has.interaction    = submodels[[1L]]$info$has.interaction,
      higherOrderLVs     = submodels[[1L]]$info$higherOrderLVs,
      indsHigherOrderLVs = submodels[[1L]]$info$indsHigherOrderLVs,

      lavOptimizerSyntaxAdditions = submodels[[1L]]$lavOptimizerSyntaxAdditions
    ),

    params = list()
  )

  # Currenlty we assume covModel has an uniform structure
  model$params$constrExprs <- getConstrExprs(parTable, model$models[[1L]]$covModel$parTable)

  if (createTheta) {
    params <- createTheta(model, parTable.in = parTable)
    model$params[names(params)] <- params

    # TODO: Remove occurences of `model$theta`, and replace them with `model$params$theta`
    model$theta <- params$theta # an ugly design decision, that was made at the very start

    model$params$bounds <- getParamBounds(model)
    model$params$gradientStruct <- getGradientStruct(model, theta = params$theta)
  }

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


omegaToParTable <- function(omegaNA, omegaEst, omegaSE, omegaLabel, parTable.in = NULL) {
  rows <- rownames(omegaEst)
  cols <- colnames(omegaEst)

  C <- \(x, y) sprintf("%s:%s", x, y)
  getIntTerm <- function(lhs, rhs) {
    xz <- C(lhs, rhs)
    zx <- C(rhs, lhs)

    if (is.null(parTable.in))         xz
    else if (xz %in% parTable.in$rhs) xz
    else                              zx
  }

  parTable <- NULL
  for (row in rows) for (col in cols) {
    if (!is.na(omegaNA[row, col]) && omegaLabel[row, col] == "") next
    eta     <- getEtaRowLabelOmega(row)
    x       <- getXiRowLabelOmega(row)
    intTerm <- getIntTerm(lhs = x, rhs = col)

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
  parTable.in   <- finalModel$parTable

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
                             matricesLabel$omegaXiXi,
                             parTable.in = parTable.in)
  parTable <- rbind(parTable, newRows)

  newRows <- omegaToParTable(matricesNA$omegaEtaXi,
                             matricesEst$omegaEtaXi,
                             matricesSE$omegaEtaXi,
                             matricesLabel$omegaEtaXi,
                             parTable.in = parTable.in)
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

  out <- NULL
  for (i in seq_len(NROW(custom))) {
    lhs <- custom[i, "lhs"]
    rhs <- custom[i, "rhs"]

    se.i <- tryCatch(se[[lhs]], error = \(e) NA)
    newRow <- data.frame(lhs = lhs, op = ":=", rhs = rhs, label = lhs,
                         group = 0L, est = coefs[[lhs]], std.error = se.i)
    out <- rbind(out, newRow)
  }

  out
}


modelToParTable <- function(model, coefs = NULL, se = NULL, method = "lms", calc.se = TRUE, group = 1L) {
  parTable <- rbind(covModelToParTable(model, method = method),
                    mainModelToParTable(model, method = method))

  colsOut <- c("lhs", "op", "rhs", "label", "group", "est", "std.error")
  parTable$group <- group
  parTable <- parTable[colsOut]

  if (!is.null(coefs) && !is.null(se) && !is.null(names(se))) {
    # this is ugly but should work...
    # due to how values are read from the matrices, std.errors are overwritten
    # by the custom parameter-values (e.g., 'X=~a*x1; a==1.2' results in a std.error of 1.2, when it should be 0)
    isLabelled <- parTable$label != ""
    labels     <- parTable[isLabelled, "label"]
    parTable[isLabelled, "std.error"] <- se[labels]
    # if the std.error of a labelled parameter is 0, it is invariant, and should be NA
    # NB: there is a very small chance that a std.error of 0 is caused by a rounding error
    zeroStdError <- parTable$std.error == 0 & !is.na(parTable$std.error)
    parTable[isLabelled & zeroStdError, "std.error"] <- NA
  }

  if (!calc.se) parTable$std.error <- NA  # when std.errors are not computed, static constraints
                                          # will get Non-NA std.errors, which is incorrect
                                          # this is naturally corrected for when calculating the
                                          # std.errors, but not when calc.se == FALSE
  parTable[!is.na(parTable$std.error) &
           parTable$std.error == -999, "std.error"] <- NA  # replace -999 with NA

  parTable
}


getConvergenceMessage <- function(converged, iterations) {
  pattern <- if (isTRUE(converged)) {
    "\nmodsem (%s) ended normally after %d iterations\n\n"
  } else {
    paste0(
           "modsem (%s) did NOT end normally after %d iterations\n",
           "** WARNING ** Estimates below are most likely unreliable\n"
    )
  }
  sprintf(pattern, PKG_INFO$version, iterations)
}


finalizeModelEstimatesDA <- function(model,
                                     theta,
                                     method = c("lms","qml"),
                                     data,
                                     logLik,
                                     iterations,
                                     converged,
                                     optimizer,
                                     # SE / FIM controls
                                     calc.se = TRUE,
                                     FIM = "observed",
                                     OFIM.hessian = FALSE,
                                     EFIM.S = 3e4,
                                     EFIM.parametric = TRUE,
                                     robust.se = FALSE,
                                     epsilon = 1e-6,
                                     cr1s = TRUE,
                                     R.max = 1e6,
                                     verbose = FALSE,
                                     # method-specific extras
                                     P = NULL,                 # only used by LMS
                                     includeStartModel = FALSE,
                                     startModel = NULL) {

  method <- match.arg(method)
  NA__ <- -999

  # coefficients and (final) filled model
  lavCoefs   <- getLavCoefs(model = model, theta = theta, method = method)

  # information matrix + SE
  typeSE <- if (!calc.se) "none" else if (robust.se) "robust" else "standard"

  fim.args <- list(model = model,
                   theta = theta,
                   method = method,
                   EFIM.S = EFIM.S,
                   hessian = OFIM.hessian,
                   calc.se = calc.se,
                   EFIM.parametric = EFIM.parametric,
                   verbose = verbose,
                   FIM = FIM,
                   robust.se = robust.se,
                   epsilon = epsilon,
                   cr1s = cr1s,
                   R.max = R.max,
                   NA__ = NA__)

  # only LMS uses the conditional P; QML ignores it
  if (!is.null(P) && method == "lms") fim.args$P <- P

  FIMo <- do.call(calcFIM_da, fim.args)

  SE <- calcSE_da(calc.se = calc.se,
                  vcov = FIMo$vcov.all,
                  rawLabels = FIMo$raw.labels,
                  NA__ = NA__)

  modelSE <- getSE_Model(model, se = SE, method = method,
                         n.additions = FIMo$n.additions)

  finalModel <- getFinalModel(model = model, theta = theta, method = method,
                              modelSE = modelSE)

  parTable <- NULL
  for (g in seq_len(model$info$n.groups)) {
    submodel <- finalModel$models[[g]]

    parTable.g <- modelToParTable(submodel,
                                  coefs = lavCoefs$all,
                                  se = SE,
                                  method = method,
                                  calc.se = calc.se,
                                  group = g)

    parTable <- rbind(parTable, parTable.g)
  }

  parTable <- rbind(parTable,
                    customParamsToParTable(model, coefs = lavCoefs$all, se = SE))

  parTable <- addZStatsParTable(parTable)
  parTable <- sortParTableDA(parTable = parTable, model = model)

  out <- list(
    model            = finalModel,
    method           = method,
    optimizer        = optimizer,
    data             = data,
    theta            = theta,
    coefs.all        = lavCoefs$all,
    coefs.free       = lavCoefs$free,
    parTable         = modsemParTable(parTable),
    originalParTable = model$parTable,
    logLik           = logLik,
    iterations       = iterations,
    convergence      = isTRUE(converged),
    convergence.msg  = getConvergenceMessage(converged, iterations),
    type.se          = typeSE,
    type.estimates   = "unstandardized",
    info.quad        = if (method == "lms") getInfoQuad(model$quad) else NULL,
    FIM              = FIMo$FIM,
    vcov.all         = FIMo$vcov.all,
    vcov.free        = FIMo$vcov.free,
    information      = FIMo$type
  )

  if (isTRUE(includeStartModel))
    out$start.model <- startModel

  out
}


addZStatsParTable <- function(parTable, se.col = "std.error", est.col = "est",
                              z.col = "z.value", p.col = "p.value",
                              ci.l = "ci.lower", ci.u = "ci.upper") {
  parTable[[z.col]] <- parTable[[est.col]] / parTable[[se.col]]
  parTable[[p.col]] <- 2 * stats::pnorm(-abs(parTable[[z.col]]))
  parTable[[ci.l]]  <- parTable[[est.col]] - CI_WIDTH * parTable[[se.col]]
  parTable[[ci.u]]  <- parTable[[est.col]] + CI_WIDTH * parTable[[se.col]]

  parTable
}


getFinalModel <- function(model, theta, method, modelSE = NULL) {
  finalModel <- fillModel(model, theta, fillPhi = method == "lms", method = method)

  # keep NA "skeletons" for printing and SE attachment
  emptyModel <- getEmptyModel(group.info = model$info$group.info,
                              cov.syntax = model$models[[1L]]$cov.syntax,
                              parTableCovModel = model$models[[1L]]$covModel$parTable,
                              mean.observed = model$info$mean.observed,
                              method = method)

  for (g in seq_along(finalModel$models)) {
    submodel <- finalModel$models[[g]]

    if (!is.null(modelSE)) {
      submodel$matricesSE <- modelSE$models[[g]]$matrices
      submodel$covModelSE <- modelSE$models[[g]]$covModel
    }

    submodel$matricesNA <- emptyModel$models[[g]]$matrices
    submodel$covModelNA <- emptyModel$models[[g]]$covModel

    finalModel$models[[g]] <- submodel
  }

  finalModel
}
