# Functions for Specifiying lms and qml model. modsem(method = c("lms", "qml"))
# Last updated: 29.05.2024


specifyModelLmsQml <- function(syntax, data = NULL, method = "lms", m = 16,
                               cov_syntax = NULL, double = FALSE) {
  parTable <- modsem::modsemify(syntax)
  structExprs <- parTable[parTable$op == "~" & 
                          parTable$rhs != "1", ]
  measrExprs <- parTable[parTable$op == "=~", ]

  # endogenous variables (etas)
  etas <- getSortedEtas(structExprs)
  numEtas <- length(etas)
  if (numEtas == 0) stop2("No etas in model")

  indsEtas <- lapplyNamed(etas, FUN = function(eta, measrExprs)
                          measrExprs[measrExprs$lhs == eta, "rhs"],
                          measrExprs = measrExprs,
                          names = etas)
  numIndsEtas <- vapply(indsEtas,
                       FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  
  # exogenouts variables (xis) and interaction terms 
  intTerms <- structExprs[grepl(":", structExprs$rhs), ] 
  varsInts <- lapplyNamed(intTerms$rhs,
                          FUN = stringr::str_split_1,
                          pattern = ":",
                          names = stringr::str_remove_all(intTerms$rhs, ":"))
  allVarsInInts <- unique(unlist(varsInts))
  
  xis <- parTable[parTable$op == "=~" &
                  !parTable$lhs %in% etas, "lhs"] |> unique()
  if (length(xis) == 0) stop2("No xis in model")

  # Sorting xis so that it is ordered with the xis in interactions first
  omegaAndSortedXis <- sortXisAndOmega(xis, varsInts, etas, intTerms,
                                       method = method, double = double)
  xis <- omegaAndSortedXis$sortedXis
  numXis <- length(xis)

  indsXis <- lapplyNamed(xis[!xis %in% etas], FUN = function(xi, measrExprs)
                         measrExprs[measrExprs$lhs == xi, "rhs"],
                         measrExprs = measrExprs,
                         names = xis[!xis %in% etas])
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis <- unlist(indsXis)
  numAllIndsXis <- length(allIndsXis)

  # measurement model x
  lambdaX <- matrix(0, nrow = numAllIndsXis, ncol = numXis,
                     dimnames = list(allIndsXis, xis))
  lastRowPreviousXi <- 0
  for (i in seq_along(xis)) {
    if (xis[[i]] %in% etas) next # if var is an eta, we should just skip
    rowIndices <- 1:numIndsXis[[xis[[i]]]] + lastRowPreviousXi
    lambdaX[rowIndices, i] <-
      c(1, rep(NA, numIndsXis[[xis[[i]]]] - 1))
    lastRowPreviousXi <- lastRowPreviousXi + numIndsXis[[xis[[i]]]]
  }

  tauX <- matrix(NA, nrow = numAllIndsXis, ncol=1,
                 dimnames = list(allIndsXis, NULL))

  thetaDelta <- matrix(0, nrow = numAllIndsXis, ncol = numAllIndsXis,
                       dimnames = list(allIndsXis, allIndsXis))
  diag(thetaDelta) <- NA

  # measurement model y
  lambdaY <- matrix(0, nrow = numAllIndsEtas, ncol = numEtas,
                    dimnames = list(allIndsEtas, etas))
  lastRowPreviousEta <- 0
  for (i in seq_along(etas)) {
    rowIndices <- 1:numIndsEtas[[i]] + lastRowPreviousEta
    lambdaY[rowIndices, i] <-
      c(1, rep(NA, numIndsEtas[[i]] - 1))
    lastRowPreviousEta <- lastRowPreviousEta + numIndsEtas[[i]]
  }
  
  tauY <- matrix(NA, nrow = numAllIndsEtas, ncol = 1,
                 dimnames = list(allIndsEtas, NULL))

  thetaEpsilon <- matrix(0, nrow = numAllIndsEtas, ncol = numAllIndsEtas,
                         dimnames = list(allIndsEtas, allIndsEtas))
  diag(thetaEpsilon) <- NA

  # structural model 
  gammaXi <- matrix(0, nrow = numEtas, ncol = numXis,
                  dimnames = list(etas, xis))
  exprsGammaXi <- structExprs[structExprs$lhs %in% etas & 
                    !grepl(":", structExprs$rhs) & 
                    structExprs$rhs %in% xis, ] 
  if (nrow(exprsGammaXi) > 0) {
    apply(exprsGammaXi, MARGIN = 1, FUN = function(row) 
          gammaXi[row[["lhs"]], row[["rhs"]]] <<- NA)
  }

  gammaEta <- matrix(0, nrow = numEtas, ncol = numEtas,
                  dimnames = list(etas, etas))
  exprsGammaEta <- structExprs[structExprs$lhs %in% etas & 
                    !grepl(":", structExprs$rhs) & 
                    structExprs$rhs %in% etas, ]
  if (nrow(exprsGammaEta) > 0) { 
    apply(exprsGammaEta, MARGIN = 1, FUN = function(row) 
          gammaEta[row[["lhs"]], row[["rhs"]]] <<- NA)
  }

  # residuals etas
  psi <- matrix(0, nrow = numEtas, ncol = numEtas,
                dimnames = list(etas, etas))
  diag(psi) <- NA
  
  # covariance matrix xis
  phi <- diag(numXis)
  colnames(phi) <- rownames(phi) <- xis 
  A <- phi
  A[lower.tri(A, diag = TRUE)] <- NA
  if (method == "qml") {
    phi <- A 
    A[TRUE] <- 0
  } 
  if (!is.null(cov_syntax)) {
    phi[TRUE] <- 0
    A[TRUE] <- 0
  }
  
  # mean etas
  alpha <- matrix(0, nrow=numEtas, ncol=1, dimnames = list(etas, "alpha"))

  # quadratic terms 
  omegaEtaXi <- omegaAndSortedXis$omegaEtaXi
  omegaXiXi <- omegaAndSortedXis$omegaXiXi

  # Select scaling variables for qml 
  selectScalingY <- !is.na(lambdaY)
  selectBetaRows <- apply(lambdaY, MARGIN = 1, 
                            FUN = function(row) any(is.na(row)))
  emptyR <- matrix(0, nrow = length(allIndsEtas) - numEtas, 
                      ncol = length(allIndsEtas))
  rownames(emptyR) <- 
    rownames(lambdaY)[selectBetaRows]
  colnames(emptyR) <- allIndsEtas

  lastRow <- 0
  lastCol <- 0
  for (i in seq_len(numEtas)) {
    nInds <- numIndsEtas[[etas[[i]]]] - 1
    # free params 
    emptyR[seq_len(nInds) + lastRow, lastCol + 1] <- NA
    emptyR[seq_len(nInds) + lastRow, 
              seq_len(nInds) + lastCol + 1] <- diag(nInds)
    lastRow <- lastRow + nInds
    lastCol <- lastCol + nInds + 1
  }

  # theta epsilon
  scalingInds <- allIndsEtas[!allIndsEtas %in% rownames(emptyR)]  
  selectThetaEpsilon <- as.logical.matrix(thetaEpsilon)
  selectThetaEpsilon[TRUE] <- FALSE
  diag(selectThetaEpsilon)[scalingInds] <- TRUE
  subThetaEpsilon <- matrix(0, nrow = length(scalingInds), 
                            ncol = length(scalingInds), 
                            dimnames = list(scalingInds, 
                                            scalingInds))
  diag(subThetaEpsilon) <- NA
  Ieta <- diag(numEtas)
  
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
  
  if (!is.null(data)) {
    # sort Data before optimizing starting params
    sortedData <- sortData(data, allIndsXis,  allIndsEtas)
    model$data <- sortedData
    completeCases <- stats::complete.cases(model$data)
    if (any(!completeCases)) {
      warning2("Removing missing values case-wise.")
      model$data <- model$data[completeCases, ]
      model$info$N <- nrow(model$data)
    }
  } else {
    model$data <- NULL
    model$info$N <- NA
  }

  model
}


isInList <- function(elems, list) {
  matches <- vapply(list,
    FUN.VALUE = vector("logical", 1L),
    FUN = function(node, elems) {
      all(elems %in% node)
    },
    elems = elems
  )
  any(matches)
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


sortXisAndOmega <- function(xis, varsInts, etas, intTerms, method = "lms",
                            double = FALSE) {
  # allVarsInInts should be sorted according to which variables 
  # occur in the most interaction terms (makes it more efficient)
  allVarsInInts <- unique(unlist(varsInts))
  freqInIntTerms <- lapply(varsInts, FUN = unique) |> unlist() |>
    table() |> oneWayTableToDataFrame()

  sortedXis <- c(allVarsInInts, xis[!xis %in% allVarsInInts])
  nonLinearXis <- character(0L)
  for (interaction in varsInts) {
    if (any(interaction %in% nonLinearXis) && !double ||
        all(interaction %in% nonLinearXis) && double) next # no need to add it again

    if (length(interaction) > 2) {
      stop2("Only interactions between two variables are allowed")
    } else if (all(interaction %in% etas)) {
      stop2("Interactions between two endogenous variables are not allowed, ",
           "see \nvignette(\"interaction_two_etas\", \"modsem\")")
    } 

    choice <- interaction[which(!interaction %in% etas)] 
    if (length(choice) > 1 && !double) {
      freq <- freqInIntTerms[choice, "freq"]
      choice <- choice[whichIsMax(freq)][[1]] # pick first if both are equal
    }

    nonLinearXis <- c(nonLinearXis, choice)
  }
  linearXis <- xis[!xis %in% nonLinearXis]
  sortedXis <- c(nonLinearXis, linearXis)

  # submatrices for omegas
  omegaXiXi <- NULL
  for (eta in etas) {
    subOmegaXiXi <- matrix(0, nrow = length(xis), ncol = length(xis),
                           dimnames = list(sortedXis, sortedXis))
    lapply(varsInts[intTerms$lhs == eta], FUN = function(row) {
       if (!all(row %in% sortedXis)) return(NULL) 
       whichRow <- which(row %in% nonLinearXis)[[1]] # if quadratic term pick first
       whichCol <- ifelse(whichRow == 1, 2, 1)
       subOmegaXiXi[row[[whichRow]], row[[whichCol]]] <<- NA
    })
    omegaXiXi <- rbind(omegaXiXi, subOmegaXiXi)
  }

  omegaEtaXi <- NULL
  for (eta in etas) {
    subOmegaEtaXi <- matrix(0, nrow = length(xis), ncol = length(etas),
                            dimnames = list(sortedXis, etas))
    lapply(varsInts[intTerms$lhs == eta], FUN = function(row) {
       if (any(row %in% etas) & !all(row %in% etas)) {
         whichXi <- which(!row %in% etas)
         whichEta <- which(row %in% etas)
         subOmegaEtaXi[row[[whichXi]], row[[whichEta]]] <<- NA
       }
    })
    omegaEtaXi <- rbind(omegaEtaXi, subOmegaEtaXi)
  }
  list(sortedXis = sortedXis, omegaXiXi = omegaXiXi, omegaEtaXi = omegaEtaXi,
       k = length(nonLinearXis))
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


# Vector for diagonal indices
diag_ind <- function(num) diag(matrix(seq_len(num^2), num))


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
