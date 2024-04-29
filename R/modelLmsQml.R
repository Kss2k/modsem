specifyLmsModel <- function(syntax, data, method = "lms", m = 16) {
  # The goal here is to create the model with its matrices
  parTable <- modsem::modsemify(syntax)
  # Some general information:
  structExprs <- parTable[parTable$op == "~" & 
                          parTable$rhs != "1", ]
  measrExprs <- parTable[parTable$op == "=~", ]

  # Etas ------------------------------------------------------------------------
  etas <- structExprs$lhs |>
    unique()
  numEtas <- length(etas)
  if (numEtas == 0) stop("No etas in model")
  indsEtas <- lapplyNamed(etas,
                    FUN = function(eta, measrExprs)
                      measrExprs[measrExprs$lhs == eta, "rhs"],
                    measrExprs = measrExprs,
                    names = etas)
  numIndsEtas <- vapply(indsEtas,
                       FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  
  # Xis and Interaction Terms --------------------------------------------------
  intTerms <- structExprs[grepl(":", structExprs$rhs), ] 
  
  # variablese in interaction terms
  varsInts <- lapplyNamed(intTerms$rhs,
                         FUN = stringr::str_split_1,
                         pattern = ":",
                         names = stringr::str_remove_all(intTerms$rhs, ":"))
  allVarsInInts <- unique(unlist(varsInts))
  
  # now etas are included as xis, but their loadings are in lambdaX is 0
  xis <- parTable[parTable$op == "=~" &
                  !parTable$lhs %in% etas, "lhs"] |> unique()
  if (length(xis) == 0) stop("No xis in model")

  # Sorting Xis so that it is ordered with the xis in interactions first
  omegaAndSortedXis <- sortXisAndOmega(xis, varsInts, etas, intTerms)
  xis <- omegaAndSortedXis$sortedXis
  numXis <- length(xis)

  indsXis <- lapplyNamed(xis[!xis %in% etas],
                    FUN = function(xi, measrExprs)
                      measrExprs[measrExprs$lhs == xi, "rhs"],
                    measrExprs = measrExprs,
                    names = xis[!xis %in% etas])
  numIndsXis <- vapply(indsXis,
                       FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis <- unlist(indsXis)
  numAllIndsXis <- length(allIndsXis)

  # Measurement model ----------------------------------------------------------
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
  # same as in univariate equations
  lambdaY <- matrix(0, nrow = numAllIndsEtas, ncol = numEtas,
                     dimnames = list(allIndsEtas, etas))
  lastRowPreviousEta <- 0
  for (i in seq_along(etas)) {
    rowIndices <- 1:numIndsEtas[[i]] + lastRowPreviousEta
    lambdaY[rowIndices, i] <-
      c(1, rep(NA, numIndsEtas[[i]] - 1))
    lastRowPreviousEta <- lastRowPreviousEta + numIndsEtas[[i]]
  }
  
  # Structural (linear) coefficients -------------------------------------------
  gammaXi <- matrix(0, nrow = numEtas, ncol = numXis,
                  dimnames = list(etas, xis))
  exprsGammaXi <- structExprs[structExprs$lhs %in% etas & 
                    !grepl(":", structExprs$rhs) & 
                    structExprs$rhs %in% xis, ] 
  if (nrow(exprsGammaXi) > 0) apply(exprsGammaXi, 
        MARGIN = 1,  
        FUN = function(row) gammaXi[row[["lhs"]], row[["rhs"]]] <<- NA)

  gammaEta <- matrix(0, nrow = numEtas, ncol = numEtas,
                  dimnames = list(etas, etas))
  exprsGammaEta <- structExprs[structExprs$lhs %in% etas & 
                    !grepl(":", structExprs$rhs) & 
                    structExprs$rhs %in% etas, ]
  if (nrow(exprsGammaEta) > 0) apply(exprsGammaEta, MARGIN = 1, 
          FUN = function(row) gammaEta[row[["lhs"]], row[["rhs"]]] <<- NA)

  # Covariance Structure residuals observed variables --------------------------
  thetaDelta <- matrix(0, nrow = numAllIndsXis, ncol = numAllIndsXis,
                        dimnames = list(allIndsXis, allIndsXis))
  diag(thetaDelta) <- NA

  thetaEpsilon <- matrix(0, nrow = numAllIndsEtas, ncol = numAllIndsEtas,
                          dimnames = list(allIndsEtas, allIndsEtas))
  diag(thetaEpsilon) <- NA

  # Covariance Structure latents -----------------------------------------------
  psi <- matrix(0, nrow = numEtas, ncol = numEtas,
                dimnames = list(etas, etas))
  diag(psi) <- NA

  phi <- diag(numXis)
  colnames(phi) <- rownames(phi) <- xis 
  A <- phi
  A[lower.tri(A, diag = TRUE)] <- NA
  if (method == "qml") {
    phi <- A 
    A[TRUE] <- 0
  }
  
  # Mean Structure -------------------------------------------------------------
  tauX <- matrix(NA, nrow = numAllIndsXis, ncol=1,
                 dimnames = list(allIndsXis, NULL))
  tauY <- matrix(NA, nrow = numAllIndsEtas, ncol = 1,
                 dimnames = list(allIndsEtas, NULL))

  alpha <- matrix(0, nrow=numEtas, ncol=1, dimnames = list(etas, "alpha"))
  # if (method == "qml") {
  #    alpha[TRUE] <- NA
  #    tauX[TRUE] <- 0 
  #    tauY[TRUE] <- 0
  # }

  # Quadratic Terms ------------------------------------------------------------
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
  # NA  1  0  0  0
  # NA  0  1  0  0
  # 0   0  0  NA 1  
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
  selectThetaEpsilon <- thetaEpsilon
  selectThetaEpsilon[TRUE] <- FALSE
  diag(selectThetaEpsilon)[scalingInds] <- TRUE
  subThetaEpsilon <- matrix(0, nrow = length(scalingInds), 
                            ncol = length(scalingInds), 
                            dimnames = list(scalingInds, 
                                            scalingInds))
  diag(subThetaEpsilon) <- NA
  # B = (Ieta - gammaEta - gammaXi %*% A)
  # B ^ -1 = (cholB %*% t(cholB)) ^ -1
  # B ^ -1 = t(cholB) ^ -1 %*% cholB ^ -1
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
                     N = nrow(data)),
                quad = quad,
                matrices = matrices,
                syntax = syntax,
                parTable = parTable)
  # sort Data before optimizing starting params
  sortedData <- sortData(data, allIndsXis,  allIndsEtas)
  model$data <- sortedData
  model$theta <- createParamVector(model)
  model$info$bounds <- getParamBounds(model)
  # Adding Expression for evaluating variances for etas in Phi
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
  theta <- c("lambdaX" = lambdaX,
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
  theta <- theta[is.na(theta)]
  if (is.null(start)) {
   theta <- vapply(theta,
      FUN.VALUE = vector("numeric", 1L),
      FUN = function(x) stats::runif(1)
    )
  }
  theta
}


fillModel <- function(model, theta, fillPhi = FALSE) {
  numXis <- model$info$numXis
  numEtas <- model$info$numEtas
  if (is.null(names(theta))) names(theta) <- names(model$theta)
  matrices <- model$matrices
  matrices$lambdaX[is.na(matrices$lambdaX)] <-
    fetch(theta, "lambdaX[0-9]*$")
  matrices$lambdaY[is.na(matrices$lambdaY)] <-
    fetch(theta, "lambdaY[0-9]*$")
  matrices$thetaDelta[is.na(matrices$thetaDelta)] <-
    fetch(theta, "thetaDelta[0-9]*")
  matrices$thetaEpsilon[is.na(matrices$thetaEpsilon)] <-
    fetch(theta, "thetaEpsilon[0-9]*")
  matrices$A[is.na(matrices$A)] <- fetch(theta, "^A[0-9]*$")
  matrices$phi <- fillMatrixPhi(theta, matrices, model)
  matrices$psi[is.na(matrices$psi)] <-
    fetch(theta, "^psi[0-9]*$")
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
  model$matrices <- matrices
  model
}


fillMatrixPhi <- function(theta, matrices, model) {
  if (!any(is.na(matrices$phi))) return(matrices$phi)
  phi <- matrices$phi
  phi[is.na(phi)] <- fetch(theta, "^phi[0-9]*")
  phi[upper.tri(phi)] <- t(phi[lower.tri(phi)])
  phi
}


fillSymmetric <- function(mat, values) {
  mat[is.na(mat)] <- values
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}


sortXisAndOmega <- function(xis, varsInts, etas, intTerms) {
  # i <= k, i < j, i = row, j = col
  allVarsInInts <- unique(unlist(varsInts))
  sortedXis <- c(allVarsInInts, xis[!xis %in% allVarsInInts])

  nonLinearXis <- character(0L)
  for (interaction in varsInts) {
    if (!any(interaction %in% nonLinearXis)) {
      if (all(interaction %in% etas)) 
        stop("Interactions between two endogenous variables are not allowed")
      choice <- which(!interaction %in% etas)[[1]]
      nonLinearXis <- c(nonLinearXis, interaction[[choice]])
    }
  }
  linearXis <- xis[!xis %in% nonLinearXis]
  sortedXis <- c(nonLinearXis, linearXis)

  # submatrices for omegas
  omegaXiXi <- NULL
  subOmegaXiXi <- matrix(0, nrow = length(xis), ncol = length(xis),
                         dimnames = list(sortedXis, sortedXis))
  for (eta in etas) {
    lapply(varsInts[intTerms$lhs == eta],
            FUN = function(row)
              if (all(row %in% sortedXis)) subOmegaXiXi[row[[1]], row[[2]]] <<- NA)
    omegaXiXi <- rbind(omegaXiXi, subOmegaXiXi)
  }

  omegaEtaXi <- NULL
  subOmegaEtaXi <- matrix(0, nrow = length(xis), ncol = length(etas),
                          dimnames = list(sortedXis, etas))
  for (eta in etas) {
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

  matricesNA$lambdaY[matricesEst$lambdaX == 1] <- NA
  matricesSE$lambdaY[matricesEst$lambdaX == 1] <- NA
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
  selectPhi <- seq_len(numXis)
  if (method == "lms") {
    phiNA <- matricesNA$A[selectPhi, selectPhi]
    phiEst <- matricesEst$phi[selectPhi, selectPhi]
    phiSE <- matricesSE$A[selectPhi, selectPhi]
  } else {
    phiNA <- matricesNA$phi[selectPhi, selectPhi]
    phiEst <- matricesEst$phi[selectPhi, selectPhi]
    phiSE <- matricesSE$phi[selectPhi, selectPhi]
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
