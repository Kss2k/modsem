# Functions for constructing matrices for LMS and QML. 
# Last updated: 06.06.2024

constructLambda <- function(lVs, indsLVs, parTable, autoConstraints = TRUE) {
  numLVs <- length(lVs) 
  indsLVs <- indsLVs[lVs] # make sure it is sorted
  numIndsLVs <- lapply(indsLVs, FUN = length)
  allIndsLVs <- unlist(indsLVs)
  numAllIndsLVs <- length(allIndsLVs)

  lambda <- matrix(0, nrow = numAllIndsLVs, ncol = numLVs,
                    dimnames = list(allIndsLVs, lVs))
  lastRowPreviousLV <- 0
  if (autoConstraints) firstVal <- 1 else firstVal <- NA
  for (i in seq_along(lVs)) {
    rowIndices <- seq_len(numIndsLVs[[i]]) + lastRowPreviousLV
    lambda[rowIndices, i] <-
      c(firstVal, rep(NA, numIndsLVs[[i]] - 1))
    lastRowPreviousLV <- lastRowPreviousLV + numIndsLVs[[i]]
  }

  constExprs <- parTable[parTable$op == "=~" & 
                         parTable$rhs %in% allIndsLVs &
                         parTable$lhs %in% lVs &
                         canBeNumeric(parTable$mod), ]
  for (i in seq_len(NROW(constExprs))) {
    lhs <- constExprs[i, "lhs"]
    rhs <- constExprs[i, "rhs"]
    lambda[rhs, lhs] <- as.numeric(constExprs[i, "mod"])
  }

  labelLambda <- as.character.matrix(lambda, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == "=~" & 
                           parTable$rhs %in% allIndsLVs &
                           parTable$lhs %in% lVs &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    rhs <- dynamicExprs[i, "rhs"]
    mod <- dynamicExprs[i, "mod"]
    lambda[rhs, lhs] <- 0
    labelLambda[rhs, lhs] <- mod
  }
  
  list(numeric = lambda, label = labelLambda)
}


constructTau <- function(lVs, indsLVs, parTable, meanStructure = TRUE) {
  indsLVs <- indsLVs[lVs] # make sure it is sorted
  numIndsLVs <- lapply(indsLVs, FUN = length)
  allIndsLVs <- unlist(indsLVs)
  numAllIndsLVs <- length(allIndsLVs)

  if (meanStructure) default <- NA else default <- 0

  tau <- matrix(default, nrow = numAllIndsLVs, ncol = 1,
                dimnames = list(allIndsLVs, "1"))
  constExprs <- parTable[parTable$op == "~" & 
                         parTable$rhs == "1" &
                         parTable$lhs %in% allIndsLVs &
                         canBeNumeric(parTable$mod,
                                      includeNA = TRUE), ]
  for (i in seq_len(NROW(constExprs))) {
    lhs <- constExprs[i, "lhs"]
    tau[lhs, 1] <- as.numeric(constExprs[i, "mod"])
  }
  
  labelTau <- as.character.matrix(tau, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == "~" & 
                           parTable$rhs == "1" &
                           parTable$lhs %in% allIndsLVs &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    mod <- dynamicExprs[i, "mod"]
    tau[lhs, 1] <- 0
    labelTau[lhs, 1] <- mod
  }
  
  list(numeric = tau, label = labelTau)
}


constructTheta <- function(lVs, indsLVs, parTable, autoConstraints = TRUE) {
  numLVs <- length(lVs) 
  indsLVs <- indsLVs[lVs] # make sure it is sorted
  numIndsLVs <- lapply(indsLVs, FUN = length)
  allIndsLVs <- unlist(indsLVs)
  numAllIndsLVs <- length(allIndsLVs)

  theta <- matrix(0, nrow = numAllIndsLVs, ncol = numAllIndsLVs,
                  dimnames = list(allIndsLVs, allIndsLVs))
  diag(theta) <- NA

  if (autoConstraints) {
    for (lV in lVs) { # set to 0 if there is only a single indicator
      if (numIndsLVs[[lV]] == 1) theta[indsLVs[[lV]], indsLVs[[lV]]] <- 0
    }
  }

  constExprs <- parTable[parTable$op == "~~" & 
                         parTable$lhs %in% allIndsLVs &
                         parTable$rhs %in% allIndsLVs &
                         canBeNumeric(parTable$mod, 
                                      includeNA = TRUE), ] 
  for (i in seq_len(NROW(constExprs))) {
    lhs <- constExprs[i, "lhs"]
    rhs <- constExprs[i, "rhs"]
    theta[lhs, rhs] <- theta[rhs, lhs] <- as.numeric(constExprs[i, "mod"])
  }
  
  labelTheta <- as.character.matrix(theta, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == "~~" & 
                           parTable$lhs %in% allIndsLVs &
                           parTable$rhs %in% allIndsLVs &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    rhs <- dynamicExprs[i, "rhs"]
    mod <- dynamicExprs[i, "mod"]
    theta[lhs, rhs] <- theta[rhs, lhs] <- 0
    labelTheta[lhs, rhs] <- labelTheta[rhs, lhs] <- mod
  }
  theta[upper.tri(theta)] <- 0

  list(numeric = theta, label = labelTheta)
}


constructGamma <- function(DVs, IVs, parTable) {
  structExprs <- parTable[parTable$op == "~" &
                          parTable$rhs != "1", ]
  numDVs <- length(DVs)
  numIVs <- length(IVs)
  gamma <- matrix(0, nrow = numDVs, ncol = numIVs,
                  dimnames = list(DVs, IVs))

  exprsGamma <- structExprs[structExprs$lhs %in% DVs & 
                            !grepl(":", structExprs$rhs) & 
                            structExprs$rhs %in% IVs & 
                            canBeNumeric(structExprs$mod,
                                         includeNA = TRUE), ] 

  for (i in seq_len(NROW(exprsGamma))) {
    lhs <- exprsGamma[i, "lhs"]
    rhs <- exprsGamma[i, "rhs"]
    gamma[lhs, rhs] <- as.numeric(exprsGamma[i, "mod"])
  }

  labelGamma <- as.character.matrix(gamma, empty = TRUE)
  dynamicExprs <- structExprs[structExprs$lhs %in% DVs & 
                              !grepl(":", structExprs$rhs) & 
                              structExprs$rhs %in% IVs & 
                              !canBeNumeric(structExprs$mod,
                                            includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    rhs <- dynamicExprs[i, "rhs"]
    mod <- dynamicExprs[i, "mod"]
    gamma[lhs, rhs] <- 0
    labelGamma[lhs, rhs] <- mod
  }

  list(numeric = gamma, label = labelGamma)
}


constructPsi <- function(etas, parTable) {
  numEtas <- length(etas)
  psi <- matrix(0, nrow = numEtas, ncol = numEtas,
                dimnames = list(etas, etas)) 
  diag(psi) <- NA

  constExprs <- parTable[parTable$op == "~~" & 
                         parTable$lhs %in% etas &
                         parTable$rhs %in% etas &
                         canBeNumeric(parTable$mod,
                                      includeNA = TRUE), ]
  for (i in seq_len(NROW(constExprs))) {
    lhs <- constExprs[i, "lhs"]
    rhs <- constExprs[i, "rhs"]
    psi[lhs, rhs] <- psi[rhs, lhs] <- as.numeric(constExprs[i, "mod"])
  }
  psi[upper.tri(psi)] <- 0 

  labelPsi <- as.character.matrix(psi, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == "~~" & 
                           parTable$lhs %in% etas &
                           parTable$rhs %in% etas &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    rhs <- dynamicExprs[i, "rhs"]
    mod <- dynamicExprs[i, "mod"]
    psi[lhs, rhs] <- psi[rhs, lhs] <- 0
    labelPsi[lhs, rhs] <- labelPsi[rhs, lhs] <- mod
  }

  list(numeric = psi, label = labelPsi)
}


constructPhi <- function(xis, method = "lms", cov_syntax = NULL, 
                         parTable) {
  numXis <- length(xis)
  phi <- matrix(0, nrow = numXis, ncol = numXis,
                dimnames = list(xis, xis))
  if (method != "lms" && is.null(cov_syntax)) {
    phi[lower.tri(phi, diag = TRUE)] <- NA
  }
  
  constExprs <- parTable[parTable$op == "~~" & 
                         parTable$lhs %in% xis &
                         parTable$rhs %in% xis &
                         canBeNumeric(parTable$mod), ]
  for (i in seq_len(NROW(constExprs))) {
    lhs <- constExprs[i, "lhs"]
    rhs <- constExprs[i, "rhs"]
    phi[lhs, rhs] <- phi[rhs, lhs] <- as.numeric(constExprs[i, "mod"])
  }
  phi[upper.tri(phi)] <- 0

  labelPhi <- as.character.matrix(phi, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == "~~" & 
                           parTable$lhs %in% xis &
                           parTable$rhs %in% xis &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    rhs <- dynamicExprs[i, "rhs"]
    mod <- dynamicExprs[i, "mod"]
    phi[lhs, rhs] <- phi[rhs, lhs] <- 0
    labelPhi[lhs, rhs] <- labelPhi[rhs, lhs] <- mod
  }

  list(numeric = phi, label = labelPhi)
}


constructA <- function(xis, method = "lms", cov_syntax = NULL,
                       parTable) {
  numXis <- length(xis)
  A <- matrix(0, nrow = numXis, ncol = numXis,
              dimnames = list(xis, xis))
  if (method == "lms" && is.null(cov_syntax)) {
    A[lower.tri(A, diag = TRUE)] <- NA
  }

  constExprs <- parTable[parTable$op == "~~" & 
                         parTable$lhs %in% xis &
                         parTable$rhs %in% xis &
                         canBeNumeric(parTable$mod), ]
  for (i in seq_len(NROW(constExprs))) {
    lhs <- constExprs[i, "lhs"]
    rhs <- constExprs[i, "rhs"]
    A[lhs, rhs] <- A[lhs, rhs] <- as.numeric(constExprs[i, "mod"])
  }
  A[upper.tri(A)] <- 0

  labelA <- as.character.matrix(A, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == "~~" & 
                           parTable$lhs %in% xis &
                           parTable$rhs %in% xis &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    rhs <- dynamicExprs[i, "rhs"]
    mod <- dynamicExprs[i, "mod"]
    A[lhs, rhs] <- A[lhs, rhs] <- 0
    labelA[lhs, rhs] <- labelA[lhs, rhs] <- mod
  }

  list(numeric = A, label = labelA)
}


constructAlpha <- function(etas, parTable, autoConstraints = TRUE,
                           meanStructure = TRUE) {
  numEtas <- length(etas)
  if (autoConstraints && meanStructure) default <- 0 else default <- NA
  alpha <- matrix(default, nrow = numEtas, ncol = 1, 
                  dimnames = list(etas, "1"))

  constExprs <- parTable[parTable$op == "~" & 
                         parTable$rhs == "1" &
                         parTable$lhs %in% etas &
                         canBeNumeric(parTable$mod, 
                                      includeNA = TRUE), ]
  for (i in seq_len(NROW(constExprs))) {
    eta <- constExprs[i, "lhs"]
    alpha[eta, 1] <- as.numeric(constExprs[i, "mod"])
  }
  
  labelAlpha <- as.character.matrix(alpha, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == "~" & 
                           parTable$rhs == "1" &
                           parTable$lhs %in% etas &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  for (i in seq_len(NROW(dynamicExprs))) {
    eta <- dynamicExprs[i, "lhs"]
    mod <- dynamicExprs[i, "mod"]
    alpha[eta, 1] <- 0
    labelAlpha[eta, 1] <- mod
  }

  list(numeric = alpha, label = labelAlpha)
}


selectScalingY <- function(lambdaY, method = "qml") {
  if (method != "qml") return(NULL)
  matrix(apply(lambdaY, MARGIN = 2, FUN = isScalingY), 
         nrow = nrow(lambdaY), ncol = ncol(lambdaY), 
         dimnames = dimnames(lambdaY))
}


selectBetaRows <- function(lambdaY, method = "qml") {
  if (method != "qml") return(NULL)
  scalingYs <- selectScalingY(lambdaY, method = "qml")
  matrix(apply(scalingYs, MARGIN = 1, FUN = function(x )!all(x)),
         nrow = nrow(lambdaY), ncol = 1, 
         dimnames = list(rownames(lambdaY), "1"))
}


constructR <- function(etas, indsEtas, lambdaY, method = "qml") {
  if (method != "qml") return(NULL)
  numEtas <- length(etas)
  numIndsEtas <- lapply(indsEtas, FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  selectBetaRows <- selectBetaRows(lambdaY, method = method)
  rowNamesR <- rownames(lambdaY)[selectBetaRows]
  
  R <- matrix(0, nrow = numAllIndsEtas - numEtas, 
              ncol = numAllIndsEtas, 
              dimnames = list(rowNamesR[seq_len(numAllIndsEtas - numEtas)],
                              allIndsEtas))

  lastRow <- lastCol <- 0
  for (i in seq_len(numEtas)) {
    nInds <- numIndsEtas[[etas[[i]]]] - 1
    if (nInds == 0) stop("Etas in QML must have at least two indicators")
    # free params 
    R[seq_len(nInds) + lastRow, lastCol + 1] <- NA
    R[seq_len(nInds) + lastRow, 
      seq_len(nInds) + lastCol + 1] <- diag(nInds)
    lastRow <- lastRow + nInds
    lastCol <- lastCol + nInds + 1
  }
  R
}


getScalingInds <- function(indsEtas, R, method = "qml") {
  if (method != "qml") return(NULL)
  allIndsEtas <- unlist(indsEtas)
  scalingInds <- allIndsEtas[!allIndsEtas %in% rownames(R)]
  scalingInds
}


selectThetaEpsilon <- function(indsEtas, thetaEpsilon, scalingInds,
                                method = "qml") {
  if (method != "qml") return(NULL)
  selectThetaEpsilon <- as.logical.matrix(thetaEpsilon)
  selectThetaEpsilon[TRUE] <- FALSE
  diag(selectThetaEpsilon)[scalingInds] <- TRUE
  selectThetaEpsilon
}


constructSubThetaEpsilon <- function(indsEtas, thetaEpsilon, scalingInds,
                                     method = "qml") {
  if (method != "qml") return(NULL)
  subThetaEpsilon <- matrix(0, nrow = length(scalingInds), 
                            ncol = length(scalingInds), 
                            dimnames = list(scalingInds, 
                                            scalingInds))
  diag(subThetaEpsilon) <- NA
  subThetaEpsilon
}


sortXisConstructOmega <- function(xis, varsInts, etas, intTerms, 
                                  method = "lms", double = FALSE) {
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

    choice <- unique(interaction[which(!interaction %in% etas)])
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
  labelOmegaXiXi <- NULL
  for (eta in etas) {
    subOmegaXiXi <- matrix(0, nrow = length(xis), ncol = length(xis),
                           dimnames = list(sortedXis, sortedXis))
    subLabelOmegaXiXi <- as.character.matrix(subOmegaXiXi, empty = TRUE)

    lapply(varsInts[intTerms$lhs == eta], FUN = function(row) {
       if (!all(row %in% sortedXis)) return(NULL) 
       whichRow <- which(row %in% nonLinearXis)[[1]] # if quadratic term pick first
       whichCol <- ifelse(whichRow == 1, 2, 1)
       
       subOmegaXiXi[row[[whichRow]], row[[whichCol]]] <<- 
         getFreeOrConsIntTerms(row, eta, intTerms)
       subLabelOmegaXiXi[row[[whichRow]], row[[whichCol]]] <<-
         getLabelIntTerms(row, eta, intTerms)
    })
    omegaXiXi <- rbind(omegaXiXi, subOmegaXiXi)
    labelOmegaXiXi <- rbind(labelOmegaXiXi, subLabelOmegaXiXi)
  }

  omegaEtaXi <- NULL
  labelOmegaEtaXi <- NULL 
  for (eta in etas) {
    subOmegaEtaXi <- matrix(0, nrow = length(xis), ncol = length(etas),
                            dimnames = list(sortedXis, etas))
    subLabelOmegaEtaXi <- as.character.matrix(subOmegaEtaXi, empty = TRUE)

    lapply(varsInts[intTerms$lhs == eta], FUN = function(row) {
       if (any(row %in% etas) & !all(row %in% etas)) {
         whichXi <- which(!row %in% etas)
         whichEta <- which(row %in% etas)

         subOmegaEtaXi[row[[whichXi]], row[[whichEta]]] <<- 
           getFreeOrConsIntTerms(row, eta, intTerms)
         subLabelOmegaEtaXi[row[[whichXi]], row[[whichEta]]] <<- 
           getLabelIntTerms(row, eta, intTerms)
       }
    })
    omegaEtaXi <- rbind(omegaEtaXi, subOmegaEtaXi)
    labelOmegaEtaXi <- rbind(labelOmegaEtaXi, subLabelOmegaEtaXi)
  }
  list(sortedXis = sortedXis, 
       omegaXiXi = list(numeric = omegaXiXi, label = labelOmegaXiXi),
       omegaEtaXi = list(numeric = omegaEtaXi, label = labelOmegaEtaXi),
       k = length(nonLinearXis))
}
