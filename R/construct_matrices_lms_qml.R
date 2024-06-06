# Functions for constructing matrices for LMS and QML. 
# Last updated: 06.06.2024

constructLambda <- function(lVs, indsLVs) {
  numLVs <- length(lVs) 
  indsLVs <- indsLVs[lVs] # make sure it is sorted
  numIndsLVs <- lapply(indsLVs, FUN = length)
  allIndsLVs <- unlist(indsLVs)
  numAllIndsLVs <- length(allIndsLVs)

  lambda <- matrix(0, nrow = numAllIndsLVs, ncol = numLVs,
                    dimnames = list(allIndsLVs, lVs))
  lastRowPreviousLV <- 0
  for (i in seq_along(lVs)) {
    rowIndices <- 1:numIndsLVs[[i]] + lastRowPreviousLV
    lambda[rowIndices, i] <-
      c(1, rep(NA, numIndsLVs[[i]] - 1))
    lastRowPreviousLV <- lastRowPreviousLV + numIndsLVs[[i]]
  }
  lambda
}


constructTau <- function(lVs, indsLVs) {
  indsLVs <- indsLVs[lVs] # make sure it is sorted
  numIndsLVs <- lapply(indsLVs, FUN = length)
  allIndsLVs <- unlist(indsLVs)
  numAllIndsLVs <- length(allIndsLVs)
  matrix(NA, nrow = numAllIndsLVs, ncol = 1,
         dimnames = list(allIndsLVs, NULL))
}


constructTheta <- function(lVs, indsLVs) {
  numLVs <- length(lVs) 
  indsLVs <- indsLVs[lVs] # make sure it is sorted
  numIndsLVs <- lapply(indsLVs, FUN = length)
  allIndsLVs <- unlist(indsLVs)
  numAllIndsLVs <- length(allIndsLVs)

  theta <- matrix(0, nrow = numAllIndsLVs, ncol = numAllIndsLVs,
                  dimnames = list(allIndsLVs, allIndsLVs))
  diag(theta) <- NA
  for (lV in lVs) { # set to 0 if there is only a single indicator
    if (numIndsLVs[[lV]] > 1) next
    theta[indsLVs[[lV]], indsLVs[[lV]]] <- 0
  }
  theta
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
                            structExprs$rhs %in% IVs, ] 
  if (NROW(exprsGamma) == 0) return(gamma)
  apply(exprsGamma, MARGIN = 1, FUN = function(row) 
        gamma[row[["lhs"]], row[["rhs"]]] <<- NA)
  gamma
}


constructPsi <- function(etas) {
  numEtas <- length(etas)
  psi <- matrix(0, nrow = numEtas, ncol = numEtas,
                dimnames = list(etas, etas)) 
  diag(psi) <- NA
  psi
}


constructPhi <- function(xis, method = "lms", cov_syntax = NULL) {
  numXis <- length(xis)
  phi <- matrix(0, nrow = numXis, ncol = numXis,
                dimnames = list(xis, xis))
  if (method != "lms" && is.null(cov_syntax)) {
    phi[lower.tri(phi, diag = TRUE)] <- NA
  }
  phi
}


constructA <- function(xis, method = "lms", cov_syntax = NULL) {
  numXis <- length(xis)
  A <- matrix(0, nrow = numXis, ncol = numXis,
              dimnames = list(xis, xis))
  if (method == "lms" && is.null(cov_syntax)) {
    A[lower.tri(A, diag = TRUE)] <- NA
  }
  A
}


constructAlpha <- function(etas) {
  numEtas <- length(etas)
  matrix(0, nrow = numEtas, ncol = 1, 
         dimnames = list(etas, "alpha"))
}


selectScalingY <- function(lambdaY, method = "qml") {
  if (method != "qml") return(NULL)
  !is.na(lambdaY)
}


selectBetaRows <- function(lambdaY, method = "qml") {
  if (method != "qml") return(NULL)
  apply(lambdaY, MARGIN = 1, FUN = function(row) any(is.na(row)))
}


constructR <- function(etas, indsEtas, lambdaY, method = "qml") {
  if (method != "qml") return(NULL)
  numEtas <- length(etas)
  numIndsEtas <- lapply(indsEtas, FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  selectBetaRows <- selectBetaRows(lambdaY, method = method)

  R <- matrix(0, nrow = numAllIndsEtas - numEtas, 
              ncol = numAllIndsEtas, 
              dimnames = list(rownames(lambdaY)[selectBetaRows], 
                              allIndsEtas))

  lastRow <- lastCol <- 0
  for (i in seq_len(numEtas)) {
    nInds <- numIndsEtas[[etas[[i]]]] - 1
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
