# Functions for constructing matrices for LMS and QML.
# Last updated: 06.06.2024
setMatrixConstraints <- function(X, parTable, op, RHS, LHS, type, nonFreeParams) {
  fillConstExprs(X, parTable = parTable, op = op, RHS = RHS, LHS = LHS,
                 type = type, nonFreeParams = nonFreeParams) |>
    fillDynExprs(parTable = parTable, op = op, RHS = RHS, LHS = LHS, type = type)
}


fillConstExprs <- function(X, parTable, op, RHS, LHS, type, nonFreeParams = TRUE) {
  constExprs <- parTable[parTable$op == op &
                         parTable$rhs %in% RHS &
                         parTable$lhs %in% LHS &
                         canBeNumeric(parTable$mod, includeNA = !nonFreeParams), ]
  constExprs[constExprs$op == "~1", "rhs"] <- "1"

  setVal <- getSetValFunc(type)
  for (i in seq_len(NROW(constExprs))) {
    lhs <- constExprs[i, "lhs"]
    rhs <- constExprs[i, "rhs"]
    val <- as.numeric(constExprs[i, "mod"])
    X   <- setVal(X = X, rhs = rhs, lhs = lhs, val = val)
  }

  if (type == "symmetric") X[upper.tri(X)] <- 0
  X
}


fillDynExprs <- function(X, parTable, op, RHS, LHS, type) {
  # dynamic exprs need a corresponding matrix of labels
  labelX <- as.character.matrix(X, empty = TRUE)
  dynamicExprs <- parTable[parTable$op == op &
                           parTable$rhs %in% RHS &
                           parTable$lhs %in% LHS &
                           !canBeNumeric(parTable$mod,
                                         includeNA = TRUE), ]
  dynamicExprs[dynamicExprs$op == "~1", "rhs"] <- "1"

  setVal <- getSetValFunc(type)
  for (i in seq_len(NROW(dynamicExprs))) {
    lhs <- dynamicExprs[i, "lhs"]
    rhs <- dynamicExprs[i, "rhs"]
    mod <- dynamicExprs[i, "mod"]

    X      <- setVal(X = X, rhs = rhs, lhs = lhs, val = 0)
    labelX <- setVal(X = labelX, rhs = rhs, lhs = lhs, val = mod)
  }

  list(numeric = X, label = labelX)
}


getSetValFunc <- function(type) {
  switch(type,
         rhs       = setValRhsFirst,
         lhs       = setValLhsFirst,
         symmetric = setValSymmetric,
         stop("Unrecognized type, this is probably a bug!"))
}


setValRhsFirst <- function(X, rhs, lhs, val) {
  X[rhs, lhs] <- val
  X
}


setValLhsFirst <- function(X, rhs, lhs, val) {
  X[lhs, rhs] <- val
  X
}


setValSymmetric <- function(X, rhs, lhs, val) {
  X[lhs, rhs]     <- X[rhs, lhs] <- val
  X
}


getEmptyPhi <- function(phi) {
  labelPhi <- as.character.matrix(phi, empty = TRUE)
  list(numeric = phi, label = labelPhi)
}


notFilledLambda <- function(ind, lambda) {
  !any(is.na(lambda[ind, ])) && all(lambda[ind, ] == 0)
}


constructLambda <- function(lVs, indsLVs, parTable, auto.constraints = TRUE) {
  numLVs        <- length(lVs)
  indsLVs       <- indsLVs[lVs] # make sure it is sorted
  allIndsLVs    <- unique(unlist(indsLVs))
  numAllIndsLVs <- length(allIndsLVs)
  firstVal      <- ifelse(auto.constraints, 1, NA)

  lambda <- matrix(0, nrow = numAllIndsLVs, ncol = numLVs,
                   dimnames = list(allIndsLVs, lVs))

  for (lV in lVs) {
    firstFilled <- FALSE
    for (ind in indsLVs[[lV]]) {
      # TODO: make this work when duplicates appear seperately in lambdaY and lambdaX
      if (!firstFilled && auto.constraints && notFilledLambda(ind, lambda)) {
        lambda[ind, lV] <- firstVal
        firstFilled <- TRUE
      } else lambda[ind, lV] <- NA
    }
  }

  setMatrixConstraints(X = lambda, parTable = parTable, op = "=~",
                       RHS = allIndsLVs, LHS = lVs, type = "rhs",
                       nonFreeParams = TRUE) # first params are by default set to 1
}


constructTau <- function(lVs, indsLVs, parTable, mean.observed = TRUE) {
  indsLVs       <- indsLVs[lVs] # make sure it is sorted
  allIndsLVs    <- unique(unlist(indsLVs))
  numAllIndsLVs <- length(allIndsLVs)
  default       <- ifelse(mean.observed, NA, 0)
  lavOptimizerSyntaxAdditions <- ""

  tau <- matrix(default, nrow = numAllIndsLVs, ncol = 1,
                dimnames = list(allIndsLVs, "1"))
  for (lV in lVs) { # set first ind to 0, if lV has meanstructure
    subPT <- parTable[parTable$lhs == lV & parTable$op == "~1", ]
    if (NROW(subPT)) {
      firstInd         <- indsLVs[[lV]][[1]]
      tau[firstInd, 1] <- 0
      lavOptimizerSyntaxAdditions <-
        getFixedInterceptSyntax(indicator = firstInd, parTable = parTable,
                                syntax = lavOptimizerSyntaxAdditions)
    }
  }

  c(setMatrixConstraints(X = tau, parTable = parTable, op = "~1",
                         RHS = "", LHS = allIndsLVs, type = "lhs",
                         nonFreeParams = FALSE),
    list(syntaxAdditions = lavOptimizerSyntaxAdditions))
}


constructTheta <- function(lVs, indsLVs, parTable, auto.constraints = TRUE) {
  numLVs        <- length(lVs)
  indsLVs       <- indsLVs[lVs] # make sure it is sorted
  numIndsLVs    <- lapply(indsLVs, FUN = length)
  allIndsLVs    <- unique(unlist(indsLVs))
  numAllIndsLVs <- length(allIndsLVs)

  theta <- matrix(0, nrow = numAllIndsLVs, ncol = numAllIndsLVs,
                  dimnames = list(allIndsLVs, allIndsLVs))
  diag(theta) <- NA

  if (auto.constraints) {
    for (lV in lVs) { # set to 0 if there is only a single indicator
      if (numIndsLVs[[lV]] != 1) next
      theta[indsLVs[[lV]], indsLVs[[lV]]] <- 0
    }
  }

  setMatrixConstraints(X = theta, parTable = parTable, op = "~~",
                       RHS = allIndsLVs, LHS = allIndsLVs, type = "symmetric",
                       nonFreeParams = FALSE)
}


constructGamma <- function(DVs, IVs, parTable) {
  exprsGamma <- parTable[parTable$op == "~" & !grepl(":", parTable$rhs), ]
  numDVs <- length(DVs)
  numIVs <- length(IVs)
  gamma  <- matrix(0, nrow = numDVs, ncol = numIVs, dimnames = list(DVs, IVs))

  setMatrixConstraints(X = gamma, parTable = exprsGamma, op = "~", RHS = IVs,
                       LHS = DVs, type = "lhs", nonFreeParams = FALSE)
}


constructPsi <- function(etas, parTable) {
  numEtas   <- length(etas)
  psi       <- matrix(0, nrow = numEtas, ncol = numEtas,
                      dimnames = list(etas, etas))
  diag(psi) <- NA

  setMatrixConstraints(X = psi, parTable = parTable, op = "~~",
                       RHS = etas, LHS = etas, type = "symmetric",
                       nonFreeParams = FALSE)
}


constructPhi <- function(xis, method = "lms", cov.syntax = NULL,
                         parTable) {
  numXis <- length(xis)
  phi    <- matrix(0, nrow = numXis, ncol = numXis,
                   dimnames = list(xis, xis))
  if (method != "lms" && is.null(cov.syntax)) {
    phi[lower.tri(phi, diag = TRUE)] <- NA
    setMatrixConstraints(X = phi, parTable = parTable, op = "~~",
                         RHS = xis, LHS = xis, type = "symmetric",
                         nonFreeParams = FALSE)
  } else getEmptyPhi(phi = phi)
}


constructA <- function(xis, method = "lms", cov.syntax = NULL,
                       parTable) {
  numXis <- length(xis)
  A      <- matrix(0, nrow = numXis, ncol = numXis,
                   dimnames = list(xis, xis))
  if (method == "lms" && is.null(cov.syntax)) {
    A[lower.tri(A, diag = TRUE)] <- NA
    setMatrixConstraints(X = A, parTable = parTable, op = "~~",
                         RHS = xis, LHS = xis, type = "symmetric",
                         nonFreeParams = FALSE)
  } else getEmptyPhi(phi = A)
}


constructAlpha <- function(etas, parTable, auto.constraints = TRUE,
                           mean.observed = TRUE) {
  numEtas <- length(etas)
  if (auto.constraints && mean.observed) default <- 0 else default <- NA
  alpha <- matrix(default, nrow = numEtas, ncol = 1,
                  dimnames = list(etas, "1"))

  setMatrixConstraints(X = alpha, parTable = parTable, op = "~1",
                       RHS = "", LHS = etas, type = "lhs",
                       nonFreeParams = FALSE)
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
  hasMultipleInds <- vapply(indsEtas, FUN = function(x) length(x) > 1,
                            FUN.VALUE = logical(1L))
  etas <- names(indsEtas)[hasMultipleInds]
  indsEtas <- indsEtas[etas]
  if (length(etas) == 0) return(NULL)

  numEtas <- length(etas)
  numIndsEtas <- lapply(indsEtas, FUN = length)
  allIndsEtas <- unique(unlist(indsEtas))
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
    if (nInds == 0) stop2("Etas in QML must have at least two indicators")
    # free params
    R[seq_len(nInds) + lastRow, lastCol + 1] <- NA
    R[seq_len(nInds) + lastRow,
      seq_len(nInds) + lastCol + 1] <- diag(nInds)
    lastRow <- lastRow + nInds
    lastCol <- lastCol + nInds + 1
  }
  R
}


getLatentEtasQml <- function(indsEtas, method = "qml") {
  if (method != "qml") return(NULL)
  hasMultipleInds <- vapply(indsEtas, FUN = function(x) length(x) > 1,
                            FUN.VALUE = logical(1L))
  etas <- names(indsEtas)[hasMultipleInds]
  if (length(etas) == 0) return(NULL)
  etas
}


getColsU <- function(etas, indsEtas, lambdaY, method = "qml") {
  numEtas         <- length(etas)
  numIndsEtas     <- lapply(indsEtas, FUN = length)
  allIndsEtas     <- unique(unlist(indsEtas))
  numAllIndsEtas  <- length(allIndsEtas)
  selectBetaRows  <- selectBetaRows(lambdaY, method = "qml")
  rowNamesR       <- rownames(lambdaY)[selectBetaRows]

  hasMultipleInds <- vapply(indsEtas, FUN = function(x) length(x) > 1,
                            FUN.VALUE = logical(1L))
  if (sum(hasMultipleInds) == 0) return(NULL)

  colsU <- rowNamesR[seq_len(numAllIndsEtas - sum(hasMultipleInds))]
  colsU[is.na(colsU)] <- paste0("__FILL__ZERO__", seq_len(sum(is.na(colsU))))

  colsU
}


constructFullU <- function(fullL2, N, etas, method = "qml") {
  if (method != "qml" || N == 0) return(NULL)

  if (is.null(fullL2)) nCols <- length(etas)
  else nCols <- NCOL(fullL2)

  matrix(0, nrow = N, ncol = nCols, dimnames = list(NULL, colnames(fullL2)))
}


constructFullR <- function(etas, indsEtas, lambdaY, method = "qml") {
  if (method != "qml") return(NULL)

  numEtas <- length(etas)
  numIndsEtas <- lapply(indsEtas, FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  selectBetaRows <- selectBetaRows(lambdaY, method = "qml")
  rowNamesR <- rownames(lambdaY)[selectBetaRows]
  hasMultipleInds <- vapply(indsEtas, FUN = function(x) length(x) > 1,
                            FUN.VALUE = logical(1L))
  if (sum(hasMultipleInds) == 0) return(NULL)

  matrix(0, nrow = numAllIndsEtas - sum(hasMultipleInds),
         ncol = numAllIndsEtas,
         dimnames = list(rowNamesR[seq_len(numAllIndsEtas - sum(hasMultipleInds))],
                                           allIndsEtas))

}


constructFullSigma2ThetaEpsilon <- function(psi, method = "qml") {
  if (method != "qml") return(NULL)
  matrix(0, nrow = nrow(psi), ncol = ncol(psi),
         dimnames = dimnames(psi))
}


getSelectSubSigma2ThetaEpsilon <- function(fullSigma2ThetaEpsilon,
                                           latentEtas, method = "qml") {
  if (method != "qml") return(NULL)
  select <- as.logical.matrix(fullSigma2ThetaEpsilon)
  select[TRUE] <- FALSE
  select[latentEtas, latentEtas] <- TRUE
  select
}


constructFullL2 <- function(colsU, etas, method = "qml") {
  if (method != "qml") return(NULL)

  if (is.null(colsU)) nCols <- length(etas)
  else nCols <- length(colsU)
  matrix(0, nrow = length(etas), ncol = nCols,
         dimnames = list(etas, colsU))
}


getSelectSubL2 <- function(fullL2, colsU, latentEtas, method = "qml") {
  if (method != "qml") return(NULL)
  select <- as.logical.matrix(fullL2)
  select[TRUE] <- FALSE
  select[latentEtas, !grepl("__FILL__ZERO__", colnames(select))] <- TRUE
  select
}


getScalingInds <- function(indsEtas, R, latentEtas, method = "qml") {
  if (method != "qml") return(NULL)
  allIndsEtas <- unique(unlist(indsEtas[latentEtas]))
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
                            dimnames = list(scalingInds, scalingInds))
  diag(subThetaEpsilon) <- NA
  subThetaEpsilon
}


getScalingLambdaY <- function(lambdaY, indsEtas, etas, method = "qml") {
  if (method != "qml") return(NULL)
  hasMultipleInds <- vapply(indsEtas, FUN = function(x) length(x) > 1,
                            FUN.VALUE = logical(1L))
  latentEtas <- names(indsEtas)[hasMultipleInds]
  indsLatentEtas <- unlist(indsEtas[latentEtas])
  lambdaY[indsLatentEtas, latentEtas]
}


sortXisConstructOmega <- function(xis, varsInts, etas, intTerms,
                                  method = "lms", double = FALSE) {
  listSortedXis  <- sortXis(xis = xis, varsInts = varsInts, etas = etas,
                            intTerms = intTerms, double = double)
  sortedXis    <- listSortedXis$sortedXis
  nonLinearXis <- listSortedXis$nonLinearXis

  omegaXiXi <- constructOmegaXiXi(xis = xis, etas = etas,
                                  sortedXis = sortedXis,
                                  nonLinearXis = nonLinearXis,
                                  varsInts = varsInts,
                                  intTerms = intTerms)
  omegaEtaXi <- constructOmegaEtaXi(xis = xis, etas = etas,
                                    sortedXis = sortedXis,
                                    nonLinearXis = nonLinearXis,
                                    varsInts = varsInts,
                                    intTerms = intTerms)

  list(sortedXis = sortedXis, omegaXiXi = omegaXiXi,
       omegaEtaXi = omegaEtaXi, k = length(nonLinearXis))
}


sortXis <- function(xis, varsInts, etas, intTerms, double) {
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

    stopif(length(interaction) > 2, "Only interactions between two variables are allowed")
    stopif(all(interaction %in% etas), "Interactions between two endogenous ",
           "variables are not allowed, see \nvignette(\"interaction_two_etas\", \"modsem\")")

    choice <- unique(interaction[which(!interaction %in% etas)])
    if (length(choice) > 1 && !double) {
      freq <- freqInIntTerms[choice, "freq"]
      choice <- choice[whichIsMax(freq)][[1]] # pick first if both are equal
    }

    nonLinearXis <- c(nonLinearXis, choice)
  }

  linearXis <- xis[!xis %in% nonLinearXis]

  list(linearXis = linearXis, sortedXis = c(nonLinearXis, linearXis),
       nonLinearXis = nonLinearXis)
}


constructOmegaEtaXi <- function(xis, etas, sortedXis, nonLinearXis,
                                varsInts, intTerms) {
  omega      <- NULL
  labelOmega <- NULL

  for (eta in etas) {
    subOmega <- matrix(0, nrow = length(xis), ncol = length(etas),
                            dimnames = list(sortedXis, etas))
    subLabelOmega <- as.character.matrix(subOmega, empty = TRUE)

    for (row in varsInts[intTerms$lhs == eta]) {
      if (!any(row %in% etas) || all(row %in% etas)) next

      whichXi  <- which(!row %in% etas)
      whichEta <- which(row %in% etas)

      subOmega[row[[whichXi]], row[[whichEta]]] <-
        getFreeOrConstIntTerms(row, eta, intTerms)
      subLabelOmega[row[[whichXi]], row[[whichEta]]] <-
        getLabelIntTerms(row, eta, intTerms)
    }

    omega      <- rbind(omega, labelRowsOmega(subOmega, eta = eta))
    labelOmega <- rbind(labelOmega, labelRowsOmega(subLabelOmega, eta = eta))
  }
  list(numeric = omega, label = labelOmega)
}


constructOmegaXiXi <- function(xis, etas, sortedXis, nonLinearXis,
                               varsInts, intTerms) {
  omega      <- NULL
  labelOmega <- NULL
  for (eta in etas) {
    subOmega <- matrix(0, nrow = length(sortedXis), ncol = length(sortedXis),
                       dimnames = list(sortedXis, sortedXis))
    subLabelOmega <- as.character.matrix(subOmega, empty = TRUE)

    for (row in varsInts[intTerms$lhs == eta]) {
       if (!all(row %in% sortedXis)) next

       whichRow <- which(row %in% nonLinearXis)[[1]] # if quadratic term pick first
       whichCol <- ifelse(whichRow == 1, 2, 1)

       subOmega[row[[whichRow]], row[[whichCol]]] <-
         getFreeOrConstIntTerms(row, eta, intTerms)
       subLabelOmega[row[[whichRow]], row[[whichCol]]] <-
         getLabelIntTerms(row, eta, intTerms)
    }

    omega      <- rbind(omega, labelRowsOmega(subOmega, eta = eta))
    labelOmega <- rbind(labelOmega, labelRowsOmega(subLabelOmega, eta = eta))
  }

  list(numeric = omega, label = labelOmega)
}


labelRowsOmega <- function(X, eta) {
  rownames(X) <- paste0(eta, "~", rownames(X))
  X
}
