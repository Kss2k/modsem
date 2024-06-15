simulateDataParTable <- function(parTable, N) {
  # endogenous variables (etas)model
  etas <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  numEtas <- length(etas)
  
  indsEtas <- getIndsLVs(parTable, etas)
  numIndsEtas <- vapply(indsEtas, FUN.VALUE = vector("integer", 1L),
                        FUN = length)
  allIndsEtas <- unlist(indsEtas)
  numAllIndsEtas <- length(allIndsEtas)
  
  # exogenouts variables (xis) and interaction terms 
  xis <- getXis(parTable, checkAny = TRUE)
  numXis <- length(xis)

  indsXis <- getIndsLVs(parTable, xis)
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis <- unlist(indsXis)
  numAllIndsXis <- length(allIndsXis)
  
  # interaction terms
  intTerms <- getIntTerms(parTable)
  intTermRows <- getIntTermRows(parTable)
  varsIntTerms <- getVarsInts(intTermRows, removeColonNames = FALSE)
  if (any(vapply(varsIntTerms, FUN.VALUE = numeric(1L), FUN = length) > 2)) {
    stop2("Cannot simulate data for interaction effects with more than two ", 
          "components, yet")
  }
  # simulate data for xis 
  phi <- getEstPhiParTable(parTable)

  dataLVs <- as.matrix(mvtnorm::rmvnorm(n = N, mean = rep(0, numXis), sigma = phi))
  colnames(dataLVs) <- xis

  subVarsIntTerms <- varsIntTerms
  for (eta in etas) {
    toBuildXZ <- vapply(subVarsIntTerms, FUN.VALUE = logical(1L),
                        FUN = function(x) all(x %in% colnames(dataLVs)))
    XZ <- mutliplyPairs(dataLVs, XZ = subVarsIntTerms[toBuildXZ])
    subVarsIntTerms <- subVarsIntTerms[!toBuildXZ]
    dataLVs <- cbind(dataLVs, XZ)
   
    structExprsEta <- parTable[parTable$lhs == eta & 
                               parTable$op == "~" & 
                               parTable$rhs != "1", , drop = FALSE]
    alpha <- parTable[parTable$lhs == eta & 
                      parTable$op == "~" & 
                      parTable$rhs == "1", "est"]
    if (NROW(alpha) == 0) alpha <- 0

    y <- rep(alpha, length = N)
    for (i in seq_len(NROW(structExprsEta))) {
      row <- structExprsEta[i, , drop = FALSE]
      y <-  y + row$est * dataLVs[ , row$rhs]
    }

    psi <- parTable[parTable$lhs == parTable$rhs &
                    parTable$lhs == eta & parTable$op == "~~", 
                    "est"]
    y <- y + stats::rnorm(N, mean = 0, sd = sqrt(psi))
    dataLVs <- cbind(dataLVs, matrix(y, nrow = N, dimnames = list(NULL, eta)))   
  }
    
  dataXZs <- dataLVs[, intTerms]
  dataLVs <- dataLVs[, c(xis, etas)]
  dataOVs <- matrix(0, nrow = N, ncol = numAllIndsXis + numAllIndsEtas,
                    dimnames = list(NULL, c(allIndsXis, allIndsEtas)))
  indsLVs <- c(indsXis, indsEtas)
  interceptVector <- rep(1, N)
  for (lV in c(xis, etas)) {
    inds <- indsLVs[[lV]] 
    tau <- parTable[parTable$lhs %in% inds & 
                    parTable$op == "~" & 
                    parTable$rhs == "1", ] |> 
      sortTau(inds)
    lambda <- parTable[parTable$lhs == lV & 
                       parTable$op == "=~" & 
                       parTable$rhs %in% inds, ] |>
      sortLambda(inds)
    residuals <- getResidualsInds(parTable, lambda$rhs)                  
    dataOVs[, inds] <- 
      interceptVector %*% t(tau$est) + 
      dataLVs[, lV] %*% t(lambda$est) + 
      mvtnorm::rmvnorm(N, rep(0, NROW(lambda)), residuals)
  }
  list(oV = dataOVs, lV = dataLVs)
}


getEstPhiParTable <- function(parTable) {
  xis <- getXis(parTable, checkAny = TRUE)
  phi <- matrix(0, nrow = length(xis), ncol = length(xis), 
                dimnames = list(xis, xis))
  vcovExpres <- parTable[parTable$lhs %in% xis &
                         parTable$op == "~~" & 
                         parTable$rhs %in% xis, ]
  for (i in seq_len(nrow(vcovExpres))) {
    lhs <- vcovExpres[i, "lhs"]
    rhs <- vcovExpres[i, "rhs"]
    est <- vcovExpres[i, "est"]
    phi[lhs, rhs] <- phi[rhs, lhs] <- est
  }
  phi
}


getResidualsInds <- function(parTable, inds) {
  vcovRes <- matrix(0, nrow = length(inds), ncol = length(inds), 
                dimnames = list(inds, inds))
  vcovExpres <- parTable[parTable$lhs %in% inds &
                         parTable$op == "~~" & 
                         parTable$rhs %in% inds, ]
  for (i in seq_len(nrow(vcovExpres))) {
    lhs <- vcovExpres[i, "lhs"]
    rhs <- vcovExpres[i, "rhs"]
    est <- vcovExpres[i, "est"]
    vcovRes[lhs, rhs] <- vcovRes[rhs, lhs] <- est
  }
  vcovRes
}


mutliplyPairs <- function(X, XZ) {
  if (!is.list(XZ)) stop("Expected xz to be a list: ", XZ)
  prods <- matrix(0, nrow = NROW(X), ncol = length(XZ), 
                  dimnames = list(NULL, names(XZ)))
  for (i in seq_len(length(XZ))) {
    col <- names(XZ)[[i]]
    xz <- XZ[[i]]
    prods[, col] <- X[ , xz[[1]]] * X[ , xz[[2]]]
  }
  prods
}


sortLambda <- function(lambda, inds) {
  lambdaSorted <- NULL 
  for (ind in inds) {
    lambdaSorted <- rbind(lambdaSorted, lambda[lambda$rhs == ind, ])
  }
  lambdaSorted
}


sortTau <- function(tau, inds) {
  tauSorted <- NULL 
  for (ind in inds) {
    tauSorted <- rbind(tauSorted, tau[tau$lhs == ind, ])
  }
  tauSorted
}
