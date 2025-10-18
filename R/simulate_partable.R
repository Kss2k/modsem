simulateDataParTable <- function(parTable, N, colsOVs = NULL, colsLVs = NULL) {
  parTable <- getMissingGroups(parTable)
  groups   <- getGroupsParTable(parTable)

  LV <- vector("list", length = length(groups))
  OV <- vector("list", length = length(groups))

  for (g in groups) {
    simg <- simulateDataParTableGroup(
      parTable = parTable[parTable$group == g, , drop = FALSE],
      N        = N,
      colsOVs  = colsOVs,
      colsLVs  = colsLVs
    )
    
    LV[[g]] <- simg$lV
    OV[[g]] <- simg$OV
  }

  list(LV = LV, OV = OV)
}


simulatedGroupsToDf <- function(sim, type = "OV") {
  simt <- sim[[toupper(type)]]

  do.call(rbind, lapply(seq_along(simt), FUN = \(g)
                        cbind(as.data.frame(simt), data.frame(group = g))))
}


simulateDataParTableGroup <- function(parTable, N, colsOVs = NULL, colsLVs = NULL) {
  # endogenous variables (etas)
  etas    <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  numEtas <- length(etas)

  indsEtas    <- getIndsLVs(parTable, etas)
  numIndsEtas <- vapply(indsEtas, FUN.VALUE = vector("integer", 1L),
                        FUN = length)
  allIndsEtas    <- unique(unlist(indsEtas))
  numAllIndsEtas <- length(allIndsEtas)

  # exogenouts variables (xis) and interaction terms
  xis    <- getXis(parTable, checkAny = TRUE)
  numXis <- length(xis)

  indsXis    <- getIndsLVs(parTable, xis)
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis    <- unique(unlist(indsXis))
  numAllIndsXis <- length(allIndsXis)

  # interaction terms
  intTerms <- getIntTerms(parTable)
  intTermRows <- getIntTermRows(parTable)
  varsIntTerms <- getVarsInts(intTermRows, removeColonNames = FALSE)

  stopif(any(vapply(varsIntTerms, FUN.VALUE = numeric(1L), FUN = length) > 2),
         "Cannot simulate data for interaction effects with more than two ",
         "components, yet")

  # simulate data for xis
  phi <- rmvnormParTable(parTable, type = "phi", N = N)
  psi <- rmvnormParTable(parTable, type = "psi", N = N)
  theta <- rmvnormParTable(parTable, type = "theta", N = N)

  dataLVs <- phi

  subVarsIntTerms <- varsIntTerms
  for (eta in etas) {
    toBuildXZ <- vapply(subVarsIntTerms, FUN.VALUE = logical(1L),
                        FUN = function(x) all(x %in% colnames(dataLVs)))
    XZ <- mutliplyPairs(dataLVs, XZ = subVarsIntTerms[toBuildXZ])
    subVarsIntTerms <- subVarsIntTerms[!toBuildXZ]
    dataLVs <- cbind(dataLVs, XZ)

    structExprsEta <- parTable[parTable$lhs == eta & parTable$op == "~", ,
                               drop = FALSE]
    alpha <- parTable[parTable$lhs == eta & parTable$op == "~1", "est"]
    if (NROW(alpha) == 0) alpha <- 0

    y <- rep(alpha, length = N)
    for (i in seq_len(NROW(structExprsEta))) {
      row <- structExprsEta[i, , drop = FALSE]
      y <-  y + row$est * dataLVs[ , row$rhs]
    }

    y <- y + psi[, eta]
    dataLVs <- cbind(dataLVs, matrix(y, nrow = N, dimnames = list(NULL, eta)))
  }

  dataXZs <- dataLVs[, intTerms]
  dataLVs <- dataLVs[, c(xis, etas)]
  dataOVs <- matrix(0, nrow = N, ncol = numAllIndsXis + numAllIndsEtas,
                    dimnames = list(NULL, c(allIndsXis, allIndsEtas)))
  indsLVs <- c(indsXis, indsEtas)
  interceptVector <- rep(1, N)

  for (lV in c(xis, etas)) {
    inds   <- indsLVs[[lV]]
    tau    <- getIntercepts(inds, parTable = parTable)
    alpha  <- getMean(lV, parTable = parTable)
    lambda <- getLambda(lV = lV, inds = inds, parTable = parTable)
    dataOVs[, inds] <-
      interceptVector %*% t(tau) +
      (alpha + dataLVs[, lV]) %*% t(lambda) +
      theta[, inds]
  }

  if (!is.null(colsOVs)) dataOVs <- dataOVs[ , colsOVs]
  if (!is.null(colsLVs)) dataLVs <- dataLVs[ , colsLVs]

  list(oV = dataOVs, lV = dataLVs)
}


rmvnormParTable <- function(parTable, type = "phi", N) {
  vars <- switch(type,
                 phi = getXis(parTable, checkAny = TRUE),
                 psi = getSortedEtas(parTable, checkAny = TRUE, isLV = TRUE),
                 theta = getInds(parTable))

  vcov <- matrix(0, nrow = length(vars), ncol = length(vars),
                dimnames = list(vars, vars))

  vcovExpres <- parTable[parTable$lhs %in% vars &
                         parTable$op == "~~" &
                         parTable$rhs %in% vars, ]

  for (i in seq_len(nrow(vcovExpres))) {
    lhs <- vcovExpres[i, "lhs"]
    rhs <- vcovExpres[i, "rhs"]
    est <- vcovExpres[i, "est"]
    vcov[lhs, rhs] <- vcov[rhs, lhs] <- est
  }

  if (type == "phi") beta0 <- getIntercepts(vars, parTable = parTable)
  else beta0 <- rep(0, length(vars))

  X <- as.matrix(mvtnorm::rmvnorm(n = N, mean = beta0, sigma = vcov))
  colnames(X) <- vars
  X
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


getLambda <- function(lV, inds, parTable) {
  lambda <- parTable[parTable$lhs == lV &
                     parTable$op == "=~" &
                     parTable$rhs %in% inds, ]
  out <- lambda$est
  names(out) <- lambda$rhs
  out[inds]
}
