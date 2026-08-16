simulateDataParTable <- function(parTable, N, colsOVs = NULL, colsLVs = NULL, seed = NULL) {
  if (!is.null(seed) && exists(".Random.seed")) .Random.seed.orig <- .Random.seed
  else                                          .Random.seed.orig <- NULL

  on.exit({
    if (!is.null(.Random.seed.orig)) .Random.seed <<- .Random.seed.orig
  })

  if (!is.null(seed))
    set.seed(seed)

  parTable <- addMissingGroups(parTable)
  groups   <- getGroupsParTable(parTable)

  LV <- vector("list", length = length(groups))
  OV <- vector("list", length = length(groups))

  N.g <- floor(N / length(groups))

  for (g in groups) {
    simg <- simulateDataParTableGroup(
      parTable = parTable[parTable$group == g, , drop = FALSE],
      N        = N.g,
      colsOVs  = colsOVs,
      colsLVs  = colsLVs
    )

    LV[[g]] <- simg$lV
    OV[[g]] <- simg$oV
  }

  list(LV = LV, OV = OV)
}


simulatedGroupsToDf <- function(sim, type = "OV") {
  simt <- sim[[toupper(type)]]

  do.call(rbind, lapply(seq_along(simt), FUN = \(g)
                        cbind(as.data.frame(simt), data.frame(group = g))))
}


simulateDataParTableGroup <- function(parTable, N, colsOVs = NULL, colsLVs = NULL) {
  # Output columns
  if (is.null(colsOVs)) colsOVs <- getOVs(parTable)
  if (is.null(colsLVs)) colsLVs <- getLVs(parTable)

  # Here we simplify the model by converting all =~ to ~
  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  # Thus we can treat indicators the same as endogenous variables
  idx <- which(op == "=~")
  parTable[idx, "lhs"] <- rhs[idx]
  parTable[idx, "op"]  <- "~"
  parTable[idx, "rhs"] <- lhs[idx]

  # endogenous variables (etas)
  etas <- getSortedEtas(parTable, isLV = FALSE, checkAny = TRUE)
  xis  <- getXis(parTable, checkAny = TRUE, isLV = FALSE)
  vars <- union(xis, etas)

  # sample disturbances
  Psi <- matrix(
    0, nrow = length(vars), ncol = length(vars),
    dimnames = list(vars, vars)
  )

  psi.idx <- which(
    parTable$lhs %in% vars & parTable$op == "~~" & parTable$rhs %in% vars
  )

  for (i in psi.idx) {
    lhs <- parTable[i, "lhs"]
    rhs <- parTable[i, "rhs"]
    est <- parTable[i, "est"]
    Psi[lhs, rhs] <- Psi[rhs, lhs] <- est
  }

  alpha <- getIntercepts(vars, parTable = parTable)

  # Keep a copy of Zeta
  Zeta <- as.matrix(mvtnorm::rmvnorm(n = N, mean = alpha, sigma = Psi))
  colnames(Zeta) <- vars
  Eta <- Zeta

  for (eta in etas) {
    idx <- which(parTable$lhs == eta & parTable$op == "~")

    y <- Eta[,eta, drop = TRUE]
    for (i in idx) {
      row <- parTable[i, , drop = FALSE]
      gamma <- row$est
      x <- getValuesFromRhs(rhs = row$rhs, mat = Eta)
      y <-  y + gamma * x
    }

    Eta[,eta] <- y
  }

  dataOVs <- Eta[,colsOVs, drop = FALSE]
  dataLVs <- Eta[,colsLVs, drop = FALSE]

  list(oV = dataOVs, lV = dataLVs)
}


getValuesFromRhs <- function(rhs, mat) {
  mod_stopif(!length(rhs), "rhs is of length zero!")

  if (!grepl(":", rhs))
    return(mat[,rhs, drop = TRUE])

  elems <- stringr::str_split_1(rhs, pattern = ":")

  if (length(elems) == 1L) # should have been caught earlier, but just in case...
    return(mat[,elems, drop = TRUE])

  mod_stopif(!length(elems),
    "zero elements in product:", rhs
  )

  mod_stopif(!all(elems %in% colnames(mat)),
    "Missing values for elements in product!",
    paste0(setdiff(elems, colnames(mat)), collapse = ", ")
  )

  xz <- rep(1, NROW(mat))

  for (x in elems)
    xz <- xz * mat[,x, drop = TRUE]

  xz
}
