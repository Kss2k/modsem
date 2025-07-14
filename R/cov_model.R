# functions for computing the constrained covariance matrix,
# based on causal relationships. This can make the lms method more flexible,
# as you can split the model into a non-linear, and linear part. allowing
# you to use (normally distributed) endogenous variables as non-normal
# as of now the mean-structure is excluded
covModel <- function(syntax = NULL, method = "lms", parTable = NULL, 
                     xis.main = NULL, parTable.main = NULL) {
  if (is.null(parTable) && !is.null(syntax)) parTable <- modsemify(syntax)
  if (is.null(parTable)) {
    return(list(matrices = NULL, freeParams = 0, info = list(etas = NULL, xis = NULL),
                theta = NULL, syntax = NULL, parTable = NULL))
  }

  parTable.full <- rbind(parTable.main, parTable)

  if (NROW(parTable) && any(parTable$op == "~")) {
    etas     <- getSortedEtas(parTable, isLV = FALSE, checkAny = TRUE)
    numEtas  <- length(etas)
    xis      <- getXis(parTable, checkAny = TRUE, isLV = FALSE)
    xis      <- unique(c(xis, xis.main[!xis.main %in% etas]))
                       
    numXis   <- length(xis)
    isSimple <- FALSE
  } else {
    etas     <- character(0L)
    numEtas  <- 0L
    xis      <- xis.main
    numXis   <- length(xis)
    isSimple <- TRUE

    parTable <- parTable.main[(parTable.main$lhs %in% xis.main &
                              parTable.main$rhs %in% xis.main &
                              parTable.main$op == "~~") |
                              parTable.main$op %in% c(":=", "=="), , drop = FALSE]
  }

  # Gamma
  listGammaXi <- constructGamma(etas, xis, parTable = parTable)
  gammaXi <- listGammaXi$numeric
  labelGammaXi <- listGammaXi$label

  listGammaEta <- constructGamma(etas, etas, parTable = parTable)
  gammaEta <- listGammaEta$numeric
  labelGammaEta <- listGammaEta$label

  # covariance matrices
  listPsi <- constructPsi(etas, parTable = parTable.full) # we need to the full parTable
                                                          # to identify pure etas
  psi <- listPsi$numeric
  labelPsi <- listPsi$label

  listPhi <- constructPhi(xis, method = "qml", parTable = parTable) # no need to treat methods differently here...
  phi <- listPhi$numeric
  labelPhi <- listPhi$label

  matrices <- list(
    gammaXi = gammaXi,
    gammaEta = gammaEta,
    psi = psi,
    phi = phi)

  labelMatrices <- list(
    gammaXi = labelGammaXi,
    gammaEta = labelGammaEta,
    psi = labelPsi,
    phi = labelPhi)

  model <- list(info =
                list(etas = etas,
                     numEtas = numEtas,
                     xis = xis,
                     numXis = numXis,
                     is.simple = isSimple),
                matrices = matrices,
                labelMatrices = labelMatrices,
                syntax = syntax,
                parTable = parTable)

  model
}


countFreeCovModel <- function(matrices) {
  vapply(matrices, FUN.VALUE = integer(1L),
         FUN = function(x) sum(is.na(x))) |> sum()
}


expectedCovModel <- function(model, method = "lms", sortedXis) {
  gammaXi <- model$matrices$gammaXi
  gammaEta <- model$matrices$gammaEta

  phi <- model$matrices$phi
  psi <- model$matrices$psi

  if (!model$info$is.simple) {
    Binv <- solve(diag(nrow(gammaEta)) - gammaEta)
    covEtaEta <- Binv %*% (gammaXi %*% phi %*% t(gammaXi) + psi) %*% t(Binv)
    covEtaXi <- Binv %*% gammaXi %*% phi
    sigma <- rbind(cbind(covEtaEta, covEtaXi),
                   cbind(t(covEtaXi), phi))

  } else sigma <- phi

  sigma <- sigma[sortedXis, sortedXis]

  if (method == "lms") {
    sigma <- tryCatch(t(chol(sigma)),
                      error = function(e) {
                        sigma[TRUE] <- NaN
                        sigma
                      })
  }
  sigma
}


covModelToParTable <- function(model, method = "lms") {
  matricesEst <- model$covModel$matrices
  matricesSE <- model$covModelSE$matrices
  matricesNA <- model$covModelNA$matrices
  matricesLabel <- model$covModel$labelMatrices

  if (is.null(matricesEst) || is.null(matricesNA)) return(NULL)
  if (is.null(matricesSE)) matricesSE <- matricesNA

  etas <- model$info$etas
  numXis <- model$info$numXis
  parTable <- NULL

  if (!model$covModel$info$is.simple) {
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
    
    newRows <- matrixToParTable(matricesNA$psi,
                                matricesEst$psi,
                                matricesSE$psi,
                                matricesLabel$psi,
                                op = "~~",
                                rowsLhs = FALSE)
    parTable <- rbind(parTable, newRows)
  }

  newRows <- matrixToParTable(matricesNA$phi,
                              matricesEst$phi,
                              matricesSE$phi,
                              matricesLabel$phi,
                              op = "~~",
                              rowsLhs = FALSE,
                              symmetric = TRUE)
  parTable <- rbind(parTable, newRows)


  parTable <- lapplyDf(parTable, FUN = function(x) replace(x, x == -999, NA))
  # return
  parTable
}
