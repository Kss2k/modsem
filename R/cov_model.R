# functions for computing the constrained covariance matrix,
# based on causal relationships. This can make the lms method more flexible,
# as you can split the model into a non-linear, and linear part. allowing
# you to use (normally distributed) endogenous variables as non-normal
# as of now the mean-structure is excluded
covModel <- function(syntax = NULL,
                     method = "lms",
                     parTable = NULL,
                     xis.main = NULL,
                     parTable.main = NULL,
                     orthogonal.x = FALSE,
                     orthogonal.y = FALSE) {
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

    parTable <- parTable.main[
      (parTable.main$lhs %in% xis.main &
       parTable.main$rhs %in% xis.main &
       parTable.main$op == "~~") |
       parTable.main$op %in% c(":=", "=="), , drop = FALSE
    ]
  }

  # Gamma, i.e., [[gammaEta, gammaXi], [0, 0]] over c(etas, xis).
  # The xi rows are all zero, as xis are exogenous by construction
  listGamma <- constructGamma(c(etas, xis), c(etas, xis), parTable = parTable)
  gamma <- listGamma$numeric
  labelGamma <- listGamma$label

  # covariance matrices
  # We need the full partable to look for all relevant covariances
  listPsi <- constructPsi(
    c(etas, xis),
    parTable = parTable.full,
    orthogonal.y = orthogonal.y,
    orthogonal.x = orthogonal.x,
    has.exo = TRUE
  )

  psi <- listPsi$numeric
  labelPsi <- listPsi$label

  matrices <- list(gamma = gamma, psi = psi)

  labelMatrices <- list(gamma = labelGamma, psi = labelPsi)

  model <- list(
    info = list(
      etas = etas,
      numEtas = numEtas,
      xis = xis,
      numXis = numXis,
      is.simple = isSimple
    ),
    matrices = matrices,
    labelMatrices = labelMatrices,
    syntax = syntax,
    parTable = parTable
  )

  model
}


countFreeCovModel <- function(matrices) {
  vapply(matrices, FUN.VALUE = integer(1L),
         FUN = function(x) sum(is.na(x))) |> sum()
}


expectedCovModel <- function(model, method = "lms", sortedXis) {
  gamma <- model$matrices$gamma
  psi <- model$matrices$psi

  # `gamma` is all-zero when `is.simple`, in which case this is just the identity
  Binv <- solve(diag(nrow(gamma)) - gamma)

  sigma <- Binv %*% psi %*% t(Binv)
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
  matricesEst   <- model$covModel$matrices
  matricesSE    <- model$covModelSE$matrices
  matricesNA    <- model$covModelNA$matrices
  matricesLabel <- model$covModel$labelMatrices

  if (is.null(matricesEst) || is.null(matricesNA)) return(NULL)
  if (is.null(matricesSE)) matricesSE <- matricesNA

  etas <- model$info$etas
  numXis <- model$info$numXis
  parTable <- NULL

  if (!model$covModel$info$is.simple) {
    # coefficients Structural Model
    newRows <- matrixToParTable(matricesNA$gamma,
                                matricesEst$gamma,
                                matricesSE$gamma,
                                matricesLabel$gamma,
                                op = "~",
                                rowsLhs = TRUE)
    parTable <- rbind(parTable, newRows)
  }

  newRows <- matrixToParTable(matricesNA$psi,
                              matricesEst$psi,
                              matricesSE$psi,
                              matricesLabel$psi,
                              op = "~~",
                              rowsLhs = FALSE,
                              symmetric = TRUE)
  parTable <- rbind(parTable, newRows)


  parTable <- lapplyDf(parTable, FUN = function(x) replace(x, x == -999, NA))
  # return
  parTable
}
