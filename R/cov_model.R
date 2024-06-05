# functions for computing the constrained covariance matrix, 
# based on causal relationships. This can make the lms method more flexible, 
# as you can split the model into a non-linear, and linear part. allowing 
# you to use (normally distributed) endogenous variables as non-normal
# as of now the mean-structure is excluded
covModel <- function(syntax, method = "lms") {
  if (is.null(syntax)) return(list(matrices = NULL, freeParams = 0, 
                                   theta = NULL, syntax = syntax))
  parTable <- modsemify(syntax)
  structExprs <- parTable[parTable$op == "~" & 
                          parTable$rhs != "1", ]

  # endogenous variables (etas)
  etas <- getSortedEtas(structExprs)
  numEtas <- length(etas)
  if (numEtas == 0) stop2("No etas in model")

  xis <- parTable[parTable$op == "~" &
                  !parTable$rhs %in% etas & 
                  parTable$rhs != "1", "rhs"] |> unique()
  numXis <- length(xis)
  if (numXis == 0) stop2("No xis in model")

  gammaXi <- matrix(0, nrow = numEtas, ncol = numXis,
                  dimnames = list(etas, xis))
  exprsGammaXi <- structExprs[structExprs$lhs %in% etas & 
                              !structExprs$rhs %in% etas, ] 

  if (nrow(exprsGammaXi) > 0) {
    apply(exprsGammaXi, MARGIN = 1, FUN = function(row) 
          gammaXi[row[["lhs"]], row[["rhs"]]] <<- NA)
  }

  gammaEta <- matrix(0, nrow = numEtas, ncol = numEtas,
                  dimnames = list(etas, etas))
  exprsGammaEta <- structExprs[structExprs$lhs %in% etas & 
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
  A <- diag(numXis)
  rownames(A) <- colnames(A) <- xis
  phi <- A
  if (method == "lms") A[lower.tri(A, diag = TRUE)] <- NA
  else if (method == "qml") phi[lower.tri(phi, diag = TRUE)] <- NA

  matrices <- list(gammaXi = gammaXi,
                   gammaEta = gammaEta,
                   A = A,
                   psi = psi,
                   phi = phi) 

  model <- list(info =
                list(etas = etas,
                     numEtas = numEtas,
                     xis = xis,
                     numXis = numXis),
                matrices = matrices,
                syntax = syntax,
                parTable = parTable, 
                freeParams = countFreeCovModel(matrices),
                theta = createParamVectorCovModel(matrices))
   
  model
}


fillCovModel <- function(covModel, covTheta, fillPhi = FALSE, method = "lms") {
  matrices <- covModel$matrices 
  if (is.null(matrices)) return(NULL)
  if (is.null(names(covTheta))) names(covTheta) <- names(covModel$theta)
  matrices <- covModel$matrices
  matrices$psi[is.na(matrices$psi)] <-
    fetch(covTheta, "^psi[0-9]*$")
  matrices$gammaEta[is.na(matrices$gammaEta)] <-
    fetch(covTheta, "gammaEta[0-9]*$")
  matrices$gammaXi[is.na(matrices$gammaXi)] <-
    fetch(covTheta, "gammaXi[0-9]*$")
  if (method == "lms") {
    matrices$A[is.na(matrices$A)] <-
      fetch(covTheta, "^A[0-9]*$")
  } else if (method == "qml") {
    matrices$phi <- fillSymmetric(matrices$phi, fetch(covTheta, "^phi[0-9]*$"))
  }
  if (fillPhi) matrices$phi <- matrices$A %*% t(matrices$A)
  covModel$matrices <- matrices 
  covModel
}


createParamVectorCovModel <- function(matrices, start = NULL) {
  set.seed(123)
  phi <- as.vector(matrices$phi)
  A <- as.vector(matrices$A)
  psi <- as.vector(matrices$psi)
  alpha <- as.vector(matrices$alpha)
  gammaXi <- as.vector(matrices$gammaXi)
  gammaEta <- as.vector(matrices$gammaEta)
  theta <- c("phi" = phi,
             "A" = A,
             "psi" = psi,
             "gammaXi" = gammaXi,
             "gammaEta" = gammaEta)
  theta <- theta[is.na(theta)]
  if (is.null(start)) {
   theta <- vapply(theta, FUN.VALUE = vector("numeric", 1L),
                   FUN = function(x) stats::runif(1))
  }
  theta
}


countFreeCovModel <- function(matrices) {
  vapply(matrices, FUN.VALUE = integer(1L), 
         FUN = function(x) sum(is.na(x))) |> sum()
}


expectedCovModel <- function(model, method = "lms", sortedXis) {
  gammaXi <- model$matrices$gammaXi
  gammaEta <- model$matrices$gammaEta

  if (method == "lms") {
    A <- model$matrices$A
    phi <- A %*% t(A)
  } else if (method == "qml") {
    phi <- model$matrices$phi
  }
  psi <- model$matrices$psi
  
  Binv <- solve(diag(nrow(gammaEta)) - gammaEta)
  covEtaEta <- Binv %*% (gammaXi %*% phi %*% t(gammaXi) + psi) %*% t(Binv)
  covEtaXi <- Binv %*% gammaXi %*% phi
  sigma <- rbind(cbind(covEtaEta, covEtaXi),
                 cbind(t(covEtaXi), phi))
  sigma <- sigma[sortedXis, sortedXis]

  if (method == "lms") return(t(chol(sigma)))
  sigma
}


covModelToParTable <- function(model, method = "lms") {
  matricesEst <- model$covModel$matrices
  matricesSE <- model$covModelSE$matrices
  matricesNA <- model$covModelNA$matrices
  if (is.null(matricesEst) || is.null(matricesNA) ||
      is.null(matricesSE)) return(NULL)
  etas <- model$info$etas
  numXis <- model$info$numXis
  parTable <- NULL

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
