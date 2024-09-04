createLavLabels <- function(matrices, subset) {
  lambdaX <- createLabelsMatrix(matrices$lambdaX, op = "~")
  lambdaY <- createLabelsMatrix(matrices$lambdaY, op = "~")
  thetaDelta <- createLabelsMatrix(matrices$thetaDelta, op = "~~")
  thetaEpsilon <- createLabelsMatrix(matrices$thetaEpsilon, op = "~~")
  phi <- createLabelsMatrix(matrices$phi, op = "~~")
  A <- createLabelsMatrix(matrices$A, op = "~~")
  psi <- createLabelsMatrix(matrices$psi, op = "~~")
  tauX <- createLabelsMatrix(matrices$tauX, op = "~", first = "rows")
  tauY <- createLabelsMatrix(matrices$tauY, op = "~", first = "rows")
  alpha <- createLabelsMatrix(matrices$alpha, op = "~", first = "rows")
  beta0 <- createLabelsMatrix(matrices$beta0, op = "~")
  gammaXi <- createLabelsMatrix(matrices$gammaXi, op = "~", first = "rows")
  gammaEta <- createLabelsMatrix(matrices$gammaEta, op = "~", first = "rows")
  omegaXiXi <- createLabelsOmega(matrices$omegaXiXi, 
                                 etas = rownames(matrices$gammaEta))
  omegaEtaXi <- createLabelsOmega(matrices$omegaEtaXi, 
                                  etas = rownames(matrices$gammaEta))
  labels <- c("lambdaX" = lambdaX,
              "lambdaY" = lambdaY,
              "tauX" = tauX,
              "tauY" = tauY,
              "thetaDelta" = thetaDelta,
              "thetaEpsilon" = thetaEpsilon,
              "phi" = phi,
              "A" = A,
              "psi" = psi,
              "alpha" = alpha,
              "beta0" = beta0,
              "gammaXi" = gammaXi,
              "gammaEta" = gammaEta,
              "omegaXiXi" = omegaXiXi,
              "omegaEtaXi" = omegaEtaXi)
  labels[subset]
}


createLavLabelsCov <- function(matrices, subset) {
  if (is.null(matrices)) return(NULL)

  phi <- createLabelsMatrix(matrices$phi, op = "~~")
  A <- createLabelsMatrix(matrices$A, op = "~~")
  psi <- createLabelsMatrix(matrices$psi, op = "~~")
  gammaXi <- createLabelsMatrix(matrices$gammaXi, op = "~", 
                                first = "rows")
  gammaEta <- createLabelsMatrix(matrices$gammaEta, op = "~",
                                 first = "rows")
  labels <- c("phi" = phi,
                "A" = A,
                "psi" = psi,
                "gammaXi" = gammaXi,
                "gammaEta" = gammaEta)
  labels[subset]
}


createLabelsMatrix <- function(X, op = "~", first = "cols") {
  labels <- character(0L)
  rows <- rownames(X)
  cols <- colnames(X)

  if (first == "cols") {
    for (i in seq_len(ncol(X))) {
      for (j in seq_len(nrow(X))) {
        labels <- c(labels, paste0(cols[[i]], op, rows[[j]])) 
      }
    } 
  } else if (first == "rows") {
    for (i in seq_len(ncol(X))) {
      for (j in seq_len(nrow(X))) {
        labels <- c(labels, paste0(rows[[j]], op, cols[[i]])) 
      }
    } 
  }
  labels
}


createLabelsOmega <- function(X, etas) {
  numEtas <- length(etas) 
  subNrow <- nrow(X) %/% numEtas
  subSeqRows <- seq_len(subNrow) 
  rows <- rownames(X)
  cols <- colnames(X)
  labels <- character(0L)

  for (i in seq_len(ncol(X))) {
    for (eta_i in seq_len(numEtas)) {
      for (j in subSeqRows + (eta_i - 1) * subNrow) {
        labels <- c(labels, paste0(etas[eta_i], "~", rows[[j]], ":", cols[[i]]))
      }
    }
  }

  labels
}
