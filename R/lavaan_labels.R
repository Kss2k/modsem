combineLavLabels <- function(lavLabelsCov, lavLabelsMain, currentLabels) {
  lavLabels <- c(lavLabelsCov, lavLabelsMain)
  finalLabels <- currentLabels
  finalLabels[finalLabels %in% names(lavLabels)] <- 
    lavLabels[names(lavLabels) %in% finalLabels]
  finalLabels
}


createLavLabels <- function(matrices, subset, etas) {
  lambdaX      <- createLabelsMatrix(matrices$lambdaX, op = "~")
  lambdaY      <- createLabelsMatrix(matrices$lambdaY, op = "~")
  thetaDelta   <- createLabelsMatrix(matrices$thetaDelta, op = "~~")
  thetaEpsilon <- createLabelsMatrix(matrices$thetaEpsilon, op = "~~")
  phi          <- createLabelsMatrix(matrices$phi, op = "~~")
  A            <- createLabelsMatrix(matrices$A, op = "~~")
  psi          <- createLabelsMatrix(matrices$psi, op = "~~")
  tauX         <- createLabelsMatrix(matrices$tauX, op = "~", first = "rows")
  tauY         <- createLabelsMatrix(matrices$tauY, op = "~", first = "rows")
  alpha        <- createLabelsMatrix(matrices$alpha, op = "~", first = "rows")
  beta0        <- createLabelsMatrix(matrices$beta0, op = "~")
  gammaXi      <- createLabelsMatrix(matrices$gammaXi, op = "~", first = "rows")
  gammaEta     <- createLabelsMatrix(matrices$gammaEta, op = "~", first = "rows")
  omegaXiXi    <- createLabelsOmega(matrices$omegaXiXi, etas = etas)
  omegaEtaXi   <- createLabelsOmega(matrices$omegaEtaXi, etas = etas)

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

  phi      <- createLabelsMatrix(matrices$phi, op = "~~")
  A        <- createLabelsMatrix(matrices$A, op = "~~")
  psi      <- createLabelsMatrix(matrices$psi, op = "~~")
  gammaXi  <- createLabelsMatrix(matrices$gammaXi, op = "~", first = "rows")
  gammaEta <- createLabelsMatrix(matrices$gammaEta, op = "~", first = "rows")

  labels <- c("phi" = phi,
              "A" = A,
              "psi" = psi,
              "gammaXi" = gammaXi,
              "gammaEta" = gammaEta)

  labels[subset]
}


createLabelsMatrix <- function(X, op = "~", first = "cols") {
  labels <- character(0L)
  rows   <- rownames(X)
  cols   <- colnames(X)

  # this is ugly, but... we have to read by cols first
  getLabel <- switch(first, cols = function(col, row) paste0(col, op, row),
                     rows = function(col, row) paste0(row, op, col))
  
  for (i in seq_len(ncol(X))) for (j in seq_len(nrow(X))) {
    labels <- c(labels, getLabel(col = cols[[i]], row = rows[[j]])) 
  }

  labels
}


createLabelsOmega <- function(X, etas) {
  rows       <- rownames(X)
  cols       <- colnames(X)
  numEtas    <- length(etas) 
  subNrow    <- nrow(X) %/% numEtas
  subSeqRows <- seq_len(subNrow) 

  labels <- character(0L)
  # this is ugly, but... we have to read by cols first, but also assign the
  # correct eta -- I don't want to do a check in the middle to see if j is 
  # in the range of a specified eta
  for (i in seq_len(ncol(X))) {
    for (eta_i in seq_len(numEtas)) {
      for (j in subSeqRows + (eta_i - 1) * subNrow) {
        labels <- c(labels, paste0(etas[eta_i], "~", rows[[j]], ":", cols[[i]]))
      }
    }
  }

  labels
}
