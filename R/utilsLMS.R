sortData <- function(data, allIndsXis, allIndsEtas) {
  if (!all(c(allIndsXis, allIndsEtas) %in% colnames(data))) 
    stop("Missing Observed Variables in Data")
  as.matrix(data[c(allIndsXis, allIndsEtas)])
}



countFreeParams <- function(model) {
  matrices <- model$matrices[c("lambdaY", 
                               "lambdaX",
                               "tauY",
                               "tauX",
                               "thetaEpsilon",
                               "thetaDelta",
                               "gammaXi",
                               "gammaEta",
                               "omegaXiXi",
                               "omegaEtaXi",
                               "psi",
                               "A",
                               "phi",
                               "alpha")]
  vapply(matrices,
         FUN.VALUE = vector("integer", 1L),
         FUN = function(x) sum(is.na(x))) |>
    sum()
}



fillModelParams <- function(model, paramVec) {
  freeParamsMatrices <- vapply(model$matrices,
    FUN.VALUE = vector("integer", 1L),
    FUN = function(x) sum(is.na(x))
  )
  for (i in seq_along(model$matrices)) {
    if (freeParamsMatrices[[i]] == 0) {
      next
    }
    model$matrices[[i]][is.na(model$matrices[[i]])] <- paramVec[1:freeParamsMatrices[[i]]]
    paramVec <- paramVec[-(1:freeParamsMatrices[[i]])]
  }
  model
}



labelParamVec <- function(model, start) {
  paramVec <- c()
  modFilled <- fillModelParams(model, start)
  for (i in names(model$matrices)) {
    # Phi should just be set to I^n, since it must be positive definite
    if (i == "phi") {
      phi <- model$matrices[[i]]
      phi[lower.tri(phi)] <- 0
      diag(phi) <- 1
      labeledPhi <- createLabeledVec(phi, "phi")[as.vector(lower.tri(phi, diag = TRUE))]
      paramVec <- c(paramVec, labeledPhi)
    } else {
      currentMnotFilled <- model$matrices[[i]]
      currentMFilled <- modFilled$matrices[[i]]
      labeledParams <- createLabeledVec(currentMFilled, i)
      labeledParams <- labeledParams[as.vector(is.na(currentMnotFilled))]
      paramVec <- c(paramVec, labeledParams)
    }
  }
  paramVec
}


createLabeledVec <- function(x, label) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  labels <- paste0(label, seq_along(x))
  structure(x, names = labels)
}


lapplyNamed <- function(X, FUN, ..., names = X) {
  structure(lapply(X, FUN, ...),
    names = names
  )
}


# Faster version of mvtnorm::dmvnorm() given that sigma is positive 
# definite. In some cases sigma will switch between beign positive 
# definite and not positive definite. dMvn gives slightly different
# results than dmvnorm, where it causes the EM algorthim to jump 
# between two different local minima. To avoid this the algorithm 
# switches to precision = TRUE, only using dmvnorm. 
dMvn <- function(X, mean, sigma, log = FALSE, precision = FALSE) {
  if (precision) return(mvtnorm::dmvnorm(X, mean, sigma, log))
  sigma <- round(sigma, 12)
  return(tryCatch(mvnfast::dmvn(X, mean, sigma, log, ncores = 2),
            error = function(e) mvtnorm::dmvnorm(X, mean, sigma, log)))
}
