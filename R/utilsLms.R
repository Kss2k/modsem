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


fetch <- function(x, pattern = ".*") {
  x[grepl(pattern, names(x))]
}


stripMatrices <- function(matrices, fill = -1) {
  lapply(matrices, function(mat) {
    mat[!is.na(mat)] <- fill
    mat
  })
}


removeInteractions <- function(model) {
  model$matrices$OmegaEtaXi[TRUE] <- 0 
  model$matrices$OmegaXiXi[TRUE] <- 0
  model
}


# Faster version of mvtnorm::dmvnorm() given that sigma is positive 
dMvn <- function(X, mean, sigma, log = FALSE) {
  sigma <- round(sigma, 12)
  return(tryCatch(mvnfast::dmvn(X, mean, sigma, log, ncores = 2),
                  error = function(e) mvtnorm::dmvnorm(X, mean, sigma, log)))
}


calcSE <- function(negHessian, names = NULL) {
  stdError <- tryCatch({
      sqrt(diag(solve(negHessian)))
    }, error=function(e) {
      NA
    }, warning=function(w) {
       if (grepl("NaN", conditionMessage(w))) {
         suppressWarnings(sqrt(diag(solve(negHessian))))
      } else{
         sqrt(diag(solve(negHessian)))
      }
    })
  if (all(is.na(stdError))) 
    warning("SE's could not be computed, negative Hessian is singular.")
  if (any(is.nan(stdError))) 
    warning("SE's for some coefficients could not be computed.") 
  
  if (!is.null(names)) names(stdError) <- names
  stdError
}
