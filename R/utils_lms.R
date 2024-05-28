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


diagPartitionedMat <- function(X, Y) {
  if (is.null(X)) return(Y) else if (is.null(Y)) return(X)
  structure(rbind(cbind(X, matrix(0, nrow = NROW(X), ncol = NCOL(Y))), 
               cbind(matrix(0, nrow = NROW(Y), ncol = NCOL(X)), Y)),
            dimnames = list(c(rownames(X), rownames(Y)), 
                            c(colnames(X), colnames(Y))))
}


format_numeric <- function(x, digits = 3) {
  if (is.numeric(x)) {
    format(round(x, digits), nsmall = digits, digits = digits, 
           trim = FALSE, justify = "right")
  } else {
    format(x, trim = FALSE, justify = "right") 
  }
}


getFirstXinY <- function(x, y) {
  # pick the value in x which appears first in y
  for (y_i in y) {
    if (y_i %in% x) return(unique(x[x == y_i]))
  }

  # if no match, just return the first elem in x
  x[[1]]
}
