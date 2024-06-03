# Utils for lms approach
# Last updated: 29.05.2024


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
    sum() + model$covModel$freeParams
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
# there are some drawbacks to using mvnfast. In particular, 
# its a little less consistent
dMvn <- function(X, mean, sigma, log = FALSE) {
  sigma <- round(sigma, 12) # makes it more reliable
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


formatNumeric <- function(x, digits = 3) {
  if (is.numeric(x)) {
    format(round(x, digits), nsmall = digits, digits = digits, 
           trim = FALSE, justify = "right")
  } else {
    format(x, trim = FALSE, justify = "right") 
  }
}


oneWayTableToDataFrame <- function(table) {
  df <- as.data.frame(table)
  out <- data.frame(freq = df$Freq)
  rownames(out) <- df$Var1
  out
}


whichIsMax <- function(x) {
  which(x == max(x))
}
