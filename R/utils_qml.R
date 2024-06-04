tr <- function(mat) sum(diag(mat))


traceOmegaXiXi <- function(omega, numEta, numXi) {
  lastRow <- 0 
  lastCol <- 0  
  trace <- numeric(numEta)
  for (i in seq_len(numEta)) {
    trace[[i]] <- tr(omega[seq_len(numXi) + (i - 1) * numXi, 
                           seq_len(numXi) + (i - 1) * numXi]) 
  }
  trace
}


diagPartitioned <- function(matrix, length) {
  out <- matrix(0, nrow = length * nrow(matrix), 
                ncol = length * ncol(matrix))
  nrows <- nrow(matrix)
  rows <- seq_len(nrows)
  ncols <- ncol(matrix)
  cols <- seq_len(ncols)
  for (i in seq_len(length)) {
    out[rows + (i - 1) * nrows, 
        cols + (i - 1) * ncols] <- matrix
  }
  out
}


repPartitionedRows <- function(matrix, length = 1) {
  if (length <= 1) return(matrix)
  out <- matrix 
  for (i in seq_len(length - 1)) {
    out <- rbind(out, matrix)
  }
  out
}


repPartitionedCols <- function(matrix, length = 1) {
  if (length <= 1) return(matrix)
  out <- matrix 
  for (i in seq_len(length - 1)) {
    out <- cbind(out, matrix)
  }
  out
}


#' @export
as.logical.matrix <- function(x, ...) {
  structure(x != 0, 
            dim = dim(x), 
            dimnames = dimnames(x))
}


# not finished yet
# calcStandardError.modsemQML <- function(object, ...) {
#   # not correct yet
#   H <- object$negHessian 
#   invH <- solve(H)
#   N <- object$object$info$N
#   gradient <- gradientLogLikQml(object$emptyModel, object$coefficients)
#   J <- outer(gradient, gradient)
#   Jstar <- (1 / N) * (invH %*% J %*% invH)
#   sqrt(diag(Jstar))
# }
