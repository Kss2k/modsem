createLabelTheta <- function(labelMatrices, start = NULL) {
  labels <- lapply(labelMatrices, FUN = function(x) {
    select <- apply(x, MARGIN = 2, FUN = function(z) 
                    !canBeNumeric(z, includeNA = TRUE)) 
    as.vector(x[select])
  }) |> unlist() |> unique()

  if (is.null(start)) {
    start <- vapply(labels, FUN.VALUE = vector("numeric", 1L),
                    FUN = function(x) stats::runif(1))
  }

  names(start) <- labels
  start
}


fillMatricesLabels <- function(matrices, labelMatrices, labelTheta, 
                               constraintExprs = NULL) {
  labels <- names(labelTheta)
  purrr::map2(.x = matrices, .y = labelMatrices, .f = function(M, L) {
                lapply(labels, FUN = function(label) {
                         M[L == label] <<- labelTheta[[label]]
                })
                M
  })
}
