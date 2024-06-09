createThetaLabel <- function(labelMatrices, labelMatricesCov, 
                             constrExprs, start = NULL) {
  matrices <- c(labelMatrices, labelMatricesCov)  
  labels <- lapply(matrices, FUN = function(x) {
    select <- apply(x, MARGIN = 2, FUN = function(z) 
                    !canBeNumeric(z, includeNA = TRUE)) 
    as.vector(x[select])
  }) |> unlist() |> unique()
  # remove fixed params
  if (!is.null(constrExprs)) {
    labels <- labels[!labels %in% constrExprs$fixedParams]
  }

  if (is.null(start)) {
    start <- vapply(labels, FUN.VALUE = vector("numeric", 1L),
                    FUN = function(x) stats::runif(1))
  }

  names(start) <- labels
  start
}


calcThetaLabel <- function(theta, constrExprs) {
  if (is.null(constrExprs)) return(theta) 
  Theta <- c(as.list(theta), constrExprs$thetaFixedParams) 

  for (i in seq_along(constrExprs$fixedParams)) {
    Theta[constrExprs$fixedParams[[i]]] <- 
      eval(constrExprs$exprs[[i]], envir = Theta)
  }  

  unlist(Theta)
}


fillMatricesLabels <- function(matrices, labelMatrices, thetaLabel, 
                               constraintExprs = NULL) {
  if (is.null(thetaLabel)) return(matrices)
  labels <- names(thetaLabel)
  purrr::map2(.x = matrices, .y = labelMatrices, .f = function(M, L) {
                lapply(labels, FUN = function(label) {
                         M[L == label] <<- thetaLabel[[label]]
                })
                M
  })
}
