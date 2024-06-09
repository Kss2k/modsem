#' Extract parameterEstimates from an estimated model
#'
#' @param object An object of class `modsem_pi`, `modsem_lms_qml`, or `modsem_mplus`
#' @param ... Additional arguments passed to other functions
#' @export 
parameter_estimates <- function(object, ...) {
  UseMethod("parameter_estimates")
}


#' Extract or modify parTable from an estimated model with estimated variances of interaction terms
#'
#' @param object An object of class `modsem_lms_qml`,  `modsem_mplus`, 
#' or a parTable of class `data.frame`
#' @param ... Additional arguments passed to other functions
#' @export
var_interactions <- function(object, ...) {
  UseMethod("var_interactions")
}


#' @export
var_interactions.data.frame <- function(object, ...) {
  parTable <- object[c("lhs", "op", "rhs", "est")]
  intTerms <- unique(parTable[grepl(":", parTable$rhs), "rhs"])

  for (i in seq_len(length(intTerms))) {
    # interaction term = XY
    XY <- stringr::str_split_fixed(intTerms[[i]], ":", 2) 
    varX <- calcCovParTable(parTable, XY[[1]], XY[[1]])
    varY <- calcCovParTable(parTable, XY[[2]], XY[[2]])
    covXY <- calcCovParTable(parTable, XY[[1]], XY[[2]])
    newRow <- data.frame(lhs = intTerms[[i]],
                         op = "~~",
                         rhs = intTerms[[i]],
                         est = varX * varY + covXY ^ 2)
    parTable <- rbind(parTable, newRow)
  }
  parTable
}


#' Get standardized estimates
#'
#' @param object An object of class `modsem_lms_qml`,  `modsem_mplus`, 
#' or a parTable of class `data.frame`
#' @param ... Additional arguments passed to other functions
#' @export
standardized_estimates <- function(object, ...) {
  UseMethod("standardized_estimates")
}


#' @export 
standardized_estimates.data.frame <- function(object, ...) {
  parTable <- object[c("lhs", "op", "rhs", "est")]
  parTable$mod <- as.character(parTable$est)
  parTable <- parTable[c("lhs", "op", "rhs", "est")] 
  etas <- getEtas(parTable)
  #for (eta in etas) {
  #    
  #}
}

