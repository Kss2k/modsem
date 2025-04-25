
#' Estimate baseline model for the LMS and QML approach
#'
#' @param object An object of class `modsem_da`.
#' @param warn_no_interaction Logical. If `TRUE`, a warning is issued if no interaction terms are found in the model.
#' @param ... Additional arguments passed to the `modsem_da` function, overriding
#'   the arguments in the original model.
#' @description Estimates a baseline model (H0) from a given model (H1) and compares the fit of both models.
#' The baseline model is estimated by removing all interaction terms from the model.
#' @rdname estimate_h0
#' @export
#' @examples
#' \dontrun{
#' m1 <- "
#'  # Outer Model
#'  X =~ x1 + x2 + x3
#'  Y =~ y1 + y2 + y3
#'  Z =~ z1 + z2 + z3
#'
#'  # Inner model
#'  Y ~ X + Z + X:Z
#' "
#'
#' est_h1 <- modsem(m1, oneInt, "lms")
#' est_h0 <- estimate_h0(est_h1, calc.se=FALSE) # std.errors are not needed
#' compare_fit(est_h0, est_h1)
#' }
#' @export
estimate_h0 <- function(object, warn_no_interaction = TRUE, ...) {
  argList    <- object$args
  parTable   <- object$originalParTable
  data       <- object$data
  method     <- object$method
  cov.syntax <- object$model$covModel$syntax

  newArgList <- list(...) 
  newArgNames <- intersect(names(argList), names(newArgList))
  argList[newArgNames] <- newArgList[newArgNames] 

  tryCatch({
    strippedParTable <- removeUnknownLabels(parTable[!grepl(":", parTable$rhs), ])
    
    if (NROW(strippedParTable) == NROW(parTable)) {
      warnif(warn_no_interaction, "No interaction terms found in the model. ",
             "The baseline model is identical to the original model, ",
             "and won't be estimated!")
      return(NULL)
    } 
    
    syntax <- parTableToSyntax(strippedParTable)
    argList <- c(
        list(model.syntax = syntax, data = data, method = method,
             cov.syntax = cov.syntax), argList
    )
    
    do.call(modsem_da, args = argList)
  },
  error = function(e) {
    warning2("Null model could not be estimated. ",
             "Error message: ", e$message)
    NULL
  })
}



