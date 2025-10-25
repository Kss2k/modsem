
#' Estimate baseline model for \code{modsem} models
#'
#' @param object An object of class \code{\link{modsem_da}} or \code{\link{modsem_pi}}.
#' @param warn_no_interaction Logical. If `TRUE`, a warning is issued if no interaction terms are found in the model.
#' @param ... Additional arguments passed to the `modsem_da` function, overriding
#'   the arguments in the original model.
#' @description Estimates a baseline model (H0) from a given model (H1).
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
#' # LMS approach
#' est_h1 <- modsem(m1, oneInt, "lms")
#' est_h0 <- estimate_h0(est_h1, calc.se=FALSE) # std.errors are not needed
#' compare_fit(est_h1 = est_h1, est_h0 = est_h0)
#'
#' # Double centering approach
#' est_h1 <- modsem(m1, oneInt, method = "dblcent")
#' est_h0 <- estimate_h0(est_h1, oneInt)
#'
#' compare_fit(est_h1 = est_h1, est_h0 = est_h0)
#'
#' # Constrained approach
#' est_h1 <- modsem(m1, oneInt, method = "ca")
#' est_h0 <- estimate_h0(est_h1, oneInt)
#'
#' compare_fit(est_h1 = est_h1, est_h0 = est_h0)
#' }
#' @export
estimate_h0 <- function(object, warn_no_interaction = TRUE, ...) {
  UseMethod("estimate_h0", object)
}


#' @describeIn estimate_h0 Estimate baseline model for \code{\link{modsem_da}} objects
#' @export
estimate_h0.modsem_da <- function(object, warn_no_interaction = TRUE, ...) {
  argList    <- object$args
  parTable   <- object$originalParTable
  data       <- object$data$data.full
  method     <- object$method
  cov.syntax <- object$model$covModel$syntax

  newArgList <- list(...)
  newArgNames <- intersect(names(argList), names(newArgList))
  argList[newArgNames] <- newArgList[newArgNames]

  tryCatch({
    # labels should probably be replaced instead...
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


#' @describeIn estimate_h0 Estimate baseline model for \code{\link{modsem_pi}} objects
#' @param reduced Should the baseline model be a reduced version of the model?
#'   If \code{TRUE}, the latent product term and its (product) indicators are kept in the model,
#'   but the interaction coefficients are constrained to zero. If \code{FALSE}, the
#'   interaction terms are removed completely from the model. Note that the models will no longer be
#'   nested, if the interaction terms are removed from the model completely.
#' @export
estimate_h0.modsem_pi <- function(object, warn_no_interaction = TRUE,
                                  reduced = TRUE, ...) {
  pars       <- parameter_estimates(object)
  input      <- object$input
  method     <- object$method
  syntax     <- input$syntax
  data       <- input$data
  parTable   <- input$parTable

  # Get arguments
  lavArgs    <- input$lavArgs
  modsemArgs <- input$modsemArgs
  lavArgs$meanstructure <- any(pars$op == "~1")
  argList    <- c(modsemArgs, lavArgs)

  newArgList <- list(...)
  newArgNames <- intersect(names(argList), names(newArgList))
  argList[newArgNames] <- newArgList[newArgNames]

  if (!is.null(lavArgs$cluster)) {
    warning2("Baseline model cannot (yet) be estimated with clustering!", .immediate. = FALSE)
    return(NULL)
  }

  tryCatch({
    if (!parTableHasInteraction(parTable)) {
      warnif(warn_no_interaction,
             "No interaction terms found in the model. ",
             "The baseline model is identical to the original model, ",
             "and won't be estimated!")
      return(NULL)

    } else if (parTableHasHigherOrderInteraction(parTable)) {
      warning2("Unable to estimate baseline model for models with higher-order interaction terms!",
               immediate. = FALSE)
      return(NULL)
    }

    isInteractionTerm <- grepl(":", parTable$rhs)
    constrainedParTable <- parTable

    # If there are any label modifiers, we must not remove them completely
    # from the model, as they might be used elsewhere in the model
    preExistingMod <- constrainedParTable[isInteractionTerm, "mod"]
    preExistingMod <- preExistingMod[preExistingMod != ""]
    preExistingLab <- unique(preExistingMod[!canBeNumeric(preExistingMod)])

    # Redine the labels as zero, if they exist
    if (length(preExistingLab)) {
      redefinedLabels <- data.frame(lhs = preExistingLab, op  = ":=", rhs = "0",
                                    mod = preExistingLab)
      constrainedParTable <- rbind(redefinedLabels, constrainedParTable)
    }

    # We don't want to remove the interaction term, since we want the baseline model
    # to have the same variables as the original model (inluding the product indicators).
    # However, we want to constrain the interaction terms to zero.

    isInteractionTerm <- grepl(":", constrainedParTable$rhs) # recompute as parTable might have changed
    if (reduced) constrainedParTable[isInteractionTerm, "mod"] <- "0"
    else         constrainedParTable <- constrainedParTable[!isInteractionTerm, ]

    syntax <- parTableToSyntax(constrainedParTable)
    argList <- c(list(model = syntax, data = data, method = method), argList)

    fit <- do.call(modsem_pi, args = argList)

    if (is.null(extract_lavaan(fit))) {
      warning2("`lavaan` failed when estimating the model!", immediate. = FALSE)
      return(NULL)
    }

    fit
  },
  error = function(e) {
    warning2("Baseline model could not be estimated: ", e$message,
             immediate. = FALSE)
    NULL
  })
}


parTableHasInteraction <- function(parTable) {
  any(grepl(":", parTable$rhs) & parTable$mod != "0")
}


parTableHasHigherOrderInteraction <- function(parTable) {
  any(grepl(":", parTable$rhs) & parTable$op == "=~")
}
