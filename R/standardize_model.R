#' Standardize the parameter vector and VCOV matrix of a \code{modsem_da} model
#'
#' \code{standardize_model()} post-processes the output of
#' [modsem_da()] (or of [modsem()] when \code{method = "lms"} /
#' \code{method = "qml"}), replacing the *unstandardized* coefficient vector
#' (\code{$coefs}) and its variance–covariance matrix (\code{$vcov}) with *fully
#' standardized* counterparts (including the measurement model).The procedure is purely
#' algebraic—**no re-estimation is carried out**—and is therefore fast and
#' deterministic.
#'
#' @param model A fitted object of class \code{modsem_da}.
#'   Passing any other object triggers an error.
#'
#' @return The same object (returned invisibly) with three slots overwritten  
#' \describe{
#'   \item{\code{$parTable}}{Parameter table whose columns \code{est} and \code{std.error}
#'         now hold standardized estimates and their (delta-method)
#'         standard errors, as produced by [standardized_estimates()].}
#'   \item{\code{$coefs}}{Named numeric vector of standardized coefficients.
#'         Intercepts (operator \code{~1}) are removed, because a standardized
#'         variable has mean 0 by definition.}
#'   \item{\code{$vcov}}{Variance–covariance matrix corresponding to the updated
#'         coefficient vector.  Rows/columns for intercepts are dropped, and
#'         the sub-matrix associated with rescaled parameters is adjusted via
#'         [rescaleVCOV()] so that its diagonal equals the squared
#'         standardized standard errors.}
#' }
#' The object keeps its class attributes, allowing seamless use by downstream
#' S3 methods such as \code{summary()}, \code{coef()}, or \code{vcov()}.
#'
#' Because the function merely transforms existing estimates, *parameter
#' constraints imposed through shared labels remain satisfied*. 
#'
#' @seealso
#' * [standardized_estimates()] – obtains the fully standardized
#'   parameter table used here.  
#' * [modsem_da()], [modsem()] for model fitting.
#'
#' @examples
#' \dontrun{
#' # Latent interaction estimated with LMS and standardized afterwards
#' syntax <- "
#'   X =~ x1 + x2 + x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'   Y ~ X + Z + X:Z
#' "
#' fit  <- modsem_da(syntax, data = oneInt, method = "lms")
#' sfit <- standardize_model(fit)
#'
#' # Compare unstandardized vs. standardized summaries
#' summary(fit)      # raw metric
#' summary(sfit)     # correlation metric
#' }
#'
#' @export
standardize_model <- function(model) {
  stopif(!inherits(model, "modsem_da"), "The model must be of class 'modsem_da'.")

  parTable <- standardized_estimates(model)

  labels <- ifelse(parTable$label == "", 
                   yes = paste0(parTable$lhs, parTable$op, parTable$rhs),
                   no = parTable$label)
  interceptLabels <- labels[parTable$op == "~1"]

  coefs  <- coef(model)
  params <- intersect(names(coefs), labels)
  scoefs <- structure(parTable$est, names = labels)

  coefs[params] <- scoefs[params]
  coefs <- coefs[!names(coefs) %in% interceptLabels]
  
  vcov    <- vcov(model)
  params  <- intersect(colnames(vcov), labels)
  std_se  <- structure(parTable$std.error, names = labels)

  vcov[params, params] <- rescaleVCOV(vcov[params, params], new_se = std_se[params])
  vcov <- vcov[!rownames(vcov) %in% interceptLabels, 
               !colnames(vcov) %in% interceptLabels]

  model$type.estimates <- "standardized"
  model$parTable       <- parTable
  model$coefs          <- coefs
  model$vcov           <- vcov

  model
}


rescaleVCOV <- function(vcov, new_se) {
  isZeroOrNA <- is.na(new_se) | new_se == 0
  S <- vcov[!isZeroOrNA, !isZeroOrNA]

  D_old_inv <- diag(1 / sqrt(diag(S)))
  D_new <- diag(new_se[!isZeroOrNA])

  R <- D_old_inv %*% S %*% D_old_inv
  Z <- D_new %*% R %*% D_new

  V <- vcov
  V[!isZeroOrNA, !isZeroOrNA] <- Z

  V
}
