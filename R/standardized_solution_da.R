#' Standardize a fitted \code{modsem_da} model
#'
#' \code{standardize_model()} post-processes the output of
#' \code{\link{modsem_da}()} (or of \code{\link{modsem}())} when \code{method = "lms"} /
#' \code{method = "qml"}), replacing the \emph{unstandardized} coefficient vector
#' (\code{$coefs}) and its variance–covariance matrix (\code{$vcov}) with \emph{fully
#' standardized} counterparts (including the measurement model).The procedure is purely
#' algebraic— \strong{no re-estimation is carried out} —and is therefore fast and
#' deterministic.
#'
#' @param model A fitted object of class \code{modsem_da}.
#'   Passing any other object triggers an error.
#'
#' @param monte.carlo Logical. If \code{TRUE}, the function will use Monte Carlo
#'   simulation to obtain the standard errors of the standardized estimates.
#'   If \code{FALSE}, the delta method is used.
#'   Default is \code{FALSE}.
#'
#' @param mc.reps Number of Monte Carlo replications. Default is 10,000.
#'   Ignored if \code{monte.carlo = FALSE}.
#'
#' @param ... Arguments passed on to other functions
#'
#' @return The same object (returned invisibly) with three slots overwritten
#' \describe{
#'   \item{\code{$parTable}}{Parameter table whose columns \code{est} and \code{std.error}
#'         now hold standardized estimates and their (delta-method)
#'         standard errors, as produced by \code{\link{standardized_estimates}()}.}
#'   \item{\code{$coefs}}{Named numeric vector of standardized coefficients.
#'         Intercepts (operator \code{~1}) are removed, because a standardized
#'         variable has mean 0 by definition.}
#'   \item{\code{$vcov}}{Variance–covariance matrix corresponding to the updated
#'         coefficient vector.  Rows/columns for intercepts are dropped, and
#'         the sub-matrix associated with rescaled parameters is adjusted
#'         so that its diagonal equals the squared standardized standard errors.}
#' }
#' The object keeps its class attributes, allowing seamless use by downstream
#' S3 methods such as \code{summary()}, \code{coef()}, or \code{vcov()}.
#'
#' Because the function merely transforms existing estimates, \emph{parameter
#' constraints imposed through shared labels remain satisfied}.
#'
#' @seealso
#' \describe{
#'    \item{\code{\link{standardized_estimates}()}}{Obtains the fully standardized
#'          parameter table used here.}
#'    \item{\code{\link{modsem}()}}{Fit model using LMS or QML approaches.}
#'    \item{\code{\link{modsem_da}()}}{Fit model using LMS or QML approaches.}
#' }
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
#' sfit <- standardize_model(fit, monte.carlo = TRUE)
#'
#' # Compare unstandardized vs. standardized summaries
#' summary(fit)  # unstandardized
#' summary(sfit) # standardized
#' }
#'
#' @export
standardize_model <- function(model, monte.carlo = FALSE, mc.reps = 10000, ...) {
  stopif(!inherits(model, "modsem_da"), "The model must be of class 'modsem_da'.")

  solution <- standardizedSolutionCOEFS(model, monte.carlo = monte.carlo, mc.reps = mc.reps, ...)

  coefs.all <- solution$coefs
  vcov.all  <- solution$vcov

  # It should be safe to use the same subset for both vcov.all and coefs.all
  # but in case something has gone wrong, we create seperate masks
  cnames.all <- names(coefs.all)
  vnames.all <- rownames(vcov.all)
  # intersect() shouldn't be neccessary, but just in case...
  cnames.free <- intersect(cnames.all, names(coef(model, type = "free")))
  vnames.free <- intersect(vnames.all, rownames(vcov(model, type = "free")))

  coefs.free <- coefs.all[cnames.free]
  vcov.free  <- vcov.all[vnames.free, vnames.free, drop = FALSE]

  parTable <- solution$parTable

  model$type.estimates <- "standardized"
  model$parTable       <- parTable

  model$coefs.all      <- coefs.all
  model$coefs.free     <- coefs.free

  model$vcov.all       <- vcov.all
  model$vcov.free      <- vcov.free

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
