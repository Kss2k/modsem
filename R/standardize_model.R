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
#'   If \code{FALSE}, the function will simply rescale the standard errors. This means that 
#'   the function will not be able to provide standard errors for fixed parameters. But ensures
#'   that p and z values remain unchanges for the standardized model.
#'   Default is \code{FALSE}.
#'
#' @param mc.reps Number of Monte Carlo replications. Default is 10,000.
#'   Ignored if \code{monte.carlo = FALSE}.
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
standardize_model <- function(model, monte.carlo = TRUE, mc.reps = 10000) {
  stopif(!inherits(model, "modsem_da"), "The model must be of class 'modsem_da'.")

  if (monte.carlo) {
    solution <- standardized_solution_mc(model, mc.reps = mc.reps)

    vcov <- solution$vcov
    coefs <- solution$coefs
    parTable <- solution$parTable

  } else {
    parTable <- standardized_estimates(model, monte.carlo = FALSE)

    labels <- getParTableLabels(parTable, labelCol="label")
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
  }

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
