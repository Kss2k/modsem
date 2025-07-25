#' compare model fit for \code{modsem} models
#'
#' @param est_h0 object of class \code{\link{modsem_da}} or \code{\link{modsem_pi}} representing the
#' null hypothesis model (without interaction terms).
#'
#' @param est_h1 object of class \code{\link{modsem_da}} or \code{\link{modsem_pi}} representing the
#' alternative hypothesis model (with interaction terms).
#'
#' @param ... additional arguments passed to the underlying comparison function. E.g.,
#' for \code{modsem_pi} models, this can be used to pass arguments to \code{lavaan::lavTestLRT}.
#' currently only used for \code{modsem_pi} models.
#'
#' @description Compare the fit of two models using the likelihood ratio test (LRT).
#' \code{est_h0} is the null hypothesis model, and \code{est_h1} the alternative hypothesis model.
#' Importantly, the function assumes that \code{est_h0} does not have more free parameters
#' (i.e., degrees of freedom) than \code{est_h1} (the alternative hypothesis model).
#' @rdname compare_fit
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
compare_fit <- function(est_h1, est_h0, ...) {
  UseMethod("compare_fit", est_h1)
}


#' @export
compare_fit.modsem_da <- function(est_h1, est_h0, ...) {
  stopif(!inherits(est_h1, "modsem_da") || !inherits(est_h0, "modsem_da"),
         "Both est_h1 and est_h0 must be of class 'modsem_da'.")

  if (is.null(est_h0) || is.null(est_h1)) return(NULL) # should never trigger

  df <- length(coef(est_h1, type = "free")) - length(coef(est_h0, type = "free"))
  D <- -2 * (est_h0$logLik - est_h1$logLik)
  p <- stats::pchisq(D, df = df, lower.tail = FALSE, log.p = FALSE)
  list(D = D, df = df, p = p, diff.loglik = est_h1$logLik - est_h0$logLik)
}


#' @export
compare_fit.modsem_pi <- function(est_h1, est_h0, comparison = lavaan::lavTestLRT, ...) {
  stopif(!inherits(est_h1, "modsem_pi") || !inherits(est_h0, "modsem_pi"),
         "Both est_h1 and est_h0 must be of class 'modsem_pi'.")

  est_h0 <- extract_lavaan(est_h0)
  est_h1 <- extract_lavaan(est_h1)
  comparison(est_h0, est_h1, ...)
}
