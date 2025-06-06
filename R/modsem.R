#' Estimate interaction effects in structural equation models (SEMs)
#'
#' @param model.syntax \code{lavaan} syntax
#'
#' @param data dataframe
#'
#' @param method method to use:
#' \code{"rca"} = residual centering approach (passed to \code{lavaan}),
#' \code{"uca"} = unconstrained approach (passed to \code{lavaan}),
#' \code{"dblcent"} = double centering approach (passed to \code{lavaan}),
#' \code{"pind"} = prod ind approach, with no constraints or centering (passed to \code{lavaan}),
#' \code{"lms"} = latent model structural equations (not passed to \code{lavaan}),
#' \code{"qml"} = quasi maximum likelihood estimation of latent model structural equations (not passed to \code{lavaan}),
#' \code{"custom"} = use parameters specified in the function call (passed to \code{lavaan}).
#'
#' @param ... arguments passed to other functions depending on the method (see \code{\link{modsem_pi}}, \code{\link{modsem_da}}, and \code{\link{modsem_mplus}})
#' @return \code{modsem} object with class \code{\link{modsem_pi}}, \code{\link{modsem_da}}, or \code{\link{modsem_mplus}}
#' @export
#' @description
#' \code{modsem()} is a function for estimating interaction effects between latent variables
#' in structural equation models (SEMs).
#' Methods for estimating interaction effects in SEMs can basically be split into
#' two frameworks:
#' 1. Product Indicator-based approaches (\code{"dblcent"}, \code{"rca"}, \code{"uca"},
#' \code{"ca"}, \code{"pind"})
#' 2. Distributionally based approaches (\code{"lms"}, \code{"qml"}).
#'
#' For the product indicator-based approaches, \code{modsem()} is essentially a fancy wrapper for \code{lavaan::sem()} which generates the
#' necessary syntax and variables for the estimation of models with latent product indicators.
#'
#' The distributionally based approaches are implemented separately and are
#' not estimated using \code{lavaan::sem()}, but rather using custom functions (largely
#' written in \code{C++} for performance reasons). For greater control, it is advised that
#' you use one of the sub-functions (\link{modsem_pi}, \link{modsem_da}, \link{modsem_mplus}) directly,
#' as passing additional arguments to them via \code{modsem()} can lead to unexpected behavior.
#'
#' @examples
#' library(modsem)
#' # For more examples, check README and/or GitHub.
#' # One interaction
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 +x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'
#'   # Inner model
#'   Y ~ X + Z + X:Z
#' '
#'
#' # Double centering approach
#' est1 <- modsem(m1, oneInt)
#' summary(est1)
#'
#' \dontrun{
#' # The Constrained Approach
#' est1_ca <- modsem(m1, oneInt, method = "ca")
#' summary(est1_ca)
#'
#' # LMS approach
#' est1_lms <- modsem(m1, oneInt, method = "lms", EFIM.S=1000)
#' summary(est1_lms)
#'
#' # QML approach
#' est1_qml <- modsem(m1, oneInt, method = "qml")
#' summary(est1_qml)
#' }
#'
#' # Theory Of Planned Behavior
#' tpb <- '
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#'   BEH ~ INT:PBC
#' '
#'
#' # Double centering approach
#' est_tpb <- modsem(tpb, data = TPB)
#' summary(est_tpb)
#'
#' \dontrun{
#' # The Constrained Approach
#' est_tpb_ca <- modsem(tpb, data = TPB, method = "ca")
#' summary(est_tpb_ca)
#'
#' # LMS approach
#' est_tpb_lms <- modsem(tpb, data = TPB, method = "lms")
#' summary(est_tpb_lms)
#'
#' # QML approach
#' est_tpb_qml <- modsem(tpb, data = TPB, method = "qml")
#' summary(est_tpb_qml)
#' }
modsem <- function(model.syntax = NULL,
                   data = NULL,
                   method = "dblcent",
                   ...) {
  if (is.null(model.syntax)) {
    stop2("No model.syntax provided")
  } else if (!is.character(model.syntax)) {
    stop2("The provided model syntax is not a string!")
  } else if (length(model.syntax) > 1) {
    stop2("The provided model syntax is not of length 1")
  }

  if (is.null(data)) {
    stop2("No data provided")
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (method %in% c("rca", "uca", "dblcent", "pind", "ca", "custom")) {
    modsem_pi(model.syntax, data = data, method = method, ...)
  } else if (method %in% c("lms", "qml")) {
    modsem_da(model.syntax, data = data, method = method, ...)
  } else if (method == "mplus") {
    modsem_mplus(model.syntax, data = data, ...)
  } else {
    stop2("Method not recognized")
  }
}
