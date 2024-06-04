#' Interaction between latent variables
#'
#' @param model_syntax lavaan syntax
#' @param data dataframe
#' @param method method to use:
#' "rca" = residual centering approach (passed to lavaan),
#' "uca" = unconstrained approach (passed to lavaan),
#' "dblcent" = double centering approach (passed to lavaan),
#' "pind" = prod ind approach, with no constraints or centering (passed to lavaan),
#' "lms" = laten model structural equations (not passed to lavaan).
#' "qml" = quasi maximum likelihood estimation of laten model structural equations (not passed to lavaan).
#' "custom" = use parameters specified in the function call (passed to lavaan)
#' @param standardize should data be scaled before fitting model
#' @param center should data be centered before fitting model
#' @param ... arguments passed to other functions depending on method (see modsem_pi, modsem_lms_qml, and modsem_mplus)
#' @return ModSEM object
#' @export 
#' @description
#' modsem is a function for estimating interaction effects between latent variables, 
#' in structural equation models (SEM's).
#' Methods for estimating interaction effects in SEM's can basically be split into 
#' two frameworks: 1. Product Indicator based approaches ("dblcent", "rca", "uca", 
#' "ca", "pind"), and 2. Distributionally based approaches ("lms", "qml").
#' For the product indicator based approaces, modsem() is essentially a just 
#' a fancy wrapper for lavaan::sem()  which generates the 
#' necessary syntax, and variables for the estimation of models with latent product indicators.
#' The distributionally based approaches are implemented in seperately, and are 
#' are not estimated using lavaan::sem(), but rather using custom functions (largely)
#' written in C++ for performance reasons. For greater control, it is advised that 
#' you use one of the sub-functions (modsem_pi, modsem_lms_qml, modsem_mplus) directly, 
#' as passing additional arguments to them via modsem() can lead to unexpected behavior.
#' @examples
#' library(modsem)
#' # For more examples check README and/or GitHub.
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
#' est1Constrained <- modsem(m1, oneInt, method = "ca")
#' summary(est1Constrained)
#'
#' # LMS approach
#' est1LMS <- modsem(m1, oneInt, method = "lms")
#' summary(est1LMS)
#'
#' # QML approach
#' est1QML <- modsem(m1, oneInt, method = "qml")
#' summary(est1QML)
#' 
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
#'   # Covariances
#'   ATT ~~ SN + PBC
#'   PBC ~~ SN 
#'   # Causal Relationsships
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC 
#'   BEH ~ INT:PBC  
#' '
#' 
#' # double centering approach
#' estTpb <- modsem(tpb, data = TPB)
#' summary(estTpb)
#'
#' \dontrun{
#' # The Constrained Approach 
#' estTpbConstrained <- modsem(tpb, data = TPB, method = "ca")
#' summary(estTpbConstrained)
#'
#' # LMS approach
#' estTpbLMS <- modsem(tpb, data = TPB, method = "lms")
#' summary(estTpbLMS)
#' }
modsem <- function(model_syntax = NULL,
                   data = NULL,
                   method = "dblcent",
                   standardize = FALSE,
                   center = FALSE,
                   ...) {
  if (is.null(model_syntax)) {
    stop("No model_syntax provided")
  } else if (!is.character(model_syntax)) {
    stop("The provided model syntax is not a string!")
  } else if (length(model_syntax) > 1) {
    stop("The provided model syntax is not of length 1")
  }

  if (is.null(data)) {
    stop("No data provided")
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (center) data <- lapplyDf(data, FUN = function(x) x - mean(x))
  if (standardize) data <- lapplyDf(data, FUN = scaleIfNumeric, 
                                    scaleFactor = FALSE)

  if (method %in% c("rca", "uca", "dblcent", "pind", "ca", "custom")) {
    return(modsem_pi(model_syntax, data = data, method = method, ...))
  } else if (method %in% c("lms", "qml")) {
    return(modsem_lms_qml(model_syntax, data = data, method = method, ...))
  } else if (method == "mplus") {
    return(modsem_mplus(model_syntax, data = data, ...))
  } else {
    stop("Method not recognized")
  }
}
