#' Interaction between latent variables using lms and qml approaches
#' @param modelSyntax lavaan syntax
#' @param data dataframe
#' @param method method to use:
#' "lms" = laten model structural equations (not passed to lavaan).
#' "qml" = quasi maximum likelihood estimation of laten model structural equations (not passed to lavaan).
#' @param verbose should estimation progress be shown
#' @param optimize should starting parameters be optimized
#' @param nodes number of quadrature nodes (points of integration) used in lms 
#' @param convergence convergence criterion
#' @param standardize should data be scaled before fitting model
#' @param center should data be centered before fitting model
#' @param ... arguments passed to other functions
#' @return modsem_lms or modsem_qml object
#' @export 
#' @description
#' modsem_lms_qml is a function for estimating interaction effects between latent variables, 
#' in structural equation models (SEMs), using product indicators.
#' Methods for estimating interaction effects in SEM's can basically be split into 
#' two frameworks: 1. Product Indicator based approaches ("dblcent", "rca", "uca", 
#' "ca", "pind"), and 2. Distributionally based approaches ("lms", "qml").
#' modsem_lms_qml() is essentially a just 
#' a fancy wrapper for lavaan::sem()  which generates the 
#' necessary syntax, and variables for the estimation of models with latent product indicators.
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
#' \dontrun{
#' # QML Approach
#' est1 <- modsem_lms_qml(m1, oneInt, method = "qml")
#' summary(est1)
#' 
#'
#' # Theory Of Planned Behavior
#' tpb <- ' 
#' # Outer Model (Based on Hagger et al., 2007)
#'   LATT =~ att1 + att2 + att3 + att4 + att5
#'   LSN =~ sn1 + sn2
#'   LPBC =~ pbc1 + pbc2 + pbc3
#'   LINT =~ int1 + int2 + int3
#'   LBEH =~ b1 + b2
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   # Covariances
#'   LATT ~~ LSN + LPBC
#'   LPBC ~~ LSN 
#'   # Causal Relationsships
#'   LINT ~ LATT + LSN + LPBC
#'   LBEH ~ LINT + LPBC 
#'   LBEH ~ LINT:LPBC  
#' '
#' 
#' # lms approach
#' estTpb <- modsem_lms_qml(tpb, data = TPB, method = lms)
#' summary(estTpb)
#' }
#'
modsem_lms_qml <- function(modelSyntax = NULL,
                           data = NULL,
                           method = "lms",
                           verbose = FALSE, 
                           optimize = TRUE,
                           nodes = 16, 
                           convergence = 1e-2,
                           center = FALSE, 
                           standardize = FALSE,
                       ...) {
  if (is.null(modelSyntax)) {
    stop("No modelSyntax provided")
  } else if (!is.character(modelSyntax)) {
    stop("The provided model syntax is not a string!")
  } else if (length(modelSyntax) > 1) {
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

  model <- specifyLmsModel(modelSyntax, data = data, 
                           method = method, m = nodes)
  if (optimize) model <- optimizeStartingParamsLms(model)
  switch(method, 
         "qml" = estQml(model, verbose = verbose, 
                        convergence = convergence, ...),
         "lms" = emLms(model, verbose = verbose, 
                       convergence = convergence, ...))
}
