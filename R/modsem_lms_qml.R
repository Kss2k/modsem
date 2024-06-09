# Estimate SEM using lms and qml
# Last updated: 29.05.2024

#' Interaction between latent variables using lms and qml approaches
#' @param model_syntax lavaan syntax
#' @param data dataframe
#' @param method method to use:
#' "lms" = laten model structural equations (not passed to lavaan).
#' "qml" = quasi maximum likelihood estimation of laten model structural equations (not passed to lavaan).
#' @param verbose should estimation progress be shown
#' @param optimize should starting parameters be optimized
#' @param nodes number of quadrature nodes (points of integration) used in lms, 
#' increased number gives better estimates but slower computation. How many is needed, depends on the complexity of the model 
#' for simple models 16 is enough, for more complex models 100 is usually enough.
#' @param convergence convergence criterion. Lower values give better estimates but slower computation.
#' @param center should data be centered before fitting model
#' @param standardize should data be scaled before fitting model
#' @param meanStructure should mean structure of the observed variables be estimated
#' @param double try to double the number of dimensions of integrations used in LMS,
#' this will be extremely slow, but should be more similar to mplus.
#' @param hessian should hessian (i.e., std.errors) be calculated
#' @param cov_syntax model syntax for implied covariance matrix (see 'vignette("interaction_two_etas", "modsem")')
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
#' # lms approach
#' estTpb <- modsem_lms_qml(tpb, data = TPB, method = lms)
#' summary(estTpb)
#' }
#'
modsem_lms_qml <- function(model_syntax = NULL,
                           data = NULL,
                           method = "lms",
                           verbose = FALSE, 
                           optimize = TRUE,
                           nodes = 16, 
                           convergence = 1e-2,
                           center = FALSE, 
                           standardize = FALSE,
                           meanStructure = TRUE,
                           double = FALSE, 
                           hessian = TRUE,
                           cov_syntax = NULL,
                           ...) {
  if (is.null(model_syntax)) {
    stop2("No model_syntax provided")
  } else if (!is.character(model_syntax)) {
    stop2("The provided model syntax is not a string!")
  } else if (length(model_syntax) > 1) {
    stop2("The provided model syntax is not of length 1")
  }

  if (is.null(data)) {
    stop2("No data provided")
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
    
  if (!method %in% c("lms", "qml")) {
    stop2("Method must be either 'lms' or 'qml'")
  }

  if (center) {
    data <- lapplyDf(data, FUN = function(x) x - mean(x))
  }
  if (standardize) {
    data <- lapplyDf(data, FUN = scaleIfNumeric, scaleFactor = FALSE)
  }

  model <- specifyModelLmsQml(model_syntax, data = data, 
                              method = method, m = nodes, 
                              cov_syntax = cov_syntax, 
                              double = double)
  if (optimize) model <- optimizeStartingParamsLms(model)
  switch(method, 
         "qml" = estQml(model, verbose = verbose, 
                        convergence = convergence, 
                        hessian = hessian, ...),
         "lms" = emLms(model, verbose = verbose, 
                       convergence = convergence,
                       hessian = hessian, ...))
}
