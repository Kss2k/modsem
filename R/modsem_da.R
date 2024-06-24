# Estimate SEM using distribution analytic (DA) approaches
# Last updated: 29.05.2024

#' Interaction between latent variables using lms and qml approaches
#'
#' @param model.syntax lavaan syntax
#'
#' @param data dataframe
#'
#' @param method method to use:
#' "lms" = laten model structural equations (not passed to lavaan).
#' "qml" = quasi maximum likelihood estimation of laten model structural equations (not passed to lavaan).
#'
#' @param verbose should estimation progress be shown
#'
#' @param optimize should starting parameters be optimized
#'
#' @param nodes number of quadrature nodes (points of integration) used in lms,
#' increased number gives better estimates but slower computation. How many is needed, depends on the complexity of the model
#' For simple models you somwhere between 16-24 should be enough, for more complex higher numbers may be needed. 
#' For models where there is an interaction effects between and endogenous and exogenous variable 
#' the number of nodes should at least be 32, but practically (e.g., ordinal/skewed data) more than 32 is recommended. In cases, 
#' where data is non-normal it might be better to use the qml approach instead.
#'
#' @param convergence convergence criterion. Lower values give better estimates but slower computation.
#'
#' @param center.data should data be centered before fitting model
#'
#' @param standardize.data should data be scaled before fitting model, will be overridden by
#' standardize if standardize is set to TRUE.
#' NOTE: It is recommended that you estimate the model normally and then standardize the output using
#' `standardized_estimates()`.
#'
#' @param standardize.out should output be standardized (note will alter the relationsships of
#' parameter constraints, since to parameters are scaled unevenly, even if they
#' have the same label). This does not alter the estimation of the model, only the
#' output.
#' NOTE: It is recommended that you estimate the model normally and then standardize the output using
#' `standardized_estimates()`.
#'
#' @param mean.observed should mean structure of the observed variables be estimated,
#' will be overridden by standardize if standardize is set to TRUE.
#' NOTE: Not recommended unless you know what you are doing.
#'
#' @param standardize will standardize the data before fitting the model, remove the mean
#' structure of the observed variables, and standardize the output. Note that standardize.data
#' mean.observed, standardize.out will be overridden by standardize if standardize is set to TRUE.
#' NOTE: It is recommended that you estimate the model normally and then standardize the output using
#' `standardized_estimates()`.
#'
#' @param double try to double the number of dimensions of integrations used in LMS,
#' this will be extremely slow, but should be more similar to mplus.
#'
#' @param cov.syntax model syntax for implied covariance matrix (see 'vignette("interaction_two_etas", "modsem")')
#'
#' @param double = NULL,
#'
#' @param calc.se should standard errros be computed, NOTE: If 'FALSE' information matrix will not be computed either
#'
#' @param FIM should fisher information matrix be calculated using observed of expected. must be either "observed" or "expected"
#'
#' @param EFIM.S if expected fisher information matrix is computed, EFIM.S selects the sample size of the generated data
#'
#' @param OFIM.hessian should observed fisher information be computed using hessian? if FALSE, it is computed using gradient
#'
#' @param EFIM.parametric should data for calculating expected fisher information matrix be 
#' simulated parametrically (simulated based on the assumptions- and implied parameters
#' from the model), or non-parametrically (stochastically sampled). If you believe that 
#' normality assumptions are violated, 'EFIM.parametric = FALSE' might be the better option.
#'
#' @param robust.se should robust standard errors be computed? Meant to be used for QML, 
#' can be unreliable with the LMS-approach.
#' @param max.iter max numebr of iterations
#' @param max.step max steps for the M-step in the EM algorithm (LMS)
#' @param start starting parameters 
#' @param epsilon finite difference for numerical derivatives
#'
#' @param ... additional arguments to be passed to the estimation function
#'
#' @return modsem_lms or modsem_qml object
#' @export
#' 
#' @description 
#' modsem_da is a function for estimating interaction effects between latent variables,
#' in structural equation models (SEMs), using distributional analytic (DA) approaches.
#' Methods for estimating interaction effects in SEM's can basically be split into
#' two frameworks: 1. Product Indicator based approaches ("dblcent", "rca", "uca",
#' "ca", "pind"), and 2. Distributionally based approaches ("lms", "qml").
#' modsem_da() handles the latter, and can estimate models using both qml and lms
#' necessary syntax, and variables for the estimation of models with latent product indicators.
#' NOTE: run 'default_settings_da()' to see default arguments.
#'
#' @examples
#' library(modsem)
#' # For more examples check README and/or GitHub.
#' # One interaction
#' m1 <- "
#'   # Outer Model
#'   X =~ x1 + x2 +x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'
#'   # Inner model
#'   Y ~ X + Z + X:Z
#' "
#'
#' \dontrun{
#' # QML Approach
#' est1 <- modsem_da(m1, oneInt, method = "qml")
#' summary(est1)
#'
#'
#' # Theory Of Planned Behavior
#' tpb <- "
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
#' "
#'
#' # lms approach
#' estTpb <- modsem_da(tpb, data = TPB, method = lms)
#' summary(estTpb)
#' }
#'
modsem_da <- function(model.syntax = NULL,
                      data = NULL,
                      method = "lms",
                      verbose = NULL,
                      optimize = NULL,
                      nodes = NULL,
                      convergence = NULL,
                      center.data = NULL,
                      standardize.data = NULL,
                      standardize.out = NULL,
                      standardize = NULL,
                      mean.observed = NULL,
                      cov.syntax = NULL,
                      double = NULL,
                      calc.se = NULL,
                      FIM = NULL,
                      EFIM.S = NULL, 
                      OFIM.hessian = NULL, 
                      EFIM.parametric = NULL,
                      robust.se = NULL,
                      max.iter = NULL, 
                      max.step = NULL,
                      start = NULL,
                      epsilon = NULL,
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

  args <-
    getMethodSettingsDA(method,
      args =
        list(
          verbose = verbose,
          optimize = optimize,
          nodes = nodes,
          convergence = convergence,
          center.data = center.data,
          standardize.data = standardize.data,
          standardize.out = standardize.out,
          standardize = standardize,
          mean.observed = mean.observed,
          double = double,
          calc.se = calc.se,
          FIM = FIM,
          EFIM.S = EFIM.S, 
          OFIM.hessian = OFIM.hessian, 
          EFIM.parametric = EFIM.parametric,
          robust.se = robust.se,
          max.iter = max.iter, 
          max.step = max.step,
          epsilon = epsilon
        )
    )

  if (!method %in% c("lms", "qml")) {
    stop2("Method must be either 'lms' or 'qml'")
  }

  if (args$center.data) {
    data <- lapplyDf(data, FUN = function(x) x - mean(x))
  }
  if (args$standardize.data) {
    data <- lapplyDf(data, FUN = scaleIfNumeric, scaleFactor = FALSE)
  }

  model <- specifyModelDA(model.syntax,
    data = data,
    method = method,
    m = args$nodes,
    cov.syntax = cov.syntax,
    mean.observed = args$mean.observed,
    double = args$double
  )

  if (args$optimize) model <- optimizeStartingParamsDA(model)
  
  if (!is.null(start)) {
    checkStartingParams(start, model = model) # throws error if somethings wrong
    model$theta <- start
  }

  est <- switch(method,
    "qml" = estQml(model,
      verbose = args$verbose,
      convergence = args$convergence,
      calc.se = args$calc.se,
      FIM = args$FIM,
      EFIM.S = args$EFIM.S,
      OFIM.hessian = args$OFIM.hessian, 
      EFIM.parametric = args$EFIM.parametric,
      robust.se = args$robust.se,
      max.iter = args$max.iter,
      epsilon = args$epsilon,
      ...
    ),
    "lms" = emLms(model,
      verbose = args$verbose,
      convergence = args$convergence,
      calc.se = args$calc.se,
      FIM = args$FIM,
      EFIM.S = args$EFIM.S,
      OFIM.hessian = args$OFIM.hessian, 
      EFIM.parametric = args$EFIM.parametric,
      robust.se = args$robust.se,
      max.iter = args$max.iter, 
      max.step = args$max.step,
      epsilon = args$epsilon,
      ...
    )
  )

  if (args$standardize.out) est$parTable <- standardized_estimates(est)

  est$args <- args
  est
}
