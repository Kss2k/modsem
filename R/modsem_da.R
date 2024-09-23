#' Interaction between latent variables using lms and qml approaches
#'
#' @param model.syntax \code{lavaan} syntax
#'
#' @param data dataframe
#'
#' @param method method to use:
#' \code{"lms"} = latent model structural equations (not passed to \code{lavaan}).
#' \code{"qml"} = quasi maximum likelihood estimation of latent model structural equations (not passed to \code{lavaan}).
#'
#' @param verbose should estimation progress be shown
#'
#' @param optimize should starting parameters be optimized
#'
#' @param nodes number of quadrature nodes (points of integration) used in \code{lms},
#' increased number gives better estimates but slower computation. How many are needed depends on the complexity of the model.
#' For simple models, somewhere between 16-24 nodes should be enough; for more complex models, higher numbers may be needed. 
#' For models where there is an interaction effect between an endogenous and exogenous variable, 
#' the number of nodes should be at least 32, but practically (e.g., ordinal/skewed data), more than 32 is recommended. In cases 
#' where data is non-normal, it might be better to use the \code{qml} approach instead. For large
#' numbers of nodes, you might want to change the \code{'quad.range'} argument.
#'
#' @param convergence convergence criterion. Lower values give better estimates but slower computation.
#' 
#' @param optimizer optimizer to use, can be either \code{"nlminb"} or \code{"L-BFGS-B"}. For LMS, \code{"nlminb"} is recommended. 
#' For QML, \code{"L-BFGS-B"} may be faster if there is a large number of iterations, but slower if there are few iterations.
#'
#' @param center.data should data be centered before fitting model
#'
#' @param standardize.data should data be scaled before fitting model, will be overridden by
#' \code{standardize} if \code{standardize} is set to \code{TRUE}.
#' 
#' \strong{NOTE}: It is recommended that you estimate the model normally and then standardize the output using
#' \code{\link{standardized_estimates}}.
#'
#' @param standardize.out should output be standardized (note will alter the relationships of
#' parameter constraints since parameters are scaled unevenly, even if they
#' have the same label). This does not alter the estimation of the model, only the
#' output.
#' 
#' \strong{NOTE}: It is recommended that you estimate the model normally and then standardize the output using
#' \code{\link{standardized_estimates}}.
#'
#' @param mean.observed should the mean structure of the observed variables be estimated?
#' This will be overridden by \code{standardize} if \code{standardize} is set to \code{TRUE}.
#' 
#' \strong{NOTE}: Not recommended unless you know what you are doing.
#'
#' @param standardize will standardize the data before fitting the model, remove the mean
#' structure of the observed variables, and standardize the output. Note that \code{standardize.data},
#' \code{mean.observed}, and \code{standardize.out} will be overridden by \code{standardize} if \code{standardize} is set to \code{TRUE}.
#' 
#' \strong{NOTE}: It is recommended that you estimate the model normally and then standardize the output using
#' \code{\link{standardized_estimates}}.
#'
#' @param double try to double the number of dimensions of integration used in LMS,
#' this will be extremely slow but should be more similar to \code{mplus}.
#'
#' @param cov.syntax model syntax for implied covariance matrix (see \code{vignette("interaction_two_etas", "modsem")})
#'
#' @param calc.se should standard errors be computed? \strong{NOTE}: If \code{FALSE}, the information matrix will not be computed either.
#'
#' @param FIM should the Fisher information matrix be calculated using the observed or expected values? Must be either \code{"observed"} or \code{"expected"}.
#'
#' @param EFIM.S if the expected Fisher information matrix is computed, \code{EFIM.S} selects the sample size of the generated data.
#'
#' @param OFIM.hessian should the observed Fisher information be computed using the Hessian? If \code{FALSE}, it is computed using the gradient.
#'
#' @param EFIM.parametric should data for calculating the expected Fisher information matrix be 
#' simulated parametrically (simulated based on the assumptions and implied parameters
#' from the model), or non-parametrically (stochastically sampled)? If you believe that 
#' normality assumptions are violated, \code{EFIM.parametric = FALSE} might be the better option.
#'
#' @param robust.se should robust standard errors be computed? Meant to be used for QML, 
#' can be unreliable with the LMS approach.
#'
#' @param max.iter maximum number of iterations.
#'
#' @param max.step maximum steps for the M-step in the EM algorithm (LMS).
#'
#' @param fix.estep if \code{TRUE}, the E-step will be fixed, and the prior probabilities will be set to the best prior probabilities, 
#' if the log-likelihood decreases for more than 30 iterations.
#'
#' @param start starting parameters.
#'
#' @param epsilon finite difference for numerical derivatives.
#'
#' @param quad.range range in z-scores to perform numerical integration in LMS using 
#' Gaussian-Hermite Quadratures. By default \code{Inf}, such that \code{f(t)} is integrated from -Inf to Inf, 
#' but this will likely be inefficient and pointless at a large number of nodes. Nodes outside 
#' \code{+/- quad.range} will be ignored.
#'
#' @param n.threads number of cores to use for parallel processing. If \code{NULL}, it will use <= 2 threads. 
#' If an integer is specified, it will use that number of threads (e.g., \code{n.threads = 4} will use 4 threads).
#' If \code{"default"}, it will use the default number of threads (2).
#' If \code{"max"}, it will use all available threads, \code{"min"} will use 1 thread.
#'
#' @param ... additional arguments to be passed to the estimation function.
#'
#' @return \code{modsem_da} object
#' @export
#' 
#' @description 
#' \code{modsem_da()} is a function for estimating interaction effects between latent variables
#' in structural equation models (SEMs) using distributional analytic (DA) approaches.
#' Methods for estimating interaction effects in SEMs can basically be split into
#' two frameworks: 
#' 1. Product Indicator-based approaches (\code{"dblcent"}, \code{"rca"}, \code{"uca"},
#' \code{"ca"}, \code{"pind"})
#' 2. Distributionally based approaches (\code{"lms"}, \code{"qml"}).
#' 
#' \code{modsem_da()} handles the latter and can estimate models using both QML and LMS,
#' necessary syntax, and variables for the estimation of models with latent product indicators.
#' 
#' \strong{NOTE}: Run \code{\link{default_settings_da}} to see default arguments.
#'
#' @examples
#' library(modsem)
#' # For more examples, check README and/or GitHub.
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
#'   # Causal Relationships
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#'   BEH ~ INT:PBC
#' "
#'
#' # LMS Approach
#' estTpb <- modsem_da(tpb, data = TPB, method = lms)
#' summary(estTpb)
#' }
modsem_da <- function(model.syntax = NULL,
                      data = NULL,
                      method = "lms",
                      verbose = NULL,
                      optimize = NULL,
                      nodes = NULL,
                      convergence = NULL,
                      optimizer = NULL,
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
                      fix.estep = NULL,
                      start = NULL,
                      epsilon = NULL,
                      quad.range = NULL,
                      n.threads = NULL,
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
          optimizer = optimizer,
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
          fix.estep = fix.estep,
          epsilon = epsilon,
          quad.range = quad.range,
          n.threads = n.threads
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
    double = args$double,
    quad.range = args$quad.range
  )

  if (args$optimize) {
    model <- tryCatch(optimizeStartingParamsDA(model), 
                      warning = function(w) {
                       warning2("warning when optimizing starting parameters:\n", w)
                       suppressWarnings(optimizeStartingParamsDA(model))
                      }, error = function(e) {
                       warning2("unable to optimize starting parameters:\n", e)
                       model
                      })
  }
  
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
      optimizer = args$optimizer,
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
      optimizer = args$optimizer,
      fix.estep = args$fix.estep,
      ...
    )
  )

  class(est) <- c("modsem_da", "modsem")
  if (args$standardize.out) {
    est$type.estimates <- "standardized"
    est$parTable <- standardized_estimates(est)
  }

  # clean up
  resetThreads()

  est$args <- args
  est
}
