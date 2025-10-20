#' Interaction between latent variables using LMS and QML approaches
#'
#' @param model.syntax \code{lavaan} syntax
#'
#' @param data A dataframe with observed variables used in the model.
#'
#' @param Character. A variable name in the data frame defining the groups in a multiple
#'   group analysis
#'
#' @param method method to use:
#' \describe{
#'   \item{\code{"lms"}}{latent moderated structural equations (not passed to \code{lavaan}).}
#'   \item{\code{"qml"}}{quasi maximum likelihood estimation (not passed to \code{lavaan}).}
#' }
#'
#' @param verbose should estimation progress be shown
#'
#' @param optimize should starting parameters be optimized
#'
#' @param nodes number of quadrature nodes (points of integration) used in \code{lms},
#'   increased number gives better estimates but slower computation. How many are needed depends on the complexity of the model.
#'   For simple models, somewhere between 16-24 nodes should be enough; for more complex models, higher numbers may be needed.
#'   For models where there is an interaction effect between an endogenous and exogenous variable,
#'   the number of nodes should be at least 32, but practically (e.g., ordinal/skewed data), more than 32 is recommended. In cases
#'   where data is non-normal, it might be better to use the \code{qml} approach instead.
#'   You can also consider setting \code{adaptive.quad = TRUE}.
#'
#' @param missing How should missing values be handled? If \code{"listwise"} (default) missing values
#'   are removed list-wise (alias: \code{"complete"} or \code{"casewise"}).
#'   If \code{impute} values are imputed using \code{Amelia::amelia}.
#'   If \code{"fiml"} (alias: \code{"ml"} or \code{"direct"}), full information maximum
#'   likelihood (FIML) is used. FIML can be (very) computationally intensive.
#'
#' @param convergence.abs Absolute convergence criterion.
#'   Lower values give better estimates but slower computation. Not relevant when
#'   using the QML approach. For the LMS approach the EM-algorithm stops whenever
#'   the relative or absolute convergence criterion is reached.
#'
#' @param convergence.rel Relative convergence criterion.
#'   Lower values give better estimates but slower computation.
#'   For the LMS approach the EM-algorithm stops whenever
#'   the relative or absolute convergence criterion is reached.
#'
#' @param optimizer optimizer to use, can be either \code{"nlminb"} or \code{"L-BFGS-B"}. For LMS, \code{"nlminb"} is recommended.
#'   For QML, \code{"L-BFGS-B"} may be faster if there is a large number of iterations, but slower if there are few iterations.
#'
#' @param center.data should data be centered before fitting model
#'
#' @param standardize.data should data be scaled before fitting model, will be overridden by
#'   \code{standardize} if \code{standardize} is set to \code{TRUE}.
#'
#' @param standardize.out should output be standardized (note will alter the relationships of
#'   parameter constraints since parameters are scaled unevenly, even if they
#'   have the same label). This does not alter the estimation of the model, only the
#'   output.
#'
#' \strong{NOTE}: It is recommended that you estimate the model normally and then standardize the output using
#' \code{\link{standardize_model}}, \code{\link{standardized_estimates}} or \code{summary(<modsem_da-object>, standardize=TRUE)}.
#'
#' @param mean.observed should the mean structure of the observed variables be estimated?
#'   This will be overridden by \code{standardize}, if \code{standardize} is set to \code{TRUE}.
#'
#' \strong{NOTE}: Not recommended unless you know what you are doing.
#'
#' @param standardize will standardize the data before fitting the model, remove the mean
#'   structure of the observed variables, and standardize the output. Note that \code{standardize.data},
#'   \code{mean.observed}, and \code{standardize.out} will be overridden by \code{standardize} if \code{standardize} is set to \code{TRUE}.
#'
#' \strong{NOTE}: It is recommended that you estimate the model normally and then standardize the output using
#'   \code{\link{standardized_estimates}}.
#'
#' @param double try to double the number of dimensions of integration used in LMS,
#' this will be extremely slow but should be more similar to \code{mplus}.
#'
#' @param cov.syntax model syntax for implied covariance matrix of exogenous latent variables
#'  (see \code{vignette("interaction_two_etas", "modsem")}).
#'
#' @param calc.se should standard errors be computed? \strong{NOTE}: If \code{FALSE}, the information matrix will not be computed either.
#'
#' @param FIM should the Fisher information matrix be calculated using the observed or expected values? Must be either \code{"observed"} or \code{"expected"}.
#'
#' @param EFIM.S if the expected Fisher information matrix is computed, \code{EFIM.S} selects the number of Monte Carlo samples. Defaults to 100.
#'   \strong{NOTE}: This number should likely be increased for better estimates (e.g., 1000), but it might drasticly increase computation time.
#'
#' @param OFIM.hessian Logical. If \code{TRUE} (default), standard errors are
#'   based on the negative Hessian (observed Fisher information).
#'   If \code{FALSE}, they come from the outer product
#'   of individual score vectors (OPG). For correctly specified models,
#'   these two matrices are asymptotically equivalent; yielding nearly identical
#'   standard errors in large samples. The Hessian usually shows smaller finite-sample
#'   variance (i.e., it's more consistent), and is therefore the default.
#'
#'   Note, that the Hessian is not always positive definite, and is more computationally
#'   expensive to calculate. The OPG should always be positive definite, and a lot
#'   faster to compute. If the model is correctly specified, and the sample size is large,
#'   then the two should yield similar results, and switching to the OPG can save a
#'   lot of time. Note, that the required sample size depends on the complexity of the model.
#'
#'   A large difference between Hessian and OPG suggests misspecification, and
#'   \code{robust.se = TRUE} should be set to obtain sandwich (robust) standard errors.
#'
#' @param EFIM.parametric should data for calculating the expected Fisher information matrix be
#'   simulated parametrically (simulated based on the assumptions and implied parameters
#'   from the model), or non-parametrically (stochastically sampled)? If you believe that
#'   normality assumptions are violated, \code{EFIM.parametric = FALSE} might be the better option.
#'
#' @param R.max Maximum population size (not sample size) used in the calculated of the expected
#'   fischer information matrix.
#'
#' @param robust.se should robust standard errors be computed, using the sandwich estimator?
#'
#' @param max.iter maximum number of iterations.
#'
#' @param max.step maximum steps for the M-step in the EM algorithm (LMS).
#'
#' @param start starting parameters.
#'
#' @param epsilon finite difference for numerical derivatives.
#'
#' @param quad.range range in z-scores to perform numerical integration in LMS using,
#'   when using quasi-adaptive Gaussian-Hermite Quadratures. By default \code{Inf}, such that \code{f(t)} is integrated from -Inf to Inf,
#'   but this will likely be inefficient and pointless at a large number of nodes. Nodes outside
#'   \code{+/- quad.range} will be ignored.
#'
#' @param adaptive.quad should a quasi adaptive quadrature be used? If \code{TRUE}, the quadrature nodes will be adapted to the data.
#'   If \code{FALSE}, the quadrature nodes will be fixed. Default is \code{FALSE}. The adaptive quadrature does not fit an adaptive
#'   quadrature to each participant, but instead tries to place more nodes where posterior distribution is highest. Compared with a
#'   fixed Gauss Hermite quadrature this usually means that less nodes are placed at the tails of the distribution.
#'
#' @param adaptive.frequency How often should the quasi-adaptive quadrature be calculated? Defaults to 3, meaning
#'   that it is recalculated every third EM-iteration.
#'
#' @param adaptive.quad.tol Relative error tolerance for quasi adaptive quadrature. Defaults to \code{1e-12}.
#'
#' @param n.threads number of threads to use for parallel processing. If \code{NULL}, it will use <= 2 threads.
#'   If an integer is specified, it will use that number of threads (e.g., \code{n.threads = 4} will use 4 threads).
#'   If \code{"default"}, it will use the default number of threads (2).
#'   If \code{"max"}, it will use all available threads, \code{"min"} will use 1 thread.
#'
#' @param algorithm algorithm to use for the EM algorithm. Can be either \code{"EM"} or \code{"EMA"}.
#'   \code{"EM"} is the standard EM algorithm. \code{"EMA"} is an
#'   accelerated EM procedure that uses Quasi-Newton and Fisher Scoring
#'   optimization steps when needed. Default is \code{"EM"}.
#'
#' @param em.control a list of control parameters for the EM algorithm. See \code{\link{default_settings_da}} for defaults.
#'
#' @param ordered Variables to be treated as ordered. Categories for ordered variables
#'   are scored, transforming them from ordinal scale to interval scale (\href{doi.org/10.1155/2014/304213}{Chen & Wang, 2014}).
#'   The underlying continous distributions
#'   are estimated analytically for indicators of exogenous variables, and using an ordered
#'   probit regression for indicators of endogenous variables. Factor scores are used as
#'   independent variables the ordered probit regressions. Interaction effects between
#'   the factor scores are included in the probit regression, if applicable.
#'   The estimates are more robust to unequal intervals in ordinal variables. I.e., the estimates
#'   should be more consistent, and less biased.
#'
#' @param ordered.probit.correction Should ordered indicators be transformed such that they
#'   reproduce their (probit) polychoric correlation matrix? This can be useful for
#'   ordered variables with only a few categories, or for linear models.
#'
#' @param cluster Clusters used to compute standard errors robust to non-indepence of observations. Must be paired with
#'   \code{robust.se = TRUE}.
#'
#' @param cr1s Logical; if \code{TRUE}, apply the CR1S small-sample correction factor
#'   to the cluster-robust variance estimator. The CR1S factor is
#'   \eqn{(G / (G - 1)) \cdot ((N - 1) / (N - q))}, where \eqn{G} is the number of
#'   clusters, \eqn{N} is the total number of observations, and \eqn{q} is the number
#'   of free parameters. This adjustment inflates standard errors to reduce the
#'   small-sample downward bias present in the basic cluster-robust (CR0) estimator,
#'   especially when \eqn{G} is small. If \code{FALSE}, the unadjusted CR0 estimator
#'   is used. Defaults to \code{TRUE}. Only relevant if \code{cluster} is specified.
#'
#' @param rcs Should latent variable indicators be replaced with reliability-corrected
#'   single item indicators instead? See \code{\link{relcorr_single_item}}.
#'
#' @param rcs.choose Which latent variables should get their indicators replaced with
#'   reliability-corrected single items? It is passed to \code{\link{relcorr_single_item}}
#'   as the \code{choose} argument.
#'
#' @param rcs.scale.corrected Should reliability-corrected items be scale-corrected? If \code{TRUE}
#'   reliability-corrected single items are corrected for differences in factor loadings between
#'   the items. Default is \code{TRUE}.
#'
#' @param orthogonal.x If \code{TRUE}, all covariances among exogenous latent variables only are set to zero.
#'  Default is \code{FALSE}.
#'
#' @param orthogonal.y If \code{TRUE}, all covariances among endogenous latent variables only are set to zero.
#'  If \code{FALSE} residual covariances are added between pure endogenous variables;
#'  those that are predicted by no other endogenous variable in the structural model.
#'  Default is \code{FALSE}.
#'
#' @param orthogonal.y If \code{TRUE}, all covariances among endogenous latent variables only are set to zero.
#'  If \code{FALSE} residual covariances are added between pure endogenous variables;
#'  those that are predicted by no other endogenous variable in the structural model.
#'  Default is \code{FALSE}.
#'
#' @param auto.fix.first If \code{TRUE} the factor loading of the first indicator, for
#'  a given latent variable is fixed to \code{1}. If \code{FALSE} no loadings are fixed
#'  (automatically). Note that that this might make it such that the model no longer is
#'  identified. Default is \code{TRUE}. \strong{NOTE} this behaviour is overridden
#'  if the first loading is labelled, where it gets treated as a free parameter instead. This
#'  differs from the default behaviour in \code{lavaan}.
#'
#' @param auto.fix.single If \code{TRUE}, the residual variance of
#'  an observed indicator is set to zero if it is the only indicator of a latent variable.
#'  If \code{FALSE} the residual variance is not fixed to zero, and treated as a free parameter
#'  of the model. Default is \code{TRUE}. \strong{NOTE} this behaviour is overridden
#'  if the first loading is labelled, where it gets treated as a free parameter instead.
#'
#' @param auto.split.syntax Should the model syntax automatically be split into a
#'   linear and non-linear part? This is done by moving the structural model for
#'   linear endogenous variables (used in interaction terms) into the \code{cov.syntax}
#'   argument. This can potentially allow interactions between two endogenous variables
#'   given that both are linear (i.e., not affected by interaction terms). This is
#'   \code{FALSE} by default for the LMS approach.
#'   When using the QML approach interation effects between exogenous and endogenous
#'   variables can in some cases be biased, if the model is not split beforehand.
#'   The default is therefore \code{TRUE} for the QML approach.
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
#' est_qml <- modsem_da(m1, oneInt, method = "qml")
#' summary(est_qml)
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
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#'   BEH ~ INT:PBC
#' "
#'
#' # LMS Approach
#' est_lms <- modsem_da(tpb, data = TPB, method = "lms")
#' summary(est_lms)
#' }
modsem_da <- function(model.syntax = NULL,
                      data = NULL,
                      group = NULL,
                      method = "lms",
                      verbose = NULL,
                      optimize = NULL,
                      nodes = NULL,
                      missing = NULL,
                      convergence.abs = NULL,
                      convergence.rel = NULL,
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
                      R.max = NULL,
                      max.iter = NULL,
                      max.step = NULL,
                      start = NULL,
                      epsilon = NULL,
                      quad.range = NULL,
                      adaptive.quad = NULL,
                      adaptive.frequency = NULL,
                      adaptive.quad.tol = NULL,
                      n.threads = NULL,
                      algorithm = NULL,
                      em.control = NULL,
                      ordered = NULL,
                      ordered.probit.correction = FALSE,
                      cluster = NULL,
                      cr1s = FALSE,
                      rcs = FALSE,
                      rcs.choose = NULL,
                      rcs.scale.corrected = TRUE,
                      orthogonal.x = NULL,
                      orthogonal.y = NULL,
                      auto.fix.first = NULL,
                      auto.fix.single = NULL,
                      auto.split.syntax = NULL,
                      ...) {
  if (is.null(model.syntax)) {
    stop2("No model.syntax provided")
  } else if (!is.character(model.syntax)) {
    stop2("The provided model syntax is not a string!")
  } else if (length(model.syntax) > 1) {
    stop2("The provided model syntax is not of length 1")
  }

  if (length(ordered) || any(sapply(data, FUN = is.ordered))) {
    out <- modsemOrderedScaleCorrection(
       model.syntax        = model.syntax,
       data                = data,
       method              = method,
       verbose             = verbose,
       optimize            = optimize,
       nodes               = nodes,
       missing             = missing,
       convergence.abs     = convergence.abs,
       convergence.rel     = convergence.rel,
       optimizer           = optimizer,
       center.data         = center.data,
       standardize.data    = standardize.data,
       standardize.out     = standardize.out,
       standardize         = standardize,
       mean.observed       = mean.observed,
       cov.syntax          = cov.syntax,
       double              = double,
       calc.se             = calc.se,
       FIM                 = FIM,
       EFIM.S              = EFIM.S,
       OFIM.hessian        = OFIM.hessian,
       EFIM.parametric     = EFIM.parametric,
       robust.se           = robust.se,
       R.max               = R.max,
       max.iter            = max.iter,
       max.step            = max.step,
       start               = start,
       epsilon             = epsilon,
       quad.range          = quad.range,
       adaptive.quad       = adaptive.quad,
       adaptive.frequency  = adaptive.frequency,
       adaptive.quad.tol   = adaptive.quad.tol,
       n.threads           = n.threads,
       algorithm           = algorithm,
       em.control          = em.control,
       ordered             = ordered,
       probit.correction   = ordered.probit.correction,
       cluster             = cluster,
       group               = group,
       cr1s                = cr1s,
       rcs                 = rcs,
       rcs.choose          = rcs.choose,
       rcs.scale.corrected = rcs.scale.corrected,
       orthogonal.x        = orthogonal.x,
       orthogonal.y        = orthogonal.y,
       auto.fix.first      = auto.fix.first,
       auto.fix.single     = auto.fix.single,
       auto.split.syntax   = auto.split.syntax,
       ...)

    return(out)
  }

  if (is.null(data)) {
    stop2("No data provided")
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (rcs) { # use reliability-correct single items?
    corrected <- relcorr_single_item(
      syntax          = model.syntax,
      data            = data,
      group           = group,
      choose          = rcs.choose,
      scale.corrected = rcs.scale.corrected,
      warn.lav        = FALSE
    )

    model.syntax <- corrected$syntax
    data         <- corrected$data
  }

  if ("convergence" %in% names(list(...))) {
    convergence.rel <- list(...)$convergence
    warning2("Argument 'convergence' is deprecated, use 'convergence.rel' instead.")
  }

  args <-
    getMethodSettingsDA(method,
      args =
        list(
          verbose            = verbose,
          optimize           = optimize,
          nodes              = nodes,
          convergence.abs    = convergence.abs,
          convergence.rel    = convergence.rel,
          optimizer          = optimizer,
          center.data        = center.data,
          standardize.data   = standardize.data,
          standardize.out    = standardize.out,
          standardize        = standardize,
          mean.observed      = mean.observed,
          double             = double,
          calc.se            = calc.se,
          FIM                = FIM,
          EFIM.S             = EFIM.S,
          OFIM.hessian       = OFIM.hessian,
          EFIM.parametric    = EFIM.parametric,
          robust.se          = robust.se,
          R.max              = R.max,
          max.iter           = max.iter,
          max.step           = max.step,
          epsilon            = epsilon,
          quad.range         = quad.range,
          adaptive.quad      = adaptive.quad,
          adaptive.frequency = adaptive.frequency,
          adaptive.quad.tol  = adaptive.quad.tol,
          n.threads          = n.threads,
          algorithm          = algorithm,
          em.control         = em.control,
          missing            = missing,
          orthogonal.x       = orthogonal.x,
          orthogonal.y       = orthogonal.y,
          auto.fix.first     = auto.fix.first,
          auto.fix.single    = auto.fix.single,
          auto.split.syntax  = auto.split.syntax,
          cr1s               = cr1s,
          group              = group
        )
    )

  cont.cols <- setdiff(colnames(data), c(cluster, group))

  if (args$center.data)
    data[cont.cols] <- lapply(data[cont.cols], FUN = centerIfNumeric, scaleFactor = FALSE)

  if (args$standardize.data)
    data[cont.cols] <- lapply(data[cont.cols], FUN = scaleIfNumeric, scaleFactor = FALSE)

  group.info <- getGroupInfo(
    model.syntax       = model.syntax,
    cov.syntax         = cov.syntax,
    data               = data,
    group              = group,
    auto.split.syntax  = args$auto.split.syntax
  )

  stopif(!method %in% c("lms", "qml"), "Method must be either 'lms' or 'qml'")

  model <- specifyModelDA(
    group.info         = group.info,
    method             = method,
    m                  = args$nodes,
    mean.observed      = args$mean.observed,
    double             = args$double,
    quad.range         = args$quad.range,
    adaptive.quad      = args$adaptive.quad,
    adaptive.frequency = args$adaptive.frequency,
    missing            = args$missing,
    orthogonal.x       = args$orthogonal.x,
    orthogonal.y       = args$orthogonal.y,
    auto.fix.first     = args$auto.fix.first,
    auto.fix.single    = args$auto.fix.single,
    cluster            = cluster
  )

  if (args$optimize) {
    model <- tryCatch({
      .optimize <- purrr::quietly(optimizeStartingParamsDA)
      #.optimize <- optimizeStartingParamsDA
      result    <- .optimize(model, args = args, group = group, engine = "sam")
      warnings  <- result$warnings

      if (length(warnings)) {
        fwarnings <- paste0(
          paste0(seq_along(warnings), ". ", warnings),
          collapse = "\n"
        )

        warning2("warning when optimizing starting parameters:\n", fwarnings)
      }

      result$result

    }, error = function(e) {
      warning2("unable to optimize starting parameters:\n", e)
      model
    })
  }

  if (!is.null(start)) {
    checkStartingParams(start, model = model) # throws error if somethings wrong
    model$theta <- start
  }

  # We want to limit the number of threads available to OpenBLAS.
  # Depending on the OpenBLAS version, it might not be compatible with
  # OpenMP. If `n.blas > 1L` you might end up getting this message:
  #> OpenBLAS Warning : Detect OpenMP Loop and this application may hang.
  #>                    Please rebuild the library with USE_OPENMP=1 option.
  # We don't want to restrict OpenBLAS in any other setttings in other settings,
  # e.g., lavaan::sem, so we reset after the model has been estimated.
  setThreads(n = args$n.threads, n.blas = 1L)
  on.exit(resetThreads()) # clean up at end of function

  est <- tryCatch(switch(method,
    qml = estQml(model,
      verbose         = args$verbose,
      convergence     = args$convergence.rel,
      calc.se         = args$calc.se,
      FIM             = args$FIM,
      EFIM.S          = args$EFIM.S,
      OFIM.hessian    = args$OFIM.hessian,
      EFIM.parametric = args$EFIM.parametric,
      robust.se       = args$robust.se,
      max.iter        = args$max.iter,
      epsilon         = args$epsilon,
      optimizer       = args$optimizer,
      R.max           = args$R.max,
      cr1s            = args$cr1s,
      ...
    ),
    lms = emLms(model,
      verbose           = args$verbose,
      convergence.abs   = args$convergence.abs,
      convergence.rel   = args$convergence.rel,
      calc.se           = args$calc.se,
      FIM               = args$FIM,
      EFIM.S            = args$EFIM.S,
      OFIM.hessian      = args$OFIM.hessian,
      EFIM.parametric   = args$EFIM.parametric,
      robust.se         = args$robust.se,
      max.iter          = args$max.iter,
      max.step          = args$max.step,
      epsilon           = args$epsilon,
      optimizer         = args$optimizer,
      R.max             = args$R.max,
      em.control        = args$em.control,
      algorithm         = args$algorithm,
      adaptive.quad     = args$adaptive.quad,
      quad.range        = args$quad.range,
      adaptive.quad.tol = args$adaptive.quad.tol,
      nodes             = args$nodes,
      cr1s              = args$cr1s,
      ...
  )),
  error = function(e) {
    if (args$verbose) cat("\n")
    message <- paste0("modsem [%s]: Model estimation failed!\n",
                      "Message: %s")
    stop2(sprintf(message, method, e$message))
  })

  # Finalize the model object
  # Expected means and covariances
  est$expected.matrices <- tryCatch(
    calcExpectedMatricesDA(
      parTable = est$parTable,
      xis  = getXisModelDA(model$models[[1L]]), # taking both the main model and cov model into account
      etas = getEtasModelDA(model$models[[1L]])  # taking both the main model and cov model into account
    ),
    error = function(e) {
      warning2("Failed to calculate expected matrices: ", e$message)
      NULL
    })

  # Arguments
  est$args <- args
  class(est) <- c("modsem_da", "modsem")

  # Check the results
  postCheckModel(est)

  # Return
  if (args$standardize.out) standardize_model(est) else est
}
