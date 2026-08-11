#' Interaction between latent variables using LMS and QML approaches
#'
#' @param model.syntax \code{lavaan} syntax
#'
#' @param data A dataframe with observed variables used in the model.
#'
#' @param group Character. A variable name in the data frame defining the groups in a multiple
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
#' @param nodes Number of points of integration used in \code{lms}. More points
#'   give better estimates but slower computation; how many are needed depends
#'   on the complexity of the model.
#'
#'   \strong{What \code{nodes} counts depends on \code{integration}}: for
#'   \code{"gh"} and \code{"rect"} it is the number of nodes \emph{per latent
#'   dimension}, so a product rule costs \code{nodes^k} points in \code{k}
#'   dimensions. For \code{"mc"} it is the \emph{total} number of draws, whose
#'   cost does not grow with dimension.
#'
#'   If left unset, the default is chosen from \code{integration} and
#'   \code{adaptive}:
#'
#'   \tabular{lrrrr}{
#'     \strong{integration} \tab \strong{backend} \tab \code{adaptive = "none"} \tab \code{"quasi"} \tab \code{"full"} \cr
#'     \code{"gh"}   \tab legacy \tab  24 \tab  24 \tab   5 \cr
#'     \code{"gh"}   \tab graph  \tab  15 \tab  15 \tab   5 \cr
#'     \code{"rect"} \tab legacy \tab  24 \tab  24 \tab   5 \cr
#'     \code{"rect"} \tab graph  \tab  15 \tab  15 \tab   5 \cr
#'     \code{"mc"}   \tab both   \tab 500 \tab   - \tab 250 \cr
#'   }
#'
#'   The legacy backend integrates only the nonlinear latent dimensions, so its
#'   product rule is small and extra nodes cost little -- and are measurably
#'   needed. On \code{oneInt}, a 15-node quasi-adaptive rule leaves the
#'   interaction 0.005 from the per-observation AGHQ reference, while 20 and
#'   above agree to four decimals. The graph backend integrates every latent
#'   variable, so its product rule grows much faster and 15 is the practical
#'   ceiling.
#'
#'   A per-observation rule (\code{adaptive = "full"}) needs far fewer points
#'   than a shared one because it is centred on each observation's own
#'   posterior. For models with an interaction between an endogenous and an
#'   exogenous variable, a shared rule may need considerably more than the
#'   default; where the data is non-normal, the \code{qml} approach may be a
#'   better choice than raising \code{nodes}.
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
#' @param quad.range range in z-scores over which to perform numerical
#'   integration in LMS. By default \code{Inf}, so that \code{f(t)} is
#'   integrated from -Inf to Inf, but this is likely inefficient at a large
#'   number of nodes. Nodes outside \code{+/- quad.range} are ignored. Applies
#'   to \code{integration = "gh"}; the rectangular rule has its own range.
#'
#' @param rect.range Half-width of the integration range for
#'   \code{integration = "rect"}, in z-scores. Defaults to 5, matching Mplus.
#'
#'   This matters as much as \code{nodes}. A trapezoid rule on a
#'   Gaussian-decaying integrand is not the \code{O(h^2)} scheme its name
#'   suggests: Poisson summation gives a discretisation error of
#'   \code{2 * exp(-2 * pi^2 / h^2)} for spacing \code{h}, so it converges
#'   exponentially. Accuracy is governed by whichever is larger, that term or
#'   the truncation error \code{2 * (1 - pnorm(rect.range))}. At 15 nodes over
#'   \code{+/- 5} the discretisation error is around \code{3e-17} while
#'   truncation is around \code{6e-7}, so the rule is truncation-limited by ten
#'   orders of magnitude and widening the range helps far more than adding
#'   nodes. Ignored by \code{"gh"} and \code{"mc"}, which use
#'   \code{quad.range}.
#'
#' @param integration Which numerical integration rule to use for \code{lms}.
#' \describe{
#'   \item{\code{"gh"}}{Gauss-Hermite product rule (the default). Highest
#'     accuracy per node for smooth posteriors.}
#'   \item{\code{"rect"}}{Equispaced rectangular product rule over \code{+/- 5}
#'     per dimension, weighted by the normal density. This is the rule Mplus
#'     uses under \code{INTEGRATION = STANDARD}, and is the closest match to
#'     Mplus output.}
#'   \item{\code{"mc"}}{Monte-Carlo integration from a fixed, seeded set of
#'     standard-normal draws. The only rule whose cost does not grow with the
#'     number of latent dimensions. Draws are held fixed across EM iterations,
#'     so the likelihood stays a deterministic function of the parameters.}
#' }
#'
#' @param adaptive How the integration rule is placed. Independent of
#'   \code{integration}: the chosen rule is used in every case, only its
#'   location changes.
#' \describe{
#'   \item{\code{"quasi"}}{(default) One rule per group, adapted to the
#'     aggregate posterior. For \code{"gh"} the rule is pruned one dimension at
#'     a time, dropping nodes that carry negligible information; for
#'     \code{"rect"} the range is shrunk onto where the mass actually lies and
#'     the same number of points respaced inside it. A shared rule permits a
#'     sufficient-statistics M-step, which is much faster than the alternative.}
#'   \item{\code{"none"}}{One fixed rule, shared by every observation.}
#'   \item{\code{"full"}}{One rule per observation, centred on that
#'     observation's posterior mode and scaled by its Hessian. Far more
#'     accurate per node -- hence the much lower default \code{nodes} -- but it
#'     rules out the sufficient-statistics M-step, so each iteration costs
#'     \code{O(N * nodes^k)}. Graph backend only.}
#' }
#'
#'   Not every combination is available:
#'
#'   \tabular{lll}{
#'     \strong{backend} \tab \strong{gh / rect} \tab \strong{mc} \cr
#'     \code{"legacy"} \tab \code{"none"}, \code{"quasi"} \tab - \cr
#'     \code{"graph"}  \tab \code{"none"}, \code{"quasi"}, \code{"full"} \tab \code{"none"}, \code{"full"} \cr
#'   }
#'
#' @param adaptive.quad \strong{Deprecated}; use \code{adaptive} instead.
#'   \code{TRUE} maps to \code{adaptive = "quasi"} and \code{FALSE} to
#'   \code{adaptive = "none"}.
#'
#' @param adaptive.frequency How often should the adaptive rule be recalculated?
#'   Defaults to 3, meaning every third EM-iteration. Applies to both
#'   \code{adaptive = "quasi"} (when the rule is re-pruned) and
#'   \code{adaptive = "full"} (when the per-observation modes are refitted).
#'
#' @param adaptive.quad.tol Relative error tolerance used when pruning or
#'   trimming a quasi-adaptive rule. Defaults to \code{1e-12}.
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
#' @param ordered Variables to be treated as ordered. Ordered indicators are handled
#'   with a Monte-Carlo correction for the LMS/QML estimator on standardized category
#'   scores. The fitted model is used to repeatedly simulate continuous indicators,
#'   ordinalize them to the observed cumulative category proportions, refit the
#'   standardized LMS/QML working model, and solve the resulting fixed-point problem.
#'   This follows the same general Monte-Carlo consistency logic as the MC-OrdPLSc
#'   algorithm described by Slupphaug, Mehmetoglu, and Mittner (2026,
#'   \doi{10.31234/osf.io/fwzj6_v1}).
#'
#'   Threshold standard errors are computed using a simple bootstrap of the category counts
#'   (see \code{ordered.boot.reps}), which does not propagate uncertainty from the
#'   model parameter estimates.
#'
#' @param ordered.mc.reps Integer. Monte-Carlo sample size used in each ordered MC
#'   correction step. Larger values reduce simulation noise but increase runtime.
#'
#' @param ordered.min.iter Integer. Minimum number of Robbins-Monro iterations for the
#'   ordered MC correction.
#'
#' @param ordered.max.iter Integer. Maximum number of Robbins-Monro iterations for the
#'   ordered MC correction.
#'
#' @param ordered.tol Convergence tolerance for the ordered MC Robbins-Monro updates.
#'
#' @param ordered.rng.seed Optional integer random seed used by the ordered MC
#'   correction.
#'
#' @param ordered.fixed.seed Logical. If \code{TRUE} and \code{ordered.rng.seed = NULL},
#'   a fixed seed is drawn once and reused throughout the ordered MC correction to
#'   improve numerical stability.
#'
#' @param ordered.polyak.juditsky Logical. Should Polyak-Juditsky averaging be used in
#'   the ordered MC Robbins-Monro solver?
#'
#' @param ordered.pj.extrapolate Logical. If \code{TRUE}, use extrapolation of the
#'   Polyak-Juditsky path to estimate the convergence point. If \code{FALSE}, the
#'   averaged iterate is used directly.
#'
#' @param ordered.se Character string selecting the ordered MC standard-error correction.
#'   \code{"delta"} (default) uses the delta method for all free parameters.
#'   \code{"penalized"} uses a conservative variance inflation
#'   based on the discrepancy between the naive and MC-corrected standardized estimates.
#'   \code{"naive"} uses the fast diagonal rescaling approximation.
#'   \code{"mixed"} uses the delta method for the structural path coefficients
#'   only, and penalized standard errors for the remaining parameters. This is
#'   considerably faster, and is a good option if you're only interested in
#'   the structural model.
#'
#' @param ordered.se.penalty Non-negative numeric multiplier used when
#'   \code{ordered.se = "penalized"}. The penalty adds
#'   \code{ordered.se.penalty * (theta_mc - theta_naive)^2} to the diagonal of the
#'   naive covariance matrix on the variance scale.
#'
#' @param ordered.delta.reps Integer. Monte-Carlo sample size used when approximating
#'   the ordered MC delta-method Jacobian. Only relevant if
#'   \code{ordered.se} is \code{"delta"} or \code{"mixed"}.
#'
#' @param ordered.delta.epsilon Finite-difference step size used for the ordered MC
#'   delta-method Jacobian. Only relevant if \code{ordered.se} is \code{"delta"}
#'   or \code{"mixed"}.
#'
#' @param ordered.boot.reps Integer. Number of bootstrap replications used to compute
#'   standard errors for the thresholds of ordered indicators. The bootstrap resamples
#'   the category counts at the observed sample size, holding the simulated
#'   (model-implied) reference distribution of the continuous indicators fixed. Note
#'   that this does not account for the uncertainty in the model parameter estimates,
#'   and should be seen as a rough approximation. Set to \code{0} to disable.
#'
#' @param ordered.standardize Logical. Should scored ordered indicators be standardized
#'   before the observed-data fit and after ordinalizing simulated indicators? This is
#'   recommended for numerical stability and is enabled by default.
#'
#' @param lms.backend LMS likelihood implementation. The experimental
#'   \code{"graph"} backend integrates one innovation per latent variable and
#'   evaluates latent variables in structural order, which lets it support
#'   n-way interactions and ordered indicators, and is what makes
#'   \code{adaptive = "full"} and \code{integration = "mc"} available. It
#'   requires conditionally independent indicators and does not support
#'   composites. The default \code{"legacy"} backend integrates only the
#'   nonlinear latent dimensions, so its product rule is much smaller.
#'   Selected automatically when the data has ordered indicators.
#'
#' @param link Link for ordered indicators in the LMS graph backend. The default
#'   is \code{"logit"}; \code{"probit"} is also available.
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
#' @param sampling.weights A variable name in the data frame containing sampling weight information.
#'   Depending on the sampling.weights.normalization argument, these weights may be rescaled (or not)
#'   so that their sum equals the number of observations (total or per group)
#'
#' @param sampling.weights.normalization If \code{"none"}, the sampling weights (if provided) will not be
#'   transformed. If \code{"total"}, the sampling weights are normalized by dividing by the total sum
#'   of the weights, and multiplying again by the total sample size. If \code{"group"}, the sampling
#'   weights are normalized per group: by dividing by the sum of the weights (in each group), and
#'   multiplying again by the group size. The default is \code{"total"}.
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
#'   Default is \code{FALSE}.
#'
#' @param orthogonal.y If \code{TRUE}, all covariances among endogenous latent variables only are set to zero.
#'   If \code{FALSE} residual covariances are added between pure endogenous variables;
#'   those that are predicted by no other endogenous variable in the structural model.
#'   Default is \code{FALSE}.
#'
#' @param orthogonal.y If \code{TRUE}, all covariances among endogenous latent variables only are set to zero.
#'   If \code{FALSE} residual covariances are added between pure endogenous variables;
#'   those that are predicted by no other endogenous variable in the structural model.
#'   Default is \code{FALSE}.
#'
#' @param auto.fix.first If \code{TRUE} the factor loading of the first indicator, for
#'   a given latent variable is fixed to \code{1}. If \code{FALSE} no loadings are fixed
#'   (automatically). Note that that this might make it such that the model no longer is
#'   identified. Default is \code{TRUE}. \strong{NOTE} this behaviour is overridden
#'   if the first loading is labelled, where it gets treated as a free parameter instead. This
#'   differs from the default behaviour in \code{lavaan}.
#'
#' @param auto.fix.single If \code{TRUE}, the residual variance of
#'   an observed indicator is set to zero if it is the only indicator of a latent variable.
#'   If \code{FALSE} the residual variance is not fixed to zero, and treated as a free parameter
#'   of the model. Default is \code{TRUE}. \strong{NOTE} this behaviour is overridden
#'   if the first loading is labelled, where it gets treated as a free parameter instead.
#'
#' @param fix.composite.var If \code{TRUE} (default) the block covariance structure
#'   of composite indiactors is fixed.
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
#' Ordered indicators use the LMS graph backend with conditionally independent
#' univariate response likelihoods. Measurement intercepts and response scales
#' are fixed for identification and ordered thresholds are estimated directly.
#'
#' \strong{NOTE}: Run \code{\link{default_settings_da}} to see default arguments.
#'
#' @references
#' Slupphaug, K., Mehmetoglu, M., and Mittner, M. (2026, March 21).
#' \emph{Consistent Estimates from Biased Estimators: Monte-Carlo Consistent Partial
#' Least Squares for Latent Interaction Models with Ordinal Indicators}. PsyArXiv.
#' \doi{10.31234/osf.io/fwzj6_v1}
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
                      rect.range = NULL,
                      integration = NULL,
                      adaptive = NULL,
                      adaptive.quad = NULL, # deprecated, use `adaptive`
                      adaptive.frequency = NULL,
                      adaptive.quad.tol = NULL,
                      n.threads = NULL,
                      algorithm = NULL,
                      em.control = NULL,
                      lms.backend = c("legacy", "graph"),
                      link = c("logit", "probit"),
                      ordered = NULL,
                      ordered.mc.reps = NULL,
                      ordered.min.iter = 20L,
                      ordered.max.iter = 250L,
                      ordered.tol = 1e-4,
                      ordered.rng.seed = NULL,
                      ordered.fixed.seed = FALSE,
                      ordered.polyak.juditsky = TRUE,
                      ordered.pj.extrapolate = TRUE,
                      ordered.se = c("delta", "penalized", "naive", "mixed"),
                      ordered.se.penalty = 0.25,
                      ordered.delta.reps = NULL,
                      ordered.delta.epsilon = 1e-2,
                      ordered.boot.reps = 1000L,
                      ordered.standardize = TRUE,
                      cluster = NULL,
                      cr1s = FALSE,
                      sampling.weights = NULL,
                      sampling.weights.normalization = NULL,
                      rcs = FALSE,
                      rcs.choose = NULL,
                      rcs.scale.corrected = TRUE,
                      orthogonal.x = NULL,
                      orthogonal.y = NULL,
                      auto.fix.first = NULL,
                      auto.fix.single = NULL,
                      auto.split.syntax = NULL,
                      fix.composite.var = NULL,
                      ...) {
  method <- tolower(method)
  lms.backend <- match.arg(lms.backend)
  link <- match.arg(link)

  if (!is.null(adaptive.quad)) {
    mapped <- if (isTRUE(adaptive.quad)) "quasi" else "none"
    mod_msg_warn_immediate(sprintf(
      "`adaptive.quad` is deprecated; use `adaptive = \"%s\"` instead.", mapped))
    if (is.null(adaptive)) adaptive <- mapped
  }

  if (is.null(model.syntax)) {
    mod_msg_stop("No model.syntax provided")
  } else if (!is.character(model.syntax)) {
    mod_msg_stop("The provided model syntax is not a string!")
  } else if (length(model.syntax) > 1) {
    mod_msg_stop("The provided model syntax is not of length 1")
  }

  ordered.se <- match.arg(ordered.se)

  has.ordered <- length(ordered) || any(vapply(data, is.ordered, logical(1L)))
  if (has.ordered && method == "lms") lms.backend <- "graph"

  if (has.ordered && method != "lms") {
    out <- modsemOrderedMCCorrection(
       model.syntax                   = model.syntax,
       data                           = data,
       method                         = method,
       verbose                        = verbose,
       optimize                       = optimize,
       nodes                          = nodes,
       missing                        = missing,
       convergence.abs                = convergence.abs,
       convergence.rel                = convergence.rel,
       optimizer                      = optimizer,
       center.data                    = center.data,
       standardize.data               = standardize.data,
       standardize.out                = standardize.out,
       standardize                    = standardize,
       mean.observed                  = mean.observed,
       cov.syntax                     = cov.syntax,
       double                         = double,
       calc.se                        = calc.se,
       FIM                            = FIM,
       EFIM.S                         = EFIM.S,
       OFIM.hessian                   = OFIM.hessian,
       EFIM.parametric                = EFIM.parametric,
       robust.se                      = robust.se,
       R.max                          = R.max,
       max.iter                       = max.iter,
       max.step                       = max.step,
       start                          = start,
       epsilon                        = epsilon,
       quad.range                     = quad.range,
       rect.range                     = rect.range,
       integration                    = integration,
       adaptive                       = adaptive,
       adaptive.frequency             = adaptive.frequency,
       adaptive.quad.tol              = adaptive.quad.tol,
       n.threads                      = n.threads,
       algorithm                      = algorithm,
       em.control                     = em.control,
       ordered                        = ordered,
       ordered.mc.reps                = ordered.mc.reps,
       ordered.min.iter               = ordered.min.iter,
       ordered.max.iter               = ordered.max.iter,
       ordered.tol                    = ordered.tol,
       ordered.rng.seed               = ordered.rng.seed,
       ordered.fixed.seed             = ordered.fixed.seed,
       ordered.polyak.juditsky        = ordered.polyak.juditsky,
       ordered.pj.extrapolate         = ordered.pj.extrapolate,
       ordered.se                     = ordered.se,
       ordered.se.penalty             = ordered.se.penalty,
       ordered.delta.reps             = ordered.delta.reps,
       ordered.delta.epsilon          = ordered.delta.epsilon,
       ordered.boot.reps              = ordered.boot.reps,
       ordered.standardize            = ordered.standardize,
       cluster                        = cluster,
       group                          = group,
       cr1s                           = cr1s,
       sampling.weights               = sampling.weights,
       sampling.weights.normalization = sampling.weights.normalization,
       rcs                            = rcs,
       rcs.choose                     = rcs.choose,
       rcs.scale.corrected            = rcs.scale.corrected,
       orthogonal.x                   = orthogonal.x,
       orthogonal.y                   = orthogonal.y,
       auto.fix.first                 = auto.fix.first,
       auto.fix.single                = auto.fix.single,
       auto.split.syntax              = auto.split.syntax,
       fix.composite.var              = fix.composite.var,
       ...)

    return(out)
  }

  if (is.null(data)) {
    mod_msg_stop("No data provided")
  } else {
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
    mod_msg_warn("Argument 'convergence' is deprecated, use 'convergence.rel' instead.")
  }

  args <-
    getMethodSettingsDA(method,
      args =
        list(
          verbose                        = verbose,
          optimize                       = optimize,
          nodes                          = nodes,
          convergence.abs                = convergence.abs,
          convergence.rel                = convergence.rel,
          optimizer                      = optimizer,
          center.data                    = center.data,
          standardize.data               = standardize.data,
          standardize.out                = standardize.out,
          standardize                    = standardize,
          mean.observed                  = mean.observed,
          double                         = double,
          calc.se                        = calc.se,
          FIM                            = FIM,
          EFIM.S                         = EFIM.S,
          OFIM.hessian                   = OFIM.hessian,
          EFIM.parametric                = EFIM.parametric,
          robust.se                      = robust.se,
          R.max                          = R.max,
          max.iter                       = max.iter,
          max.step                       = max.step,
          epsilon                        = epsilon,
          quad.range                     = quad.range,
          rect.range                     = rect.range,
          lms.backend                    = lms.backend,
          integration                    = integration,
          adaptive                       = adaptive,
          adaptive.frequency             = adaptive.frequency,
          adaptive.quad.tol              = adaptive.quad.tol,
          n.threads                      = n.threads,
          algorithm                      = algorithm,
          em.control                     = em.control,
          missing                        = missing,
          orthogonal.x                   = orthogonal.x,
          orthogonal.y                   = orthogonal.y,
          auto.fix.first                 = auto.fix.first,
          auto.fix.single                = auto.fix.single,
          auto.split.syntax              = auto.split.syntax,
          cr1s                           = cr1s,
          group                          = group,
          sampling.weights               = sampling.weights,
          sampling.weights.normalization = sampling.weights.normalization,
          fix.composite.var              = fix.composite.var
        )
    )
  quad.spec <- quadSpec(
    integration        = args$integration,
    adaptive           = args$adaptive,
    nodes              = args$nodes,
    quad.range         = args$quad.range,
    rect.range         = args$rect.range,
    adaptive.frequency = args$adaptive.frequency,
    adaptive.quad.tol  = args$adaptive.quad.tol
  )
  if (method == "lms") checkQuadAvailable(quad.spec, backend = lms.backend)

  # Only keys that are also `modsem_da()` arguments may live on `args`: it is
  # replayed through `do.call(modsem_da, ...)` by `estimate_h0()`, so anything
  # else surfaces there as an unused-argument error.
  args$ordered <- if (has.ordered) lmsGraphOrderedColumns(data, model.syntax, ordered) else NULL

  cont.cols <- setdiff(colnames(data), c(cluster, group, sampling.weights))

  if (args$center.data)
    data[cont.cols] <- lapply(data[cont.cols], FUN = centerIfNumeric, scaleFactor = FALSE)

  if (args$standardize.data)
    data[cont.cols] <- lapply(data[cont.cols], FUN = scaleIfNumeric, scaleFactor = FALSE)

  group.info <- parseModelArgumentsByGroupDA(
    model.syntax       = model.syntax,
    cov.syntax         = cov.syntax,
    method             = method,
    data               = data,
    group              = group,
    auto.split.syntax  = args$auto.split.syntax,
    sampling.weights   = sampling.weights,
    sampling.weights.normalization = args$sampling.weights.normalization
  )

  mod_stopif(!method %in% c("lms", "qml"), "Method must be either 'lms' or 'qml'")

  model <- specifyModelDA(
    group.info         = group.info,
    method             = method,
    m                  = args$nodes,
    mean.observed      = args$mean.observed,
    double             = args$double,
    quad               = quad.spec,
    missing            = args$missing,
    orthogonal.x       = args$orthogonal.x,
    orthogonal.y       = args$orthogonal.y,
    auto.fix.first     = args$auto.fix.first,
    auto.fix.single    = args$auto.fix.single,
    fix.composite.var  = args$fix.composite.var,
    cluster            = cluster,
    sampling.weights   = sampling.weights
  )

  if (method == "lms" && lms.backend == "graph")
    model <- lmsGraphPrepareProducts(model)

  if (method == "lms" && lms.backend == "graph" && has.ordered) {
    model <- lmsGraphPrepareOrdered(
      model = model, data = data, model.syntax = model.syntax,
      ordered = ordered, group = group, link = link
    )
  }

  # The graph backend integrates every latent variable, not just the nonlinear
  # ones, so its product rule is far larger than the legacy backend's.
  if (method == "lms" && lms.backend == "graph") {
    integration.dimension <- model$models[[1L]]$info$numXis +
      model$models[[1L]]$info$numEtas
    mod_warnif(integration.dimension > 5L && args$integration != "mc",
      paste0("The LMS graph model integrates over ", integration.dimension,
             " latent variables, so a product rule needs ", args$nodes, "^",
             integration.dimension, " points. Consider `adaptive = \"full\"` ",
             "or `integration = \"mc\"`; see `?modsem_da`."))
  }

  if (args$optimize) {
    model <- tryCatch({
      .optimize <- purrr::quietly(optimizeStartingParamsDA)
      #.optimize <- \(...) list(result = optimizeStartingParamsDA(...))

      ops <- c(group.info$parTable$op, group.info$parTableCov$op)
      engine <- if (any(ops %in% c("<", ">", "=="))) "pi" else "sam"

      result <- .optimize(
        model            = model,
        args             = args,
        group            = group,
        sampling.weights = sampling.weights,
        engine           = engine
      )

      warnings  <- result$warnings

      if (length(warnings)) {
        fwarnings <- paste0(
          paste0(seq_along(warnings), ". ", warnings),
          collapse = "\n"
        )

        mod_msg_warn_immediate(
          paste0("warning when optimizing starting parameters:\n", fwarnings)
        )
      }

      result$result

    }, error = function(e) {
      mod_msg_warn_immediate(
        paste0("unable to optimize starting parameters:\n", e)
      )

      if (is.null(max.step) && args$max.step <= 1 && method == "lms") {
        # When we don't have optimized starting parameters with LMS, the algorithm
        # will likely not converge with the default max.step
        mod_msg_note("Increasing max step in EM-algorithm to 50")
        args$max.step <<- 50
      }

      model
    })
  }

  # Starting-value optimization may rebuild the DA model, so reinstate the
  # requested integration spec afterwards. The two backends integrate over
  # different dimensions: the graph backend takes one innovation per latent
  # variable, the legacy backend only the nonlinear ones.
  if (method == "lms") for (g in seq_along(model$models)) {
    info <- model$models[[g]]$info
    k <- if (lms.backend == "graph") info$numXis + info$numEtas else
      model$models[[g]]$quad$k
    model$models[[g]]$quad <- utils::modifyList(quad.spec, list(k = k))
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
      quad.range        = args$quad.range,
      adaptive.quad.tol = args$adaptive.quad.tol,
      nodes             = args$nodes,
      cr1s              = args$cr1s,
      lms.backend       = lms.backend,
      link              = link,
      ...
  )),
  error = function(e) {
    if (args$verbose) cat("\n")
    message <- paste0("modsem [%s]: Model estimation failed!\n",
                      "Message: %s")
    mod_msg_stop(sprintf(message, method, e$message))
  })

  # Finalize the model object
  # Expected means and covariances
  expected.matrices <- tryCatch(
    calcExpectedMatricesDA(
      parTable = est$parTable,
      xis  = getXisModelDA(model$models[[1L]]), # taking both the main model and cov model into account
      etas = getEtasModelDA(model$models[[1L]])  # taking both the main model and cov model into account
    ),
    error = function(e) {
      mod_msg_warn(paste0("Failed to calculate expected matrices: ", e$message))
      NULL
    }
  )

  if (!is.null(expected.matrices)) {
    est$expected.matrices <- expected.matrices

    for (g in seq_along(est$model$models)) # attach to submodels
      est$model$models[[g]]$expected.matrices <- expected.matrices[[g]]
  }


  # Arguments
  est$args <- args
  attr(est, "isRCS_Model") <- rcs
  class(est) <- c("modsem_da", "modsem")

  # Check the results
  postCheckModel(est)

  # Return
  if (args$standardize.out) standardize_model(est) else est
}
