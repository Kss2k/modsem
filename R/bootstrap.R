#' Bootstrap a modsem Model
#'
#' A generic interface for parametric and non‑parametric bootstrap procedures
#' for structural equation models estimated with the \pkg{modsem} ecosystem.
#' The function dispatches on the class of \code{model}; currently dedicated
#' methods exist for \code{\link{modsem_pi}} (product‑indicator approach) and
#' \code{\link{modsem_da}} (distributional‑analytic approach).
#'
#' @param model A fitted \code{modsem} object.
#' @param FUN  A function that returns the statistic of interest when applied to
#'   a fitted model.  The function must accept a single argument, the model
#'   object, and should ideally return a numeric vector; see Value.
#' @param ...  Additional arguments forwarded to \code{lavaan::bootstrapLavaan} 
#'   for \code{\link{modsem_pi}} objects, or \code{\link{modsem_da}} for 
#'   \code{\link{modsem_da}} objects.
#'
#' @return Depending on the return type of \code{FUN} either
#'   \describe{
#'     \item{numeric}{A matrix with \code{R} rows (bootstrap replicates) and as
#'       many columns as \code{length(FUN(model))}.}
#'     \item{other}{A list of length \code{R}; each element is the raw output of
#'       \code{FUN}. \strong{NOTE}: Only applies for \code{\link{modsem_da}} objects}
#'   }
#' @seealso \code{\link[lavaan]{bootstrapLavaan}}, \code{\link{modsem_pi}},
#'   \code{\link{modsem_da}}
#' @export
bootstrap_modsem <- function(model, FUN, ...) {
  UseMethod("bootstrap_modsem")
}


#' @describeIn bootstrap_modsem Bootstrap a \code{modsem_pi} model by delegating
#'   to \code{\link[lavaan]{bootstrapLavaan}}.
#'
#' @details A thin wrapper around \code{lavaan::bootstrapLavaan()} that performs the
#'   necessary book‑keeping so that \code{FUN} receives a fully‑featured
#'   \code{modsem_pi} object—rather than a bare \code{lavaan} fit—at every
#'   iteration.
#'
#' @examples
#'
#' m1 <- '
#'   X =~ x1 + x2
#'   Z =~ z1 + z2
#'   Y =~ y1 + y2
#' 
#'   Y ~ X + Z + X:Z
#' '
#' 
#' fit_pi <- modsem(m1, oneInt)
#' bootstrap_modsem(fit_pi, FUN = coef, R = 10L)
#'
#' @export
bootstrap_modsem.modsem_pi <- function(model, FUN, ...) {
  wrapLavFit <- function(lavfit) {
    model$lavaan       <- lavfit
    model$coefParTable <- lavaan::parameterEstimates(lavfit)
    model$data         <- lavaan::lavInspect(lavfit, what = "data")
    model
  }

  fun.mod <- \(fit) FUN(wrapLavFit(fit))

  lavaan::bootstrapLavaan(extract_lavaan(model), FUN = fun.mod, ...)
}


#' @describeIn bootstrap_modsem Parametric or non‑parametric bootstrap for
#'   \code{modsem_da} models.
#'
#' @param R number of bootstrap replicates.
#' @param P.max ceiling for the simulated population size.
#' @param type Bootstrap flavour, see Details.
#' @param verbose Should progress information be printed to the console?
#' @param calc.se Should standard errors for each replicate. Defaults to \code{FALSE}.
#' @param optimize Should starting values be re-optimized for each replicate. Defaults to \code{FALSE}.
#'
#' @details The function internally resamples the observed data (non‑parametric
#'   case) or simulates from the estimated parameter table (parametric case),
#'   feeds the sample to \code{\link{modsem_da}}, evaluates \code{FUN} on the
#'   refitted object and finally collates the results.
#'
#' @examples
#'
#' m1 <- '
#'   X =~ x1 + x2
#'   Z =~ z1 + z2
#'   Y =~ y1 + y2
#' 
#'   Y ~ X + Z + X:Z
#' '
#' 
#' \dontrun{
#' fit_lms <- modsem(m1, oneInt, method = "lms")
#' bootstrap_modsem(fit_lms, FUN = coef, R = 10L)
#' }
#' @export
bootstrap_modsem.modsem_da <- function(model, 
                                       FUN = "coef", 
                                       R = 1000L, 
                                       P.max = 1e5,
                                       type = c("nonparametric", "parameteric"), 
                                       verbose = interactive(),
                                       calc.se = FALSE, 
                                       optimize = FALSE,
                                       ...) {
  type <- tolower(type)
  type <- match.arg(type)

  data     <- as.data.frame(modsem_inspect(model, what = "data"))
  ovs      <- colnames(data)
  N        <- nrow(data)
  P.ceil   <- P.max < N * R
  P        <- min(P.max, N * R)
  parTable <- parameter_estimates(model)

  warnif(P.max <= N, "`P.max` is less than `N`!")

  population <- switch(type,
    parameteric    = simulateDataParTable(parTable, N = P, colsOVs = ovs)$oV,
    nonparametric = data[sample(nrow(data), P, replace = TRUE), ],
    stop2("Unrecognized type!\n")
  )

  argList              <- model$args
  argList$calc.se      <- calc.se
  argList$verbose      <- verbose
  argList$model.syntax <- model$model$syntax
  argList$method       <- model$method

  if (!optimize) {
    argList$start        <- model$theta
    argList$optimize     <- FALSE
  } else {
    argList$start        <- NULL
    argList$optimize     <- TRUE 
  }

  f0 <- FUN(model)
  FUN.VEC <- is.vector(f0) || inherits(f0, "ModsemVector")

  if (FUN.VEC) {
    out <- matrix(NA, nrow = R, ncol = length(f0), 
                  dimnames = list(NULL, names(f0)))

  } else out <- vector("list", length = R)

  for (i in seq_len(R)) {
    printedLines <- utils::capture.output(split = TRUE, {
      if (verbose) printf("Bootstrap %d/%d...\n", i, R)

      sample_i  <- population[sample(P, N), ]
      argList_i <- c(argList, list(data = sample_i))

      fit_i <- tryCatch(
        do.call(modsem_da, args = argList_i, ...),
        error = \(e) {
          print(e)
          NULL
        }
      )

      if (is.null(fit_i)) next

      fi <- tryCatch(do.call(FUN, args = list(fit_i)),
                     error = \(e) NA)

      if (FUN.VEC) out[i, ] <- fi else out[[i]] <- fi
    })

    nprinted <- length(printedLines)
    if (i < R) eraseConsoleLines(nprinted)
  }

  out
}

