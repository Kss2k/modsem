#' Estimate interaction effects in structural equation models (SEMs) using a twostep procedure
#'
#' @param model.syntax \code{lavaan} syntax
#'
#' @param data dataframe
#'
#' @param method method to use:
#' \describe{
#'   \item{\code{"dblcent"}}{double centering approach (passed to \code{lavaan}).}
#'   \item{\code{"ca"}}{constrained approach (passed to \code{lavaan}).}
#'   \item{\code{"rca"}}{residual centering approach (passed to \code{lavaan}).}
#'   \item{\code{"uca"}}{unconstrained approach (passed to \code{lavaan}).}
#'   \item{\code{"pind"}}{prod ind approach, with no constraints or centering (passed to \code{lavaan}).}
#'   \item{\code{"lms"}}{latent moderated structural equations (not passed to \code{lavaan}).}
#'   \item{\code{"qml"}}{quasi maximum likelihood estimation (not passed to \code{lavaan}).}
#'   \item{\code{"custom"}}{use parameters specified in the function call (passed to \code{lavaan}).}
#' }
#'
#' @param ... arguments passed to other functions depending on the method (see \code{\link{modsem_pi}} and \code{\link{modsem_da}})
#'
#' @return \code{modsem} object with class \code{\link{modsem_pi}} or \code{\link{modsem_da}}.
#'
#' @export
#' @description
#' Estimate an interaction model using a twostep procedure. For the PI approaches, the \code{lavaan::sam} function
#'   is used to optimize the models, instead of \code{lavaan::sem}. Note that the product indicators are still used,
#'   and not the newly developed SAM approach to estimate latent interactions. For the DA approaches (LMS and QML)
#'   the measurement model is estimated using a CFA (\code{lavaan::cfa}). The structural model is estimated using
#'   \code{\link{modsem_da}}, where the estimates in the measurement model are fixed, based on the CFA estimates.
#'   Note that standard errors are uncorrected (i.e., naive), and do not account for the uncertainty in the CFA estimates.
#'   \strong{NOTE, this is an experimental feature!}
#'
#' @examples
#' library(modsem)
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
#' est_dblcent <- twostep(m1, oneInt, method = "dblcent")
#' summary(est_dblcent)
#'
#' \dontrun{
#' est_lms <- twostep(m1, oneInt, method = "lms")
#' summary(est_lms)
#'
#' est_qml <- twostep(m1, oneInt, method = "qml")
#' summary(est_qml)
#' }
#'
#' tpb_uk <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'  ATT =~ att3 + att2 + att1 + att4
#'  SN =~ sn4 + sn2 + sn3 + sn1
#'  PBC =~ pbc2 + pbc1 + pbc3 + pbc4
#'  INT =~ int2 + int1 + int3 + int4
#'  BEH =~ beh3 + beh2 + beh1 + beh4
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'  # Causal Relationsships
#'  INT ~ ATT + SN + PBC
#'  BEH ~ INT + PBC
#'  BEH ~ INT:PBC
#' "
#'
#' uk_dblcent <- twostep(tpb_uk, TPB_UK, method = "dblcent")
#' summary(uk_dblcent)
#'
#' \dontrun{
#'   uk_qml <- twostep(tpb_uk, TPB_UK, method = "qml")
#'
#'   uk_lms <- twostep(tpb_uk, TPB_UK, method = "lms", nodes = 32, adaptive.quad = TRUE)
#'   summary(uk_lms)
#' }
#' @export
twostep <- function(model.syntax, data, method = "lms", ...) {
  if (method %in% PI_METHODS)      twostepfun <- twostepPI
  else if (method %in% DA_METHODS) twostepfun <- twostepDA
  else stop2("Unsupported method: ", method)

  twostepfun(
     model.syntax = model.syntax,
     data         = data,
     method       = method,
     ...
  )
}


twostepPI <- function(..., LAVFUN = NULL) { # capture LAVFUN arg
  modsem_pi(..., LAVFUN = lavaan::sam)
}


twostepDA <- function(model.syntax, data, method = "lms", zero.tol = 1e-12,
                      fix.cov.xis = TRUE, ...) {
  stopif(!method %in% DA_METHODS, "Unsupported method: ", method)

  data <- as.data.frame(data)
  parTable      <- modsemify(model.syntax)
  parTableOuter <- parTable[parTable$op == "=~", ]
  parTableInner <- parTable[parTable$op != "=~", ]

  etas <- getEtas(parTable, isLV = FALSE)
  lvs  <- getLVs(parTable)
  syntaxCFA <- parTableToSyntax(parTableOuter)
  cfa       <- lavaan::cfa(syntaxCFA, data = data, meanstructure = TRUE)

  parTableCFA <- lavaan::parameterEstimates(cfa)

  fixCovVars <- if (fix.cov.xis) etas else lvs
  # Remove (co-) variances for etas (and optionally xis)
  parTableCFA <- parTableCFA[!(parTableCFA$op == "~~" &
                               parTableCFA$lhs %in% fixCovVars |
                               parTableCFA$rhs %in% fixCovVars), ]

  # Remove latent intercepts which are zero
  parTableCFA <- parTableCFA[!(parTableCFA$op == "~1" &
                               parTableCFA$lhs %in% lvs &
                               abs(parTableCFA$est) < zero.tol), ]

  parTableOuterFixed <- data.frame(
    lhs = parTableCFA$lhs,
    op = parTableCFA$op,
    rhs = parTableCFA$rhs,
    mod = as.character(parTableCFA$est)
  )

  if ("label" %in% names(parTableCFA)) {
    labelledParams <- parTableCFA[parTableCFA$label != "", ]
    parTableCustom <- data.frame(
      lhs = labelledParams$label,
      op = ":=",
      rhs = as.character(labelledParams$est),
      mod = ""
    )
  } else parTableCustom <- NULL

  isNonLinEta <- vapply(
    X = etas,
    FUN.VALUE = logical(1L),
    FUN = intTermsAffectLV,
    parTable = parTable
  )

  if (any(isNonLinEta)) {
    etaIntercepts <- data.frame(
      lhs = etas[isNonLinEta],
      op = "~1",
      rhs = "",
      mod = ""
    )
  } else etaIntercepts <- NULL

  parTableTwoStep <- rbind(
    parTableCustom,
    parTableOuterFixed,
    parTableInner,
    etaIntercepts
  )

  syntaxTwoStep <- parTableToSyntax(parTableTwoStep)

  twostep.spec <- list(
    syntax        = syntaxTwoStep,
    cfa           = cfa,
    parTable      = parTableTwoStep,
    parTableOuter = parTableOuterFixed,
    parTableInner = parTableInner
  )

  fit <- modsem_da(
    model.syntax = syntaxTwoStep,
    data = data,
    method = method,
    ...
  )

  # prep partable from CFA
  parTableCFA <- rename(parTableCFA, se = "se.cfa")
  parTableCFA <- parTableCFA[c("lhs", "op", "rhs", "se.cfa")]
  parTableCFA$se.cfa[parTableCFA$se.cfa <= zero.tol] <- NA
  estParTable <- parameter_estimates(fit)

  # Add reversed order covariances, to ensure matching
  covrows <- parTableCFA$op == "~~" & parTableCFA$rhs != parTableCFA$lhs
  parTableCFACov     <- parTableCFA
  parTableCFACov$rhs <- parTableCFA$lhs
  parTableCFACov$lhs <- parTableCFA$rhs
  parTableCFA <- rbind(
    parTableCFA, parTableCFACov[covrows, ]
  )

  parTable.out <- leftJoin(left  = estParTable, right = parTableCFA)

  parTable.out$std.error <- ifelse(
    is.na(parTable.out$std.error) & !is.na(parTable.out$se.cfa),
    yes = parTable.out$se.cfa,
    no  = parTable.out$std.error
  )

  parTable.out$z.value  <- parTable.out$est / parTable.out$std.error
  parTable.out$p.value  <- 2 * stats::pnorm(-abs(parTable.out$z.value))
  parTable.out$ci.lower <- parTable.out$est - CI_WIDTH * parTable.out$std.error
  parTable.out$ci.upper <- parTable.out$est + CI_WIDTH * parTable.out$std.error

  fit$parTable <- parTable.out

  parTabelLabelled <- getMissingLabels(parTable.out)

  # Get vcov
  vcov.11 <- lavaan::vcov(cfa)
  vcov.22 <- vcov(fit, type = "all")

  diff.pars <- setdiff(colnames(vcov.11), colnames(vcov.22))
  vcov.11   <- vcov.11[diff.pars, diff.pars]
  vcov.all  <- diagPartitionedMat(X = vcov.11, Y = vcov.22)

  # Get coef
  coef.1 <- lavaan::coef(cfa)
  coef.2 <- coef(fit, type = "all")

  diff.pars <- setdiff(names(coef.1), names(coef.2))
  coef.1    <- coef.1[diff.pars]
  coef.all  <- c(coef.1, coef.2)

  # Add to fit
  fit$vcov.all  <- vcov.all
  fit$coefs.all <- coef.all

  fit$twostep.spec <- twostep.spec

  fit
}
