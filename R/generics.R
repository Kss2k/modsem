#' Extract parameterEstimates from an estimated model
#'
#' @param object An object of class \code{\link{modsem_pi}}, \code{\link{modsem_da}}, or \code{\link{modsem_mplus}}
#' @param ... Additional arguments passed to other functions
#'
#' @examples
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'
#'   # Inner Model
#'   Y ~ X + Z + X:Z
#' '
#' # Double centering approach
#' est_dca <- modsem(m1, oneInt)
#'
#' pars <- parameter_estimates(est_dca) # no correction
#' 
#' # Pretty summary
#' summarize_partable(pars)
#'
#' # Only print the data.frame
#' pars
#' @export
parameter_estimates <- function(object, ...) {
  UseMethod("parameter_estimates")
}


#' Extract or modify parTable from an estimated model with estimated variances of interaction terms
#'
#' @param object An object of class \code{\link{modsem_da}},  \code{\link{modsem_mplus}},
#' or a parTable of class \code{\link{data.frame}}
#' @param ... Additional arguments passed to other functions
#' @export
var_interactions <- function(object, ...) {
  UseMethod("var_interactions")
}


#' @export
var_interactions.data.frame <- function(object, ignore.means = FALSE, ...) {
  parTable <- removeInteractionVariances(fillColsParTable(object))

  intTerms <- unique(parTable[grepl(":", parTable$rhs) &
                     parTable$op == "~", "rhs"])

  getLabel <- \(x, y) sprintf("%s~~%s", x, y)

  for (intTerm in intTerms) {
    # interaction term = XZ
    # TODO:
    #   I should also add covariances between X:Z and the other exogenous
    #   variables... (only relevant when mu(X) or mu(Z) != 0)
    XZ   <- stringr::str_split_fixed(intTerm, ":", 2)
    X    <- XZ[[1]]
    Z    <- XZ[[2]]

    muX  <- if (!ignore.means) getMean(X, parTable) else 0
    muZ  <- if (!ignore.means) getMean(Z, parTable) else 0

    lhs  <- c(X, Z, X)
    rhs  <- c(X, Z, Z)

    covs  <- calcCovParTable(lhs, rhs, parTable)
    varX  <- covs[[getLabel(X, X)]]
    varZ  <- covs[[getLabel(Z, Z)]]
    covXZ <- covs[[getLabel(X, Z)]]

    varXZ <- varX * muZ ^ 2 + varZ * muX ^ 2 +
      2 * muX * muZ * covXZ + varX * varZ + covXZ ^ 2
    # needed if mu != 0, but not sure if this is completely correct
    # when X or Z is endogenous
    covX_XZ <- varX * muZ + muX * covXZ
    covZ_XZ <- varZ * muX + muZ * covXZ
    newRows <- data.frame(lhs = c(intTerm, XZ[[1]], XZ[[2]]),
                          op = rep("~~", 3),
                          rhs = rep(intTerm, 3),
                          label = "",
                          est = c(varXZ, covX_XZ, covZ_XZ),
                          std.error = NA, z.value = NA, p.value = NA,
                          ci.lower = NA, ci.upper = NA)

    if (XZ[[1]] == XZ[[2]]) newRows <- newRows[seq_len(2), ]

    parTable <- rbind(parTable, newRows)
  }
  

  for (intTermXX in intTerms) {
    # If we have both X:X and X:Z in the model, we must include the 
    # covariance between X:X and X:Z
    elemsXX <- stringr::str_split_fixed(intTermXX, ":", 2)
    X1      <- elemsXX[[1]]
    Z1      <- elemsXX[[2]]
     
    if (X1 != Z1) next

    for (intTermXZ in intTerms) {
      elemsXZ <- stringr::str_split_fixed(intTermXZ, ":", 2)
      X2      <- elemsXZ[[1]]
      Z2      <- elemsXZ[[2]]

      if (X2 == Z2 || !any(elemsXX %in% elemsXZ)) next

      X <- X1 
      Z <- elemsXZ[elemsXZ != X]

      muX  <- getMean(X, parTable)
      muZ  <- getMean(Z, parTable)

      lhs  <- c(X, Z)
      rhs  <- c(X, X)

      covs   <- calcCovParTable(lhs, rhs, parTable)
      varX   <- covs[[1]]
      covX_Z <- covs[[2]]

      covXZ_XX <- 2 * covX_Z * (muX^2 + varX) + 2 * muX * muZ * varX
      
      newRow <- data.frame(lhs = intTermXZ,
                           op = "~~",
                           rhs = intTermXX,
                           label = "",
                           est = covXZ_XX,
                           std.error = NA, z.value = NA, p.value = NA,
                           ci.lower = NA, ci.upper = NA)
      parTable <- rbind(parTable, newRow)
    }
  }

  parTable
}


#' Get Standardized Estimates
#'
#' Computes standardized estimates of model parameters for various types of \code{\link{modsem}} objects.
#'
#' @param object An object of class \code{\link{modsem_da}}, \code{\link{modsem_mplus}},
#' \code{\link{modsem_pi}}, or a parameter table (\code{parTable}) of class \code{data.frame}.
#' @param ... Additional arguments passed to underlying methods. See specific method
#' documentation for supported arguments, including:
#' \describe{
#'   \item{\code{correction}}{(Logical) Applies only to \code{modsem_pi} objects. Whether to correct
#'   standardized estimates for the interaction term by computing \code{var(xz)} based on
#'   variances and covariance of \code{x} and \code{z}. Default is \code{FALSE}.}
#'   \item{\code{std.errors}}{(Character) Specifies the method for computing standard errors when
#'   \code{correction = TRUE}. Options are \code{"rescale"}, \code{"delta"}, and \code{"monte.carlo"}.
#'   See method \code{standardized_estimates.modsem_pi()} for details.}
#' }
#'
#' @details
#' For \code{modsem_da} and \code{modsem_mplus} objects, the interaction term is not a formal
#' variable in the model and therefore lacks a defined variance. Under assumptions of normality
#' and zero-mean variables, the interaction variance is estimated as:
#' \deqn{var(xz) = var(x) * var(z) + cov(x, z)^2}
#' This means the standardized estimate for the interaction differs from approaches like
#' \code{lavaan}, which treats the interaction as a latent variable with unit variance.
#'
#' For \code{modsem_pi} objects, the interaction term is standardized by default assuming
#' \code{var(xz) = 1}, but this can be overridden using the \code{correction} argument.
#'
#' \strong{NOTE:} Standardized estimates are always placed in the \strong{\code{est}} column,
#' not \code{est.std}, regardless of model type.
#'
#' @return A \code{data.frame} with standardized estimates in the \code{est} column.
#'
#' @examples
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'
#'   # Inner Model
#'   Y ~ X + Z + X:Z
#' '
#' # Double centering approach
#' est_dca <- modsem(m1, oneInt)
#'
#' std1 <- standardized_estimates(est_dca) # no correction
#' summarize_partable(std1)
#'
#' std2 <- standardized_estimates(est_dca, correction = TRUE) # apply correction
#' summarize_partable(std2)
#'
#' \dontrun{
#' est_lms <- modsem(m1, oneInt, method = "lms")
#' standardized_estimates(est_lms) # correction not relevant for lms
#' }
#'
#' @export
standardized_estimates <- function(object, ...) {
  UseMethod("standardized_estimates")
}


#' @export
standardized_estimates.data.frame <- function(object, intercepts = FALSE, ...) {
  parTable <- object[c("lhs", "op", "rhs", "label", "est", "std.error")]
  parTable <- centerInteractions(parTable) # re-estimate path-coefficients 
                                          # when intercepts are zero

  if (!intercepts) {
    parTable <- parTable[parTable$op != "~1", ]
  } else {
    parTable[parTable$op == "~1", "est"]       <- 0
    parTable[parTable$op == "~1", "std.error"] <- NA
  }

  parTable <- var_interactions(parTable)
  lVs      <- getLVs(parTable)
  intTerms <- getIntTerms(parTable)
  etas     <- getSortedEtas(parTable, isLV = TRUE)
  xis      <- getXis(parTable, etas = etas, isLV = TRUE)
  indsLVs  <- getIndsLVs(parTable, lVs)
  allInds  <- unique(unlist(indsLVs))

  # get variances
  variancesInds <- calcVarParTable(allInds, parTable, measurement.model = TRUE)
  variancesLVs  <- calcVarParTable(lVs, parTable, measurement.model = FALSE)
  variancesIntTerms <- structure(numeric(length(intTerms)), names = intTerms)

  for (xz in intTerms) {
    variancesIntTerms[[xz]] <- parTable[parTable$lhs == xz & 
                                        parTable$rhs == xz &
                                        parTable$op == "~~", "est"]
  }

  variances <- c(variancesInds, variancesLVs, variancesIntTerms)

  # Factor Loadings
  lambda     <- NULL
  selectRows <- NULL
  selectCols <- c("est", "std.error")

  for (lV in lVs) {
    for (ind in indsLVs[[lV]]) {
      selectRows  <- parTable$lhs == lV & parTable$op == "=~" & parTable$rhs == ind
      scalingCoef <- sqrt(variances[[lV]]) / sqrt(variances[[ind]])
      lambda      <- parTable[selectRows, selectCols]

      parTable[selectRows, selectCols] <- lambda * scalingCoef
    }
  }

  # Structural Coefficients
  gamma               <- NULL
  selectStrucExprsEta <- NULL
  structExprsEta      <- NULL
  selectStrucExprs    <- parTable$op == "~" & parTable$lhs %in% etas

  for (eta in etas) {
    selectStrucExprsEta <- selectStrucExprs & parTable$lhs == eta
    structExprsEta      <- parTable[selectStrucExprsEta, ]

    for (xi in structExprsEta$rhs) { # xi can be an interaction term
      selectRows  <- selectStrucExprsEta & parTable$rhs == xi
      scalingCoef <- sqrt(variances[[xi]]) / sqrt(variances[[eta]])
      gamma       <- parTable[selectRows, selectCols]

      parTable[selectRows, selectCols] <- gamma * scalingCoef
    }
  }

  # (Co-) Variances of xis
  selectCovXis <- parTable$op == "~~" & parTable$lhs %in% xis
  selectRows   <- NULL
  combosXis    <- getUniqueCombos(xis, match = TRUE)

  for (i in seq_len(nrow(combosXis))) {
    xis         <- combosXis[i, , drop = TRUE]
    selectRows  <- selectCovXis & parTable$lhs %in% xis & parTable$rhs %in% xis
    scalingCoef <- sqrt(variances[[xis[[1]]]]) * sqrt(variances[[xis[[2]]]])

    if (xis[[1]] != xis[[2]]) {
      selectRows <- selectRows & parTable$lhs != parTable$rhs
    }

    covXis <- parTable[selectRows, selectCols]
    parTable[selectRows, selectCols] <- covXis / scalingCoef
  }

  # Residual Variances etas
  selectRows <- NULL
  residual   <- NULL

  for (eta in etas) {
    selectRows <- parTable$lhs == eta & parTable$op == "~~" & parTable$rhs == eta
    residual   <- parTable[selectRows, selectCols]

    parTable[selectRows, selectCols] <- residual / variances[[eta]]
  }

  # residual variances inds
  for (ind in allInds) {
    selectRows <- parTable$lhs == ind & parTable$op == "~~" & parTable$rhs == ind
    residual   <- parTable[selectRows, selectCols]

    parTable[selectRows, selectCols] <- residual / variances[[ind]]
  }

  # recalculate variance of interaction terms
  # and rescale coefficients for interaction terms
  parTable <- var_interactions(parTable)
  for (xz in intTerms) {
    selectRows <- parTable$rhs == xz & parTable$op == "~"
    varXZ      <- parTable[parTable$lhs == xz & parTable$op == "~~" &
                           parTable$rhs == xz, "est"]
    gamma      <- parTable[selectRows, selectCols]

    parTable[selectRows, selectCols] <- gamma / sqrt(varXZ) # unstandardizing, since varXZ != 1 | cov(X, Z) != 0
  }

  # recalculate custom parameters
  constrExprs <- sortConstrExprsFinalPt(parTable)
  parTable <- parTable[parTable$op != ":=", ]

  for (i in seq_len(NROW(constrExprs))) {
    row <- constrExprs[i, , drop=FALSE]

    expr      <- parse(text=constrExprs[i, "rhs"])
    labelList <- parTableLabelsToList(parTable) # must be updated for each iteration
  
    oldVal <- row$est
    newVal <- eval(expr, envir = labelList)

    ratio  <- newVal / oldVal
    values <- row[ , selectCols]
   
    row[, selectCols] <- values * ratio

    parTable <- rbind(parTable, row)
  }

  parTable$z.value  <- parTable$est / parTable$std.error
  parTable$p.value  <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - CI_WIDTH * parTable$std.error
  parTable$ci.upper <- parTable$est + CI_WIDTH * parTable$std.error

  parTable
}


#' Inspect model information
#'
#' @param object fittet model to inspect
#' @param what what to inspect
#' @param ... Additional arguments passed to other functions
#' @description function used to inspect fittet object. Similar to \code{lavaan::lavInspect}
#' argument \code{what} decides what to inspect
#' @details For \code{\link{modsem_pi}} objects, it is just a wrapper for \code{lavaan::lavInspect}.
#'  For \code{\link{modsem_da}} objects an internal function is called, which takes different 
#'  keywords for the \code{what} argument.
#' @export
modsem_inspect <- function(object, what = NULL, ...) {
  UseMethod("modsem_inspect")
}



#' Wrapper for vcov
#'
#' @param object fittet model to inspect
#' @param ... additional arguments
#' @description wrapper for vcov, to be used with modsem::modsem_vcov, since
#' vcov is not in the namespace of modsem, but stats
#' since vcov is not in the namespace of modsem, but stats
#' @export
modsem_vcov <- function(object, ...) {
  vcov(object, ...)
}


#' Wrapper for coef
#'
#' @param object fittet model to inspect
#' @param ... additional arguments
#' @description wrapper for coef, to be used with modsem::modsem_coef, since
#' coef is not in the namespace of modsem, but stats
#' since coef is not in the namespace of modsem, but stats
#' @export
modsem_coef <- function(object, ...) {
  coef(object, ...)
}


#' Wrapper for nobs
#'
#' @param object fittet model to inspect
#' @param ... additional arguments
#' @description wrapper for nobs, to be used with modsem::modsem_nobs, since
#' nobs is not in the namespace of modsem, but stats
#' since nobs is not in the namespace of modsem, but stats
#' @export
modsem_nobs <- function(object, ...) {
  nobs(object, ...)
}


#' Predict From \code{modsem} Models
#'
#' A generic function (and corresponding methods) that produces predicted
#' values or factor scores from \code{\link{modsem}} models.
#'
#' @param object An object of class \code{modsem_pi} or \code{modsem_da},
#'   respectively.
#' @param ... Further arguments passed to \code{lavaan::predict};
#'   currently ignored by the \code{\link{modsem_da}} method.
#'
#' @return
#' * For \code{\link{modsem_pi}}: whatever \code{lavaan::predict()}, which usually 
#'   returns a matrix of factor scores.
#' * For \code{\link{modsem_da}}: a numeric matrix \eqn{n \times p}, where \eqn{n} is the number of
#'   (complete) observations in the dataset, and \eqn{p} the number of latent variables. Each
#'   column contains either raw or standardised factor scores, depending on the
#'   \code{standardized} argument.
#'
#' @examples
#' m1 <- '
#' # Outer Model
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#' 
#' # Inner Model
#'   Y ~ X + Z + X:Z
#' '
#'
#' est_dca <- modsem(m1, oneInt, method = "dblcent")
#' modsem_predict(est_dca)
#'
#' \dontrun{
#' est_lms <- modsem(m1, oneInt, method = "lms")
#' modsem_predict(est_lms)
#' }
#'
#' @export
modsem_predict <- function(object, ...) {
  UseMethod("modsem_predict")
}
