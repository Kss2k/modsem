#' Extract parameterEstimates from an estimated model
#'
#' @param object An object of class \code{\link{modsem_pi}}, \code{\link{modsem_da}}, or \code{\link{modsem_mplus}}
#' @param ... Additional arguments passed to other functions
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
var_interactions.data.frame <- function(object, ...) {
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

    muX  <- getMean(X, parTable)
    muZ  <- getMean(Z, parTable)

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


#' Get standardized estimates
#'
#' @param object An object of class \code{\link{modsem_da}}, \code{\link{modsem_mplus}},
#' or a \code{parTable} of class \code{data.frame}
#' @param ... Additional arguments passed to other functions
#' @details For \code{modsem_da}, and \code{modsem_mplus} objects,
#' the interaction term is not standardized such that \code{var(xz) = 1}.
#' The interaction term is not an actual variable in the model, meaning that it does not
#' have a variance. It must therefore be calculated from the other parameters in the model.
#'
#' Assuming normality and zero-means, the variance is calculated as
#' \code{var(xz) = var(x) * var(z) + cov(x, z)^2}. Thus setting the variance of the interaction
#' term to 1 would only be 'correct' if the correlation between \code{x} and \code{z} is zero.
#' This means that the standardized estimates for the interaction term will
#' be different from those using \code{lavaan}, since there the interaction term is an
#' actual latent variable in the model, with a standardized variance of 1. 
#'
#' In \code{\link{modsem_pi}} the interaction term is standardized such that \code{var(xz) = 1}.
#' It is possible to apply a correction to the standardized estimates from 
#' \code{\link{modsem_pi}}, by passing \code{correction = TRUE}.
#'
#' \strong{NOTE} that the standardized coefficent will be placed in the \strong{\code{est}} column
#' (\strong{not \code{est.std}}) for all models, inluding those from \code{\link{modsem_pi}}. 
#' This is different from the results from \code{lavaan::standardizedSolution} where the 
#' standardized estimates are placed in the \code{est.std} column.
#'
#' @return
#' A \code{data.frame} with the standardized estimates of the model parameters, in the 
#' \code{est} column.
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
#' standardized_estimates(est_dca) # no correction
#' standardized_estimates(est_dca, correction = TRUE) # apply correction
#'
#' \dontrun{
#' est_lms <- modsem(m1, oneInt, method = "lms")
#' standardized_estimates(est_lms) # correction not relevant for lms
#' }
#' @export
standardized_estimates <- function(object, ...) {
  UseMethod("standardized_estimates")
}


#' @export
standardized_estimates.data.frame <- function(object, intercepts = FALSE, ...) {
  parTable <- object[c("lhs", "op", "rhs", "label", "est", "std.error")]
  parTable <- centerInteraction(parTable) # re-estimate path-coefficients 
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
#' @description function used to inspect fittet object. similar to `lavInspect()`
#' argument 'what' decides what to inspect
#' @details for `modsem_da`, and `modsem_lavaan`
#' for `modsem_lavaan`, it is just a wrapper for `lavInspect()`
#' for `modsem_da` and `` what can either be "all", "matrices", "optim",
#' or just the name of what to extract.
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
