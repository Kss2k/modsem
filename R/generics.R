#' Extract parameterEstimates from an estimated model
#'
#' @param object An object of class `modsem_pi`, `modsem_da`, or `modsem_mplus`
#' @param ... Additional arguments passed to other functions
#' @export 
parameter_estimates <- function(object, ...) {
  UseMethod("parameter_estimates")
}


#' Extract or modify parTable from an estimated model with estimated variances of interaction terms
#'
#' @param object An object of class `modsem_da`,  `modsem_mplus`, 
#' or a parTable of class `data.frame`
#' @param ... Additional arguments passed to other functions
#' @export
var_interactions <- function(object, ...) {
  UseMethod("var_interactions")
}


#' @export
var_interactions.data.frame <- function(object, ...) {
  parTable <- fillColsParTable(object)
  intTermVarRows <- parTable$lhs == parTable$rhs &
    grepl(":", parTable$lhs) & parTable$op == "~~"
  parTable <- parTable[!intTermVarRows, ]

  intTerms <- unique(parTable[grepl(":", parTable$rhs) & 
                     parTable$op == "~", "rhs"])

  for (i in seq_len(length(intTerms))) {
    # interaction term = XY
    XY <- stringr::str_split_fixed(intTerms[[i]], ":", 2) 
    muX <- getMean(XY[[1]], parTable)
    muY <- getMean(XY[[2]], parTable)
    varX <- calcCovParTable(XY[[1]], XY[[1]], parTable)
    varY <- calcCovParTable(XY[[2]], XY[[2]], parTable)
    covXY <- calcCovParTable(XY[[1]], XY[[2]], parTable)
    varXZ <- muX ^ 2 * muY ^ 2 + varX * muY ^ 2 + varY * muX ^ 2 +
      2 * muX * muY * covXY + varX * varY + covXY ^ 2
    newRow <- data.frame(lhs = intTerms[[i]],
                         op = "~~",
                         rhs = intTerms[[i]],
                         label = "",
                         est = varXZ,
                         std.error = NA, z.value = NA, p.value = NA,
                         ci.lower = NA, ci.upper = NA)
    parTable <- rbind(parTable, newRow)
  }
  parTable
}


#' Get standardized estimates
#'
#' @param object An object of class `modsem_da`,  `modsem_mplus`, 
#' or a parTable of class `data.frame`
#' @param ... Additional arguments passed to other functions
#' @details for `modsem_lms`, `modsem_qml` and `modsem_mplus` objects, 
#' the interaction term is not standardized such that var(xz) = 1. 
#' The interaction term is not an actual variable in the model, meaning that it does not 
#' have a variance. It must therefore be calculated from the other parameters in the model.
#' Assuming normality and zero-means the variance is calculated as 
#' `var(xz) = var(x) * var(z) + cov(x, z)^2`. Thus setting the variance of the interaction 
#' term to 1, would only be 'correct' if the correlation between x and z is zero.
#' This means that the standardized estimates for the interaction term will 
#' be different from those using lavaan, since there the interaction term is an 
#' actual latent variable in the model, with a standardized variance of 1.
#' @export
standardized_estimates <- function(object, ...) {
  UseMethod("standardized_estimates")
}


#' @export 
standardized_estimates.data.frame <- function(object, intercepts = FALSE, ...) {
  parTable <- object[c("lhs", "op", "rhs", "label", "est", "std.error")]
  if (!intercepts) { # remove intercepts
    parTable <- centerInteraction(parTable)
    parTable <- parTable[!(parTable$rhs == "1" & parTable$op == "~"), ]
  }
  parTable <- var_interactions(parTable)

  lVs <- getLVs(parTable)
  intTerms <- getIntTerms(parTable)
  etas <- getSortedEtas(parTable, isLV = TRUE)
  xis <- getXis(parTable, etas = etas, isLV = TRUE)
  indsLVs <- getIndsLVs(parTable, lVs)
  allInds <- unique(unlist(indsLVs))

  variances <- vector("list", length = length(allInds) + length(lVs) +
                      length(intTerms))
  names(variances) <- c(allInds, lVs, intTerms)

  # get variances 
  for (x in allInds) {
    variances[[x]] <- calcCovParTable(x, x, parTable, 
                                      measurement.model = TRUE)
  } 
  for (lV in lVs) {
    variances[[lV]] <- calcCovParTable(lV, lV, parTable, 
                                       measurement.model = FALSE)
  } 
  for (xz in intTerms) {
    variances[[xz]] <- parTable[parTable$lhs == xz & 
                                parTable$rhs == xz & 
                                parTable$op == "~~", "est"]

  }

  # factor loadings
  lambda <- NULL 
  selectRows <- NULL 
  selectCols <- c("est", "std.error")
  for (lV in lVs) {
    for (ind in indsLVs[[lV]]) {
      selectRows <- parTable$lhs == lV & parTable$op == "=~" & 
        parTable$rhs == ind
      lambda <- parTable[selectRows, selectCols]
      parTable[selectRows, selectCols] <- 
        lambda * (sqrt(variances[[lV]]) / sqrt(variances[[ind]]))
    }
  }

  # structural coefficients
  gamma <- NULL
  selectStrucExprsEta <- NULL
  structExprsEta <- NULL
  selectStrucExprs <- parTable$op == "~" & parTable$rhs != "1" &
                         parTable$lhs %in% etas
  for (eta in etas) {
    selectStrucExprsEta <- selectStrucExprs & parTable$lhs == eta
    structExprsEta <- parTable[selectStrucExprsEta, ]

    for (xi in structExprsEta$rhs) {
      selectRows <- selectStrucExprsEta & parTable$rhs == xi
      gamma <- parTable[selectRows, selectCols]
      parTable[selectRows, selectCols] <- 
        gamma * (sqrt(variances[[xi]]) / sqrt(variances[[eta]]))
    }
  }
  
  # variances / covariances of xis 
  selectCovXis <- parTable$op == "~~" & parTable$lhs %in% xis
  selectRows <- NULL
  combosXis <- getUniqueCombos(xis, match = TRUE)
  for (i in seq_len(nrow(combosXis))) {
    xis <- combosXis[i, , drop = TRUE]

    selectRows <- selectCovXis & parTable$lhs %in% xis & 
      parTable$rhs %in% xis 
    if (xis[[1]] != xis[[2]]) {
      selectRows <- selectRows & parTable$lhs != parTable$rhs
    }

    covXis <- parTable[selectRows, selectCols]
    parTable[selectRows, selectCols] <- 
      covXis / (sqrt(variances[[xis[[1]]]]) * sqrt(variances[[xis[[2]]]]))
  }

  # residual variances etas
  selectRows <- NULL 
  residual <- NULL
  for (eta in etas) {
    selectRows <- parTable$lhs == eta & parTable$op == "~~" & 
      parTable$rhs == eta
    residual <- parTable[selectRows, selectCols]
    projected <- calcCovParTable(eta, eta, parTable) - residual 
    parTable[selectRows, selectCols] <- residual / variances[[eta]]
  }

  # residual variances inds
  for (ind in allInds) {
    selectRows <- parTable$lhs == ind & parTable$op == "~~" & 
      parTable$rhs == ind
    residual <- parTable[selectRows, selectCols]
    
    parTable[selectRows, selectCols] <- residual / variances[[ind]]
  }
  
  # recalculate variance of interaction terms
  # and rescale coefficients for interaction terms
  parTable <- var_interactions(parTable)
  for (xz in intTerms) {
    selectRows <- parTable$rhs == xz & parTable$op == "~"
    varXZ <- parTable[parTable$lhs == xz & 
                      parTable$op == "~~" & 
                      parTable$rhs == xz, "est"]

    gamma <- parTable[selectRows, selectCols]
    parTable[selectRows, selectCols] <- gamma / sqrt(varXZ)
  }
  
  parTable$z.value <- parTable$est / parTable$std.error
  parTable$p.value <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - 1.96 * parTable$std.error
  parTable$ci.upper <- parTable$est + 1.96 * parTable$std.error

  parTable
}


#' Inspect model information
#'
#' @param object fittet model to inspect
#' @param what what to inspect
#' @param ... Additional arguments passed to other functions
#' @description function used to inspect fittet object. similar to `lavInspect()`
#' argument 'what' decides what to inspect
#' @details for `modsem_lms`, `modsem_qml` and `modsem_lavaan` 
#' for `modsem_lavaan`, it is just a wrapper for `lavInspect()`
#' for `modsem_lms` and `modsem_qml` what can either be "all", "matrices", "optim", 
#' or just the name of what to extract. 
#' @export
modsem_inspect <- function(object, what = NULL, ...) {
  UseMethod("modsem_inspect")
}
