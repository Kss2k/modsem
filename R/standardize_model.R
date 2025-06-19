#' Standardize a fitted \code{modsem_da} model
#'
#' \code{standardize_model()} post-processes the output of
#' \code{\link{modsem_da}()} (or of \code{\link{modsem}())} when \code{method = "lms"} /
#' \code{method = "qml"}), replacing the \emph{unstandardized} coefficient vector
#' (\code{$coefs}) and its variance–covariance matrix (\code{$vcov}) with \emph{fully
#' standardized} counterparts (including the measurement model).The procedure is purely
#' algebraic— \strong{no re-estimation is carried out} —and is therefore fast and
#' deterministic.
#'
#' @param model A fitted object of class \code{modsem_da}.
#'   Passing any other object triggers an error.
#'
#' @param monte.carlo Logical. If \code{TRUE}, the function will use Monte Carlo
#'   simulation to obtain the standard errors of the standardized estimates.
#'   If \code{FALSE}, the delta method is used.
#'   Default is \code{FALSE}.
#'
#' @param mc.reps Number of Monte Carlo replications. Default is 10,000.
#'   Ignored if \code{monte.carlo = FALSE}.
#' 
#' @param ... Arguments passed on to other functions
#'
#' @return The same object (returned invisibly) with three slots overwritten  
#' \describe{
#'   \item{\code{$parTable}}{Parameter table whose columns \code{est} and \code{std.error}
#'         now hold standardized estimates and their (delta-method)
#'         standard errors, as produced by \code{\link{standardized_estimates}()}.}
#'   \item{\code{$coefs}}{Named numeric vector of standardized coefficients.
#'         Intercepts (operator \code{~1}) are removed, because a standardized
#'         variable has mean 0 by definition.}
#'   \item{\code{$vcov}}{Variance–covariance matrix corresponding to the updated
#'         coefficient vector.  Rows/columns for intercepts are dropped, and
#'         the sub-matrix associated with rescaled parameters is adjusted 
#'         so that its diagonal equals the squared standardized standard errors.}
#' }
#' The object keeps its class attributes, allowing seamless use by downstream
#' S3 methods such as \code{summary()}, \code{coef()}, or \code{vcov()}.
#'
#' Because the function merely transforms existing estimates, \emph{parameter
#' constraints imposed through shared labels remain satisfied}. 
#'
#' @seealso
#' \describe{
#'    \item{\code{\link{standardized_estimates}()}}{Obtains the fully standardized
#'          parameter table used here.}
#'    \item{\code{\link{modsem}()}}{Fit model using LMS or QML approaches.}
#'    \item{\code{\link{modsem_da}()}}{Fit model using LMS or QML approaches.}
#' }
#' @examples
#' \dontrun{
#' # Latent interaction estimated with LMS and standardized afterwards
#' syntax <- "
#'   X =~ x1 + x2 + x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'   Y ~ X + Z + X:Z
#' "
#' fit  <- modsem_da(syntax, data = oneInt, method = "lms")
#' sfit <- standardize_model(fit, monte.carlo = TRUE)
#'
#' # Compare unstandardized vs. standardized summaries
#' summary(fit)  # unstandardized
#' summary(sfit) # standardized 
#' }
#'
#' @export
standardize_model <- function(model, monte.carlo = FALSE, mc.reps = 10000, ...) {
  stopif(!inherits(model, "modsem_da"), "The model must be of class 'modsem_da'.")

  solution <- standardized_solution_COEFS(model, monte.carlo = monte.carlo, mc.reps = mc.reps, ...)

  vcov <- solution$vcov
  coefs <- solution$coefs
  parTable <- solution$parTable

  model$type.estimates <- "standardized"
  model$parTable       <- parTable
  model$coefs          <- coefs
  model$vcov           <- vcov

  model
}


rescaleVCOV <- function(vcov, new_se) {
  isZeroOrNA <- is.na(new_se) | new_se == 0
  S <- vcov[!isZeroOrNA, !isZeroOrNA]

  D_old_inv <- diag(1 / sqrt(diag(S)))
  D_new <- diag(new_se[!isZeroOrNA])

  R <- D_old_inv %*% S %*% D_old_inv
  Z <- D_new %*% R %*% D_new

  V <- vcov
  V[!isZeroOrNA, !isZeroOrNA] <- Z

  V
}


standardized_solution_COEFS <- function(object, monte.carlo = FALSE, mc.reps = 10000, tolerance.zero = 1e-10, seed = 123, 
                                        delta.epsilon = 1e-8, ...) {
  set.seed(seed)

  stopif(!inherits(object, "modsem_da"), "The object must be of class 'modsem_da'.")

  parTable <- parameter_estimates(object)
  parTable <- parTable[c("lhs", "op", "rhs", "label", "est", "std.error")]
  parTable <- centerInteraction(parTable) # re-estimate path-coefficients 
                                          # when intercepts are zero
  parTable <- var_interactions(removeInteractionVariances(parTable))

  lVs      <- getLVs(parTable)
  intTerms <- getIntTerms(parTable)
  etas     <- getSortedEtas(parTable, isLV = TRUE)
  xis      <- getXis(parTable, etas = etas, isLV = TRUE)
  indsLVs  <- getIndsLVs(parTable, lVs)
  allInds  <- unique(unlist(indsLVs))

  originalLabels <- parTable$label
  labels         <- getParTableLabels(parTable, labelCol="label")
  parTable$label <- labels
  labels         <- unique(labels)

  # parTable$est <- NULL # no longer needed
  parTable$std.error <- NA
  
  V     <- vcov(object)
  coefs <- structure(parTable$est, names = labels)
  V     <- expandVCOV(V, labels=labels)
  coefs <- coefs[labels]

  legalNames <- stringr::str_replace_all(labels, OP_REPLACEMENTS)

  if (monte.carlo) {
    COEFS <- as.data.frame(mvtnorm::rmvnorm(mc.reps, mean = coefs, sigma = V))
  } else { # delta method
    k <- length(coefs)
    COEFS <- matrix(coefs, nrow=k, ncol=k, byrow=TRUE)
    DeltaEpsilon <- diag(delta.epsilon, k)
    isFixed <- abs(diag(V)) < tolerance.zero
    DeltaEpsilon[isFixed, isFixed] <- 0
    COEFS_p <- as.data.frame(COEFS + DeltaEpsilon)
    COEFS_m <- as.data.frame(COEFS - DeltaEpsilon)
    COEFS <- rbind(COEFS_p, COEFS_m)
  }

  # Get legal parameter names
  names(COEFS) <- legalNames
  names(coefs) <- legalNames
  COEFS <- rbind(as.data.frame(as.list(coefs)), COEFS) # first row is the original values

  colnames(V) <- legalNames
  rownames(V) <- legalNames

  parTable$label <- stringr::str_replace_all(parTable$label, OP_REPLACEMENTS)
  parTable       <- parTable[c("lhs", "op", "rhs", "label")]

  # Unstandardized copies
  parTable.ustd  <- parTable
  COEFS.ustd     <- COEFS

  # get variances
  varianceEquations <- list()
  variances <- list()

  for (x in allInds) {
    # get the variance equation for each variable
    eqVarX <- getCovEqExpr(x=x, y=x, parTable=parTable, measurement.model=TRUE)
    varianceEquations[[x]] <- eqVarX  
    variances[[x]] <- eval(eqVarX, envir = COEFS)
  }
  
  for (x in lVs) {
    # get the variance equation for each variable
    eqVarX <- getCovEqExpr(x=x, y=x, parTable=parTable)
    varianceEquations[[x]] <- eqVarX  
    variances[[x]] <- eval(eqVarX, envir = COEFS)
  }

  for (xz in intTerms) {
    labelVarXZ <- getLabelVarXZ(xz)
    eqVarXZ <- parse(text=labelVarXZ)
    varianceEquations[[xz]] <- eqVarXZ
    variances[[xz]] <- eval(eqVarXZ, envir = COEFS)
  }

  # Factor Loadings
  lambda     <- NULL
  selectRows <- NULL

  for (lV in lVs) {
    for (ind in indsLVs[[lV]]) {
      selectRows  <- parTable$lhs == lV & parTable$op == "=~" & parTable$rhs == ind
      label <- parTable[selectRows, "label"]

      # est in parTable
      scalingCoef <- sqrt(variances[[lV]]) / sqrt(variances[[ind]])
      lambda      <- COEFS[[label]] * scalingCoef

      COEFS[[label]] <- lambda
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

    for (xi in structExprsEta$rhs) {
      selectRows  <- selectStrucExprsEta & parTable$rhs == xi
      scalingCoef <- sqrt(variances[[xi]]) / sqrt(variances[[eta]])
      label       <- parTable[selectRows, "label"]
      gamma       <- COEFS[[label]] * scalingCoef
      
      COEFS[[label]] <- gamma
    }
  }

  # (Co-) Variances of xis
  selectCovXis <- parTable$op == "~~" & 
    (parTable$lhs %in% c(xis, intTerms) | parTable$lhs %in% c(xis, intTerms))

  covRowsXis <- parTable[selectCovXis, , drop = FALSE]
  for (i in seq_len(nrow(covRowsXis))) {
    lhs         <- covRowsXis$lhs[[i]]
    rhs         <- covRowsXis$rhs[[i]]
    xis         <- c(lhs, rhs)
    selectRows  <- selectCovXis & parTable$lhs %in% xis & parTable$rhs %in% xis
    scalingCoef <- sqrt(variances[[lhs]]) * sqrt(variances[[rhs]])

    if (lhs != rhs) selectRows <- selectRows & parTable$lhs != parTable$rhs

    label <- parTable[selectRows, "label"]
    covs <- COEFS[[label]] / scalingCoef

    COEFS[[label]] <- covs
  }

  # Residual Variances etas
  selectRows <- NULL
  residual   <- NULL
  for (eta in etas) {
    selectRows <- parTable$lhs == eta & parTable$op == "~~" & parTable$rhs == eta
    label <- parTable[selectRows, "label"]
    residual <- COEFS[[label]] / variances[[eta]]

    COEFS[[label]] <- residual
  }

  # residual variances inds
  for (ind in allInds) {
    selectRows <- parTable$lhs == ind & parTable$op == "~~" & parTable$rhs == ind
    label <- parTable[selectRows, "label"]
    residual <- COEFS[[label]] / variances[[ind]]

    COEFS[[label]] <- residual
  }

  # recalculate custom parameters
  constrExprs <- sortConstrExprsFinalPt(parTable)
  for (i in seq_len(NROW(constrExprs))) {
    row <- constrExprs[i, , drop=FALSE]
    label <- row$label
    expr <- parse(text=constrExprs[i, "rhs"])
    newVals <- eval(expr, envir = COEFS)
  
    COEFS[[label]] <- newVals
  }

  nzeros <- round(log10(1 / tolerance.zero), 0)
  coefs <- unlist(COEFS[1, ])
  COEFS <- round(COEFS[2:(nrow(COEFS)), ],  nzeros) # skip first row

  if (monte.carlo) {
    vcov <- stats::cov(COEFS)
  } else {
    # delta method
    p <- 1:k
    m <- (k+1):(k + k)

    COEFS_p <- t(COEFS[p, ])
    COEFS_m <- t(COEFS[m, ])

    J <- (COEFS_p - COEFS_m) / (2 * delta.epsilon)
    vcov <- J %*% V %*% t(J)
  }

  # fill parTable
  std.errors <- suppressWarnings(sqrt(diag(vcov)))
  warnFunc <- function(type, row) {
    warning2("Unable to calculate standardized ", type, " for: ", 
             paste0(row$lhs, row$op, row$rhs))
  }

  verboseLabels <- stringr::str_replace_all(parTable$label, OP_REPLACEMENTS)
  for (i in seq_len(nrow(parTable))) {
    row <- parTable[i, , drop = FALSE]
    label <- verboseLabels[[i]]

    if (!label %in% names(coefs)) {
      warnFunc("coefficient", row)
      est <- NA
    } else est <- coefs[[label]]

    if (!label %in% names(std.errors)) {
      warnFunc("std.error", row)
      se <- NA 
    } else se <- std.errors[[label]]

    parTable[i, "est"] <- est
    parTable[i, "std.error"] <- se
  }

  # Fill in NA on zero-standard errors, and create vcov and coefs
  parTable[!is.na(parTable$std.error) & 
           abs(parTable$std.error) < tolerance.zero, "std.error"] <- NA

  isFree <- !is.na(parTable$std.error)
  coefs <- structure(parTable$est[isFree], names = parTable$label[isFree])
  params <- intersect(names(coefs), rownames(vcov))

  coefs <- coefs[params]
  vcov <- vcov[params, params]

  cleanedParamLabels <- stringr::str_replace_all(params, OP_REPLACEMENTS_INV)
  names(coefs) <- cleanedParamLabels 
  colnames(vcov) <- rownames(vcov) <- cleanedParamLabels

  # Finalize parTable
  parTable[!parTable$label %in% originalLabels, "label"] <- "" 
  parTable$z.value  <- parTable$est / parTable$std.error
  parTable$p.value  <- 2 * stats::pnorm(-abs(parTable$z.value))
  parTable$ci.lower <- parTable$est - CI_WIDTH * parTable$std.error
  parTable$ci.upper <- parTable$est + CI_WIDTH * parTable$std.error


  # Correct Scale of interaction term
  intTermsList <- structure(stringr::str_split(intTerms, ":"),
                            names = intTerms)
  parTable <- correctStdSolutionCOEFS(
    parTable = parTable.ustd, 
    parTable.std = parTable, 
    intTerms = intTermsList, 
    cols = c("est", "std.error")
  )

  # TODO: Scale vcov and coefs too!!
  warning("vcov and coefs must be scaled!")
  warning("Calculation of custom parameters after standardization is incorrect!")

  list(parTable  = parTable,
       coefs     = coefs, 
       vcov      = vcov,
       COEFS     = COEFS)
}


correctStdSolutionCOEFS <- function(parTable, 
                                    COEFS.std, 
                                    COEFS.ustd, 
                                    variances,
                                    intTerms, 
                                    cols = c("est", "se")) {
  getLabel <- \(lhs, op, rhs) paste0(lhs, OP_REPLACEMENTS[[op]], rhs)

  for (XZ in intTerms) {
    elems <- stringr::str_split_fixed(intTerms[[XZ]], ":", 2)
    X     <- elems[[1]]
    Z     <- elems[[2]]

    vars  <- variances[c(X, Z)]
    sds   <- sqrt(vars)

    rowsXZ <- parTable[parTable$rhs == XZ & parTable$op == "~", , drop = FALSE]
    Y <- rowsXZ$lhs[[1]]

    if (!length(Y)) {
      warning2("No endogenous variable found for interaction term '", XZ, "'.",
               immediate. = FALSE)
      next
    }

    # Find the relevant exogenous variables
    struct <- parTable[parTable$lhs == y & parTable$op == "~", , drop = FALSE]
    xis    <- struct$rhs

    # Unstandardized terms
    varY  <- variances[[Y]]
    sdY   <- sqrt(varY)
    B3    <- getCOEFS(y = Y, x = XZ, COEFS = COEFS.std)

    # Correlations
    combosXis <- getUniqueCombos(xis)
    combosXis <- combosXis[combosXis[[1]] == XZ |
                           combosXis[[2]] == XZ, , drop = FALSE]

    lxis <- combosXis[[1]]
    rxis <- combosXis[[2]]

    corrs <- matrix(NA, nrow = NROW(COEFS), ncol = length(lxis))
    for (i in seq_len(ncol(corrs))) {
      eqCorr <- getCovEqExpr(x = lxis[i], y = rxis[i], parTable = parTable)
      corrs[, i] <- eval(eqCorr, envir = COEFS.std)
    }

    # Incorrectly standardized terms
    lcoefsIncorrect <- getCOEFS(x = lxis, y = y, COEFS = COEFS.std)
    rcoefsIncorrect <- getCOEFS(x = rxis, y = y, COEFS = COEFS.std)

    corrtermsIncorrect <- sum(2 * lcoefsIncorrect * rcoefsIncorrect * corrs)
    
    b3Incorrect <- getCoefs(y = y, x = XZ, parTable = parTable.std)
    projVarY_XZ <- b3Incorrect ^ 2 + corrtermsIncorrect # this should be the same 
                                                        # for both the correctly, and
                                                        # incorrectly standardized terms
                                                        # and is the identity which makes it
                                                        # possible to standardize the terms correctly

    # Correctly standardized terms
    b3Correct <- B3 * abs(prod(sds) / sdy) # in case some variances are negative
                                           # we want to make sure we don't flip
                                           # the sign...
    
    lcoefsCorrect <- lcoefsIncorrect
    rcoefsCorrect <- rcoefsIncorrect
    lcoefsCorrect[lxis == XZ] <- b3Correct
    rcoefsCorrect[rxis == XZ] <- b3Correct
    corrtermsCorrect <- sum(2 * lcoefsCorrect * rcoefsCorrect * corrs)

    # Calculate the correct standard deviation of the interaction term
    # Using the identity:
    #   projVarY_XZ = b3Correct ^ 2 * sd(xz) ^ 2 + sd(xz) * corrterms
    # Solve for sd(xz) using the quadratic formula:
    #   sd(xz) = (- corrterms +/- sqrt(corrterms^2 + 4 * b3Correct ^ 2 * projVarY_XZ)) / 2 * b3Correct ^ 2
    numerator <- -corrtermsCorrect + sign(projVarY_XZ) * 
      sqrt(corrtermsCorrect^2 + 4 * (b3Correct ^ 2) * projVarY_XZ)
    denominator <- 2 * (b3Correct ^ 2)
    sdXZ <- numerator / denominator # correctly standardized sd(xz)

    # Correct parTable
    # We apply scaling factors to the parameters and the standard errors.
    # In theory we should use the delta method to calculate the new standard errors
    # but in practice it shouldn't really matter...
    lequalY  <- parTable.std$lhs == y
    requalXZ <- parTable.std$rhs == XZ
    lequalXZ <- parTable.std$lhs == XZ
    isCov    <- parTable.std$op == "~~"
    isCoef   <- parTable.std$op == "~"

    isCovTerm  <- xor(lequalXZ, requalXZ) & isCov
    isCoefTerm <- lequalY & requalXZ & isCoef
    isVarTerm  <- lequalXZ & requalXZ & isCov

    # Scaling factors
    scalefVar  <- (sdXZ^2) / 1
    scalefCoef <- b3Correct / b3Incorrect
    scalefCov  <- sdXZ

    # Apply
    parTable.std[isCovTerm, cols]  <- parTable.std[isCovTerm, cols]  * scalefCov
    parTable.std[isCoefTerm, cols] <- parTable.std[isCoefTerm, cols] * scalefCoef
    parTable.std[isVarTerm, cols]  <- parTable.std[isVarTerm, cols]  * scalefVar
  }

  etas <- unique(parTable.std[parTable.std$op == "~", "lhs"])
  varEtas <- calcVarParTable(etas, parTable.std)

  warnif(
    any(abs(varEtas - 1) > 1e-10), 
    "Some variances are not equal to 1! ",
    "This indicates that the solution was not standardized correctly!",
    immediate. = FALSE
  )

  warning("Custom parameters must also be corrected in `correctStdSolution()`!")
  attr(parTable.std, "var.etas") <- varEtas
  parTable.std
}


correctStdSolutionPI <- function(object, parTable.std) {
  parTable <- parameter_estimates(object)
  intTerms <- object$elementsInProdNames

  correctStdSolution(parTable = parTable, parTable.std = parTable.std,
                     cols = c("est", "se"))
}


correctStdSolution <- function(parTable, parTable.std, intTerms, cols = c("est", "se")) {
  # Function for correcting the standardized solution for interaction terms
  # lavaan standardizes the interaction term (X:Z) to have a variance of 1,
  # which is not correct if cov(X, Z) != 0.
  for (xz in names(intTerms)) {
    elems <- intTerms[[xz]]
    vars  <- calcVarParTable(elems, parTable)
    sds   <- sqrt(vars)

    rowsXZ <- parTable[parTable$rhs == xz & parTable$op == "~", , drop = FALSE]
    y <- rowsXZ$lhs[[1]]

    if (!length(y)) {
      warning2("No endogenous variable found for interaction term '", xz, "'.",
               immediate. = FALSE)
      next
    }

    # Find the relevant exogenous variables
    struct <- parTable[parTable$lhs == y & parTable$op == "~", , drop = FALSE]
    xis    <- struct$rhs

    # Unstandardized terms
    vary  <- calcVarParTable(y, parTable)
    sdy   <- sqrt(vary)
    B3    <- getCoefs(y = y, x = xz, parTable = parTable)

    # Correlations
    combosXis <- getUniqueCombos(xis)
    combosXis <- combosXis[combosXis[[1]] == xz |
                           combosXis[[2]] == xz, , drop = FALSE]

    lxis <- combosXis[[1]]
    rxis <- combosXis[[2]]

    corrs <- calcCovParTable(x = lxis, y = rxis, parTable = parTable.std)

    # Incorrectly standardized terms
    lcoefsIncorrect <- getCoefs(x = lxis, y = y, parTable = parTable.std)
    rcoefsIncorrect <- getCoefs(x = rxis, y = y, parTable = parTable.std)

    corrtermsIncorrect <- sum(2 * lcoefsIncorrect * rcoefsIncorrect * corrs)
    
    b3Incorrect <- getCoefs(y = y, x = xz, parTable = parTable.std)
    projVarY_XZ <- b3Incorrect ^ 2 + corrtermsIncorrect # this should be the same 
                                                        # for both the correctly, and
                                                        # incorrectly standardized terms
                                                        # and is the identity which makes it
                                                        # possible to standardize the terms correctly

    # Correctly standardized terms
    b3Correct <- B3 * abs(prod(sds) / sdy) # in case some variances are negative
                                           # we want to make sure we don't flip
                                           # the sign...
    
    lcoefsCorrect <- lcoefsIncorrect
    rcoefsCorrect <- rcoefsIncorrect
    lcoefsCorrect[lxis == xz] <- b3Correct
    rcoefsCorrect[rxis == xz] <- b3Correct
    corrtermsCorrect <- sum(2 * lcoefsCorrect * rcoefsCorrect * corrs)

    # Calculate the correct standard deviation of the interaction term
    # Using the identity:
    #   projVarY_XZ = b3Correct ^ 2 * sd(xz) ^ 2 + sd(xz) * corrterms
    # Solve for sd(xz) using the quadratic formula:
    #   sd(xz) = (- corrterms +/- sqrt(corrterms^2 + 4 * b3Correct ^ 2 * projVarY_XZ)) / 2 * b3Correct ^ 2
    numerator <- -corrtermsCorrect + sign(projVarY_XZ) * 
      sqrt(corrtermsCorrect^2 + 4 * (b3Correct ^ 2) * projVarY_XZ)
    denominator <- 2 * (b3Correct ^ 2)
    sdXZ <- numerator / denominator # correctly standardized sd(xz)

    # Correct parTable
    # We apply scaling factors to the parameters and the standard errors.
    # In theory we should use the delta method to calculate the new standard errors
    # but in practice it shouldn't really matter...
    lequalY  <- parTable.std$lhs == y
    requalXZ <- parTable.std$rhs == xz
    lequalXZ <- parTable.std$lhs == xz
    isCov    <- parTable.std$op == "~~"
    isCoef   <- parTable.std$op == "~"

    isCovTerm  <- xor(lequalXZ, requalXZ) & isCov
    isCoefTerm <- lequalY & requalXZ & isCoef
    isVarTerm  <- lequalXZ & requalXZ & isCov

    # Scaling factors
    scalefVar  <- (sdXZ^2) / 1
    scalefCoef <- b3Correct / b3Incorrect
    scalefCov  <- sdXZ

    # Apply
    parTable.std[isCovTerm, cols]  <- parTable.std[isCovTerm, cols]  * scalefCov
    parTable.std[isCoefTerm, cols] <- parTable.std[isCoefTerm, cols] * scalefCoef
    parTable.std[isVarTerm, cols]  <- parTable.std[isVarTerm, cols]  * scalefVar
  }

  etas <- unique(parTable.std[parTable.std$op == "~", "lhs"])
  varEtas <- calcVarParTable(etas, parTable.std)

  warnif(
    any(abs(varEtas - 1) > 1e-10), 
    "Some variances are not equal to 1! ",
    "This indicates that the solution was not standardized correctly!",
    immediate. = FALSE
  )

  warning("Custom parameters must also be corrected in `correctStdSolution()`!")
  attr(parTable.std, "var.etas") <- varEtas
  parTable.std
}


