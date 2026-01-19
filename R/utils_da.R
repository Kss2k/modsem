OP_REPLACEMENTS <- c("~~"  = "___COVARIANCE___",
                     "=~"  = "___MEASUREMENT___",
                     ":="  = "___CUSTOM___",
                     "~"   = "___REGRESSION___",
                     ":"   = "___INTERACTION___",
                     "<->" = "___MPLUS_COVARIANCE___",
                     "<-"  = "___MPLUS_REGRESSION___",
                     "\\|" = "___THRESHOLD___")

OP_OV_INT_1 <- "."
OP_OV_INT_2 <- OP_REPLACEMENTS[[":"]]
OP_REPLACEMENTS_INV <- structure(names(OP_REPLACEMENTS), names = OP_REPLACEMENTS)


CONSTRAINT_OPS <- c("==", ">", "<", ":=")
BOUNDUARY_OPS <- c(">", "<")


TEMP_OV_PREFIX <- ".TEMP_OV__"


getFreeParams <- function(model) {
  model$freeParams
}


fetch <- function(x, pattern = ".*") {
  x[grepl(pattern, names(x))]
}


stripMatrices <- function(matrices, fill = -1) {
  lapply(matrices, function(mat) {
    mat[!is.na(mat)] <- fill
    mat
  })
}


removeInteractions <- function(model) {
  model$matrices$OmegaEtaXi[TRUE] <- 0
  model$matrices$OmegaXiXi[TRUE] <- 0
  model
}


# Faster version of mvtnorm::dmvnorm() given that sigma is positive
# there are some drawbacks to using mvnfast. In particular,
# its a little less consistent
dmvn <- function(X, mean, sigma, log = FALSE) {
  tryCatch({
    csigma <- suppressWarnings(chol(sigma))
    mvnfast::dmvn(X, mu=mean, sigma=csigma, log=log,
                  ncores = ThreadEnv$n.threads, isChol=TRUE)
  }, error = function(e) mvtnorm::dmvnorm(X, mean, sigma, log))
}


totalLogDmvn <- function(mu, sigma, nu, S, n, d) {
  # X: n x d matrix of observations (rows = samples, cols = features)
  # mu: d-dimensional mean vector
  # sigma: d x d covariance matrix
  if (any(diag(sigma) < 0))
    return(NaN)

  tryCatch({
    sigma_inv <- chol2inv(chol(sigma))

    log_det_sigma <- determinant(sigma, logarithm = TRUE)$modulus

    trace_term <- sum(sigma_inv * S)  # Efficient trace of product
    mean_diff <- matrix(nu - mu, nrow = 1)
    mahalanobis_term <- n * (mean_diff %*% sigma_inv %*% t(mean_diff))

    log_likelihood <- -0.5 * (n * d * log(2 * pi) + n * log_det_sigma + trace_term + mahalanobis_term)
    as.numeric(log_likelihood)

  }, error = function(e) NA)
}


totalLogDvmnW <- function(X, mu, sigma, nu, S, tgamma, n, d) {
  # X: n x d matrix of data
  # mu: d-dimensional mean vector
  # sigma: d x d covariance matrix
  # gamma: n-dimensional vector of weights
  # nu: weighted observed means
  # S: weighted observed covariance matrix
  if (any(diag(sigma) < 0))
    return(NaN)

  tryCatch({
    sigma_inv <- chol2inv(chol(sigma))

    log_det_sigma <- determinant(sigma, logarithm = TRUE)$modulus

    trace_term <- sum(sigma_inv * S)
    mean_diff <- matrix(nu - mu, nrow = 1)
    mahalanobis_term <- tgamma * (mean_diff %*% sigma_inv %*% t(mean_diff))

    log_likelihood <- -0.5 * (tgamma * d * log(2 * pi) + tgamma * log_det_sigma + trace_term + mahalanobis_term)
    as.numeric(log_likelihood)

  },  error = function(e) NA)
}


diagPartitionedMat <- function(X, Y) {
  isNULL <- \(x) is.null(x) || length(x) == 0

  if (isNULL(X)) return(Y) else if (isNULL(Y)) return(X)

  if (is.null(dimnames(X)) && is.null(dimnames(Y))) {
    rownames <- NULL
    colnames <- NULL

  } else if (is.null(dimnames(X))) {
    rownames <- c(rep("", NROW(X)), rownames(Y))
    colnames <- c(rep("", NCOL(X)), colnames(Y))

  } else if (is.null(dimnames(Y))) {
    rownames <- c(rownames(X), rep("", NROW(Y)))
    colnames <- c(colnames(X), rep("", NCOL(Y)))

  } else {
    rownames <- c(rownames(X), rownames(Y))
    colnames <- c(colnames(X), colnames(Y))
  }

  structure(rbind(cbind(X, matrix(0, nrow = NROW(X), ncol = NCOL(Y))),
                  cbind(matrix(0, nrow = NROW(Y), ncol = NCOL(X)), Y)),
            dimnames = list(rownames, colnames))
}


formatNumeric <- function(x, digits = 3, scientific = FALSE,
                          justify = "right", width = NULL) {
  digits_fmt <- if (is.finite(digits)) max(0L, as.integer(digits)) else 3L
  digits_fmt_fmt <- max(1L, digits_fmt)
  if (is.numeric(x)) {
    x_round <- round(x, digits_fmt)
    format(x_round, nsmall = digits_fmt, digits = digits_fmt_fmt,
           trim = FALSE, justify = justify, scientific = scientific,
           width = width)
  } else {
    format(x, trim = FALSE, justify = justify, scientific = scientific,
           width = width)
  }
}


oneWayTableToDataFrame <- function(table) {
  df <- as.data.frame(table)
  out <- data.frame(freq = df$Freq)
  rownames(out) <- df$Var1
  out
}


whichIsMax <- function(x) {
  which(x == max(x))
}


getK_NA <- function(omega, labelOmega) {
  na    <- apply(omega, MARGIN = 1, FUN = function(x) any(is.na(x)))
  label <- apply(labelOmega, MARGIN = 1, FUN = function(x) any(x != ""))
  sum(na | label)
}


sortData <- function(data, allIndsXis, allIndsEtas) {
  inds    <- unique(c(allIndsXis, allIndsEtas))
  ovs     <- colnames(data)
  missing <- inds[!inds %in% ovs]

  stopif(!all(inds %in% ovs), "Missing observed variables in data:\n  ",
         paste(missing, collapse = ", "))

  data[unique(c(allIndsXis, allIndsEtas))]
}


anyAllNA <- function(data) {
  any(vapply(data, FUN.VALUE = logical(1L), function(x) all(is.na(x))))
}


castDataNumericMatrix <- function(data) {
  force(data) # evaluate to check for errors
  data <- tryCatch({
    numericData <- lapplyDf(data, FUN = as.numeric)
  },
  warning = function(w) {
    warning2("Warning in converting data to numeric: \n", w)
    numericData <- suppressWarnings(lapplyDf(data, FUN = as.numeric))
    stopif(anyAllNA(numericData), "Unable to convert data to type numeric")
    numeric
  },
  error = function(e) {
    stop2("Unable to convert data to type numeric")
  })
  as.matrix(data)
}


patternizeMissingDataFIML <- function(data) {
  # if we are not using fiml, the missing data should already have been removed...
  CLUSTER <- attr(data, "cluster")
  WEIGHTS <- attr(data, "weights")

  Y   <- as.matrix(data)
  obs <- !is.na(Y)

  rowMissingAll <- apply(obs, MARGIN = 1, FUN = \(x) !any(x))
  colMissingAll <- apply(obs, MARGIN = 2, FUN = \(x) !any(x))
  stopif(any(colMissingAll),
         "Please remove variables with only missing values:\n  ",
         paste0(colnames(obs)[colMissingAll], collapse = ", "))

  patterns <- unique(obs, MARING = 2)

  if (any(rowMissingAll)) { # remove patterns where all are missing
    # This shouldn't really happen, as it should be handled already in
    # `handleMissingData()`. Regardless, we should handle it if it ever happens
    return(patternizeMissingDataFIML(data[!rowMissingAll, , drop = FALSE]))
  }

  ids <- seq_len(NROW(patterns))
  n   <- NROW(Y)
  k   <- NCOL(Y)

  rowidx <- vector("list", length = NROW(ids))
  colidx <- vector("list", length = NROW(ids))
  data.split  <- vector("list", length = NROW(ids))
  n.pattern  <- numeric(NROW(ids))
  d.pattern  <- numeric(NROW(ids))

  for (id in ids) {
    mask  <- matrix(patterns[id, ], nrow = n, ncol = k, byrow = TRUE)
    match <- apply(obs == mask, MARGIN = 1, FUN = all)
    ridx  <- which(match)
    cidx  <- which(patterns[id, ])

    rowidx[[id]] <- ridx
    colidx[[id]] <- cidx
    data.split[[id]] <- Y[ridx, cidx, drop = FALSE]
    n.pattern[[id]] <- sum(match)
    d.pattern[[id]] <- length(cidx)
  }

  list(
    ids        = ids,
    rowidx     = rowidx,
    colidx     = colidx,
    colidx0    = lapply(colidx, FUN = \(idx) idx - 1),
    patterns   = patterns,
    data.split = data.split,
    n.pattern  = n.pattern,
    d.pattern  = d.pattern,
    n          = NROW(data),
    k          = NCOL(data),
    p          = length(ids),
    data.full  = data,
    is.fiml    = length(ids) > 1L,
    cluster    = CLUSTER,
    weights    = WEIGHTS
  )
}


handleMissingData <- function(data, missing = "listwise", CLUSTER = NULL, WEIGHTS = NULL) {
  missing       <- tolower(missing)
  completeCases <- stats::complete.cases(data)
  anyMissing    <- any(!completeCases)
  allMissing    <- all(!completeCases)

  if (!anyMissing){
    attr(data, "cluster") <- CLUSTER
    attr(data, "weights") <- WEIGHTS
    return(data)

  } else if (allMissing) {
    missingAllCol <- apply(data, MARGIN = 2, FUN = \(x) all(is.na(x)))
    colsMissing   <- colnames(data)[missingAllCol]

    stop2("Please remove variables with only missing values:\n  ",
          paste0(colsMissing, collapse = ", "))
  }

  if (missing %in% c("listwise", "casewise", "complete")) {
    warning2("Removing missing values list-wise!\n",
             "Consider using `missing=\"fiml\"`, `missing=\"impute\"`, ",
             "or the `modsem_mimpute()` function!\n")

    out <- data[completeCases, ]
    attr(out, "cluster") <- CLUSTER[completeCases]
    attr(out, "weights") <- WEIGHTS[completeCases]

    return(out)

  } else if (missing == "impute") {
    message("Imputing missing values. Consider using the `modsem_mimpute()` function!")

    imp  <- Amelia::amelia(data, m = 1, p2s = 0)

    imp1 <- as.matrix(as.data.frame(imp$imputations[[1]]))
    attr(imp1, "cluster") <- CLUSTER

    return(imp1)

  } else if (missing %in% c("fiml", "ml", "direct")) {
    attr(data, "cluster") <- CLUSTER
    attr(data, "weights") <- WEIGHTS

    rowMissingAll <- apply(data, MARGIN = 1, FUN = \(x) all(is.na(x)))
    data          <- data[!rowMissingAll, , drop = FALSE] # we've already know that
                                                          # all(rowMissingAll) != TRUE
    return(data)

  } else {
    stop2(sprintf("Unrecognized value for `missing`: `%s`", missing))
  }
}


prepDataModsemDA <- function(data, allIndsXis, allIndsEtas, missing = "listwise",
                             cluster = NULL, sampling.weights = NULL) {

  if (is.null(data) || !NROW(data))
    return(list(data.full = NULL, n = 0, k = 0, p = 0, cluster = NULL))

  if (!is.null(cluster)) {
    stopif(length(cluster) > 1L, "`cluster` must be a single variable!")

    CLUSTER <- as.factor(data[, cluster])

  } else CLUSTER <- NULL

  if (!is.null(sampling.weights)) {
    stopif(length(sampling.weights) > 1L, "`sampling.weights` must be a single variable!")
    WEIGHTS <- data[, sampling.weights]

  } else WEIGHTS <- NULL

  sortData(data, allIndsXis,  allIndsEtas) |>
    castDataNumericMatrix() |>
    handleMissingData(missing = missing, CLUSTER = CLUSTER, WEIGHTS = WEIGHTS) |>
    patternizeMissingDataFIML()
}


canBeNumeric <- function(x, includeNA = FALSE) {
  if (includeNA) x[x == ""] <- 0
  !is.na(suppressWarnings(as.numeric(x)))
}


createDoubleIntTerms <- function(x, z = NULL, sep = ":") {
  if (is.null(z)) {
    z <- x[[2]]
    x <- x[[1]]
  }
  c(paste0(x, sep, z), paste0(z, sep, x))
}


getFreeOrConstIntTerms <- function(varsInInt, eta, intTerms) {
  expr <- intTerms[intTerms$lhs == eta & intTerms$rhs %in%
                   createDoubleIntTerms(varsInInt), "mod"]
  if (canBeNumeric(expr, includeNA = TRUE)) as.numeric(expr) else 0
}


getLabelIntTerms <- function(varsInInt, eta, intTerms) {
  expr <- intTerms[intTerms$lhs == eta & intTerms$rhs %in%
                   createDoubleIntTerms(varsInInt), "mod"]
  if (!canBeNumeric(expr)) return(expr)
  ""
}


getEmptyModel <- function(group.info, cov.syntax, parTableCovModel,
                          mean.observed = TRUE, method = "lms") {
  group.info$parTable$mod <- ""
  group.info$parTable <- removeConstraintExpressions(group.info$parTable)
  group.info["data"]  <- list(data = NULL)

  if (!is.null(cov.syntax) && is.null(parTableCovModel)) {
    parTableCovModelOrig <- modsemify(cov.syntax)

    parTableCovModel <- expandParTableByGroup(
      parTableCovModelOrig, group.levels = group.info$group.levels
    )

    group.info$parTableCov.orig <- parTableCovModelOrig
    group.info$parTableCov      <- parTableCovModel

  } else if (NROW(parTableCovModel)) {
    parTableCovModel$mod <- ""
    parTableCovModel <- removeConstraintExpressions(parTableCovModel)

    group.info$parTableCov <- parTableCovModel
  }

  specifyModelDA(
    group.info       = group.info,
    method           = method,
    cov.syntax       = cov.syntax,
    mean.observed    = mean.observed,
    auto.fix.first   = FALSE,
    auto.fix.single  = FALSE,
    createTheta      = FALSE,
    checkModel       = FALSE
  )
}


#' @export
as.character.matrix <- function(x, empty = TRUE, ...) {
  if (empty) x[TRUE] <- ""
  matrix(as.character(x), nrow = NROW(x), ncol = NCOL(x),
          dimnames = dimnames(x))
}


replaceNonNaModelMatrices <- function(model, value = -999) {
  .fillna <- function(x) {
    x[!is.na(x)] <- value
    x
  }

  for (g in seq_along(model$models))
    model$models[[g]]$matrices <- lapply(model$models[[g]]$matrices, FUN = .fillna)

  model
}


removeUnknownLabels <- function(parTable) {
  ops    <- c("==", ">", "<", ":=")
  labels <- unique(parTable$mod[parTable$mod != ""])
  parTable[!parTable$op %in% ops |
           (parTable$op %in% ops &
            parTable$lhs %in% parTable$mod &
            !constraintsContainUnmatchedLabels(parTable, labels)), ]
}


getLabeledParamsLavaan <- function(parTable, fixedParams = NULL) {
  if (is.null(parTable$label)) return(NULL)
  labelRows <- parTable[parTable$label != "" &
                        !parTable$label %in% fixedParams,
                        c("est", "label"), drop = FALSE] |>
    uniqueByVar("label")

  theta <- as.numeric(labelRows$est)
  names(theta) <- labelRows$label
  theta
}


uniqueByVar <- function(df, var) {
  vals <- unique(df[[var]])
  udf  <- NULL

  for (v in vals) {
    match <- df[df[[var]] == v, , drop = FALSE][1, , drop = FALSE]
    udf <- rbind(udf, match)
  }

  udf
}


removeInteractionVariances <- function(parTable) {
  parTable[!(parTable$op == "~~" &
             (grepl(":", parTable$lhs) | grepl(":", parTable$rhs))), ]
}


tr <- function(mat) sum(diag(mat))


traceOmegaXiXi <- function(omega, numEta, numXi) {
  lastRow <- 0
  lastCol <- 0
  trace <- numeric(numEta)
  for (i in seq_len(numEta)) {
    trace[[i]] <- tr(omega[seq_len(numXi) + (i - 1) * numXi,
                           seq_len(numXi) + (i - 1) * numXi])
  }
  trace
}


diagPartitioned <- function(matrix, length) {
  out <- matrix(0, nrow = length * nrow(matrix),
                ncol = length * ncol(matrix))
  nrows <- nrow(matrix)
  rows <- seq_len(nrows)
  ncols <- ncol(matrix)
  cols <- seq_len(ncols)
  for (i in seq_len(length)) {
    out[rows + (i - 1) * nrows,
        cols + (i - 1) * ncols] <- matrix
  }
  out
}


repPartitionedRows <- function(matrix, length = 1) {
  if (length <= 1) return(matrix)
  out <- matrix
  for (i in seq_len(length - 1)) {
    out <- rbind(out, matrix)
  }
  out
}


repPartitionedCols <- function(matrix, length = 1) {
  if (length <= 1) return(matrix)
  out <- matrix
  for (i in seq_len(length - 1)) {
    out <- cbind(out, matrix)
  }
  out
}


diagBindSquareMatrices <- function(X, Y) {
  XY <- matrix(0, nrow = NROW(X), ncol = NCOL(Y),
               dimnames = list(rownames(X), colnames(Y)))
  rbind(cbind(X, XY), cbind(t(XY), Y))
}


#' @export
as.logical.matrix <- function(x, ...) {
  structure(x != 0,
            dim = dim(x),
            dimnames = dimnames(x))
}


isScalingY <- function(x) {
  (seq_along(x) %in% which(x == 0)) | seq_along(x) %in% which.max(x == 1)
}


runningAverage <- function(x, n = 10) {
  if (length(x) < n) return(0)
  last <- length(x)
  if (last < n) first <- 1 else first <- last - n + 1
  mean(x[first:last])
}


nNegativeLast <- function(x, n = 10) {
  if (length(x) < n) return(0)
  last <- length(x)
  if (last < n) first <- 1 else first <- last - n + 1
  sum(x[first:last] < 0)
}


getDegreesOfFreedom <- function(p, coef, mean.structure = TRUE) {
  if (!length(p)) return(NA_real_)

  momentsPerGroup <- function(pp) {
    if (mean.structure) {
      pp * (pp + 3) / 2
    } else {
      pp * (pp + 1) / 2
    }
  }

  m_total <- sum(vapply(p, momentsPerGroup, numeric(1L)))
  m_total - length(coef)
}


nFreeInterceptsDA <- function(model) {
  parTable   <- getMissingLabels(parameter_estimates(model))
  intercepts <- parTable$label[parTable$op == "~1"]
  coefs <- coef(model, type = "free")

  sum(intercepts %in% names(coefs))
}


getInfoQuad <- function(quad) {
  list(dim = quad$k, nodes.dim = quad$m, nodes.total = quad$m ^ quad$k)
}


getFixedInterceptSyntax <- function(indicator, syntax, parTable) {
  if (is.null(indicator) || is.null(syntax) ||
      NROW(parTable[parTable$lhs == indicator &
           parTable$op == "~1", ])) return(syntax)
  else addition <-  paste0("\n", indicator, " ~ 0 * 1")
  paste0(syntax, addition)
}


getEtaRowLabelOmega <- function(label) {
  stringr::str_split_1(label, "~")[[1]]
}


getXiRowLabelOmega <- function(label) {
  stringr::str_split_1(label, "~")[[2]]
}


expandVCOV <- function(vcov, labels) {
  labels.vcov <- colnames(vcov)
  labels.vv <- intersect(labels, labels.vcov)
  labels.zz <- setdiff(labels, labels.vcov)

  m <- length(labels.vv)
  k <- length(labels.zz)

  Vvv <- vcov[labels.vv, labels.vv]
  Vzz <- matrix(0, nrow = k, ncol = k, dimnames = list(labels.zz, labels.zz))
  Vvz <- matrix(0, nrow = m, ncol = k, dimnames = list(labels.vv, labels.zz))

  V <- rbind(cbind(Vvv, Vvz), cbind(t(Vvz), Vzz))
  V[labels, labels] # sort
}


expandCoef <- function(coef, labels) {
  labels.coef <- names(coef)
  labels.c    <- intersect(labels, labels.coef)
  labels.z    <- setdiff(labels, labels.coef)

  m <- length(labels.c)
  k <- length(labels.z)

  coef.c <- coef[labels.c]
  coef.z <- stats::setNames(numeric(length(labels.z)), nm = labels.z)

  c(coef.c, coef.z)
}


getLabelVarXZ <- function(intTerm) {
  XZ   <- stringr::str_split_fixed(intTerm, ":", 2)
  X    <- XZ[[1]]
  Z    <- XZ[[2]]

  labelXZ <- paste0(X, OP_REPLACEMENTS[[":"]], Z)
  paste0(labelXZ, OP_REPLACEMENTS[["~~"]], labelXZ)
}


intTermsAffectLV <- function(lV, parTable, etas = NULL) {
  if (grepl(":", lV)) return(TRUE)
  else if (!lV %in% parTable[parTable$op == "~", "lhs"]) return(FALSE)

  gamma <- parTable[parTable$lhs == lV & parTable$op == "~", , drop = FALSE]

  if (NROW(gamma) == 0) return(FALSE)

  for (i in seq_len(NROW(gamma))) {
    rhs <- gamma$rhs[i]
    if (intTermsAffectLV(rhs, parTable)) return(TRUE)
  }

  FALSE
}


getLambdaParTable <- function(parTable, rows = NULL, cols = NULL, fill.missing = FALSE) {
  lVs <- getLVs(parTable)

  indsLVs <- getIndsLVs(parTable, lVs = lVs)
  allInds <- unique(unlist(indsLVs))

  lambda <- matrix(0, nrow = length(allInds), ncol = length(lVs),
                   dimnames = list(allInds, lVs))
  for (lV in lVs) for (ind in indsLVs[[lV]]) {
    lambda[ind, lV] <- parTable[parTable$lhs == lV & parTable$rhs == ind, "est"]
  }

  if (is.null(rows)) rows <- rownames(lambda)
  if (is.null(cols)) cols <- colnames(lambda)

  if (fill.missing) {
    missingCols <- setdiff(cols, colnames(lambda))
    missingRows <- setdiff(rows, rownames(lambda))

    if (length(missingCols)) {
      newCols <- matrix(0, nrow = NROW(lambda), ncol = length(missingCols),
                        dimnames = list(rownames(lambda), missingCols))
      lambda <- cbind(lambda, newCols)
    }

    if (length(missingRows) && fill.missing) {
      newRows <- matrix(0, nrow = length(missingRows), ncol = NCOL(lambda),
                        dimnames = list(missingRows, colnames(lambda)))
      lambda <- rbind(lambda, newRows)
    }
  }

  lambda <- lambda[rows, cols, drop = FALSE]
}


getLevelsParTable <- function(parTable) {
  parTable$lhs <- stringr::str_remove_all(parTable$lhs, pattern = ":")
  parTable$rhs <- stringr::str_remove_all(parTable$rhs, pattern = ":")

  xis <- getXis(parTable, isLV = FALSE)
  etas <- getEtas(parTable, isLV = FALSE)
  vars <- c(xis, etas)

  MAPPED <- c()

  getLevelVarParTable <- function(x) {
    if      (x %in% names(MAPPED)) return(MAPPED[[x]])
    else if (x %in% xis) {
      MAPPED[[x]] <<- 1L
      return(MAPPED[[x]])
    }

    if (!x %in% etas) stop("Unexpected type of variable: ", x)

    predictors <- parTable[parTable$op == "~" & parTable$lhs == x, , drop = FALSE]

    if (!NROW(predictors)) {
      MAPPED[[x]] <<- 1L
      return(MAPPED[[x]])
    }

    maxlevel <- 0L
    for (i in seq_len(NROW(predictors))) {
      predictor <- predictors$rhs[i]
      predlevel <- getLevelVarParTable(predictor)

      if (predlevel > maxlevel) maxlevel <- predlevel
    }

    MAPPED[[x]] <<- maxlevel + 1L
    MAPPED[[x]]
  }

  for (var in vars) {
    if (var %in% names(MAPPED)) next
    getLevelVarParTable(x = var)
  }

  MAPPED
}


isPureEta <- function(eta, parTable) {
  predictors <- unique(parTable[parTable$op == "~", "rhs"])
  indicators <- unique(parTable[parTable$op == "=~", "rhs"])
  !eta %in% c(predictors, indicators)
}


getCoefMatricesDA <- function(parTable,
                              xis = NULL,
                              etas = NULL,
                              intTerms = NULL,
                              centered = TRUE) {

  parTable <- removeInteractionVariances(parTable)

  if (centered) {
    parTable <- centerInteractions(parTable, center.means = FALSE) |>
      var_interactions(ignore.means = TRUE) |>
      meanInteractions(ignore.means = TRUE)
  } else {
    parTable <- var_interactions(parTable, ignore.means = FALSE) |>
      meanInteractions(ignore.means = FALSE)
  }

  if (is.null(intTerms))
    intTerms <- unique(parTable[grepl(":", parTable$rhs), "rhs"])
  if (is.null(etas))
    etas <- getEtas(parTable, isLV = TRUE)
  if (is.null(xis))
    xis  <- getXis(parTable, etas = etas, isLV = TRUE)

  xis <- unique(c(xis, intTerms)) # treat int-terms as normal xis
  lVs  <- c(xis, etas)

  indsLV <- getIndsLVs(parTable, lVs = lVs)
  inds <- unique(unlist(indsLV))

  # Create lambda
  lambda <- getLambdaParTable(parTable, rows = inds, cols = lVs, fill.missing = TRUE)

  # Create Gamma
  gammaXi <- matrix(0, nrow = length(etas), ncol = length(xis),
                    dimnames = list(etas, xis))
  gammaEta <- matrix(0, nrow = length(etas), ncol = length(etas),
                     dimnames = list(etas, etas))

  for (eta in etas) {
    reg <- parTable[parTable$lhs == eta & parTable$op == "~", , drop = FALSE]

    for (i in seq_len(NROW(reg))) {
      predictor <- reg[i, "rhs"]
      est       <- reg[i, "est"]

      if      (predictor %in% xis)  gammaXi[eta, predictor] <- est
      else if (predictor %in% etas) gammaEta[eta, predictor] <- est
      else warning("Unexpected type of predictor: ", predictor)
    }
  }

  createCov <- function(vars) {
    cov <- matrix(0, nrow = length(vars), ncol = length(vars),
                  dimnames = list(vars, vars))
    covRows <- parTable[parTable$op == "~~" & parTable$lhs %in% vars &
                        parTable$rhs %in% vars, , drop = FALSE]
    for (i in seq_len(NROW(covRows))) {
      lhs <- covRows$lhs[i]
      rhs <- covRows$rhs[i]
      est <- covRows$est[i]

      if (lhs %in% vars && rhs %in% vars) {
        cov[lhs, rhs] <- cov[rhs, lhs] <- est
      } else warning("Unexpected type of variable in covariance: ", lhs, " and ", rhs)
    }

    cov
  }

  psi   <- createCov(etas)
  phi   <- createCov(xis)
  theta <- createCov(inds)

  createBeta <- function(var) {
    beta <- matrix(0, nrow = length(var), ncol = 1,
                   dimnames = list(var, "~1"))

    betaRows <- parTable[parTable$op == "~1" & parTable$lhs %in% var, , drop = FALSE]
    for (i in seq_len(NROW(betaRows))) {
      lhs <- betaRows$lhs[i]
      est <- betaRows$est[i]

      if (lhs %in% var) {
        beta[lhs, "~1"] <- est
      } else warning("Unexpected type of variable in beta: ", lhs)
    }

    beta
  }

  alpha <- createBeta(etas)
  beta0 <- createBeta(xis)
  tau   <- createBeta(inds)

  Binv <- solve(diag(nrow(gammaEta)) - gammaEta)

  list(gammaXi = gammaXi, gammaEta = gammaEta, Binv = Binv, psi = psi,
       phi = phi, theta = theta, alpha = alpha, beta0 = beta0, tau = tau,
       lambda = lambda, inds = inds, xis = xis, etas = etas, lVs = lVs)
}


calcExpectedMatricesDA <- function(parTable, xis = NULL, etas = NULL, intTerms = NULL) {
  lapply(
    X = getGroupsParTable(parTable),
    FUN = function(g) {
      parTable.g <- parTable[parTable$group == g, , drop = FALSE]
      calcExpectedMatricesDA_Group(parTable.g, xis = xis, etas = etas, intTerms = intTerms)
    }
  )
}


calcExpectedMatricesDA_Group <- function(parTable, xis = NULL, etas = NULL, intTerms = NULL) {
  matricesCentered <- getCoefMatricesDA(parTable, xis = xis, etas = etas,
                                        intTerms = intTerms, centered = TRUE)
  matricesNonCentered <- getCoefMatricesDA(parTable, xis = xis, etas = etas,
                                           intTerms = intTerms, centered = FALSE)

  lVs  <- matricesCentered$lVs
  xis  <- matricesCentered$xis
  etas <- matricesCentered$etas
  inds <- matricesCentered$inds

  # Sigma ----------------------------------------------------------------------
  # Uses centered solution
  gammaXi  <- matricesCentered$gammaXi
  gammaEta <- matricesCentered$gammaEta
  phi      <- matricesCentered$phi
  psi      <- matricesCentered$psi
  Binv     <- matricesCentered$Binv
  tau      <- matricesCentered$tau
  lambda   <- matricesCentered$lambda
  alpha    <- matricesCentered$alpha
  beta0    <- matricesCentered$beta0
  theta    <- matricesCentered$theta

  covEtaEta <- Binv %*% (gammaXi %*% phi %*% t(gammaXi) + psi) %*% t(Binv)
  covEtaXi <- Binv %*% gammaXi %*% phi
  sigma.lv <- rbind(cbind(phi, t(covEtaXi)),
                    cbind(covEtaXi, covEtaEta))
  sigma.ov <- lambda %*% sigma.lv %*% t(lambda) + theta

  # lower left corner cov-lv-ov
  sigma.lv.ov <- lambda %*% sigma.lv
  sigma.ov.lv <- t(sigma.lv.ov)

  sigma.all <- rbind(cbind(sigma.lv, sigma.ov.lv),
                     cbind(sigma.lv.ov, sigma.ov))

  # Residuals and R^2 ----------------------------------------------------------
  # Uses centered solution
  eta.all <- c(etas, inds)
  var.eta.all <- diag(sigma.all[eta.all, eta.all, drop = FALSE])
  res.eta.all <- c(diag(psi), diag(theta))

  r2.all <- (var.eta.all - res.eta.all) / var.eta.all
  r2.lv  <- r2.all[etas]
  r2.ov  <- r2.all[inds]

  res.all <- 1 - r2.all
  res.lv  <- res.all[etas]
  res.ov  <- res.all[inds]

  # Mu -------------------------------------------------------------------------
  # Uses uncentered solution
  gammaXiNc  <- matricesNonCentered$gammaXi
  gammaEtaNc <- matricesNonCentered$gammaEta
  phiNc      <- matricesNonCentered$phi
  psiNc      <- matricesNonCentered$psi
  BinvNc     <- matricesNonCentered$Binv
  tauNc      <- matricesNonCentered$tau
  lambdaNc   <- matricesNonCentered$lambda
  alphaNc    <- matricesNonCentered$alpha
  beta0Nc    <- matricesNonCentered$beta0

  mu.eta <- BinvNc %*% (alphaNc + gammaXiNc %*% beta0Nc)
  mu.lv  <- rbind(beta0Nc, mu.eta)
  mu.ov  <- tauNc + lambdaNc %*% mu.lv
  mu.all <- rbind(mu.lv, mu.ov)

  list(
    sigma.all = sigma.all,
    sigma.lv  = sigma.lv,
    sigma.ov  = sigma.ov,
    mu.all    = mu.all,
    mu.lv     = mu.lv,
    mu.ov     = mu.ov,
    r2.all    = r2.all,
    r2.lv     = r2.lv,
    r2.ov     = r2.ov,
    res.all   = res.all,
    res.lv    = res.lv,
    res.ov    = res.ov,
    lambda    = lambda,
    gammaXi   = gammaXi,
    gammaEta  = gammaEta,
    psi       = psi,
    phi       = phi,
    theta     = theta
  )
}


getInternalCoefsDA <- function(model) {
  model$theta
}


getXisModelDA <- function(model) {
  covModelEtas <- model$covModel$info$etas
  allXis <- unique(c(model$covModel$info$xis, model$info$xis))
  allXis[!allXis %in% covModelEtas]
}


getEtasModelDA <- function(model) {
  unique(c(model$covModel$info$etas, model$info$etas))
}


splitParTableEtas <- function(parTable, parTableCov = NULL, splitEtas, allEtas,
                              nonLinearEtas = NULL) {
  if (!length(splitEtas) || !NROW(parTable))
    return(list(parTable = parTable, parTableCov = parTableCov))

  for (eta in splitEtas) {
    isReg <- parTable$lhs == eta & parTable$op == "~"
    isCov <- (parTable$lhs == eta | parTable$rhs == eta) & parTable$op == "~~"
    isLinear <- !parTable$lhs %in% nonLinearEtas &
                !parTable$rhs %in% nonLinearEtas

    parTableEta <- parTable[isReg | isCov & isLinear, ]

    vars <- unique(c(parTableEta$rhs, parTableEta$lhs))
    downstreamEtas <- vars[vars != eta & vars %in% allEtas]

    parTable    <- parTable[!(isReg | isCov & isLinear), ]
    parTableCov <- rbind(parTableCov, parTableEta)

    split <- splitParTableEtas(parTable = parTable,
                               parTableCov = parTableCov,
                               splitEtas = downstreamEtas,
                               allEtas = allEtas)

    parTable    <- split$parTable
    parTableCov <- split$parTableCov
  }

  list(parTable = parTable, parTableCov = parTableCov)
}


splitParTable <- function(parTable) {
  intTerms <- getIntTerms(parTable)
  empty <- list(parTable = parTable, parTableCov = NULL)
  if (!length(intTerms))
    return(empty)

  intVars  <- unique(unlist(stringr::str_split(intTerms, pattern = ":")))
  etas    <- getEtas(parTable)

  if (!length(intTerms))
    return(empty)

  isNonLinear <- vapply(
    X = etas,
    FUN.VALUE = logical(1L),
    FUN = intTermsAffectLV,
    parTable = parTable
  )

  nonLinearEtas <- etas[isNonLinear]
  linearEtas    <- etas[!isNonLinear]
  etasCovModel  <- intVars[intVars %in% linearEtas]

  if (!length(etasCovModel))
    return(empty)

  splitParTableEtas(
    parTable = parTable, allEtas = etas,
    splitEtas = etasCovModel,
    nonLinearEtas = nonLinearEtas
  )
}


sortParTableDA <- function(parTable, model) {
  parTable.input <- rbind(
    model$models[[1L]]$covModel$parTable, # by definition lower order etas
    model$models[[1L]]$parTable
  )

  etas.input <- getEtas(parTable.input)
  etas.model <- getEtasModelDA(model$models[[1L]]) # sorted to be lower-triangular in gammaEta
                                                   # and includes etas from cov-model
  etas <- unique(c(etas.input, etas.model))

  # xis <- getXisModelDA(model)
  # we don't want xis as they are sorted in the model
  # xis in the model are sorted to fit `OmegaXiXi`
  # and `OmegaEtaXi`, which doesn't match the sorting
  # in the model.syntax input. Instead we use getXis() on
  # parTable.input

  xis            <- getXis(parTable.input, checkAny = FALSE, etas = etas)
  indsXis        <- model$info$allIndsXis
  indsEtas       <- model$info$allIndsEtas
  higherOrderLVs <- model$info$higherOrderLVs

  isHigherOrderXi  <- xis  %in% higherOrderLVs
  isHigherOrderEta <- etas %in% higherOrderLVs

  xisLow   <- xis[!isHigherOrderXi]
  xisHigh  <- xis[isHigherOrderXi]
  etasLow  <- etas[!isHigherOrderEta]
  etasHigh <- etas[isHigherOrderEta]

  opOrder <- c("=~", "~", "~1", "~~", "|", ":=")
  varOrder <- unique(c(indsXis, indsEtas, xisLow, etasLow, xisHigh, xisLow))
  groupOrder <- c(getGroupsParTable(parTable), 0)

  getScore <- function(x, order.by) {
    order.by <- unique(c(order.by, x)) # ensure that all of x is in order.by
    mapping  <- stats::setNames(seq_along(order.by), nm = as.character(order.by))
    score    <- mapping[as.character(x)]

    if (length(score) != length(x)) {
      warning2("Sorting of parameter estimates failed!\n",
               immediate. = FALSE)

      return(seq_along(x))
    }

    score
  }

  scoreGroup <- getScore(x = parTable$group, order.by = groupOrder)
  scoreLhs   <- getScore(x = parTable$lhs, order.by = varOrder)
  scoreOp    <- getScore(x = parTable$op,  order.by = opOrder)
  scoreRhs   <- getScore(x = parTable$rhs, order.by = varOrder)

  out <- parTable[order(scoreGroup, scoreOp, scoreLhs, scoreRhs), , drop = FALSE]
  rownames(out) <- NULL

  out
}


updateStatusLog <- function(iterations, mode, logLikNew, deltaLL, relDeltaLL, verbose = FALSE) {
  if (verbose) {
    clearConsoleLine()
    printf("\rIter=%d Mode=%s LogLik=%.2f \u0394LL=%.2g rel\u0394LL=%.2g",
           iterations, mode, logLikNew, deltaLL, relDeltaLL)
  }
}


higherOrderStruct2Measr <- function(parTable, indsHigherOrderLVs) {
  if (is.null(indsHigherOrderLVs) || !length(indsHigherOrderLVs))
    return(parTable)

  higherOrderLVs <- names(indsHigherOrderLVs)

  for (lVh in higherOrderLVs) {
    inds <- indsHigherOrderLVs[[lVh]]
    mask <- parTable$rhs == lVh & parTable$op == "~" & parTable$lhs %in% inds

    if (!any(mask)) next

    lhs <- parTable[mask, "lhs"]
    rhs <- parTable[mask, "rhs"]

    parTable[mask, "lhs"] <- rhs
    parTable[mask, "op"]  <- "=~"
    parTable[mask, "rhs"] <- lhs
  }

  parTable
}


higherOrderMeasr2Struct <- function(parTable) {
  higherOrderLVs <- getHigherOrderLVs(parTable)

  if (is.null(higherOrderLVs) || !length(higherOrderLVs))
    return(parTable)

  for (lVh in higherOrderLVs) {
    mask <- parTable$lhs == lVh & parTable$op == "=~"

    if (!any(mask)) next

    lhs <- parTable[mask, "lhs"]
    rhs <- parTable[mask, "rhs"]

    parTable[mask, "lhs"] <- rhs
    parTable[mask, "op"]  <- "~"
    parTable[mask, "rhs"] <- lhs
  }

  parTable
}


recalcInterceptsY <- function(parTable.nlin, parTable.lin) {
  out <- NULL

  parTable.nlin <- addMissingGroups(parTable.nlin)
  parTable.lin  <- addMissingGroups(parTable.lin)

  for (g in getGroupsParTable(parTable.nlin)) {
    parTable.nlin.g <- parTable.nlin[parTable.nlin$group == g, , drop = FALSE]
    parTable.lin.g  <- parTable.lin[parTable.lin$group == g, , drop = FALSE]
    out <- rbind(out, recalcInterceptsY_Group(parTable.nlin = parTable.nlin.g,
                                              parTable.lin  = parTable.lin.g))
  }

  rbind(out, getZeroGroupParTable(parTable.nlin))
}


recalcInterceptsY_Group <- function(parTable.nlin, parTable.lin) {
  # fix intercept for indicators of endogenous variables, based on means
  # of interaction terms
  # intercepts are from a linear (H0) model, combined with a non-linear SAM
  # structural model. We want the mean structure to be coherent with those
  # from a full non-linear model
  parTable      <- parTable.nlin
  parTable.nlin <- meanInteractions(parTable.nlin, ignore.means = FALSE)

  for (eta in getEtas(parTable)) {
    meta.lin  <- getMean(eta, parTable = parTable.lin)
    meta.nlin <- getMean(eta, parTable = parTable.nlin)

    inds <- parTable[parTable$lhs == eta & parTable$op == "=~", "rhs"]
    inds <- unique(inds)

    ieta <- getIntercept(eta, parTable)

    if (ieta != 0) {
      cond <- parTable$op == "~1" & parTable$lhs == eta
      parTable[cond, "est"] <- parTable[cond, "est"] - (meta.nlin - meta.lin) # correct for the difference
      meta.lin  <- getMean(eta, parTable = parTable.nlin) # update for the rest of the code
    }

    for (ind in inds) {
      cond <- parTable$lhs == ind & parTable$op == "~1"
      lambda <- parTable[parTable$lhs == eta &
                         parTable$op == "=~" &
                         parTable$rhs == ind, "est"]

      if (!length(lambda) || !any(cond))
        next

      intercept <- parTable[cond, "est"]
      mu.lin  <- intercept + lambda * meta.lin
      mu.nlin <- intercept + lambda * meta.nlin
      parTable[cond, "est"] <- intercept - (mu.nlin - mu.lin) # correct for the difference
    }
  }

  parTable
}


getGroupsParTable <- function(parTable) {
  sort(unique(parTable$group[parTable$group > 0L]))
}


getZeroGroupParTable <- function(parTable) {
  parTable[parTable$group <= 0L, , drop = FALSE]
}


prepareDataGroupDA <- function(group, data, sampling.weights, sampling.weights.normalization) {
  if (!is.null(sampling.weights)) {
    stopif(length(sampling.weights) != 1L, "sampling.weights variable must be of length 1!")
    stopif(!sampling.weights %in% names(data),
           sprintf("sampling.weights variable '%s' not found in `data`.", sampling.weights))

    weights <- data[[sampling.weights]]
    stopif(any(is.na(weights)), "`sampling.weights` cannot have missing values!")
    stopif(any(weights < 0), "`sampling.weights` cannot have negative values!")
  }

  if ((!is.null(sampling.weights) && tolower(sampling.weights.normalization) == "total")) {
    weights <- data[[sampling.weights]]
    data[[sampling.weights]] <- NROW(data) * weights / sum(weights)
  }

  if (is.null(group)) {
    if ((!is.null(sampling.weights) && tolower(sampling.weights.normalization) == "group")) {
      weights <- data[[sampling.weights]]
      data[[sampling.weights]] <- NROW(data) * weights / sum(weights)
    }

    return(list(
      has.groups = FALSE,
      group = NULL,
      group.raw = NULL,
      group.var = NULL,
      levels = NULL,
      n.groups = 1L,
      indices = list(seq_len(NROW(data))),
      data = data,
      data.raw = data
    ))
  }

  data.raw <- data

  stopif(!is.data.frame(data), "`data` must be a data.frame when grouping is used.")
  n <- NROW(data)
  stopif(n == 0L, "Grouping requires non-empty data.")

  group.raw <- NULL
  group.var <- NULL

  stopif(length(group) != 1L, "Grouping variable must be of length 1!")
  stopif(!group %in% names(data),
         sprintf("Grouping variable '%s' not found in `data`.", group))

  group.raw <- data[[group]]
  data <- data[, setdiff(names(data), group), drop = FALSE]
  group.var <- group

  if ((!is.null(sampling.weights) && tolower(sampling.weights.normalization) == "group")) {
    weights <- data[[sampling.weights]]

    for (g in group.raw) {
      cond      <- group.raw == g
      weights.g <- data[cond, sampling.weights]
      data[cond, sampling.weights] <- sum(cond) * weights.g / sum(weights.g)
    }
  }

  stopif(length(group.raw) != n, "Length of `group` must match the number of rows in `data`!")
  stopif(any(is.na(group.raw)), "`group` cannot contain missing values!")

  group.raw <- as.character(group.raw) # match lavaan behaviour, ignore factor levels
  levels_order <- unique(group.raw)
  group_factor <- factor(group.raw, levels = levels_order)

  group.levels <- levels(group_factor)
  n.groups <- length(group.levels)
  indices <- split(seq_len(n), group_factor)

  list(
    has.groups = n.groups > 1L,
    group = group_factor,
    group.raw = group.raw,
    group.var = group.var,
    levels = group.levels,
    n.groups = n.groups,
    indices = indices,
    data = data,
    data.raw = data.raw
  )
}


expandGroupModifier <- function(mod, n.groups) {
  if (n.groups <= 1L)
    return(rep(mod, n.groups))

  if (length(mod) == 0L || mod == "")
    return(rep("", n.groups))

  if (is.na(mod))
    return(rep(NA_character_, n.groups))

  mod_trim <- trimws(mod)
  if (!nzchar(mod_trim))
    return(rep("", n.groups))

  is_c_call <- grepl("^c\\s*\\(", mod_trim) && grepl("\\)$", mod_trim)
  if (!is_c_call)
    return(rep(mod_trim, n.groups))

  inside <- substr(mod_trim,
                   start = regexpr("\\(", mod_trim, perl = TRUE) + 1L,
                   stop = nchar(mod_trim) - 1L)
  tokens <- strsplit(inside, ",", fixed = FALSE)[[1]]
  tokens <- trimws(tokens)

  if (!length(tokens)) {
    tokens <- rep("", n.groups)
  } else if (length(tokens) == 1L) {
    tokens <- rep(tokens, n.groups)
  } else {
    stopif(length(tokens) != n.groups,
           sprintf("Found %d modifiers but expected %d groups.", length(tokens), n.groups))
  }

  tokens
}


expandGroupModifiers <- function(mod, n.groups) {
  do.call(rbind, lapply(mod, FUN = expandGroupModifier, n.groups = n.groups))
}


expandParTableByGroup <- function(parTable, group.levels) {
  if (!NROW(parTable))
    return(parTable) # don't return NULL, as a zero-dim parTable has a special meaning
                     # for cov models

  n.groups <- length(group.levels)

  if (n.groups <= 1L) {
    if (!"mod" %in% names(parTable)) parTable$mod <- ""
    parTable$group <- 1L
    return(parTable[, c("lhs", "op", "rhs", "group", "mod")])
  }

  constraints <- parTable[parTable$op %in% CONSTRAINT_OPS, , drop = FALSE]
  baseTable   <- parTable[!parTable$op %in% CONSTRAINT_OPS, , drop = FALSE]

  if (!"mod" %in% names(baseTable)) baseTable$mod <- ""

  MOD <- expandGroupModifiers(mod = baseTable$mod, n.groups = n.groups)

  colsOut <- c("lhs", "op", "rhs", "group", "mod")
  parTableFull <- NULL
  for (i in seq_len(n.groups)) {
    parTable.g       <- baseTable
    parTable.g$mod   <- MOD[, i]
    parTable.g$group <- i

    parTable.g <- parTable.g[colsOut]
    parTableFull <- rbind(parTableFull, parTable.g)
  }

  if (NROW(constraints)) {
    if (!"mod" %in% names(constraints)) constraints$mod <- ""
    constraints$group <- 0L
    constraints <- constraints[colsOut]
  } else {
    constraints <- data.frame(lhs = character(),
                              op = character(),
                              rhs = character(),
                              group = integer(),
                              mod = character())
  }

  rbind(parTableFull, constraints)
}


parseModelArgumentsByGroupDA <- function(model.syntax, cov.syntax,
                                         method, data, group,
                                         auto.split.syntax = FALSE,
                                         sampling.weights = NULL,
                                         sampling.weights.normalization = "total",
                                         rng.indicators = 1) {
  parTable    <- modsemify(model.syntax)
  parTableCov <- modsemify(cov.syntax)

  # Check for observed (structural) variables
  structovs <- getStructOVs(rbind(parTable, parTableCov))
  ovs       <- getOVs(rbind(parTable, parTableCov))
  lvs       <- getLVs(parTable)
  all.vars  <- unique(c(ovs, lvs, colnames(data)))

  missing <- setdiff(ovs, colnames(data))
  stopif(length(missing), "Missing observed variables in data:\n  ",
         paste(missing, collapse = ", "))

  varsInts  <- getVarsInts(getIntTermRows(parTable), removeColonNames = FALSE)
  isOV_Int  <- vapply(varsInts, FUN.VALUE = logical(1L), FUN = \(x) all(x %in% ovs))
  ovIntTerms <- names(varsInts)[isOV_Int]
  ovIntTermMappings <- character(0L)

  for (ovInt in ovIntTerms) {
    vars <- varsInts[[ovInt]]

    ovIntNew1 <- stringr::str_replace_all(
      string = ovInt, pattern = ":",
      replacement = OP_OV_INT_1
    )

    ovIntNew2 <- stringr::str_replace_all(
      string = ovInt, pattern = ":",
      replacement = OP_OV_INT_2
    )

    if (!ovIntNew1 %in% all.vars) ovIntNew <- ovIntNew1
    else                          ovIntNew <- ovIntNew2

    warnif(ovIntNew %in% all.vars,
           sprintf("Overwriting %s variable in model...", ovIntNew))

    parTable[parTable$lhs == ovInt, "lhs"] <- ovIntNew
    parTable[parTable$rhs == ovInt, "rhs"] <- ovIntNew
    data[[ovIntNew]] <- apply(data[vars], MARGIN = 1L, FUN = prod)

    structovs <- c(structovs, ovIntNew)
    ovIntTermMappings[[ovIntNew]] <- ovInt
  }

  for (ov in structovs) {
    tmp.ov <- paste0(TEMP_OV_PREFIX, ov)

    # Replace ov in measurement model with tmp.ov
    # parTable[parTable$op == "=~" & parTable$rhs == ov, "rhs"] <- tmp.ov

    data[[tmp.ov]] <- data[[ov]]
    parTable <- rbind(
      parTable,
      data.frame(lhs = ov, op = "=~", rhs = tmp.ov, mod = "1"),
      data.frame(lhs = ov, op = "~1", rhs = "", mod = ""),
      data.frame(lhs = tmp.ov, op = "~1", rhs = "", mod = "0")
    )
  }

  if (length(structovs))
    model.syntax <- parTableToSyntax(parTable)

  group.info <- prepareDataGroupDA(group = group, data = data, sampling.weights = sampling.weights,
                                   sampling.weights.normalization = sampling.weights.normalization)

  if (auto.split.syntax && is.null(parTableCov) && is.null(cov.syntax)) {
    split <- splitParTable(parTable)

    parTable    <- split$parTable
    parTableCov <- split$parTableCov

    syntax     <- parTableToSyntax(parTable)
    cov.syntax <- parTableToSyntax(parTableCov)
  }

  group.info$syntax     <- model.syntax
  group.info$cov.syntax <- cov.syntax

  group.info$parTable.orig     <- parTable
  group.info$parTableCov.orig  <- parTableCov
  group.info$ovIntTerms        <- names(ovIntTermMappings)
  group.info$ovIntTermMappings <- ovIntTermMappings
  group.info$structovs         <- structovs

  group.levels <- group.info$levels
  if (is.null(group.levels)) group.levels <- ""

  group.info$group.levels <- group.levels
  group.info$parTable    <- expandParTableByGroup(parTable, group.levels = group.levels)
  group.info$parTableCov <- expandParTableByGroup(parTableCov, group.levels = group.levels)

  group.info
}


removeTempOV_RowsParTable <- function(parTable) {
  tmp <- startsWith(parTable$lhs, TEMP_OV_PREFIX) | startsWith(parTable$rhs, TEMP_OV_PREFIX)
  parTable[!tmp, , drop = FALSE]
}
