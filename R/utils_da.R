OP_REPLACEMENTS <- c("~~" = "___COVARIANCE___",
                     "=~" = "___MEASUREMENT___",
                     ":=" = "___CUSTOM___",
                     "~"  = "___REGRESSION___",
                     ":"  = "___INTERACTION___")
OP_REPLACEMENTS_INV <- structure(names(OP_REPLACEMENTS), names = OP_REPLACEMENTS)


CONSTRAINT_OPS <- c("==", ">", "<", ":=")
BOUNDUARY_OPS <- c(">", "<")


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
  if (is.null(X)) return(Y) else if (is.null(Y)) return(X)

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


formatNumeric <- function(x, digits = 3) {
  if (is.numeric(x)) {
    format(round(x, digits), nsmall = digits, digits = digits,
           trim = FALSE, justify = "right")
  } else {
    format(x, trim = FALSE, justify = "right")
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
  inds    <- c(allIndsXis, allIndsEtas)
  ovs     <- colnames(data)
  missing <- inds[!inds %in% ovs]

  stopif(!all(inds %in% ovs), "Missing observed variables in data:\n  ",
         missing)

  data[c(allIndsXis, allIndsEtas)]
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


handleMissingData <- function(data, impute.na = FALSE, seed = 123) {
  completeCases <- stats::complete.cases(data)
  anyMissing    <- any(!completeCases)

  if (!anyMissing) return(data)

  if (!impute.na) {
    warning2("Removing missing values case-wise!\nConsider using `impute.na = TRUE`, or impute yourself!")

    return(data[completeCases, ])

  } else {
    message("Imputing missing values. Consider imputing yourself!")
    set.seed(seed)

    imp  <- Amelia::amelia(data, m = 1, p2s = 0)
    imp1 <- as.matrix(as.data.frame(imp[[1]]))
    colnames(imp1) <- rownames(imp$mu)
  
    return(imp1)
  }
}


cleanAndSortData <- function(data, allIndsXis, allIndsEtas, impute.na = FALSE) {
  if (is.null(data)) return(NULL)
  # sort Data before optimizing starting params
  sortData(data, allIndsXis,  allIndsEtas) |>
    castDataNumericMatrix() |> 
    handleMissingData(impute.na = impute.na)
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
  if (canBeNumeric(expr, includeNA = TRUE)) return(as.numeric(expr))
  0
}


getLabelIntTerms <- function(varsInInt, eta, intTerms) {
  expr <- intTerms[intTerms$lhs == eta & intTerms$rhs %in%
                   createDoubleIntTerms(varsInInt), "mod"]
  if (!canBeNumeric(expr)) return(expr)
  ""
}


getEmptyModel <- function(parTable, cov.syntax, parTableCovModel,
                          mean.observed = TRUE, method = "lms") {
  parTable$mod <- ""
  parTable <- removeConstraintExpressions(parTable)

  if (NROW(parTableCovModel)) {
    parTableCovModel$mod <- ""
    parTableCovModel <- removeConstraintExpressions(parTableCovModel)
  }

  specifyModelDA(
    parTable         = parTable, method = method,
    cov.syntax       = cov.syntax,
    parTableCovModel = parTableCovModel,
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
  model$matrices <- lapply(model$matrices, function(x) {
    x[!is.na(x)] <- value
    x
  })
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


getDegreesOfFreedom <- function(m, coef) {
  t <- (m * (m + 1)) / 2
  df <- t - length(coef)
  nMeans <- sum(grepl("tau|alpha", names(coef)))
  df + nMeans
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
  labels_vcov <- colnames(vcov)
  labels_vv <- intersect(labels, labels_vcov) 
  labels_zz <- setdiff(labels, labels_vcov)

  m <- length(labels_vv)
  k <- length(labels_zz)

  Vvv <- vcov[labels_vv, labels_vv]
  Vzz <- matrix(0, nrow = k, ncol = k, dimnames = list(labels_zz, labels_zz))
  Vvz <- matrix(0, nrow = m, ncol = k, dimnames = list(labels_vv, labels_zz))

  V <- rbind(cbind(Vvv, Vvz), cbind(t(Vvz), Vzz))
  V[labels, labels] # sort
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


getLambdaParTable <- function(parTable, rows = NULL, cols = NULL) {
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
  !eta %in% predictors
}
