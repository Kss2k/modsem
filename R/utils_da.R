# Utils for lms approach
# Last updated: 31.07.2024


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
  return(tryCatch(mvnfast::dmvn(X, mean, sigma, log, ncores = 2), #ThreadEnv$n.threads),
                  error = function(e) mvtnorm::dmvnorm(X, mean, sigma, log)))
}


diagPartitionedMat <- function(X, Y) {
  if (is.null(X)) return(Y) else if (is.null(Y)) return(X)
  structure(rbind(cbind(X, matrix(0, nrow = NROW(X), ncol = NCOL(Y))), 
               cbind(matrix(0, nrow = NROW(Y), ncol = NCOL(X)), Y)),
            dimnames = list(c(rownames(X), rownames(Y)), 
                            c(colnames(X), colnames(Y))))
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
  if (!all(c(allIndsXis, allIndsEtas) %in% colnames(data))) 
    stop2("Missing Observed Variables in Data")
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
    if (anyAllNA(numericData)) stop2("Unable to conver data to type numeric") 
    numeric
  }, 
  error = function(e) {
    stop2("Unable to convert data to type numeric")
  })
  as.matrix(data)
}


filterData <- function(data) {
  completeCases <- stats::complete.cases(data)
  if (any(!completeCases)) warning2("Removing missing values case-wise.")
  data[completeCases, ]
}


cleanAndSortData <- function(data, allIndsXis, allIndsEtas) {
  if (is.null(data)) return(NULL)
  # sort Data before optimizing starting params
  sortData(data, allIndsXis,  allIndsEtas) |> 
    castDataNumericMatrix() |> filterData()
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
                          method = "lms") {
  parTable$mod <- ""
  parTable <- removeConstraintExpressions(parTable)

  if (!is.null(parTableCovModel)) {
    parTableCovModel$mod <- ""
    parTableCovModel <- removeConstraintExpressions(parTableCovModel)
  }

  specifyModelDA(parTable = parTable, method = method,
                 cov.syntax = cov.syntax,
                 parTableCovModel = parTableCovModel,
                 auto.constraints = FALSE, createTheta = FALSE,
                 checkModel = FALSE)
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


checkModel <- function(model, covModel = NULL, method = "lms") {
  modelXis <- model$info$xis
  if (!is.null(covModel$info)) {
    covModelEtas <- covModel$info$etas
    covModelXis <- covModel$info$xis 
    if (!all(c(covModelXis, covModelEtas) %in% modelXis)) {
      stop2("All latent variables in the cov-model must be an exogenous variable in the main model")
    }
    if (!all(modelXis %in% c(covModelXis, covModelEtas))) {
      stop2("All exogenous variables in main model must be part of the cov-model")
    }
  }

  checkNodesLms(parTable = rbind(model$parTable, covModel$parTable),
                nodes = model$quad$m, method = method)
  NULL
}


checkNodesLms <- function(parTable,
                          nodes, 
                          method = "lms",
                          minNodesXiXi = 16,
                          minNodesXiEta = 32, 
                          minNodesEtaEta = 48) {
  if (method == "qml") return(NULL)

  etas <- getEtas(parTable, isLV = TRUE)
  xis <- getXis(parTable, etas = etas, isLV = TRUE)
  varsInts <- getVarsInts(getIntTermRows(parTable))

  nodesXiXi_ok <- TRUE
  nodesXiEta_ok <- TRUE
  nodesEtaEta_ok <- TRUE
  
  lapply(varsInts, FUN = function(x) {
    if (all(x %in% xis)) nodesXiXi_ok <<- nodes >= minNodesXiXi
    else if (all(x %in% etas)) nodesEtaEta_ok <<- nodes >= minNodesEtaEta
    else if (any(x %in% etas)) nodesXiEta_ok <<- nodes >= minNodesXiEta
    else warning2("Unable to classify latent variables in interaction terms")
  }) 

  if (!nodesXiXi_ok) {
    warning2("It is recommended that you have at least ", minNodesXiXi,  " nodes ",
             "for interaction effects between exogenous variables in the lms approach ", 
             "'nodes = ", nodes, "'")
  } 
  if (!nodesXiEta_ok) {
    warning2("It is recommended that you have at least ", minNodesXiEta, " nodes ", 
             "for interaction effects between exogenous and endogenous variables in the lms approach ", 
             "'nodes = ", nodes, "'")
  }   
  if (!nodesEtaEta_ok) {
    warning2("It is recommended that you have at least ", minNodesEtaEta, " nodes ",
             "for interaction effects between endogenous variables in the lms approach ", 
             "'nodes = ", nodes, "'")
  }
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
