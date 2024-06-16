# Utils for lms approach
# Last updated: 29.05.2024


sortData <- function(data, allIndsXis, allIndsEtas) {
  if (!all(c(allIndsXis, allIndsEtas) %in% colnames(data))) 
    stop2("Missing Observed Variables in Data")
  as.matrix(data[c(allIndsXis, allIndsEtas)])
}


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
dMvn <- function(X, mean, sigma, log = FALSE) {
  sigma <- round(sigma, 12) # makes it more reliable
  return(tryCatch(mvnfast::dmvn(X, mean, sigma, log, ncores = 2),
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


fillSymmetric <- function(mat, values) {
  mat[is.na(mat)] <- values
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}


getK_NA <- function(omegaEta) {
  sum(apply(omegaEta, 1, function(x) any(is.na(x))))
}


cleanAndSortData <- function(data, allIndsXis, allIndsEtas) {
  if (is.null(data)) return(NULL)
  # sort Data before optimizing starting params
  data <- sortData(data, allIndsXis,  allIndsEtas)
  completeCases <- stats::complete.cases(data)
  if (any(!completeCases)) warning2("Removing missing values case-wise.")
  data[completeCases, ]
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


getFreeOrConsIntTerms <- function(varsInInt, eta, intTerms) {
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
                     auto.constraints = FALSE, create.theta = FALSE)
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
  fixedParams <- unique(parTable[parTable$op %in% c("==", ">", "<"), ]$lhs)
  parTable[!parTable$lhs %in% fixedParams &
           !parTable$op %in% c("==", ">", "<") &
           !parTable$lhs %in% parTable$mod, ]
}


getLabeledParamsLavaan <- function(parTable, fixedParams = NULL) {
  if (is.null(parTable$label)) return(NULL)
  labelRows <- parTable[parTable$label != "" &
                        !parTable$label %in% fixedParams,
                        c("est", "label"), drop = FALSE] |> 
    lapply(unique)

  theta <- as.numeric(labelRows$est)
  names(theta) <- labelRows$label
  theta
}


checkModel <- function(model, covModel = NULL) {
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

  NULL
}


removeInteractionVariances <- function(parTable) {
  parTable[!(parTable$op == "~~" & grepl(":", parTable$lhs) & 
             grepl(":", parTable$rhs) & parTable$rhs == parTable$lhs), ]
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


#' @export
as.logical.matrix <- function(x, ...) {
  structure(x != 0, 
            dim = dim(x), 
            dimnames = dimnames(x))
}


isScalingY <- function(x) {
  seq_along(x) %in% which(x == 1 | x == 0)
}