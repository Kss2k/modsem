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


calcSE <- function(hessian, names = NULL, NA__ = -999) {
  stdError <- tryCatch({
      sqrt(diag(solve(hessian)))
    }, error=function(e) {
      rep(NA, nrow(hessian))
    }, warning=function(w) {
       if (grepl("NaN", conditionMessage(w))) {
         suppressWarnings(sqrt(diag(solve(hessian))))
      } else{
         sqrt(diag(solve(hessian)))
      }
    })
  if (all(is.na(stdError))) 
    warning2("SE's could not be computed, negative Hessian is singular.")
  if (any(is.nan(stdError))) 
    warning2("SE's for some coefficients could not be computed.") 
  
  if (!is.null(names)) names(stdError) <- names
  stdError[is.na(stdError)] <- NA__
  stdError
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
  if (!is.null(parTableCovModel)) parTableCovModel$mod <- ""
  specifyModelLmsQml(parTable = parTable, method = method,
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


getConstrExprs <- function(parTable, parTableCov) {
  parTable <- rbind(parTable, parTableCov)
  rows <- parTable[parTable$op %in% c("==", ">", "<"), ]
  if (NROW(rows) == 0) return(NULL)

  fixedParams <- unique(rows$lhs)
  thetaFixedParams <- vector("list", length(fixedParams))
  names(thetaFixedParams) <- fixedParams

  exprs <- lapply(rows$rhs, function(expr) parse(text = expr))
  list(fixedParams = fixedParams, thetaFixedParams = thetaFixedParams,
       exprs = exprs)
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

