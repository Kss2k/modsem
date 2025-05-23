# Functions for tracing paths in SEMs, used both for calculating (co-)variances,
# as well as calculating formulas for (co-)variances.
prepParTable <- function(parTable, addCovPt = TRUE, maxlen = 100, paramCol = "mod") {
  # Remove any potential ':' from the model
  parTable <- lapplyDf(parTable, strRemovIfString, pattern = ":")
  parTable <- parTable[c("lhs", "op", "rhs", paramCol)]

  # get relevant variables
  structuralVars <- parTable[parTable$op == "~", c("lhs", "rhs")] |>
    unlist() |> c(getLVs(parTable)) |> unique()

  parTable <- parTable[(parTable$lhs %in% structuralVars & 
                       parTable$rhs %in% structuralVars) |
                       parTable$op == "~~", ]
  # redefine higher order (i.e.,`=~` -> `~`)
  if (any(parTable$op == "=~")) parTable <- redefineMeasurementModel(parTable)

  # add missing labels
  if (paramCol %in% c("mod", "label")) {
    labels <- apply(parTable[parTable[[paramCol]] == "", c("lhs", "op", "rhs")],
                    MARGIN = 1, FUN = stringr::str_c, collapse = "")
    parTable[parTable[[paramCol]] == "", paramCol] <- labels
  }

  if (addCovPt) {
    reversedCovPaths <- parTable[parTable$lhs != parTable$rhs & parTable$op == "~~", ]
    colnames(reversedCovPaths) <- c("rhs", "op", "lhs", paramCol)
    parTable <- rbind(parTable, reversedCovPaths)
  }

  parTable
}


tracePathsRecursively = function(x, y, pt, maxlen, currentPath = c(),
                                 covCount = 0, from = "lhs", depth = 1, 
                                 paramCol = "mod") {
  if (depth > maxlen) {
    warning2("Encountered a non-recursive model (infinite loop) when tracing paths")
    return(NULL)
  }

  if (covCount > 1) return(NULL)
  if (x == y && from == "rhs") return(list(currentPath))

  branches <- pt[pt[[from]] == x, ]
  if (NROW(branches) == 0) return(NULL)

  paths <- list()
  for (i in seq_len(nrow(branches))) {
    branch       <- branches[i, ]
    nextFrom     <- from
    nextCovCount <- covCount
    nextVar <- ifelse(from == "lhs", yes = branch$rhs, no = branch$lhs)

    # reverse direction of travel if covariance (but not before selecting next variable)
    if (branch$op == "~~") {
      nextCovCount <- covCount + 1
      nextFrom <- ifelse(from == "lhs", yes = "rhs", no = "lhs")
    }

    newPaths <- tracePathsRecursively(x = nextVar, y = y, pt = pt, maxlen = maxlen,
                                      currentPath = c(currentPath, branch[[paramCol]]),
                                      covCount = nextCovCount, from = nextFrom, 
                                      depth = depth + 1, paramCol = paramCol)
    paths <- c(paths, newPaths)

  }

  paths
}


cleanTracedPaths <- function(paths) {
  squaredExprs <- purrr::map_chr(paths, function(x) {
     reps <- as.data.frame(table(sort(x)))
     dplyr::if_else(reps[[2]] > 1,
                    true = paste(reps[[1]], reps[[2]], sep = " ^ "),
                    false = reps[[1]]) |>
     stringr::str_c(collapse = " * ")
  })
  reps <- as.data.frame(table(squaredExprs))
  dplyr::if_else(reps[[2]] > 1,
                 true = paste(reps[[2]], reps[[1]], sep = " * "),
                 false = reps[[1]]) |>
  stringr::str_c(collapse = " + ")
}


generateSyntax <- function(x, y, pt, maxlen = 100, parenthesis = TRUE, paramCol = "mod",...) {
  pt    <- prepParTable(pt, paramCol = paramCol, ...)
  paths <- tracePathsRecursively(x = x, y = y, pt = pt, maxlen = maxlen, paramCol = paramCol)

  if (!length(paths)) return(NA)
  cleaned <- cleanTracedPaths(paths)
  if (parenthesis) cleaned <- paste0("(", cleaned, ")")
  cleaned
}


addMissingCovariances <- function(pt, paramCol = "mod") {
  pt <- pt[!pt$op %in% c("=~", "~1"), ]

  xis  <- getXis(pt, checkAny = FALSE, isLV = FALSE)
  xis  <- xis[!grepl(":", xis)]

  if (length(xis)) {
    cxis <- getUniqueCombos(xis, match = TRUE)
    phi  <- data.frame(lhs = cxis[[1]], op = "~~", rhs = cxis[[2]], mod = "")
  } else phi <- NULL

  etas <- getEtas(pt, checkAny = FALSE, isLV = FALSE)
  if (length(etas)) {
    psi  <- data.frame(lhs = etas, op = "~~", rhs = etas, mod = "")
  } else psi <- NULL

  if (is.null(psi) && is.null(phi)) return(pt)

  covs <- rbind(psi, phi)
  for (i in seq_len(nrow(covs))) {
    row <- covs[i, , drop = FALSE]
    if (isRowInParTable(row = row, pt = pt, ignore = paramCol)) next
    pt <- rbind(pt, row)
  }

  pt
}


#' Estimate formulas for (co-)variance paths using Wright's path tracing rules
#'
#' @param pt A data frame with columns \code{lhs}, \code{op}, \code{rhs}, and \code{mod}, from \link{modsemify}
#' @param x Source variable
#' @param y Destination variable
#' @param parenthesis If \code{TRUE}, the output will be enclosed in parenthesis
#' @param missing.cov If \code{TRUE}, covariances missing from the model syntax will be added
#' @param measurement.model If \code{TRUE}, the function will use the measurement model
#' @param maxlen Maximum length of a path before aborting
#' @param paramCol The column name in \code{pt} that contains the parameter labels
#' @param ... Additional arguments passed to \link{trace_path}
#'
#' @return A string with the estimated path (simplified if possible)
#' @export
#' @description
#' This function estimates the path from \code{x} to \code{y} using the path tracing rules.
#' Note that it only works with structural parameters, so \code{"=~"} are ignored, unless
#' \code{measurement.model = TRUE}.
#' If you want to use the measurement model,
#' \code{"~"} should be in the \code{mod} column of \code{pt}.
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
#' pt <- modsemify(m1)
#' trace_path(pt, x = "Y", y = "Y", missing.cov = TRUE) # variance of Y
trace_path <- function(pt, x, y, parenthesis = TRUE, missing.cov = FALSE,
                       measurement.model = FALSE, maxlen = 100, 
                       paramCol = "mod", ...) {
  if (measurement.model) pt <- redefineMeasurementModel(pt)
  if (missing.cov) pt <- addMissingCovariances(pt, paramCol = paramCol)
  generateSyntax(x = x, y = y, pt = pt, maxlen = maxlen, 
                 parenthesis = parenthesis, paramCol = paramCol, ...)
}


getCovEqExpr <- function(x, y, parTable, paramCol = "label", ...) {
  parse(text=trace_path(x = x, y = y, pt = parTable, paramCol = paramCol, ...))
}
