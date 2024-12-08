warning2 <- function(...) {
  warning(..., call. = FALSE, immediate. = TRUE)
}


stop2 <- function(...) {
  stop(..., call. = FALSE)
}


stopif <- function(cond, ...) {
  if (cond) stop2(...)
}


warnif <- function(cond, ...) {
  if (cond) warning2(...)
}


# utils for all methods
calcCovParTable <- function(x, y, parTable, measurement.model = FALSE) {
  parTable$mod <- paste0("(", as.character(parTable$est), ")")
  parTable <- parTable[c("lhs", "op", "rhs", "mod")]
  eval(parse(text = trace_path(parTable, x, y, measurement.model = measurement.model)))
}


reverseIntTerm <- function(xz) {
  stopif(length(xz) > 1, "xz must be a single string")
  stringr::str_c(rev(stringr::str_split_1(xz, ":")), collapse = ":")
}


getEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  cond <- parTable$op == "~"
  if (isLV) {
    lVs <- unique(parTable[parTable$op == "=~", "lhs"])
    cond <- cond & parTable$lhs %in% lVs
  }

  etas <- unique(parTable[cond, "lhs"])
  stopif(checkAny && !length(etas), "No etas found")
  etas
}


getSortedEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  structExprs  <- parTable[parTable$op == "~", ]
  unsortedEtas <- getEtas(parTable, isLV = isLV, checkAny = checkAny)
  sortedEtas   <- character(0L)

  while (length(sortedEtas) < length(unsortedEtas) && nrow(structExprs) > 0) {
    stopif(all(unique(structExprs$lhs) %in% structExprs$rhs), "Model is non-recursive")

    for (i in seq_len(nrow(structExprs))) {
      if ((eta <- structExprs[i, "lhs"]) %in% structExprs$rhs) next

      sortedEtas  <- c(eta, sortedEtas)
      structExprs <- structExprs[!grepl(eta, structExprs$lhs), ]
      break
    }
  }

  if (!all(sortedEtas %in% unsortedEtas) &&
      length(sortedEtas) != length(unsortedEtas)) {
      warning("unable to sort etas")
      return(unsortedEtas)
  }

  sortedEtas
}


getXis <- function(parTable, etas = NULL, isLV = TRUE, checkAny = TRUE) {
  if (is.null(etas)) etas <- getEtas(parTable, isLV = isLV)
  # add all LVs which are not etas
  xis <- unique(parTable[parTable$op == "=~" & !parTable$lhs %in% etas, "lhs"])

  if (!isLV) { # add any other variabels found in structural expressions
    xis <- unique(c(xis, parTable[parTable$op == "~" &
                                  !parTable$rhs %in% etas, "rhs"]))
  }

  xis <- xis[!grepl(":", xis)] # remove interaction terms

  stopif(checkAny && !length(xis), "No xis found")
  xis
}
   

getIndicators <- function(parTable, observed=TRUE) {
  indicators <- unique(parTable[!grepl(":", parTable$rhs) & 
                                parTable$op == "=~", "rhs"])

  if (observed) indicators <- indicators[!indicators %in% getLVs(parTable)]
  indicators
}


getProdNames <- function(parTable) {
  unique(parTable[grepl(":", parTable$rhs) & 
         parTable$op %in% c("~", "=~"), "rhs"])
}


getLVs <- function(parTable) {
  unique(parTable[parTable$op == "=~", "lhs"])
}


getOVs <- function(parTable = NULL, model.syntax = NULL) {
  if (!is.null(model.syntax)) parTable <- modsemify(model.syntax)
  stopif(is.null(parTable), "Missing parTable")

  lVs    <- getLVs(parTable)
  select <- parTable$op %in% c("=~", "~", "~~")
  vars   <- unique(c(parTable$lhs[select], parTable$rhs[select]))

  vars[!vars %in% lVs & !grepl(":", vars)]
}


getHigherOrderLVs <- function(parTable) {
  lVs                  <- getLVs(parTable)
  isHigherOrder        <- logical(length(lVs))
  names(isHigherOrder) <- lVs

  for (lV in lVs) {
    inds <- parTable[parTable$lhs == lV & parTable$op == "=~", "rhs"] |>
      stringr::str_split(pattern = ":") |> unlist()
     
    if (any(inds %in% lVs)) isHigherOrder[[lV]] <- TRUE
  }

  lVs[isHigherOrder]
}


getIndsLVs <- function(parTable, lVs) {
  measrExprs <- parTable[parTable$op == "=~" & parTable$lhs %in% lVs, ]
  stopif(!NROW(measrExprs), "No measurement expressions found, for", lVs)
  lapplyNamed(lVs, FUN = function(lV) measrExprs[measrExprs$lhs == lV, "rhs"],
              names = lVs)
}


getInds <- function(parTable) {
  unique(unlist(getIndsLVs(parTable, lVs = getLVs(parTable))))
}


getIntTermRows <- function(parTable) {
  structExprs <- parTable[parTable$op == "~", ]
  structExprs[grepl(":", structExprs$rhs), ]
}


getIntTerms <- function(parTable) {
  structExprs <- parTable[parTable$op == "~", ]
  unique(structExprs[grepl(":", structExprs$rhs), "rhs"])
}


getVarsInts <- function(intTerms, removeColonNames = TRUE) {
  if (removeColonNames) names <- stringr::str_remove_all(intTerms$rhs, ":")
  else names <- intTerms$rhs
  lapplyNamed(intTerms$rhs, FUN = stringr::str_split_1, pattern = ":",
              names = names)
}


maxchar <- function(x) {
  max(nchar(x), na.rm = TRUE)
}


fillColsParTable <- function(parTable) {
  colNames <- c("lhs", "op", "rhs", "label", "est",
                "std.error", "z.value", "p.value", "ci.lower", "ci.upper")
  parTable[colNames[!colNames %in% colnames(parTable)]] <- NA
  parTable[colNames]
}


#  function for getting unique combinations of two values in x
getUniqueCombos <- function(x, match = FALSE) {
  # Base case, x is 1 length long and there are no unique combos
  if (length(x) <= 1) return(NULL)

  rest <- getUniqueCombos(x[-1], match = FALSE)
  combos <- data.frame(V1 = rep(x[[1]], length(x) - 1),
                       V2 = x[-1])
  if (match) combos <- rbind(data.frame(V1 = x, V2 = x), combos)
  rbind(combos, rest)
}


lastRow <- function(x) {
  x[nrow(x), ]
}


lapplyMatrix <- function(X, FUN, FUN.VALUE, ...) {
  matrix(vapply(X, FUN.VALUE = FUN.VALUE, FUN = FUN, ...),
         nrow = length(FUN.VALUE), ncol = length(X), ...)
}


# Wrapper of lapply where elements are names based on names argument, by default names
# are based on X
lapplyNamed <- function(X, FUN, ..., names = X) {
  structure(lapply(X, FUN, ...),
            names = names)
}


lapplyDf <- function(df, FUN, ...) {
  structure(lapply(df, FUN, ...),
            names = names(df),
            row.names = seq_len(nrow(df)),
            class = "data.frame")
}


isModsemObject <- function(x) {
  inherits(x, c("modsem_pi", "modsem_da", "modsem_mplus"))
}


getIntercept <- function(x, parTable) {
  if (length(x) > 1) stop2("x must be a single string")

  intercept <- parTable[parTable$lhs == x & parTable$op == "~1", "est"]

  if (length(intercept) == 0) return(0)
  intercept
}


getIntercepts <- function(x, parTable) {
  out <- vapply(x, FUN.VALUE = numeric(1L), FUN = function(x_i)
                getIntercept(x_i, parTable = parTable))
  names(out) <- x
  out
}


getMean <- function(x, parTable) {
  stopif(length(x) > 1, "x must be a single string")

  meanY <- getIntercept(x, parTable = parTable)
  gamma <- parTable[parTable$lhs == x & parTable$op == "~", , drop = FALSE]

  if (NROW(gamma) == 0) return(meanY)
  for (i in NROW(gamma)) {
    meanX <- getMean(gamma[i, "rhs"], parTable = parTable)
    meanY <- meanY + gamma[i, "est"] * meanX
  }

  meanY
}


getMeans <- function(x, parTable) {
  out <- vapply(x, FUN.VALUE = numeric(1L), FUN = function(x_i)
                getMean(x_i, parTable = parTable))
  names(out) <- x
  out
}


centerInteraction <- function(parTable) {
  rows <- getIntTermRows(parTable)
  for (i in NROW(rows)) {
    Y <- rows[i, "lhs"]
    XZ <- unlist(stringr::str_split(rows[i, "rhs"], ":"))
    X <- XZ[[1]]
    Z <- XZ[[2]]

    meanX <- getMean(X, parTable)
    meanZ <- getMean(Z, parTable)

    gammaXZ <- rows[i, "est"]
    gamma <- parTable[parTable$lhs == Y & parTable$op == "~", , drop = FALSE]
    gammaX <- gamma[gamma$rhs == X, "est"] + gammaXZ * meanZ

    gammaZ <- gamma[gamma$rhs == Z, "est"] + gammaXZ * meanX

    parTable[parTable$lhs == Y & parTable$op == "~" &
             parTable$rhs == X, "est"] <- gammaX

    parTable[parTable$lhs == Y & parTable$op == "~" &
             parTable$rhs == Z, "est"] <- gammaZ
  }

  parTable
}


getWarningWrapper <- function(silent = FALSE) { # function factory
  if (silent) return(suppressWarnings)
  function(x) x
}


isRowInParTable <- function(row, pt, ignore = NULL) {
  if (!is.null(ignore)) {
    row <- row[!names(row) %in% ignore]
    pt  <- pt[!colnames(pt) %in% ignore]
  }

  for (i in seq_len(nrow(pt))) {
    x <- unlist(row)
    y <- unlist(pt[i, ])
    if (all(x == y)) return(TRUE)
  }

  return(FALSE)
}


rename <- function(X, ...) {
  newNames <- list(...)
  oldNames <- names(newNames)

  for (old in oldNames) {
    names(X)[names(X) == old] <- newNames[[old]]
  }

  X
}


printf <- function(...) {
  cat(sprintf(...))
  utils::flush.console()
}
  

clearConsoleLine <- function() {
  printf(paste0("\r", strrep(" ", getOption("width", default=0L)), "\r"))
}


getDiffTwoMax <- function(x) {
  if (length(x) < 2) return(NA)
  y <- sort(x, decreasing = TRUE)
  y[[1]] - y[[2]]
}


stripColonsParTable <- function(parTable) {
  parTable$lhs <- stringr::str_remove_all(parTable$lhs, ":")
  parTable$rhs <- stringr::str_remove_all(parTable$rhs, ":")
  parTable
}


getMissingLabels <- function(parTable) {
  if (!"label" %in% colnames(parTable)) parTable$label <- ""
  missing <- is.na(parTable$label) | parTable$label == ""
  labels <- sprintf("%s%s%s", parTable$lhs, parTable$op, parTable$rhs)
  parTable[missing, "label"] <- labels[missing]
  parTable
}
