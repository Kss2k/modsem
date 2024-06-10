# utils for all methods
calcCovParTable <- function(x, y, parTable, measurement.model = FALSE) {
  parTable$mod <- as.character(parTable$est)
  parTable <- parTable[c("lhs", "op", "rhs", "mod")]
  eval(parse(text = tracePath(parTable, x, y, 
                              measurement.model = measurement.model)))
}


reverseIntTerm <- function(xz) {
  if (length(xz) > 1) stop2("xz must be a single string")
  stringr::str_c(rev(stringr::str_split_1(xz, ":")), collapse = ":")
}


getEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  cond <- parTable$op == "~" & parTable$lhs != "1"
  if (isLV) {
    lVs <- unique(parTable[parTable$op == "=~", "lhs"])
    cond <- cond & parTable$lhs %in% lVs
  }
  etas <- unique(parTable[cond, "lhs"])
  if (checkAny && length(etas) == 0) stop("No etas found")
  etas
}


getSortedEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  structExprs <- parTable[parTable$op == "~" & parTable$rhs != "1", ]
  unsortedEtas <- getEtas(parTable, isLV = isLV, checkAny = checkAny)
  sortedEtas <- c()
  while (length(sortedEtas) < length(unsortedEtas) && nrow(structExprs) > 0) {
    if (all(unique(structExprs$lhs) %in% structExprs$rhs)) {
      stop("Model is non-recursive")
    }
    for (i in seq_len(nrow(structExprs))) {
      eta <- structExprs[i, "lhs"]
      if (eta %in% structExprs$rhs) next
      sortedEtas <- c(sortedEtas, eta)
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
  if (!isLV) {
    xis <- unique(parTable[parTable$rhs != "1" &
                  parTable$lhs %in% etas, "rhs"])
    return(xis)
  } 
  xis <- unique(parTable[parTable$op == "=~" &
                !parTable$lhs %in% etas, "lhs"])
  if (checkAny && length(xis) == 0) stop("No xis found")
  xis
}


getLVs <- function(parTable) {
  unique(parTable[parTable$op == "=~", "lhs"])
}


getIndsLVs <- function(parTable, lVs) {
  measrExprs <- parTable[parTable$op == "=~" & 
                         parTable$lhs %in% lVs, ]
  if (NROW(measrExprs) == 0) stop("No measurement expressions found, for", lVs)
  lapplyNamed(lVs, FUN = function(lV) measrExprs[measrExprs$lhs == lV, "rhs"],
              names = lVs)
}


getIntTerms <- function(parTable) {
  structExprs <- parTable[parTable$op == "~" & parTable$rhs != "1", ]
  structExprs[grepl(":", structExprs$rhs), ] 
}


getVarsInts <- function(intTerms) {
  lapplyNamed(intTerms$rhs, FUN = stringr::str_split_1, pattern = ":", 
              names = stringr::str_remove_all(intTerms$rhs, ":"))
}


maxchar <- function(x) {
  max(nchar(x), na.rm = TRUE)
}


fillColsParTable <- function(parTable) {
  colNames <- c("lhs", "op", "est", "std.error", "z.value", "p.value", 
                "ci.lower", "ci.upper")
  parTable[colNames[!colNames %in% colnames(parTable)]] <- NA
  parTable
}


#  function for getting unique combinations of two values in x
getUniqueCombos <- function(x, match = FALSE) {
  # Base case, x is 1 length long and there are no unique combos
  if (length(x) <= 1) {
    return(NULL)
  }
  rest <- getUniqueCombos(x[-1], match = FALSE)
  combos <- data.frame(V1 = rep(x[[1]], length(x) - 1),
                       V2 = x[-1])
  if (match) combos <- rbind(data.frame(V1 = x, V2 = x), combos)
  rbind(combos, rest)
}
