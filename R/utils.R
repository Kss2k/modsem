warning2 <- function(..., call. = FALSE, immediate. = TRUE) {
  warning(..., call. = call., immediate. = immediate.)
}


stop2 <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}


stopif <- function(cond, ...) {
  if (cond) stop2(...)
}


warnif <- function(cond, ..., .newline = FALSE) {
  if (cond && .newline) cat("\n")
  if (cond) warning2(...)
}


calcCovParTable <- function(x, y, parTable, measurement.model = FALSE, maxlen = 100) {
  stopif(length(x) != length(y), "x and y must be the same length")

  if (measurement.model) parTable <- redefineMeasurementModel(parTable)
  parTable <- parTable[!parTable$op %in% c(":=", "~1", "=~"), ]
  parTable <- prepParTable(parTable[c("lhs", "op", "rhs", "est")], paramCol = "est")

  tracePathsNumericCpp(x = x, y = y, parTable = parTable)
}


calcVarParTable <- function(x, parTable, measurement.model = FALSE, maxlen = 100) {
  var <- calcCovParTable(x = x, y = x, parTable = parTable,
                         measurement.model = measurement.model, maxlen = maxlen)
  structure(var, names = x)
}


reverseIntTerm <- function(xz) {
  stopif(length(xz) > 1, "xz must be a single string")
  stringr::str_c(rev(stringr::str_split_1(xz, ":")), collapse = ":")
}


getEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  lVs <- unique(parTable[parTable$op == "=~", "lhs"])
  cond.lhs <- parTable$op == "~"
  cond.rhs <- parTable$op == "=~" & parTable$rhs %in% lVs

  if (isLV) cond.lhs <- cond.lhs & parTable$lhs %in% lVs

  etas.lhs <- parTable[cond.lhs, "lhs"]
  etas.rhs <- parTable[cond.rhs, "rhs"]

  etas <- unique(c(etas.lhs, etas.rhs))
  stopif(checkAny && !length(etas), "No etas found")

  etas
}


getSortedEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  unsortedEtas <- getEtas(parTable, isLV = isLV, checkAny = checkAny)

  cond1 <- parTable$op == "~"
  cond2 <- parTable$op == "=~" & parTable$rhs %in% unsortedEtas

  structExprs <- parTable[cond1, , drop = FALSE]
  measrExprs  <- parTable[cond2, , drop = FALSE]

  if (NROW(measrExprs)) {
    measr2struct <- measrExprs
    measr2struct$lhs <- measrExprs$rhs
    measr2struct$op  <- "~"
    measr2struct$rhs <- measrExprs$lhs

    structExprs <- rbind(structExprs, measr2struct)
  }

  sortedEtas  <- character(0L)

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

  if (!any(isHigherOrder)) NULL else lVs[isHigherOrder]
}


isHigherOrderParTable <- function(parTable) {
  length(getHigherOrderLVs(parTable)) > 0
}


isClustered <- function(object) {
  length(modsem_inspect(object, what = "cluster")) > 0
}


getIndsLVs <- function(parTable, lVs, isOV = FALSE, ovs = NULL) {
  if (!length(lVs)) return(NULL)

  measr <- parTable[parTable$op == "=~" & parTable$lhs %in% lVs, ]
  stopif(!NROW(measr), "No measurement expressions found, for", lVs)

  if (isOV) .f <- \(lV) measr[measr$lhs == lV & measr$rhs %in% ovs, "rhs"]
  else      .f <- \(lV) measr[measr$lhs == lV, "rhs"]

  lapplyNamed(lVs, FUN = .f, names = lVs)
}


getInds <- function(parTable, ...) {
  unique(unlist(getIndsLVs(parTable, lVs = getLVs(parTable)), ...))
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
  parTable #[colNames]
}


#  function for getting unique combinations of two values in x
getUniqueCombos <- function(x, match = FALSE) {
  # Base case, x is 1 length long and there are no unique combos
  if (length(x) <= 1) return(NULL)

  rest <- getUniqueCombos(x[-1], match = FALSE) # this is not a mistake!
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


getIntercept <- function(x, parTable, col = "est") {
  if (length(x) > 1) stop2("x must be a single string")

  intercept <- parTable[parTable$lhs == x & parTable$op == "~1", col]

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
  for (i in seq_len(NROW(gamma))) {
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


centerInteractions <- function(parTable, center.means = TRUE) {
  parTable <- meanInteractions(parTable)
  rows <- getIntTermRows(parTable)

  for (i in seq_len(NROW(rows))) {
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

  if (center.means) {
    innerVars <- unique(unlist(parTable[parTable$op == "~", c("rhs", "lhs")]))
    parTable[parTable$lhs %in% innerVars & parTable$op == "~1", "est"] <- 0
  }

  parTable
}


meanInteractions <- function(parTable, ignore.means = FALSE) {
  intTerms <- unique(parTable[grepl(":", parTable$rhs), "rhs"])

  # remove existing
  parTable <- parTable[!(parTable$op == "~1" & parTable$lhs %in% intTerms),
                       , drop = FALSE]

  present  <- colnames(parTable)
  newcols  <- c("lhs", "op", "rhs", "label", "est")
  newcols  <- intersect(newcols, present)
  fillcols <- setdiff(present, newcols)

  for (intTerm in intTerms) {
    XZ <- unlist(stringr::str_split(intTerm, ":"))
    X <- XZ[[1]]
    Z <- XZ[[2]]

    meanX <- if (!ignore.means) getMean(X, parTable) else 0
    meanZ <- if (!ignore.means) getMean(Z, parTable) else 0

    covXZ <- calcCovParTable(x = X, y = Z, parTable = parTable)
    meanXZ <- meanX * meanZ + covXZ

    newRow <- data.frame(lhs = intTerm, op = "~1", rhs = "",
                         label = "", est = meanXZ)

    # Get correct cols
    newRow <- newRow[newcols]
    newRow[fillcols] <- NA

    parTable <- rbind(parTable, newRow)
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


rename <- function(.X, ...) {
  newNames <- list(...)
  oldNames <- names(newNames)

  for (old in oldNames) {
    names(.X)[names(.X) == old] <- newNames[[old]]
  }

  .X
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
  abs.diff <- y[[1]] - y[[2]]
  rel.diff  <- (y[[1]] - y[[2]]) / y[[2]]

  list(rel = rel.diff, abs = abs.diff)
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


strRemovIfString <- function(x, pattern) {
  if (is.character(x)) stringr::str_remove_all(x, pattern=pattern) else x
}


getParTableLabels <- function(parTable, labelCol="label", replace.dup = FALSE) {
  if (replace.dup) {
		labels <- unique(parTable[[labelCol]][parTable[[labelCol]] != ""])

		for (label in labels) {
			match <- parTable[[labelCol]] == label
			first <- which.max(match)
			parTable[match, labelCol] <- ""
			parTable[first, labelCol] <- label
		}
  }

  custom <- parTable$op == ":="
  parTable[custom, c("op", "rhs")] <- ""

  ifelse(parTable[[labelCol]] == "",
         yes = paste0(parTable$lhs, parTable$op, parTable$rhs),
         no = parTable[[labelCol]])
}


CI_WIDTH <- qnorm(.975)


padCharMatrix <- function(X, n=1) {
  pad <- strrep(" ", n)
  matrix(paste0(pad, X), nrow = nrow(X), ncol = ncol(X),
         dimnames = dimnames(X))
}


timeExpr <- function(expr) {
  start <- Sys.time()
  out <- expr
  end <- Sys.time()
  print(end - start)
  invisible(out)
}


formatParameters <- function(...) {
  pad <- "  "

  args <- as.character(substitute(list(...)))[-1]
  names <- format(args, justify = "left")
  values <- format(c(...), justify = "right")

  paste0(sprintf("%s%s = %s\n", pad, names, values),
         collapse = "")
}


padString <- function(x, pad) {
  s <- stringr::str_replace_all(x, pattern = "\n", replacement = paste0("\n", pad))
  paste0(pad, s, "\n")
}


getCoefsRows <- function(y, x, parTable) {
  row <- parTable[parTable$lhs %in% y &
                  parTable$rhs %in% x &
                  parTable$op == "~", , drop = FALSE]
}


getCoefs <- function(y, x, parTable, col = "est") {
  coefRows <- getCoefsRows(y, x, parTable)
  coefs <- structure(coefRows[[col]], names = coefRows$rhs)
  coefs[x]
}


getCOEFS <- function(y, x, COEFS, parTable = NULL) {
  if (is.null(parTable)) {
    y <- stringr::str_replace_all(y, ":", OP_REPLACEMENTS[[":"]])
    x <- stringr::str_replace_all(x, ":", OP_REPLACEMENTS[[":"]])
    labels <- paste0(y, OP_REPLACEMENTS[["~"]], x)

  } else {
    rows   <- getCoefsRows(y, x, parTable)
    names  <- paste0(rows$lhs, rows$op, rows$rhs)
    labels <- structure(rows$label, names = names)
    labels <- labels[paste0(y, "~", x)]
  }

  if (length(labels) == 1) COEFS[[labels]] else COEFS[labels]
}


parTableToEnv <- function(parTable, prefix = "") {
  lhs <- parTable$lhs
  rhs <- parTable$rhs
  op <- stringr::str_replace_all(parTable$op, OP_REPLACEMENTS)

  names <- paste0(prefix, lhs, op, rhs)
  values <- structure(parTable$est, names = names)

  list2env(values)
}


getConstrExprs <- function(parTable, parTableCov) {
  parTable <- unique(rbind(parTable, parTableCov))
  rows <- sortConstrExprs(parTable)
  if (is.null(rows)) return(NULL)

  fixedParams <- unique(rows$lhs)
  thetaFixedParams <- vector("list", length(fixedParams))
  names(thetaFixedParams) <- fixedParams

  exprs <- lapply(rows$rhs, function(expr) parse(text = expr))
  list(fixedParams = fixedParams, thetaFixedParams = thetaFixedParams,
       exprs = exprs)
}


sortConstrExprs <- function(parTable) {
  rows <- parTable[parTable$op %in% CONSTRAINT_OPS, ]
  if (NROW(rows) == 0) return(NULL)

  labelled <- unique(parTable$mod[parTable$mod != ""])
  isConst  <- canBeNumeric(rows$rhs)

  rows <- rows[!(isConst & rows$op %in% BOUNDUARY_OPS), ] # not relevant

  if (!all(rows$lhs %in% labelled)) {
    stop2("Unknown labels in constraints: ", rows$lhs[!rows$lhs %in% labelled])

  } else if (length(unique(rows$lhs)) != length(rows$lhs)) {
    stop2("Duplicated labels in constraints:\n", capturePrint(table(rows$lhs)))

  } else if (any(rows$op %in% BOUNDUARY_OPS)) {
    stop2("Dynamic constraints with ('<', '>') are not implemented yet!")
  }

  definedLabels <- labelled[!labelled %in% rows$lhs]
  subRows <- rows
  sortedRows <- data.frame(lhs = NULL, op = NULL, lhs = NULL, mod = NULL)
  while (NROW(subRows) > 0) {
    matchedAny <- FALSE
    for (i in seq_len(nrow(subRows))) {
      labels_i <- getVarsExpr(subRows[i, "rhs"])
      if (length(labels_i) == 0 || all(labels_i %in% definedLabels)) {
        matchedAny <- TRUE
        sortedRows <- rbind(sortedRows, subRows[i, ])
        definedLabels <- c(definedLabels, subRows[i, "lhs"])
        subRows <- subRows[-i, ]
        break
      }
    }

    stopif(!matchedAny, "Unkown labels in constraints: ",
           labels_i[!labels_i %in% definedLabels])
  }

  if (NROW(sortedRows) != NROW(rows)) {
    warning2("Something went wrong when sorting and parsing constraint-expressions ",
             "attempting to estimate model with unsorted expressions")
    return(rows)
  }

  sortedRows
}


sortConstrExprsFinalPt <- function(parTable) {
  if (!NROW(parTable)) return(NULL)

  # wrapr for sortConstrExprs meant to be used on the final output parTable
  # not for the input parTable (e.g., mod -> label)
  constrExprs <- sortConstrExprs(rename(parTable, label="mod"))

  if (is.null(constrExprs)) NULL else rename(constrExprs, mod="label")
}


subsetByGrouping <- function(df, grouping = NULL) {
  if (is.null(grouping)) return(df)
  grouping <- grouping[intersect(names(grouping), colnames(df))]
  df[apply(df[names(grouping)], MARGIN = 1, FUN = \(x) all(x == grouping)), ]
}


GINV <- function(X) {
  structure(MASS::ginv(X), dimnames = dimnames(X))
}


cov2cor <- function(vcov) {
  if (is.null(vcov))
    return(NULL)

  sd <- sqrt(abs(diag(vcov))) # use `abs()`, in case some variances are negative

  D <- diag(1 / sd)
  structure(D %*% vcov %*% D, dimnames = dimnames(vcov))
}


leftJoin <- function(left, right, by = intersect(colnames(left), colnames(right))) {
  ol <- "__orig_order_right__"
  or <- "__orig_order_left__"

  left[[ol]] <- seq_len(nrow(left))
  right[[or]] <- seq_len(nrow(right))

  # left join
  joined <- merge(left, right, by = by, all.x = TRUE)

  # order
  ordered <- joined[order(joined[[ol]], joined[[or]]), ]

  # return
  ordered[!colnames(ordered) %in% c(ol, or)]
}


eraseConsoleLines <- function(n = 1L) {
  if (n < 1L || !interactive() || !.isOnUnix())
    return(invisible()) # if we're not in interactive mode, or on UNIX
                        # there is not point

  seq <- paste0(rep("\033[1A\033[2K", n), collapse = "")
  # bring cursor to start of the line we landed on
  cat(seq, "\033[1G", sep = "")
}


tryFormatC <- function(..., .onfail = "NA") {
  tryCatch(formatC(...), error = \(e) .onfail)
}


is.diag <- function(M, zero.tol = .Machine$double.eps) {
  diag(NROW(M)) > zero.tol
}


is.invertible <- function(M) {
  tryCatch({
    solve(M)
    TRUE
  }, error = \(e) FALSE)
}


.isOnUnix <- function(.onFail = FALSE) {
  tryCatch(tolower(.Platform$OS.type) == "unix",
           error = \(e) .onFail)
}
