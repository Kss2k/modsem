colsOut <- c("lhs", "op", "rhs", "est", "std.error",
              "z.value", "p.value", "ci.lower", "ci.upper")
header <- c("Variable", "op", "Variable", "Estimate",
            "Std.Error", "z.value", "P(>|z|)", "ci.lower", "ci.upper")


shorten <- \(x, n) abbreviate(x, minlength = n, strict = TRUE)


formatParTable <- function(parTable, 
                           digits = 3, 
                           scientific = FALSE,
                           ci = FALSE, 
                           width.min = 14, 
                           pad.res = TRUE, 
                           extra.cols = NULL, 
                           shorten.lhs.header = FALSE) {
  parTable    <- fillColsParTable(parTable)
  parTable.in <- parTable # unformatted partable

  isCustom <- parTable$op == ":="
  isVar    <- parTable$op == "~~" & parTable$lhs == parTable$rhs
  isCov    <- parTable$op == "~~" & parTable$lhs != parTable$rhs
  isReg    <- parTable$op == "~"
  isIntr   <- parTable$op == "~1" # Intercept
  isMeasr  <- parTable$op == "=~"

  parTable[isCustom, "label"] <- ""

  isDoubleCol <- isReg  | isCov | isMeasr # Should output have double columns?
  isSingleCol <- isIntr | isVar | isCustom


  if (shorten.lhs.header) {
    lhsDouble <- paste(shorten(parTable$lhs[isDoubleCol], n = width - 2L),
                       parTable$op[isDoubleCol])
  } else {
    lhsDouble <- paste(parTable$lhs[isDoubleCol], parTable$op[isDoubleCol])
  }

  width <- max(width.min, maxchar(lhsDouble))

  lhsSingle <- pasteLabels(parTable$lhs[isSingleCol],
                           parTable$label[isSingleCol],
                           width = width - pad.res)
  rhsDouble <- pasteLabels(parTable$rhs[isDoubleCol],
                           parTable$label[isDoubleCol],
                           width = width - pad.res)

  parTable$p.value <- formatPval(parTable$p.value, scientific = scientific)

  
  ftext <- \(text, n = width - pad.res) format(text, width = n, justify = "left")

  parTable[isDoubleCol, "lhs"] <- ftext(lhsDouble)
  parTable[isDoubleCol, "rhs"] <- ftext(rhsDouble)
  parTable[isSingleCol, "lhs"] <- ftext(lhsSingle)
  parTable[isSingleCol, "rhs"] <- ""
  
  if (pad.res) {
    etas    <- getEtas(parTable.in, isLV = FALSE)
    inds    <- getInds(parTable.in)
    resvars <- c(etas, inds)

    isResVarCovOrInt <- (
      (isVar | isCov | isIntr) & 
      (parTable.in$rhs %in% resvars | parTable.in$lhs %in% resvars)
    )

    padVarCov <- ifelse(isResVarCovOrInt, yes = ".", no = " ")

    parTable$lhs <- paste0(padVarCov, parTable$lhs) 
    parTable$rhs <- paste0(padVarCov, parTable$rhs) 
  }

  if (!ci) {
    ci.pattern <- "^ci\\.(lower|upper)$"
    header  <- header[!grepl(ci.pattern, header)]
    colsOut <- colsOut[!grepl(ci.pattern, colsOut)]
  }

  if (!is.null(extra.cols)) {
    header  <- c(header, extra.cols)
    colsOut <- c(colsOut, extra.cols)
  }

  parTable      <- parTable[colsOut]
  names(header) <- colsOut

  valueCols <- setdiff(colnames(parTable), c("lhs", "op", "rhs"))

  for (col in valueCols) { # skip lhs, op and rhs
    if (is.numeric(parTable[[col]])) 
      parTable[[col]] <- round(parTable[[col]], digits)

    maxWidth        <- maxchar(c(header[[col]], parTable[[col]]))
    parTable[[col]] <- format(parTable[[col]], width = maxWidth,
                              digits = digits, justify = "right", 
                              scientific=FALSE)

    parTable[[col]] <- stringr::str_replace_all(parTable[[col]], "NA", "  ")
    header[[col]]   <- format(header[[col]], width = maxWidth, justify = "right")
  }

  list(parTable = parTable, header = header)
}


printParTable <- function(parTable,
                          scientific = FALSE,
                          ci = FALSE, 
                          digits = 3,
                          loadings = TRUE,
                          regressions = TRUE,
                          covariances = TRUE,
                          intercepts = TRUE,
                          variances = TRUE,
                          custom = TRUE,
                          extra.cols = NULL,
                          padWidth = 2,
                          padWidthLhs = 2,
                          spacing = 2) {
  formatted <- formatParTable(parTable, digits = digits,
                              ci = ci, scientific = scientific,
                              extra.cols = extra.cols)
  fParTable <- formatted$parTable
  header    <- formatted$header
  lhs       <- unique(fParTable$lhs)

  fullPadWidth <- padWidth + padWidthLhs + maxchar(fParTable$rhs) - 1
  pad          <- stringr::str_dup(" ", fullPadWidth)
  space        <- stringr::str_dup(" ", spacing)

  formattedHeader <- paste0(
    pad, stringr::str_c(header[-(1:3)], collapse = space), "\n"
  )

  # Measurement model
  parTableLoadings <- fParTable[fParTable$op == "=~", ]
  if (loadings && NROW(parTableLoadings) > 0) {
    cat("Latent Variables:\n", formattedHeader)

    printParTableDouble(
      parTableLoadings, 
      padWidth = padWidth, 
      padWidthLhs = padWidthLhs,
      spacing = spacing
    )
  }

  # Regressions
  parTableRegressions <- fParTable[parTable$op == "~", ]
  if (regressions && NROW(parTableRegressions) > 0) {
    cat("\nRegressions:\n", formattedHeader)

    printParTableDouble(
      parTableRegressions, 
      padWidth = padWidth, 
      padWidthLhs = padWidthLhs,
      spacing = spacing
    )
  }

  # Intercepts
  parTableIntercepts <- fParTable[parTable$op == "~1", ]
  if (intercepts && NROW(parTableIntercepts) > 0) {
    cat("\nIntercepts:\n", formattedHeader)

    printParTableSingle(
      parTableIntercepts,
      padWidth = padWidth, 
      padWidthLhs = padWidthLhs,
      spacing = spacing
    )
  }

  # Covariances
  parTableCovariances <- fParTable[parTable$op == "~~" & parTable$lhs != parTable$rhs, ]
  if (covariances && NROW(parTableCovariances) > 0) {
    cat("\nCovariances:\n", formattedHeader)

    printParTableDouble(
      parTableCovariances,
      padWidth = padWidth,
      padWidthLhs = padWidthLhs,
      spacing = spacing
    )
  }

  # Variances
  parTableVariances <- fParTable[parTable$op == "~~" & parTable$lhs == parTable$rhs, ]
  if (variances && NROW(parTableVariances) > 0) {
    cat("\nVariances:\n", formattedHeader)

    printParTableSingle(
      parTableVariances,
      padWidth = padWidth,
      padWidthLhs = padWidthLhs,
      spacing = spacing
    )
  }

  # Defined parameters
  parTableCustom <- fParTable[parTable$op == ":=", ]
  if (custom && NROW(parTableCustom) > 0) {
    cat("\nDefined Parameters:\n", formattedHeader)

    printParTableSingle(
      parTableCustom,
      padWidth = padWidth,
      padWidthLhs = padWidthLhs,
      spacing = spacing
    )
  }
  cat("\n")
}


printParTableDouble <- function(parTable, padWidth = 2, padWidthLhs = 2,
                                spacing = 2) {
  lhs <- unique(parTable$lhs)
  pad <- stringr::str_dup(" ", padWidth)
  for (l in lhs) {
    cat(paste0(pad, l), "\n")
    printRowsParTable(lhs = parTable[parTable$lhs == l, "rhs", drop = FALSE],
                      rhs = parTable[parTable$lhs == l, -(1:3), drop = FALSE],
                      padWidthLhs = padWidthLhs + padWidth,
                      spacing = spacing)
  }
}


printParTableSingle <- function(parTable, padWidth = 2, padWidthLhs = 2,
                                spacing = 2) {
  lhs <- unique(parTable$lhs)
  pad <- stringr::str_dup(" ", padWidth)

  printRowsParTable(lhs = parTable[ , "lhs", drop = FALSE],
                    rhs = parTable[ , -(1:3), drop = FALSE],
                    padWidthLhs = padWidthLhs + padWidth,
                    spacing = spacing)
}


printRowsParTable <- function(lhs, rhs, padWidthLhs = 2,
                              spacing = 2) {
  rhs    <- lapplyDf(rhs, FUN = format)
  padLhs <- stringr::str_dup(" ", padWidthLhs)
  space  <- stringr::str_dup(" ", spacing)

  out <- ""
  for (i in seq_len(nrow(rhs))) {
    fStrRhs <- stringr::str_c(rhs[i, ], collapse = space)
    out <- paste0(out, padLhs, lhs[i, 1], fStrRhs, "\n")
  }
  cat(out)
}


pasteLabels <- function(vars, labels, width = 14) {
  initWidthVar   <- round(width / 2)
  widthLabel     <- min(maxchar(labels), width - initWidthVar - 3L) # space + ()
  widthVar       <- width - widthLabel - 3L

  warnif(widthVar + widthLabel + 3L != width, "Mismatching widhts!")

  pasted <- paste0(vars, " (", labels, ")")
  widths <- nchar(pasted)

  wide         <- widths > width 
  vars[wide]   <- shorten(vars[wide], n = widthVar)
  labels[wide] <- shorten(labels[wide], n = widthLabel)

  labels[labels != ""] <- paste0("(", labels[labels != ""], ")")

  for (i in seq_along(vars)) {
    ncharVar   <- nchar(vars[[i]])
    ncharLabel <- nchar(labels[[i]])
    sep        <- stringr::str_dup(" ", width - ncharVar - ncharLabel)

    vars[[i]] <- paste0(vars[[i]], sep, labels[[i]])
  }
  vars
}


allignLhsRhs <- function(lhs, rhs, pad = "", width.out = 50) {
  if (length(lhs) != length(rhs)) {
    warning("lhs and rhs must have the same length")
    if (length(lhs) > length(rhs)) lhs <- rhs[seq_along(lhs)]
    else rhs <- lhs[seq_along(rhs)]
  }

  out <- ""
  width.out <- width.out - nchar(pad)
  for (i in seq_along(lhs)) {
    ncharLhs <- nchar(lhs[[i]])
    ncharRhs <- nchar(rhs[[i]])

    sep  <- stringr::str_dup(" ", max(0, width.out - ncharLhs - ncharRhs))
    line <- paste0(pad, lhs[[i]], sep, rhs[[i]], "\n")
    out  <- paste0(out, line)
  }
  out
}


# this is really ugly, but it is the easiest way to get the width of the
# printed table without splitting the function into multiple functions
# in a messy way
getWidthPrintedParTable <- function(parTable,
                                    scientific  = FALSE,
                                    ci          = FALSE,
                                    digits      = 3,
                                    loadings    = TRUE,
                                    regressions = TRUE,
                                    covariances = TRUE,
                                    intercepts  = TRUE,
                                    variances   = TRUE,
                                    padWidth    = 2,
                                    padWidthLhs = 2,
                                    spacing     = 2,
                                    extra.cols = NULL) {
  formatted <- formatParTable(
    parTable, 
    digits = digits,
    ci = ci, 
    scientific = scientific,
    extra.cols = extra.cols,
  )

  fParTable <- formatted$parTable
  header    <- formatted$header
  lhs       <- unique(fParTable$lhs)

  pad   <- stringr::str_dup(" ", padWidth + padWidthLhs +
                            maxchar(fParTable$rhs) - 1)
  space <- stringr::str_dup(" ", spacing)

  fStrHeader <- stringr::str_c(header[-(1:3)], collapse = space)
  formattedHeader <- paste0(pad, fStrHeader, "\n")
  nchar(formattedHeader)
}


formatPval <- function(p, scientific = TRUE) {
  if (scientific) return(format.pval(p))
  format(round(p, digits = 3), nsmall = 3)
}
