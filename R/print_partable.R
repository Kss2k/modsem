colsOut <- c("lhs", "op", "rhs", "est", "std.error", 
              "z.value", "p.value", "ci.lower", "ci.upper")
header <- c("Variable", "op", "Variable", "Estimate",
            "Std.Error", "z.value", "Pr(>|z|)", "CI.Lower", "CI.Upper")


formatParTable <- function(parTable, digits = 3, scientific = FALSE,
                           ci = FALSE, width = 14) {
  parTable <- fillColsParTable(parTable)

  isStructOrMeasure <- parTable$op %in% c("~", "=~", "~~") & 
    parTable$rhs != "1" & parTable$lhs != parTable$rhs
  parTable$lhs[isStructOrMeasure] <- 
    paste(parTable$lhs[isStructOrMeasure], parTable$op[isStructOrMeasure])


  isResVar <- parTable$op == "~~" & parTable$lhs == parTable$rhs
  parTable$lhs[parTable$rhs == "1" | isResVar] <- 
    pasteLabels(parTable$lhs[parTable$rhs == "1" | isResVar], 
                parTable$label[parTable$rhs == "1" | isResVar], 
                width = width)
  parTable$rhs[parTable$rhs != "1"] <- 
    pasteLabels(parTable$rhs[parTable$rhs != "1"], 
                parTable$label[parTable$rhs != "1"], width = width)

  parTable$lhs[!isStructOrMeasure] <- 
    format(parTable$lhs[!isStructOrMeasure], width = width, justify = "left")
  parTable$rhs <- format(parTable$rhs, width = width, justify = "left")
  parTable$p.value <- formatPval(parTable$p.value, scientific = scientific)

  if (!ci) {
    header <- header[!grepl("CI", header)]
    colsOut <- colsOut[!grepl("ci", colsOut)]
  }
  parTable <- parTable[colsOut]

  for (i in seq_len(length(colsOut) - 3) + 3) { # skip first 3 (lhs, op, rhs)
    if (is.numeric(parTable[[i]])) parTable[[i]] <- round(parTable[[i]], digits)
    maxWidth <- maxchar(c(header[[i]], parTable[[i]]))
    parTable[[i]] <- format(parTable[[i]], width = maxWidth, 
                            digits = digits, justify = "right")
    parTable[[i]] <- stringr::str_replace_all(parTable[[i]], "NA", "  ")
    header[[i]] <- format(header[[i]], width = maxWidth, justify = "right")
  }

  list(parTable = parTable, header = header)
}


printParTable <- function(parTable, 
                          scientific = FALSE, 
                          ci = FALSE, digits = 3, 
                          loadings = TRUE,
                          regressions = TRUE,
                          covariances = TRUE,
                          intercepts = TRUE,
                          variances = TRUE,
                          padWidth = 2, 
                          padWidthLhs = 2, 
                          spacing = 2) {
  formatted <- formatParTable(parTable, digits = digits, 
                              ci = ci, scientific = scientific)
  fParTable <- formatted$parTable 
  header <- formatted$header
  lhs <- unique(fParTable$lhs)
  
  pad <- stringr::str_dup(" ", padWidth + padWidthLhs + 
                          maxchar(fParTable$rhs) - 1)
  space <- stringr::str_dup(" ", spacing)

  formattedHeader <- 
    paste0(pad, stringr::str_c(header[-(1:3)], collapse = space), "\n")

  # Measurement model 
  parTableLoadings <- fParTable[fParTable$op == "=~", ]
  if (loadings && NROW(parTableLoadings) > 0) {
    cat("Latent Variables:\n", formattedHeader)
    printParTableDouble(parTableLoadings, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        spacing = spacing)
  }

  # Regressions 
  parTableRegressions <- fParTable[parTable$op == "~" & parTable$rhs != "1", ]
  if (regressions && NROW(parTableRegressions) > 0) {
    cat("\nRegressions:\n", formattedHeader)
    printParTableDouble(parTableRegressions, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        spacing = spacing)
  }

  # Intercepts 
  parTableIntercepts <- fParTable[parTable$op == "~" & parTable$rhs == "1", ]
  if (intercepts && NROW(parTableIntercepts) > 0) {
    cat("\nIntercepts:\n", formattedHeader) 
    printParTableSingle(parTableIntercepts, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        spacing = spacing)
  }
  
  # Covariances 
  parTableCovariances <- fParTable[parTable$op == "~~" & parTable$lhs != parTable$rhs, ]
  if (covariances && NROW(parTableCovariances) > 0) {
    cat("\nCovariances:\n", formattedHeader)
    printParTableDouble(parTableCovariances, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        spacing = spacing)   
  }

  # Variances 
  parTableVariances <- fParTable[parTable$op == "~~" & parTable$lhs == parTable$rhs, ]
  if (variances && NROW(parTableVariances) > 0) {
    cat("\nVariances:\n", formattedHeader)
    printParTableSingle(parTableVariances, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        spacing = spacing)

  }
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
  rhs <- lapplyDf(rhs, FUN = format) 
  padLhs <- stringr::str_dup(" ", padWidthLhs)
  space <- stringr::str_dup(" ", spacing)

  out <- ""
  for (i in seq_len(nrow(rhs))) {
    out <- paste0(out, padLhs, lhs[i, 1],  
                  stringr::str_c(rhs[i, ], collapse = space), "\n")
  }
  cat(out)
}


pasteLabels <- function(vars, labels, width = 14, widthVar = 7, widthLabel = 4) {
  pasted <- paste0(vars, " (", labels, ")")
  widths <- nchar(pasted)
  vars[widths > width] <- 
    abbreviate(vars[widths > width], minlength = widthVar)
  labels[widths > width] <- 
    abbreviate(labels[widths > width], minlength = widthLabel)
  labels[labels != ""] <- paste0("(", labels[labels != ""], ")")

  for (i in seq_along(vars)) {
    ncharVar <- nchar(vars[[i]]) 
    ncharLabel <- nchar(labels[[i]])
    sep <- stringr::str_dup(" ", width - ncharVar - ncharLabel)
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
    sep <- stringr::str_dup(" ", max(0, width.out - ncharLhs - ncharRhs))
    line <- paste0(pad, lhs[[i]], sep, rhs[[i]], "\n")
    out <- paste0(out, line)
  }
  out
}


# this is really ugly, but it is the easiest way to get the width of the 
# printed table without splitting the function into multiple functions 
# in a messy way
getWidthPrintedParTable <- function(parTable, 
                                    scientific = FALSE, 
                                    ci = FALSE, digits = 3, 
                                    loadings = TRUE,
                                    regressions = TRUE,
                                    covariances = TRUE,
                                    intercepts = TRUE,
                                    variances = TRUE,
                                    padWidth = 2, 
                                    padWidthLhs = 2, 
                                    spacing = 2) {
  formatted <- formatParTable(parTable, digits = digits, 
                              ci = ci, scientific = scientific)
  fParTable <- formatted$parTable 
  header <- formatted$header
  lhs <- unique(fParTable$lhs)
  
  pad <- stringr::str_dup(" ", padWidth + padWidthLhs + 
                          maxchar(fParTable$rhs) - 1)
  space <- stringr::str_dup(" ", spacing)

  formattedHeader <- 
    paste0(pad, stringr::str_c(header[-(1:3)], collapse = space), "\n")
  nchar(formattedHeader)
}


formatPval <- function(p, scientific = TRUE) {
  if (scientific) return(format.pval(p))
  format(round(p, digits = 3), nsmall = 3)
}
