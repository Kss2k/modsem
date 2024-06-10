colNames <- c("lhs", "op", "rhs", "est", "std.error", 
              "z.value", "p.value", "ci.lower", "ci.upper")
header <- c("Variable", "op", "Variable", "Estimate",
            "Std.Error", "z.value", "Pr(>|z|)", "CI.Lower", "CI.Upper")


formatParTable <- function(parTable, digits = 3, scientific = FALSE,
                           ci = FALSE) {
  if (!ci) {
    header <- header[!grepl("CI", header)]
    colNames <- colNames[!grepl("ci", colNames)]
  }
  parTable <- parTable[colNames]
  
  isStructOrMeasure <- parTable$op %in% c("~", "=~", "~~") & 
    parTable$rhs != "1" & parTable$lhs != parTable$rhs
  parTable$lhs[isStructOrMeasure] <- 
    paste(parTable$lhs[isStructOrMeasure], parTable$op[isStructOrMeasure])

  allVars <- c(parTable$rhs[grepl("[[:alpha:]]", parTable$rhs)], parTable$lhs)
  maxWidth <- maxchar(allVars)

  parTable$rhs <- format(parTable$rhs, width = maxWidth, justify = "left")
  parTable$lhs <- format(parTable$lhs, width = maxWidth, justify = "left")
  if (scientific) {
    parTable$p.value <- format.pval(parTable$p.value)
  } else {
    parTable$p.value <- format(round(parTable$p.value, digits = 3), 
                               nsmall = 3)
  }

  for (i in seq_len(length(colNames) - 3) + 3) { # skip first 3 (lhs, op, rhs)
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
                          padWidth = 2, padWidthLhs = 2, 
                          padWidthRhs = 4, spacing = 2) {
  formatted <- formatParTable(parTable, digits = digits, 
                              ci = ci, scientific = scientific)
  fParTable <- formatted$parTable 
  header <- formatted$header
  lhs <- unique(fParTable$lhs)
  
  pad <- stringr::str_dup(" ", 2 * padWidth + padWidthLhs + 
                          padWidthRhs + maxchar(fParTable$rhs) - 1)
  space <- stringr::str_dup(" ", spacing)

  formattedHeader <- 
    paste0(pad, stringr::str_c(header[-(1:3)], collapse = space), "\n")

  # Measurement model 
  parTableLoadings <- fParTable[fParTable$op == "=~", ]
  if (loadings && NROW(parTableLoadings) > 0) {
    cat("Latent Variables:\n", formattedHeader)
    printParTableDouble(parTableLoadings, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        padWidthRhs = padWidthRhs, spacing = spacing)
  }

  # Regressions 
  parTableRegressions <- fParTable[parTable$op == "~" & parTable$rhs != "1", ]
  if (regressions && NROW(parTableRegressions) > 0) {
    cat("\nRegressions:\n", formattedHeader)
    printParTableDouble(parTableRegressions, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        padWidthRhs = padWidthRhs, spacing = spacing)
  }

  # Intercepts 
  parTableIntercepts <- fParTable[parTable$op == "~" & parTable$rhs == "1", ]
  if (intercepts && NROW(parTableIntercepts) > 0) {
    cat("\nIntercepts:\n", formattedHeader) 
    printParTableSingle(parTableIntercepts, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        padWidthRhs = padWidthRhs, spacing = spacing)
  }
  
  # Covariances 
  parTableCovariances <- fParTable[parTable$op == "~~" & parTable$lhs != parTable$rhs, ]
  if (covariances && NROW(parTableCovariances) > 0) {
    cat("\nCovariances:\n", formattedHeader)
    printParTableDouble(parTableCovariances, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        padWidthRhs = padWidthRhs, spacing = spacing)   
  }

  # Variances 
  parTableVariances <- fParTable[parTable$op == "~~" & parTable$lhs == parTable$rhs, ]
  if (variances && NROW(parTableVariances) > 0) {
    cat("\nVariances:\n", formattedHeader)
    printParTableSingle(parTableVariances, padWidth = padWidth, padWidthLhs = padWidthLhs, 
                        padWidthRhs = padWidthRhs, spacing = spacing)

  }
}


printParTableDouble <- function(parTable, padWidth = 2, padWidthLhs = 2, 
                                padWidthRhs = 4, spacing = 2) {
  lhs <- unique(parTable$lhs)
  pad <- stringr::str_dup(" ", padWidth) 

  for (l in lhs) {
    cat(paste0(pad, l), "\n")
    printRowsParTable(lhs = parTable[parTable$lhs == l, "rhs", drop = FALSE],
                      rhs = parTable[parTable$lhs == l, -(1:3), drop = FALSE],
                      padWithLhs = padWidthLhs + padWidth, 
                      padWithRhs = padWidthRhs + padWidth,
                      spacing = spacing)
  }
}


printParTableSingle <- function(parTable, padWidth = 2, padWidthLhs = 2, 
                                padWidthRhs = 4, spacing = 2) {
  lhs <- unique(parTable$lhs)
  pad <- stringr::str_dup(" ", padWidth) 

  printRowsParTable(lhs = parTable[ , "lhs", drop = FALSE],
                    rhs = parTable[ , -(1:3), drop = FALSE],
                    padWithLhs = padWidthLhs + padWidth, 
                    padWithRhs = padWidthRhs + padWidth,
                    spacing = spacing)
}


printRowsParTable <- function(lhs, rhs, padWithLhs = 2, padWithRhs = 4,
                              spacing = 2) {
  rhs <- lapplyDf(rhs, FUN = format) 
  padLhs <- stringr::str_dup(" ", padWithLhs)
  padRhs <- stringr::str_dup(" ", padWithRhs)
  space <- stringr::str_dup(" ", spacing)

  out <- ""
  for (i in seq_len(nrow(rhs))) {
    out <- paste0(out, padLhs, lhs[i, 1], padRhs, 
                  stringr::str_c(rhs[i, ], collapse = space), "\n")
  }
  cat(out)
}
