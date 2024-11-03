# function for selecting rows in a dataframe matching values on a specified column (type chr_vec)
selectRowsByCol <- function(value, df, column) {
  column <- as.character(column)
  out <- df[df[[column]] == value,]

  if (nrow(out) <= 0) return(NULL)
  else return(out)
}


# function for selecting values in the column of a dataframe matching values on a specified column (type chr_vec)
selectValuesByCol <- function(value, df, column.value, column.match) {
  column.value <- as.character(column.value)
  column.match <- as.character(column.match)

  # this should return a vector anyways, but here i force it
  out <- as.vector(df[[column.value]][df[[column.match]] == value])

  if (length(out) <= 0) return(NULL)
  else return(out)

}


# scale numeric vector is numeric
scaleIfNumeric <- function(x, scaleFactor = TRUE) {
  if (is.null(x)) {
    warning2("x in scaleIfNumeric was NULL")
    return(NULL)
  }
  if (scaleFactor == TRUE & is.factor(x)) {
    x <- as.numeric(x)
  }
  if (is.numeric(x)) {
    (x - mean(x, na.rm = TRUE))/stats::sd(x, na.rm = TRUE)
  } else x
}


# A fancy version of purrr::list_cbind()
  # This function will remove duplicate columns, before combining the
    # dataframes
  # note this function might be slow for large df's, I might create a C++
    # version of it
combineListDf <- function(listDf) {
  # This function should work recursively
  if (is.null(listDf) || length(listDf) < 1) {
    return(NULL)

  } else if (length(listDf) == 1)  {
    return(listDf[[1]])

    # This shouldnt really be necessary, but in the case that the function
      # recieves a dataframe, we should just return it
  } else if (is.data.frame(listDf))  {
    return(listDf)
  }

  # Basecase combine the first two columns of the df
    # check if there are matching colnames (the first df has priority)
  matchingColnames <- colnames(listDf[[2]]) %in% colnames(listDf[[1]])

  if (sum(as.integer(matchingColnames)) > 0) {
    duplicates <- stringr::str_c(colnames(listDf[[2]])[matchingColnames],
                                 collapse = ", ")
    warning2("There were some duplicate product indicators, was this intended?\n",
             "The duplicates of these product indicators were removed: \n", duplicates, "\n")
  }

  combinedDf <- cbind.data.frame(listDf[[1]], listDf[[2]][,!matchingColnames, drop = FALSE])

  combineListDf(c(list(combinedDf), listDf[-(1:2)]))
}


maxDepth <- function(list, max = 2, depth = 1) {

  if (is.null(list) | !is.list(list)) {
    return(depth)
  }

  stopif(depth > max, "Incorrectly nested syntax")
  deepest <- 1
  for (i in seq_along(list)) {
    branchDepth <- maxDepth(list[[i]], max = max, depth + 1)
    if (branchDepth > deepest) {
      deepest <- branchDepth
    }
  }
  deepest

}


capturePrint <- function(x, ...) {
  paste(utils::capture.output(print(x, ...)), collapse = "\n")
}


rbindParTable <- function(parTable, newRows) {
  # Merges rows of two partables, and replaces duplicates in the lhs partable
  # check for duplicates in lhs, op & rhs
  if (is.null(newRows) || NROW(newRows) == 0 || NCOL(newRows) == 0) {
    return(parTable)
  }
  newParTableRows <- apply(newRows[c("lhs", "op", "rhs")],
                        MARGIN = 1,
                        FUN = function(row)
                          list(row)) |>
    purrr::list_flatten()

  duplicateRows <- apply(parTable[c("lhs", "op", "rhs")],
        MARGIN = 1,
        FUN = function(row, parTableRows)
          list(row) %in% parTableRows,
        parTableRows = newParTableRows)

  warnif(sum(as.integer(duplicateRows)) > 0,
          "Some duplicates in the parTable was removed, have you accidentally ",
          "specified some of these in your syntax? \n",
          capturePrint(parTable[duplicateRows, ]))

  rbind(parTable[!duplicateRows, ], newRows)
}


greplRowDf <- function(col, df) {
  if (is.null(df)) return(FALSE)
  any(apply(df, MARGIN = 2, FUN = function(dfCol) all(sort(col) == sort(dfCol))))
}


`forceRowNames<-` <- function(df, value) {
  attr(df, "row.names") <- value
  df
}


# see functions for generating formulas for
# constrained approach
`%in_paired%` <- function(x, y) {
  out <- vector("logical", length(x))
  for (i in seq_along(x)) {
    matches <- y == x[[i]]
    if (any(matches)) {
      y <- y[!matches]
      out[[i]] <- TRUE
    } else {
      out[[i]] <- FALSE
    }
  }
  out
}


defineUndefinedLabels <- function(parTable.x, parTable.y) {
  # parTable.x = original parTable without new constraints
  # parTable.y = altered parTable with new constraints (and new labels)
  # goal: make sure user-specified labels are not overwritten

  parTable.o <- parTable.y # parTable out
  parTable.x <- rename(parTable.x, mod="mod.x")
  parTable.y <- rename(parTable.y, mod="mod.y")

  # means that we ha replaced a label
  parTable.z <- merge(parTable.y, parTable.x)
  parTable.z <- parTable.z[parTable.z$mod.y != "" &  parTable.z$mod.x != "" &
                           parTable.z$mod.x != parTable.z$mod.y &
                           parTable.z$op %in% c("=~", "~1", "~", "~~"), ]

  redefinedLabels <- getRedefinedLabels(parTable.z=parTable.z)

  rbind(redefinedLabels, parTable.o) # place redefinitions at the start, as to not confuse lavaan
}


getRedefinedLabels <- function(parTable.z) {
  redefinedLabels <- NULL
  for (i in seq_len(nrow(parTable.z))) {
    row <- parTable.z[i, , drop = FALSE]

    if (isLavLabelFunction(label = row$mod.x)) { # is function? then: skip and warn
      warnReplacingLabel(old = row$mod.x, new = row$mod.y, parTable.row = row)
      next
    }

    if (canBeNumeric(row$mod.x)) { # is numeric constant? then: constrain
      redef <- data.frame(lhs = row$mod.y, op = "==", rhs = row$mod.x, mod = "")
    } else { # is label? then: redefine using `:=`
      redef <- data.frame(lhs = row$mod.x, op = ":=", rhs = row$mod.y, mod = "")
    }

    redefinedLabels <- rbind(redefinedLabels, redef)
  }

  redefinedLabels
}


getVarsPI <- function(parTable) {
  unique(c(parTable$rhs[parTable$op %in% c("~", "=~")],
           parTable$lhs[parTable$op == "~"])) |>
  stringr::str_split(pattern = ":", simplify = FALSE) |>
  unlist() |> unique()
}


isLavLabelFunction <- function(label, context, warning = FALSE) {
  grepl("^(start|equal)\\(.*\\)$", label)
}


warnReplacingLabel <- function(old, new, parTable.row) {
  context <- paste(parTable.row$lhs, parTable.row$op,
                   paste0(old, "*", parTable.row$rhs))
  warning2("Replacing `", old, "` with new label `", new, "` in: `",
           context, "`")
}


checkHigherOrderInteractions <- function(elementsInProds, parTable) {
  higherOrderLVs <- getHigherOrderLVs(parTable)
  for (xz in elementsInProds) {
    stopif(any(xz %in% higherOrderLVs), "The ':' operator is not allowed ",
           "for higher order latent variables, please redefine the interaction term",
           "as a higher order latent variable using the '=~' operator.\n",
           "Run 'vignette(\"higher_order_interactions\")' for more information")
  }
}
