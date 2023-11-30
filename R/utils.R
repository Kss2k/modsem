
# wrapper around stop
stop2 <- function(...) {
  stop(..., call. = FALSE)
}



# Wrapper around warning() with .call = FALSE and .immediate = TRUE
warning2 <- function(...) {
  warning(..., call.= FALSE, immediate. = TRUE)
}



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



# Wrapper of lapply where elements are names based on names argument, by default names
# are based on X
lapplyNamed <- function(X, FUN, ..., names = X) {

  structure(lapply(X, FUN, ...),
            names = names)
}




#' Wrapper of lapply for dataframes
#'
#' @param df a dataframe
#' @param FUN a function to apply to variables in df
#' @param ...
#'
#' @return
#' @export
#'
#' @examples lapplyDf(iris[1:4], FUN = sqrt)
lapplyDf <- function(df, FUN, ...) {
  structure(lapply(df, FUN, ...),
            names = names(df),
            row.names = 1:nrow(df),
            class = "data.frame")
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
    scale(x)

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
    warning2(
      "There were some duplicate product indicators, was this intended?\n",
      "The duplicates of these product indicators were removed: \n",
      duplicates, "\n")
  }

  combinedDf <- cbind.data.frame(
    listDf[[1]], listDf[[2]][,!matchingColnames, drop = FALSE]
    )

  combineListDf(c(list(combinedDf), listDf[-(1:2)]))

}



maxDepth <- function(list, max = 2, depth = 1) {

  if (is.null(list) | !is.list(list)) {
    return(depth)
  }

  if (depth > max) {
    stop("Incorrectly nested syntax")
  }
  deepest <- 1
  for (i in seq_along(list)) {
    branchDepth <- maxDepth(list[[i]], max = max, depth + 1)
    if (branchDepth > deepest) {
      deepest <- branchDepth
    }
  }
  deepest

}



capturePrint <- function(x) {
  paste(capture.output(print(x)), collapse = "\n")
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
  if (sum(as.integer(duplicateRows)) > 0) {
    warning("Some duplicates in the parTable was removed, have you accidentally ",
            "specified some of these in your syntax? \n",
            capturePrint(parTable[duplicateRows, ]))
  }
  rbind(parTable[!duplicateRows, ], newRows)
}
