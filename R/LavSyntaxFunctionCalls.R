

# function for finding index for matching ) to a function call
findFunctionEnd <- function(listTokens, i = 1) {
  if (i > length(listTokens)) {
    stop("No matching parenthesis for function call")

  } else if (listTokens[[i]] == ")") {
    return(i)

  } else {
    findFunctionEnd(listTokens, i + 1)
  }
}


evalLavFunction <- function(listTokens) {

  do.call(listTokens[[1]], listTokens[-c(1, length(listTokens))])
}



`mean(` <- function(...) {
  varNames <- c(...)[c(...) != ","]

  if (length(varNames) <= 1) {
    stop2("Attempted to make a parcel out of a single variable!\n")
  }
  if (!is.data.frame(LavDataToBeModified)) {
    stop2("Data for parceling should be a dataframe\n")
  }
  parcelName <- stringr::str_c(c("MEAN", varNames), collapse = "_")
  # Modify dataset outside the scope of the function
  LavDataToBeModified[[parcelName]] <<- rowMeans(LavDataToBeModified[varNames])
  parcelName
}



`sum(` <- function(...) {
  varNames <- c(...)[c(...) != ","]

  if (length(varNames) <= 1) {
    stop2("Attempted to make a parcel out of a single variable!\n")
  }
  if (!is.data.frame(LavDataToBeModified)) {
    stop2("Data for parceling should be a dataframe\n")
  }
  parcelName <- stringr::str_c(c("SUM", varNames), collapse = "_")
  # Modify dataset outside the scope of the function
  LavDataToBeModified[[parcelName]] <<- rowSums(LavDataToBeModified[varNames])
  parcelName
}



`equal(` <- function(...) {
  paste0("equal(", ..., ")")
}



`start(` <- function(...) {
  paste0("start(", ..., ")")
}
