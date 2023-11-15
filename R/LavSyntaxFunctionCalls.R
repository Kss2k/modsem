
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
  functionCall <- stringr::str_c(unlist(listTokens), collapse = "")
  functionCall
  eval(rlang::parse_expr(functionCall), envir = modsemEnvironment)
}



LavMean <- function(...) {
  data <- data.frame(...)

  if (ncol(data) <= 1) {
    stop2("Attempted to make a parcel out of a single variable!\n")
  }
  if (!is.data.frame(data)) {
    stop2("Data for parceling should be a dataframe\n")
  }
  parcelName <- stringr::str_c(c("MEAN", colnames(data)), collapse = "_")
  # Modify dataset outside the scope of the function
  modsemEnvironment$data[[parcelName]] <- rowMeans(data)
  parcelName
}



LavSum <- function(...) {
  data <- data.frame(...)

  if (ncol(data) <= 1) {
    stop2("Attempted to make a parcel out of a single variable!\n")
  }
  if (!is.data.frame(data)) {
    stop2("Data for parceling should be a dataframe\n")
  }
  parcelName <- stringr::str_c(c("SUM", colnames(data)), collapse = "_")
  # Modify dataset outside the scope of the function
  modsemEnvironment$data[[parcelName]] <- rowSums(data)
  parcelName
}



LavEqual <- function(string) {
  paste0("equal(\"", string, "\")")
}



LavStart <- function(number) {
  paste0("start(", number, ")")
}



modsemEnvironment <- rlang::env(
  LavDataToBeModified = NULL,
  mean = LavMean,
  sum = LavSum,
  equal = LavEqual,
  start = LavStart
)
