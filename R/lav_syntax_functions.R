# function for finding index for matching ) to a function call
findFunctionEnd <- function(listTokens, i = 1) {
  if (i > length(listTokens)) {
    stop2("No matching parenthesis for function call")

  } else if (listTokens[[i]] == ")") {
    return(i)

  } else {
    findFunctionEnd(listTokens, i + 1)
  }
}


LavEqual <- function(string) {
  if (!is.character(string)) {
    stop2("Expected argument in equal() to be string, got: ", string)
  } else if (length(string) > 1) {
    stop2("Expected a single string in equal(), got: ", string)
  }
  paste0("equal(\"", string, "\")")
}


LavStart <- function(number) {
  if (!is.numeric(number)) {
    stop2("Expected argument in start() to be string, got: ", number)
  } else if (length(number) > 1) {
    stop2("Expected a single number in start(), got: ", number)
  }
  paste0("start(", number, ")")
}


LavConcat <- function(...) {
  as.character(substitute(expression(c(...))))[[2]]
}


modEnv <- rlang::env(
  equal = LavEqual,
  start = LavStart,
  c = LavConcat
)
