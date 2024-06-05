# utils for all methods
calcCovParTable <- function(parTable, x, y) {
  parTable$mod <- as.character(parTable$est)
  parTable <- parTable[c("lhs", "op", "rhs", "mod")]
  eval(parse(text = tracePath(parTable, x, y)))
}


reverseIntTerm <- function(xz) {
  if (length(xz) > 1) stop2("xz must be a single string")
  stringr::str_c(rev(stringr::str_split_1(xz, ":")), collapse = ":")
}
