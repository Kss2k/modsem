# utils for all methods
calcCovParTable <- function(parTable, x, y) {
  parTable$mod <- as.character(parTable$est)
  parTable <- parTable[c("lhs", "op", "rhs", "mod")]
  eval(parse(text = tracePath(parTable, x, y)))
}
