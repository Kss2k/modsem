getExpectedMean <- function(Y, parTable) {
  alpha <- getIntercept(Y, parTable)
  X <- parTable[parTable$lhs == Y & 
                parTable$op == "~" &
                parTable$rhs != "1", "rhs"] |>
    unique()
  interceptsX <- vapply(X, FUN.VALUE = numeric(1L),
                        FUN = getIntercept, parTable = parTable) 

  gammaX <- vapply(X, FUN.VALUE = numeric(1L),
                   FUN = getGamma, X = Y, parTable = parTable) 

  alpha + sum(interceptsX * gammaX)
}


getIntercept <- function(X, parTable) {
  intercept <- parTable[parTable$lhs == X & 
                        parTable$op == "~" &
                        parTable$rhs == "1", "est"] 
  if (NROW(intercept) == 0) intercept <- 0
  intercept
}


getGamma <- function(Y, X, parTable) {
  gamma <- parTable[parTable$lhs == Y & 
                    parTable$op == "~" &
                    parTable$rhs == X, "est"]
  if (NROW(gamma) == 0) gamma <- 0
  gamma
}
