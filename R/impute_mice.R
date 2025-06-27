imputeMissingData <- function(data) {
  Amelia::amelia(data, m = 1)
}
