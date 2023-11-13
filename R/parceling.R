


getParcelInfo <- function(syntax) {
  # Gather information about parcels in syntax
  types = c("mean", "sum")
  info <- lapplyNamed(types,
                      FUN = extractParcels,
                      syntax = syntax,
                      names = types)
  if (length(info$mean$parcelExpressions) <= 0 &
      length(info$sum$parcelExpressions) <= 0) {
    return(NULL)
  }

  info
}



computeParcels <- function(parcelInfo, data) {
  sums <- parcelInfo$sum
  means <- parcelInfo$mean

  sumParcels <- structure(lapply(sums$variablesInParcels,
                                      FUN = computeRowSums,
                                      data = data),
                          names = sums$parcelNames,
                          class = "data.frame",
                          row.names = row.names(data))
  meanParcels <- structure(lapply(means$variablesInParcels,
                                  FUN = computeRowMeans,
                                  data = data),
                          names = means$parcelNames,
                          class = "data.frame",
                          row.names = row.names(data))

  combineListDf(list(data, sumParcels, meanParcels))

}



fixParcelSyntax <- function(syntax, parcelInfo) {
  parcelExpressions <- c(parcelInfo$sum$parcelExpressions,
                         parcelInfo$mean$parcelExpressions)

  replacements <- c(parcelInfo$sum$parcelNames,
                    parcelInfo$mean$parcelNames)

  for (i in seq_along(parcelExpressions)) {
    syntax <-  stringr::str_replace_all(syntax,
                                        pattern = parcelExpressions[[i]],
                                        replacement = replacements[[i]])
  }

  syntax
}



extractParcels <- function(type = "mean", syntax) {
  parcelExpressions <- stringr::str_extract_all(
    syntax, paste0(type, "\\(([^)]+)\\)")) |>
    unlist()
  regexParcelExpressions <- stringr::str_replace_all(parcelExpressions,
                                                     pattern = c("\\(" = "\\\\(",
                                                                 "\\)" = "\\\\)"))

  # parcel names
  parcelNames <- stringr::str_remove_all(parcelExpressions, "\\s+|\\(|\\)") |>
    stringr::str_replace_all(type, paste0(stringr::str_to_upper(type), "_")) |>
    stringr::str_replace_all(",", "_")
  # split each parcel into its component
  parcelComponents <-
    stringr::str_remove_all(parcelExpressions, paste0(type, "\\(|\\)|\\+|\\s+"))
  variablesInParcels <- lapplyNamed(parcelComponents,
                                    FUN = stringr::str_split_1,
                                    pattern = ",",
                                    names = parcelNames)

  list(parcelExpressions = regexParcelExpressions,
       parcelNames = parcelNames,
       variablesInParcels = variablesInParcels)
}


computeRowMeans <- function(varNames, data) {
  if (length(varNames) <= 1) {
    stop2("Attempted to make a parcel out of a single variable!\n")
  }
  if (!is.data.frame(data)) {
    stop2("Data for parceling should be a dataframe\n")
  }

  rowMeans(data[varNames])
}



computeRowSums <- function(varNames, data) {
  if (length(varNames) <= 1) {
    stop2("Attempted to make a parcel out of a single variable!\n")
  }
  if (!is.data.frame(data)) {
    stop2("Data for parceling should be a dataframe\n")
  }

  rowSums(data[varNames])
}
