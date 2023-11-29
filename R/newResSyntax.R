newGetParTableResCov <- function(relDf) {
  if (ncol(relDf) <= 1) {
    return(NULL)
  }

  prodNames <- colnames(relDf)

  sharedMatrix <- matrix("", nrow = length(prodNames), ncol = length(prodNames),
                         dimnames = list(prodNames, prodNames))

  # Now we want to specify the covariance based on shared inds
  for (i in prodNames) {
    for (j in prodNames) {
      sharedIndicators <- relDf[[i]][relDf[[i]] %in% relDf[[j]]]
      sharedMatrix[i, j] <- stringr::str_c(sharedIndicators,
                                                  collapse = "_")
    }
  }
  labelMatrix <- sharedMatrix
  labelMatrix <- ifelse(labelMatrix == "", "0", paste0("share_", labelMatrix))
  labelMatrix[upper.tri(labelMatrix, diag = TRUE)] <- ""

  print(labelMatrix)
  uniqueCombos <- getUniqueCombinations(prodNames)
  uniqueCombos[["labels"]] <- vector("character", length = nrow(uniqueCombos))
  for (i in 1:nrow(uniqueCombos)) {
    uniqueCombos[["labels"]][[i]] <- labelMatrix[uniqueCombos[i, "V2"],
                                     uniqueCombos[i, "V1"]]
  }

  apply(uniqueCombos, MARGIN = 1,
        FUN = function(x)
          createParTableRow(x[c("V1", "V2")], op = "~~", mod = x[["labels"]])) |>
    purrr::list_rbind()
}
