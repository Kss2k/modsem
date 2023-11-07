

modsem.LMS <- function(modelSpecification,
                       data,
                       qml = FALSE,
                       standardizeData = TRUE,
                       verbose = TRUE) {
  info <- modelSpecification$nlsem
  nlsemModel <- nlsem::lav2nlsem(info$modelSyntax)

  # Sorted colnames for xi variables
  sortedIndicatorsXi <- unlist(info$indicatorsXi)
  # Sorted colnames for etas
  sortedIndicatorsEta <- unlist(info$indicatorsEta)
  # Sort/subset data in sorted format (needed for nlsem::em() to work)
  nlsemData <- data[c(sortedIndicatorsXi, sortedIndicatorsEta)]

  if (standardizeData == TRUE) {
    nlsemData <- lapplyDf(nlsemData,
                          FUN = scale)
  } else if (standardizeData != TRUE) {
    warning("modsem recommends standardizing data before running an LMS ",
            "in some cases nlsem::em() wont be able to estimate logLikelyhood \n")
  }

  startingValues <- runif(nlsem::count_free_parameters(nlsemModel))
  nlsemEstimate <- nlsem::em(nlsemModel,
                             nlsemData,
                             start = startingValues,
                             qml = qml,
                             verbose = TRUE)
  nlsemEstimate
}
