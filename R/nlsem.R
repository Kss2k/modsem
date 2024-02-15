

modsem.LMS <- function(modelSpec,
                       data,
                       qml = FALSE,
                       standardizeData = TRUE,
                       verbose = TRUE, 
                       startingValues = NULL,
                       ...) {
  info <- modelSpec$nlsem
  nlsemModel <- nlsem::lav2nlsem(info$modelSyntax)
  # Sorted colnames for xi variables
  sortedIndsXi <- unlist(info$indsXi)
  # Sorted colnames for etas
  sortedIndsEta <- unlist(info$indsEta)

  # Sort/subset data in sorted format (needed for nlsem::em() to work)
  nlsemData <- data[c(sortedIndsXi, sortedIndsEta)]

  if (standardizeData) {
    nlsemData <- lapplyDf(nlsemData,
                          FUN = scale)
  } else if (!standardizeData) {
    warning("modsem recommends standardizing data before running an LMS ",
            "in some cases nlsem::em() wont be able to estimate logLikelyhood \n")
  }
  if (is.null(startingValues)) {
    startingValues <- runif(nlsem::count_free_parameters(nlsemModel))
  }
  nlsemEstimate <- nlsem::em(nlsemModel,
                             nlsemData,
                             start = startingValues,
                             qml = qml,
                             verbose = TRUE,
                             ...)
  nlsemEstimate
}
