

modsem.LMS <- function(modelSpec,
                       data,
                       qml = TRUE,
                       centerData = TRUE,
                       verbose = TRUE, 
                       start = NULL,
                       optimize = TRUE,
                       m = 16, 
                       ...) {
  browser()
  info <- modelSpec$nlsem
  nlsemModel <- nlsem::lav2nlsem(info$modelSyntax)
  parameters <- as.data.frame(nlsemModel)
  parameters$class1[grepl("^nu", parameters$label)] <- NA
  parameters$class1[grepl("^tau", parameters$label)] <- 0
  # Sorted colnames for xi variables
  sortedIndsXi <- unlist(info$indsXi)
  # Sorted colnames for etas
  sortedIndsEta <- unlist(info$indsEta)
  if (centerData || optimize) {
    if (!centerData) data <- lapplyDf(data, function(x) x - mean(x))
    parameters$class1[grepl("^nu\\.|^tau", parameters$label)] <- 0
  }

  nlsemModel <- nlsem::create_sem(parameters)
  # Sort/subset data in sorted format (needed for nlsem::em() to work)
  nlsemData <- data[c(sortedIndsXi, sortedIndsEta)]

  if (optimize) {
    modsemEst <- modsem(info$modelSyntax, data, "dblcent")$coefParTable
    startDf <- parameters[is.na(parameters$class1), ]
    startDf$class1[grepl("^Lambda.x", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% info$xiNames &
                    modsemEst$op == "=~" & 
                    !is.na(modsemEst$z)]
    startDf$class1[grepl("^Lambda.y", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% info$etaNames &
                    modsemEst$op == "=~" & 
                    !is.na(modsemEst$z)]
    startDf$class1[grepl("^Gamma", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% info$etaNames &
                    modsemEst$op == "~" &
                    modsemEst$rhs %in% info$xiNames]
    startDf$class1[grepl("^Omega", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% info$etaNames &
                    modsemEst$op == "~" &
                    modsemEst$rhs %in% info$prodTerms]
    startDf$class1[grepl("Theta.d", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% sortedIndsXi &
                    modsemEst$op == "~~" &
                    modsemEst$rhs %in% sortedIndsXi &
                    modsemEst$lhs == modsemEst$rhs]
    startDf$class1[grepl("Theta.e", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% sortedIndsEta &
                    modsemEst$op == "~~" &
                    modsemEst$rhs %in% sortedIndsEta &
                    modsemEst$lhs == modsemEst$rhs]
    startDf$class1[grepl("Psi", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% info$etaNames &
                    modsemEst$op == "~~" &
                    modsemEst$rhs %in% info$etaNames]
    startDf$class1[grepl("Phi", startDf$label)] <- 
      c(modsemEst$est[modsemEst$lhs %in% info$xiNames &
                    modsemEst$op == "~~" &
                    modsemEst$rhs %in% info$xiNames &
                    modsemEst$rhs == modsemEst$lhs], 
        modsemEst$est[modsemEst$lhs %in% info$xiNames &
                    modsemEst$op == "~~" &
                    modsemEst$rhs %in% info$xiNames &
                    modsemEst$rhs != modsemEst$lhs]
      )
    startDf$class1[grepl("alpha", startDf$label)] <- 
      modsemEst$est[modsemEst$lhs %in% info$xiNames &
                    modsemEst$op == "~~" &
                    modsemEst$rhs %in% info$xiNames &
                    modsemEst$rhs != modsemEst$lhs]

    start <- startDf$class1
    start[is.na(start)] <- runif(sum(is.na(start)))
  }
  if (is.null(start)) {
    start <- runif(nlsem::count_free_parameters(nlsemModel))
  }
  nlsemEstimate <- nlsem::em(nlsemModel,
                             nlsemData,
                             start = start, 
                             qml = qml,
                             verbose = TRUE,
                             ...)
  nlsemEstimate
}
