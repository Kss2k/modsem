checkModel <- function(model, covModel = NULL, method = "lms") {
  checkCovModelVariables(covModel = covModel, modelXis = model$info$xis)

  checkZeroVariances(model = model, method = method)

  checkNodesLms(parTableMain = model$parTable,
                parTableCov  = covModel$parTable,
                nodes = model$quad$m, method = method)

  checkOVsInStructuralModel(parTableMain = model$parTable,
                            parTableCov  = covModel$parTable)

  checkOverlappingIndicators(allIndsXis = model$info$allIndsXis,
                             allIndsEtas = model$info$allIndsEtas)
}


checkCovModelVariables <- function(covModel, modelXis) {
  if (is.null(covModel$info)) return(NULL) # nothing to check
  covModelEtas <- covModel$info$etas
  covModelXis  <- covModel$info$xis

  stopif(!all(c(covModelXis, covModelEtas) %in% modelXis),
         "All latent variables in the cov-model must be an ",
         "exogenous variable in the main model")
  stopif(!all(modelXis %in% c(covModelXis, covModelEtas)),
         "All exogenous variables in main model must be ",
         "part of the cov-model")
}


checkZeroVariances <- function(model, method = "lms") {
  if (method != "lms") return(NULL)
  
  nonLinearXis <- model$info$nonLinearXis
  inds <- model$info$indsXis[nonLinearXis]

  thetaDelta <- model$matrices$thetaDelta 

  message <- paste(
      "The variance of a moderating variable of integration",
      "has an indicator with zero residual variance!",
      "\nThis will likely not work with the LMS approach, see:",
      "\n   `vignette('observed_lms_qml', 'modsem')` for more information.",
      "\n\nThe following indicators have zero residual variance:"
  )

  m1 <- \(i) sprintf("  -> %s", i)

  width <- options()$width
  width <- if (is.null(width) || width == Inf) 30 else (width - 11) / 2

  error <- FALSE
  for (lv in nonLinearXis) for (ind in inds[[lv]]) {
    est <- thetaDelta[ind, ind]

    if (!is.na(est) && est == 0) {
      error <- TRUE
      message <- paste(message, m1(ind), sep = "\n")
    }
  }

  if (error) stop2(message)
}


checkNodesLms <- function(parTableMain,
                          parTableCov,
                          nodes,
                          method = "lms",
                          minNodesXiXi = 16,
                          minNodesXiEta = 32,
                          minNodesEtaEta = 48) {
  if (method != "lms") return(NULL)

  parTable <- rbind(parTableMain, parTableCov)

  etas     <- getEtas(parTable, isLV = TRUE)
  xis      <- getXis(parTable, etas = etas, isLV = TRUE)
  varsInts <- getVarsInts(getIntTermRows(parTable))

  nodesXiXi_ok   <- TRUE
  nodesXiEta_ok  <- TRUE
  nodesEtaEta_ok <- TRUE

  lapply(varsInts, FUN = function(x) {
    if      (all(x %in% xis))  nodesXiXi_ok   <<- nodes >= minNodesXiXi
    else if (all(x %in% etas)) nodesEtaEta_ok <<- nodes >= minNodesEtaEta
    else if (any(x %in% etas)) nodesXiEta_ok  <<- nodes >= minNodesXiEta
    else warning2("Unable to classify latent variables in interaction terms")
  })

  warnif(!nodesXiXi_ok, "It is recommended that you have at least ",
         minNodesXiXi,  " nodes for interaction effects between ",
         "exogenous variables in the lms approach 'nodes = ", nodes, "'")
  warnif(!nodesXiEta_ok, "It is recommended that you have at least ",
         minNodesXiEta, " nodes for interaction effects between exogenous ",
         "and endogenous variables in the lms approach 'nodes = ", nodes, "'")
  warnif(!nodesEtaEta_ok, "It is recommended that you have at least ",
         minNodesEtaEta, " nodes for interaction effects between endogenous ",
         "variables in the lms approach 'nodes = ", nodes, "'")
}


checkOVsInStructuralModel <- function(parTableMain, parTableCov) {
  parTable <- rbind(parTableMain, parTableCov)
  xisLVs   <- getXis(parTable, isLV = TRUE)
  xisAll   <- getXis(parTable, isLV = FALSE)

  stopif(length(xisAll) != length(xisLVs) || !all(xisLVs %in% xisAll),
         "Observed variables are not allowed in the structural model in LMS/QML directly. ",
         "Please redefine them as latent.\nSee:\n",
         "  vignette(\"observed_lms_qml\", \"modsem\")")
}


checkOverlappingIndicators <- function(allIndsXis, allIndsEtas) {
  stopif(any(allIndsXis %in% allIndsEtas),
         "The same indicator cannot be used for both an exogenous ",
         "and endogenous variable, in the same model: ",
         paste(allIndsXis[allIndsXis %in% allIndsEtas], collapse = ", "))
}


checkParTableDA <- function(parTable) {
  stopif(length(getHigherOrderLVs(parTable)) > 0,
         "Higher-order latent variables are not supported in the lms and qml approaches.")
}


checkVarsIntsDA <- function(varsInts, lVs) {
  for (xz in varsInts) {
    stopif(!all(xz %in% lVs), "Element in product term is not a latent variable: `",
           xz[!xz %in% lVs][[1]], "`!\n",
           "If it is an observed variable, please redefine it as a latent variable.\n",
           "See:\n  vignette(\"observed_lms_qml\", \"modsem\")")
  }
}
