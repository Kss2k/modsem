checkModel <- function(model, covModel = NULL, method = "lms") {
  checkCovModelVariables(covModel = covModel, modelXis = model$info$xis)

  checkNodesLms(parTableMain = model$parTable, 
                parTableCov  = covModel$parTable,
                nodes = model$quad$m, method = method)

  checkOVsInStructuralModel(parTableMain = model$parTable,
                            parTableCov  = covModel$parTable)
  NULL
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
