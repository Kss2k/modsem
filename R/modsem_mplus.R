#' Estimation latent interactions through mplus 
#'
#' @param model_syntax lavaan/modsem syntax 
#' @param data dataset
#' @param estimator estimator argument passed to mplus
#' @param type type argument passed to mplus
#' @param algorithm algorithm argument passed to mplus
#' @param process process argument passed to mplus
#' @param ... arguments passed to other functions
#'
#' @return modsem_mplus object
#' @export
#'
#' @examples
#' # Theory Of Planned Behavior
#' tpb <- ' 
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   # Covariances
#'   ATT ~~ SN + PBC
#'   PBC ~~ SN 
#'   # Causal Relationsships
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC 
#'   BEH ~ INT:PBC  
#' '
#' 
#' \dontrun{
#' estTpbMplus <- modsem_mplus(tpb, data = TPB)
#' summary(estTpbLMS)
#' }
#' 
modsem_mplus <- function(model_syntax, 
                         data, 
                         estimator = "ml", 
                         type = "random", 
                         algorithm = "integration", 
                         process = "8", 
                         ...) {
  parTable <- modsemify(model_syntax)
  indicators <- parTable[parTable$op == "=~", "rhs", drop = TRUE] |>
    unique()
  model <- MplusAutomation::mplusObject(
    TITLE = "Running Model via Mplus",
    usevariables = indicators,
    ANALYSIS =
      paste(paste("estimator =", estimator), 
            paste("type =", type), 
            paste("algorithm =", algorithm), 
            paste("process =", process, ";\n"),
            sep = ";\n"), 
    MODEL = parTableToMplusModel(parTable, ...),
    rdata = data[indicators],
  )
  results <- MplusAutomation::mplusModeler(model, 
                                           modelout = "mplusResults.inp", 
                                           run = 1L)
  coefs <- MplusAutomation::extract.mplus.model(results)
  coefsTable <- data.frame(lhsOpRhs = coefs@coef.names,
                           est = coefs@coef,
                           se = coefs@se,
                           pvalue = coefs@pvalues)
  # Measurement Model
  indicatorsCaps <- stringr::str_to_upper(indicators)
  patternMeas <-
    paste0("(", stringr::str_c(indicatorsCaps, collapse = "|"), ")") |>
    paste0("<-(?!>|Intercept)")
  measCoefNames <- grepl(patternMeas, coefsTable$lhsOpRhs, perl = TRUE)
  # Mplus has lhs/rhs in reversed order for the measurement model,
    # compared to lavaan,
  measRhs <- stringr::str_split_i(coefsTable$lhsOpRhs[measCoefNames],
                                  "<-", i = 1)
  measLhs <- stringr::str_split_i(coefsTable$lhsOpRhs[measCoefNames],
                                  "<-", i = 2)
  measModel <- data.frame(lhs = measLhs, op = "=~", rhs = measRhs) |>
    cbind(coefsTable[measCoefNames, c("est", "se", "pvalue")])

  # Structural Model
  measrRemoved <- coefsTable[!measCoefNames, ]
  patternStruct <- "<-(?!>|Intercept)"
  structCoefNames <- grepl(patternStruct, measrRemoved$lhsOpRhs, perl = TRUE)

  structLhs <- stringr::str_split_i(measrRemoved$lhsOpRhs[structCoefNames],
                                  "<-", i = 1)
  structRhs <- stringr::str_split_i(measrRemoved$lhsOpRhs[structCoefNames],
                                  "<-", i = 2)
  structModel <- data.frame(lhs = structLhs, op = "~", rhs = structRhs) |>
    cbind(measrRemoved[structCoefNames, c("est", "se", "pvalue")])

  # Variances and Covariances
  structMeasrRemoved <- measrRemoved[!structCoefNames, ]
  patternCovVar <- "<->"
  covVarCoefNames <- grepl(patternCovVar, structMeasrRemoved$lhsOpRhs, perl = TRUE)
  covVarLhs <- stringr::str_split_i(structMeasrRemoved$lhsOpRhs[covVarCoefNames],
                                  "<->", i = 1)
  covVarRhs <- stringr::str_split_i(structMeasrRemoved$lhsOpRhs[covVarCoefNames],
                                  "<->", i = 2)
  covVarModel <- data.frame(lhs = covVarLhs, op = "~~", rhs = covVarRhs) |>
    cbind(structMeasrRemoved[covVarCoefNames, c("est", "se", "pvalue")])

  # Intercepts
  covStructMeasrRemoved <- structMeasrRemoved[!covVarCoefNames, ]
  patternIntercept <- "<-Intercept"
  interceptNames <- grepl(patternIntercept, covStructMeasrRemoved$lhsOpRhs, perl = TRUE)
  interceptLhs <- stringr::str_split_i(covStructMeasrRemoved$lhsOpRhs[interceptNames],
                                  "<-", i = 1)
  interceptModel <- data.frame(lhs = interceptLhs, op = "~", rhs = 1) |>
    cbind(covStructMeasrRemoved[interceptNames, c("est", "se", "pvalue")])

  mplusParTable <- rbind(measModel, structModel, covVarModel, interceptModel)
  mplusParTable [c("lhs", "rhs")] <- lapplyDf(mplusParTable[c("lhs", "rhs")],
                                   function(x)
                                    stringr::str_remove_all(x, " "))
  mplusParTable$ci.lower <- mplusParTable$est - 1.96*mplusParTable$se
  mplusParTable$ci.upper <- mplusParTable$est + 1.96*mplusParTable$se
  mplusParTable$pvalue[mplusParTable$pvalue == 999] <- NA
  mplusParTable$label <- NA
  mplusParTable$z <- NA

  modelSpec <- list(parTable = parTable,
                    coefParTable = mplusParTable,
                    model = results,
                    coefs = coefs,
                    data = data)
  structure(modelSpec,
            class = "modsem_mplus",
            method = "Mplus")
}


parTableToMplusModel <- function(parTable, ignoreLabels = TRUE) {
  # INTERACTIONEXPRESSIOns
  interactions <- parTable[grepl(":", parTable$rhs), "rhs"]
  elemsInInts <- stringr::str_split(interactions, ":")
  newRows <- lapply(elemsInInts,
                    function(x) {
                      if (length(x) != 2) {
                        stop2("Number of variables in interaction must be two")
                      }
                      lhs <- paste0(x[[1]], x[[2]])
                      rhs <- paste(x[[1]], "XWITH", x[[2]])
                      createParTableRow(c(lhs, rhs), op = ":")
                      }) |>
    purrr::list_rbind()
  parTable <- rbind(parTable, newRows)

  parTable$op <- replaceLavOpWithMplus(parTable$op)
  out <- ""
  if (ignoreLabels) parTable[["mod"]] <- ""
  for (i in 1:nrow(parTable)) {
    if (parTable[["mod"]][i] != "") {
      warning2("Using labels in Mplus, was this intended?")
      modifier <- paste0("* (", parTable[["mod"]][[i]],")")

    } else {
      modifier <- ""
    }
    line <- paste0(parTable[["lhs"]][[i]], " ",
                   parTable[["op"]][[i]], " ",
                   parTable[["rhs"]][[i]],
                   modifier, ";" ,"\n")
    out <- paste0(out, line)
  }
  stringr::str_remove_all(out, ":")
}


replaceLavOpWithMplus <- function(op) {
  vapply(op,
         FUN = switchLavOpToMplus,
         FUN.VALUE = character(1L))
}


switchLavOpToMplus <- function(op) {
  switch(op,
         "=~" = "BY",
         "~" = "ON",
         "~~" = "WITH",
         ":" = "|",
         stop2("Operator not supported for use in Mplus: ", op, "\n"))
}
