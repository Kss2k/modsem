#' Estimation latent interactions through \code{Mplus}
#'
#' @param model.syntax lavaan/modsem syntax
#' @param data dataset
#' @param estimator estimator argument passed to \code{Mplus}
#' @param type type argument passed to \code{Mplus}
#' @param algorithm algorithm argument passed to \code{Mplus}
#' @param process process argument passed to \code{Mplus}
#' @param integration integration argument passed to \code{Mplus}
#' @param ... arguments passed to other functions
#'
#' @return modsem_mplus object
#' @export
#'
#' @examples
#' # Theory Of Planned Behavior
#' m1 <- '
#' # Outer Model
#'   X =~ x1 + x2
#'   Z =~ z1 + z2
#'   Y =~ y1 + y2 
#' 
#' # Inner model
#'   Y ~ X + Z + X:Z
#' '
#'
#' \dontrun{
#' # Check if Mplus is installed
#' run <- tryCatch({MplusAutomation::detectMplus(); TRUE},
#'                 error = \(e) FALSE)
#' 
#' if (run) {
#'   est_mplus <- modsem_mplus(m1, data = oneInt)
#'   summary(est_mplus)
#' }
#' }
#'
modsem_mplus <- function(model.syntax,
                         data,
                         estimator = "ml",
                         type = "random",
                         algorithm = "integration",
                         process = 8,
                         integration = 15,
                         ...) {
  parTable <- modsemify(model.syntax)
  indicators <- unique(parTable[parTable$op == "=~", "rhs", drop = TRUE])
  intTerms <- unique(getIntTermRows(parTable)$rhs)
  intTermsMplus <- stringr::str_remove_all(intTerms, ":") |>
    stringr::str_to_upper()

  model <- MplusAutomation::mplusObject(
    TITLE = "Running Model via Mplus",
    usevariables = indicators,
    ANALYSIS = paste0(
      paste(paste("estimator =", estimator),
            paste("type =", type),
            paste("algorithm =", algorithm),
            paste("process =", process),
            paste("integration = ", integration),
            sep = ";\n"), ";\n"), # add final ";"
    MODEL = parTableToMplusModel(parTable, ...),
    OUTPUT = "TECH3;",
    rdata = data[indicators],
  )
  results <- MplusAutomation::mplusModeler(model,
                                           modelout = "mplusResults.inp",
                                           run = 1L)
  coefs <- MplusAutomation::extract.mplus.model(results)
  coefsTable <- data.frame(label = stringr::str_remove_all(coefs@coef.names, pattern = " "),
                           est = coefs@coef,
                           std.error = coefs@se,
                           p.value = coefs@pvalues)

  # Measurement Model
  indicatorsCaps <- stringr::str_to_upper(indicators)
  patternMeas <-
    paste0("(", stringr::str_c(indicatorsCaps, collapse = "|"), ")") |>
    paste0("<-(?!>|Intercept)")
  measCoefNames <- grepl(patternMeas, coefsTable$label, perl = TRUE)

  # Mplus has lhs/rhs in reversed order for the measurement model,
    # compared to lavaan,
  measRhs <- stringr::str_split_i(coefsTable$label[measCoefNames],
                                  "<-", i = 1)
  measLhs <- stringr::str_split_i(coefsTable$label[measCoefNames],
                                  "<-", i = 2)
  measModel <- data.frame(lhs = measLhs, op = "=~", rhs = measRhs) |>
    cbind(coefsTable[measCoefNames, ])

  # Structural Model
  measrRemoved <- coefsTable[!measCoefNames, ]
  patternStruct <- "<-(?!>|Intercept)"
  structCoefNames <- grepl(patternStruct, measrRemoved$label, perl = TRUE)

  structLhs <- stringr::str_split_i(measrRemoved$label[structCoefNames],
                                  "<-", i = 1)
  structRhs <- stringr::str_split_i(measrRemoved$label[structCoefNames],
                                  "<-", i = 2)
  structModel <- data.frame(lhs = structLhs, op = "~", rhs = structRhs) |>
    cbind(measrRemoved[structCoefNames, ])

  for (i in seq_along(intTerms)) {
    xzMplus <- intTermsMplus[[i]]
    xzModsem <- intTerms[[i]]
    structModel[structModel$rhs == xzMplus, "rhs"] <- xzModsem
    structModel[structModel$lhs == xzMplus, "lhs"] <- xzModsem
  }

  # Variances and Covariances
  structMeasrRemoved <- measrRemoved[!structCoefNames, ]
  patternCovVar <- "<->"
  covVarCoefNames <- grepl(patternCovVar, structMeasrRemoved$label, perl = TRUE)
  covVarLhs <- stringr::str_split_i(structMeasrRemoved$label[covVarCoefNames],
                                  "<->", i = 1)
  covVarRhs <- stringr::str_split_i(structMeasrRemoved$label[covVarCoefNames],
                                  "<->", i = 2)
  covVarModel <- data.frame(lhs = covVarLhs, op = "~~", rhs = covVarRhs) |>
    cbind(structMeasrRemoved[covVarCoefNames, ])

  # Intercepts
  covStructMeasrRemoved <- structMeasrRemoved[!covVarCoefNames, ]
  patternIntercept <- "<-Intercept"
  interceptNames <- grepl(patternIntercept, covStructMeasrRemoved$label, perl = TRUE)
  interceptLhs <- stringr::str_split_i(covStructMeasrRemoved$label[interceptNames],
                                  "<-", i = 1)
  interceptModel <- data.frame(lhs = interceptLhs, op = "~1", rhs = "") |>
    cbind(covStructMeasrRemoved[interceptNames, ])

  mplusParTable <- rbind(measModel, structModel, covVarModel, interceptModel)
  mplusParTable [c("lhs", "rhs")] <-
    lapplyDf(mplusParTable[c("lhs", "rhs")], function(x)
             stringr::str_remove_all(x, " "))

  mplusParTable$ci.lower <- mplusParTable$est - CI_WIDTH * mplusParTable$std.error
  mplusParTable$ci.upper <- mplusParTable$est + CI_WIDTH * mplusParTable$std.error
  mplusParTable$p.value[mplusParTable$p.value == 999] <- NA
  mplusParTable$z.value <- mplusParTable$est / mplusParTable$std.error
  mplusParTable$z.value[is.infinite(mplusParTable$z.value)] <- NA

  # coef and vcov
  isfree <- coefsTable$p.value != 999
  labels <- coefsTable$label[isfree]
  coef <- structure(coefsTable$est[isfree], names = labels)

  vcov <- tryCatch({
    # vcov <- MplusAutomation::get_tech3(results)$paramCov
    # `MplusAutomation::get_tech3` seems very unstable...
    # I often get this error:
    # `Error in get_results(x = results, element = "tech3"): could not find function "get_results"`
    # `MplusAutomation::get_results` seems to work better...
    vcov <- MplusAutomation::get_results(results, element = "tech3")$paramCov
    vcov[upper.tri(vcov)] <- t(vcov)[upper.tri(vcov)]
    dimnames(vcov) <- list(labels, labels)
    vcov
  }, error = function(e) {
    warning2("Unable to retrive `tech3` from `Mplus` results\n",
             "Message: ", e, immediate. = FALSE)
    vcov <- matrix(NA, nrow = sum(isfree), ncol = sum(isfree), 
                   dimnames = list(labels, labels))
  })

  modelSpec <- list(parTable = mplusParTable,
                    model = results,
                    coefs = coefs,
                    data = data, 
                    coef = coef,
                    vcov = vcov)

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


#' @export
vcov.modsem_mplus <- function(object, ...) {
  object$vcov
}


#' @export
coef.modsem_mplus <- function(object, ...) {
  object$coef
}
