#' Estimation latent interactions through \code{Mplus}
#'
#' @param model.syntax lavaan/modsem syntax
#' @param data dataset
#' @param estimator estimator argument passed to \code{Mplus}
#' @param type type argument passed to \code{Mplus}
#' @param algorithm algorithm argument passed to \code{Mplus}
#' @param processors processors argument passed to \code{Mplus}
#' @param integration integration argument passed to \code{Mplus}
#'
#' @param rcs Should latent variable indicators be replaced with reliability-corrected
#'   single item indicators instead? See \code{\link{relcorr_single_item}}.
#'
#' @param rcs.choose Which latent variables should get their indicators replaced with
#'   reliability-corrected single items? It is passed to \code{\link{relcorr_single_item}}
#'   as the \code{choose} argument.
#'   
#' @param rcs.scale.corrected Should reliability-corrected items be scale-corrected? If \code{TRUE}
#'   reliability-corrected single items are corrected for differences in factor loadings between
#'   the items. Default is \code{TRUE}.
#'
#' @param output.std Should \code{STANDARDIZED} be added to \code{OUTPUT}?
#'
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
                         processors = 2,
                         integration = 15,
                         rcs = FALSE,
                         rcs.choose = NULL,
                         rcs.scale.corrected = TRUE,
                         output.std = TRUE,
                         ...) {
  if (rcs) { # use reliability-correct single items?
    corrected <- relcorr_single_item(
      syntax          = model.syntax, 
      data            = data,
      choose          = rcs.choose,
      scale.corrected = rcs.scale.corrected,
      warn.lav        = FALSE
    )

    model.syntax <- corrected$syntax
    data         <- corrected$data
  }

  parTable <- modsemify(model.syntax)

  # Abbreviate variable names
  names <- unique(c(parTable$rhs, parTable$lhs))
  names <- names[!grepl(":", names)]
  abbreviated <- abbreviate(names, minlength = 8L, strict = TRUE)
 
  # Abbreviate names in intTerms
  intTerms <- unique(parTable$rhs[grepl(":", parTable$rhs)])
  abbrevIntTerm <- function(xz)
    paste0(abbreviated[stringr::str_split(xz, pattern = ":")[[1]]], collapse = ":")

  newIntTerms <- vapply(intTerms, FUN.VALUE = character(1L), FUN = abbrevIntTerm)
  names(newIntTerms) <- intTerms
  abbreviated <- c(abbreviated, newIntTerms)

  # Fix names in parTable
  rmask <- parTable$rhs != ""
  lmask <- parTable$lhs != ""
  parTable$rhs[rmask] <- abbreviated[parTable$rhs[rmask]]
  parTable$lhs[lmask] <- abbreviated[parTable$lhs[lmask]]
  
  # Fix names in data
  dmask <- colnames(data) %in% names(abbreviated)
  colnames(data)[dmask] <- abbreviated[colnames(data)[dmask]]

  indicators <- unique(parTable[parTable$op == "=~", "rhs", drop = TRUE])
  intTerms <- unique(getIntTermRows(parTable)$rhs)
  intTermsMplus <- stringr::str_remove_all(intTerms, ":") |>
    stringr::str_to_upper()

  # Set OUTPUT
  OUTPUT <- "TECH3;"

  getIntTermLength <- \(xz) length(stringr::str_split(xz, pattern = ":")[[1]])
  kway <- max(0, vapply(intTerms, FUN.VALUE = integer(1L), FUN = getIntTermLength))
  
  if (output.std && kway <= 2)
    OUTPUT <- paste(OUTPUT, "STANDARDIZED;", sep = "\n")

  # Estimate model
  model <- MplusAutomation::mplusObject(
    TITLE = "Running Model via Mplus",
    usevariables = indicators,
    ANALYSIS = paste0(
      paste(paste("estimator =", estimator),
            paste("type =", type),
            paste("algorithm =", algorithm),
            paste("processors =", processors),
            paste("integration = ", integration),
            sep = ";\n"), ";\n"), # add final ";"
    MODEL = parTableToMplusModel(parTable, ...),
    OUTPUT = OUTPUT,
    rdata = data[indicators],
  )

  results <- MplusAutomation::mplusModeler(model,
                                           modelout = "mplusResults.inp",
                                           run = 1L)
  coefsTable    <- coef(results)
  mplusParTable <- mplusTableToParTable(coefsTable, 
                                        intTerms = intTerms, 
                                        intTermsMplus = intTermsMplus,
                                        indicators = indicators)

  # coef and vcov
  isfree <- coefsTable$pval != 999
  labels <- stringr::str_remove_all(coefsTable$Label[isfree], pattern = " ")
  coefs  <- structure(coefsTable$est[isfree], names = labels)

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
  
  std <- tryCatch({
    MplusAutomation::get_results(results, element = "standardized")
  }, error = function(e) {
    warning2("Unable to retrive `standardized` from `Mplus` results\n",
             "Message: ", e, immediate. = FALSE)
    std <- NULL
  })

  modelSpec <- list(
    parTable = mplusParTable,
    model    = results,
    coefs    = coefs,
    data     = data, 
    coef     = coef,
    vcov     = vcov,
    std      = std,
    info     = list(indicators    = indicators,
                    intTerms      = intTerms,
                    intTermsMplus = intTermsMplus)
  )

  structure(modelSpec,
            class = "modsem_mplus",
            method = "Mplus")
}


xwith <- function(elems) {
  if (length(elems) <= 1) return(NULL)
  lhs <- paste0(elems[[1]], elems[[2]])
  rhs <- paste(elems[[1]], "XWITH", elems[[2]])

  elems2 <- c(lhs, elems[-c(1, 2)])
  
  unique(rbind(
   createParTableRow(c(lhs, rhs), op = ":"),
   xwith(elems2)
  ))
}


parTableToMplusModel <- function(parTable) {
  # INTERACTIONEXPRESSIOns
  interactions <- parTable[grepl(":", parTable$rhs), "rhs"]
  elemsInInts <- stringr::str_split(interactions, ":")
  newRows <- lapply(elemsInInts, FUN = xwith) |>
    purrr::list_rbind()
  parTable <- rbind(parTable, newRows)

  # Identify variances
  isVar <- parTable$op == "~~" & parTable$lhs == parTable$rhs

  parTable$op <- replaceLavOpWithMplus(parTable$op)
  parTable[isVar, c("op", "lhs")] <- ""

  out <- ""
  for (i in seq_len(nrow(parTable))) {
    mod   <- parTable[["mod"]][i]
    empty <- mod == ""
    islab <- !canBeNumeric(mod)

    if (!empty && islab)       modifier <- paste0(" (", mod,")")
    else if (!empty && !islab) modifier <- paste0("@", mod)
    else                       modifier <- ""

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


mplusTableToParTable <- function(coefsTable, 
                                 indicators, 
                                 intTerms, 
                                 intTermsMplus) {
  coefsTable <- rename(coefsTable, Label = "label", se = "std.error",
                       pval = "p.value")
  coefsTable$label <- stringr::str_remove_all(coefsTable$label, pattern = " ")

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

  nchars <- maxchar(c(structModel$rhs, structModel$lhs))
  for (i in seq_along(intTerms)) {
    xzMplus  <- substr(intTermsMplus[[i]], start = 1, stop = nchars)
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

  modsemParTable(mplusParTable)
}
