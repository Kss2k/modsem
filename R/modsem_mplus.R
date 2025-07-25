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
  OUTPUT <- "TECH1;\nTECH3;"

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
  TECH1 <- MplusAutomation::get_results(results, element = "tech1")
  pars.tech1 <- getOrderedParameterLabelsMplus(mplusParTable, 
                                               TECH1 = TECH1, 
                                               intTerms = intTerms, 
                                               intTermsMplus = intTermsMplus)
  labels.tech1 <- names(pars.tech1)

  isfree <- coefsTable$pval != 999
  labels <- stringr::str_remove_all(coefsTable$Label[isfree], pattern = " ")
  coefs  <- structure(coefsTable$est[isfree], names = labels)
  coefs  <- coefs[labels.tech1[labels.tech1 %in% labels]]

  vcov <- tryCatch({
    # vcov <- MplusAutomation::get_tech3(results)$paramCov
    # `MplusAutomation::get_tech3` seems very unstable...
    # I often get this error:
    # `Error in get_results(x = results, element = "tech3"): could not find function "get_results"`
    # `MplusAutomation::get_results` seems to work better...
    vcov <- MplusAutomation::get_results(results, element = "tech3")$paramCov
    vcov[upper.tri(vcov)] <- t(vcov)[upper.tri(vcov)]
    dimnames(vcov) <- list(labels.tech1, labels.tech1)
    vcov

  }, error = function(e) {
    warning2("Unable to retrive `tech3` from `Mplus` results\n",
             "Message: ", e, immediate. = FALSE)

    k <- length(pars.tech1)
    vcov <- matrix(NA, nrow = k, ncol = k, 
                   dimnames = list(labels.tech1, labels.tech1))
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
    data     = data, 
    coefs    = modsemVector(coefs),
    vcov     = modsemMatrix(vcov, symmetric = TRUE),
    std      = std,
    info     = list(indicators    = indicators,
                    intTerms      = intTerms,
                    intTermsMplus = intTermsMplus,
                    pars.tech1    = pars.tech1)
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
  coefsTable <- rename(coefsTable, Label = "label", 
                       se = "std.error", pval = "p.value")
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

  maxCharMplus <- maxchar(c(structModel$rhs, structModel$lhs))
  for (i in seq_along(intTerms)) {
    xzMplus  <- substr(intTermsMplus[[i]], start = 1, stop = maxCharMplus)
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


getOrderedParameterLabelsMplus <- function(parTable, TECH1, intTerms, intTermsMplus) {
  spec   <- TECH1$parameterSpecification
  nu     <- spec$nu
  alpha  <- spec$alpha
  lambda <- spec$lambda
  beta   <- spec$beta
  psi    <- spec$psi
  theta  <- spec$theta

  setLabel <- function(out, label, id) {
    if (!length(label)) {
      warning2("Unable to find label for parameter ", id, "!",
               immediate. = FALSE)
      out[as.character(id)] <- id

    } else out[label] <- id

    out
  }

  getReg <- function(M, row.lhs = TRUE, op = "~") {
    out <- c()

    rows <- rownames(M)
    cols <- colnames(M)
    for (i in seq_len(NROW(M))) for (j in seq_len(NCOL(M))) {
      id <- as.integer(M[i, j])
      
      if (id <= 0L)
        next

      if (row.lhs) {
        lhs <- rows[i]
        rhs <- cols[j]
      } else {
        lhs <- cols[j]
        rhs <- rows[i]
      }

      label <- parTable[parTable$op == op &
                        parTable$lhs == lhs &
                        parTable$rhs == rhs, "label"]

      out <- setLabel(out = out, label = label, id = id)
    }

    out
  }
  
  getIntercept <- function(T) {
    vars <- colnames(T)
    T <- c(T)
    out <- c()

    for (i in seq_along(T)) {
      id <- as.integer(T[i])
      
      if (id <= 0L)
        next
      
      lhs <- vars[i]

      label <- parTable[parTable$op == "~1" &
                        parTable$lhs == lhs, "label"]

      out <- setLabel(out = out, label = label, id = id)
    }

    out
  }


  getCovariance <- function(M, op = "~~") {
    out <- c()

    rows <- rownames(M)
    cols <- colnames(M)

    for (i in seq_len(NROW(M))) for (j in seq_len(i)) {
      id <- as.integer(M[i, j])
      
      if (id <= 0L)
        next

      lhs <- cols[i]
      rhs <- rows[j]

      label <- parTable[parTable$op == op &
                        ((parTable$lhs == lhs & parTable$rhs == rhs) |
                         (parTable$lhs == rhs & parTable$rhs == lhs)), "label"]

      out <- setLabel(out = out, label = label, id = id)
    }
    
    out
  }

  for (i in seq_along(intTerms)) {
    xzMplus  <- substr(intTermsMplus[[i]], start = 1L, stop = 8L)
    xzModsem <- intTerms[[i]]

    parTable[parTable$lhs == xzModsem, "lhs"] <-  xzMplus
    parTable[parTable$rhs == xzModsem, "rhs"] <-  xzMplus
  }

  out <- c(
    getIntercept(nu),
    getIntercept(alpha),
    getReg(lambda, row.lhs = FALSE, op = "=~"),
    getReg(beta, row.lhs = TRUE, op = "~"),
    getCovariance(psi),
    getCovariance(theta)
  )

  out <- out[!duplicated(out)] # unique() removes labels
  
  sort(out)
}

