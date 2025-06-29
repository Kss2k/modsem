#' Reliability‑Corrected Single‑Item SEM
#'
#' Replace (some of) the first‑order latent variables in a lavaan measurement
#' model by **single composite indicators whose error variances are fixed from
#' Cronbach's \eqn{\alpha}. The function returns a modified lavaan model
#' syntax together with an augmented data set that contains the newly created
#' composite variables, so that you can fit the full SEM in a single step.
#'
#' The resulting object can be fed directly into \code{modsem} or \code{lavaan::sem}
#' by supplying \code{syntax = $syntax} and \code{data = $data}.
#'
#' @param syntax A character string containing lavaan model syntax.  Must at
#'   least include the measurement relations (\code{=~}).
#' @param data A \code{data.frame}, \code{tibble} or object coercible to a data frame that
#'   holds the raw observed indicators.
#' @param choose *Optional.* Character vector with the names of the latent
#'   variables you wish to replace by single indicators.  Defaults to **all**
#'   first‑order latent variables in `syntax`.
#'
#' @return An object of S3 class `modsem_relcorr` (a named list) with elements:
#' \describe{
#'   \item{`syntax`}{Modified lavaan syntax string.}
#'   \item{`data`}{Data frame with additional composite indicator columns.}
#'   \item{`parTable`}{Parameter table corresponding to `syntax`.}
#'   \item{`reliability`}{Named numeric vector of reliabilities (one per latent
#'     variable).}
#'   \item{`AVE`}{Named numeric vector with Average Variance Extracted values.}
#'   \item{`lVs`}{Character vector of latent variables that were corrected.}
#'   \item{`single.items`}{Character vector with names for the generated reliability corrected items}
#' }
#'
#' @examples
#' \dontrun{
#' tpb_uk <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'  ATT =~ att3 + att2 + att1 + att4
#'  SN =~ sn4 + sn2 + sn3 + sn1
#'  PBC =~ pbc2 + pbc1 + pbc3 + pbc4
#'  INT =~ int2 + int1 + int3 + int4
#'  BEH =~ beh3 + beh2 + beh1 + beh4
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'  INT ~ ATT + SN + PBC
#'  BEH ~ INT + PBC
#'  BEH ~ INT:PBC
#' "
#' 
#' corrected <- relcorr_single_item(syntax = tpb_uk, data = TPB_UK)
#' print(corrected)
#'
#' syntax <- corrected$syntax
#' data   <- corrected$data
#' 
#' est_dca <- modsem(syntax, data = data, method = "dblcent")
#' est_lms <- modsem(syntax, data = data, method="lms", nodes=32)
#' }
#'
#' @export
relcorr_single_item <- function(syntax, data, choose = NULL) {
  data <- as.data.frame(data)

  parTable      <- modsemify(syntax)
  parTableOuter <- parTable[parTable$op == "=~", ]
  parTableInner <- parTable[parTable$op != "=~", ]

  higherOrderLVs <- getHigherOrderLVs(parTable)
  lVs            <- getLVs(parTable)
  
  choose <- if (is.null(choose)) lVs else choose
  ignore <- setdiff(lVs, choose)

  if (length(higherOrderLVs)) {
    warning2("Higher order latent variables will be ignored when ",
             "creating reliablity-corrected single items!")
    ignore <- c(ignore, higherOrderLVs)
  }

  stopif(any(!choose %in% lVs), "Could not find latent variables:\n", choose)

  ignoreRows <- parTableOuter$lhs %in% ignore

  parTableInner <- rbind(parTableInner, parTableOuter[ignoreRows, ])
  parTableOuter <- parTableOuter[!ignoreRows, ]

  stopif(!NROW(parTableOuter), "Cannot ignore all latent variables when creating ",
         "reliability corrected items!")
 
  warnif(any(parTableOuter$mod != ""), 
         "Labels and modifiers in the measurement model used to create reliability ",
         "corrected single items will be ignored!\nThis may cause the model constraints ",
         "in the structural model to break!")


  lVs <- setdiff(lVs, ignore)
  inds <- getInds(parTableOuter)
  indsLVs <- getIndsLVs(parTableOuter, lVs = lVs)
  R <- stats::cov(data[inds])
  k <- length(lVs)
  
  cfaSyntax <- parTableToSyntax(parTableOuter)
  cfa       <- lavaan::cfa(cfaSyntax, data = data)

  stdSolution <- lavaan::lavInspect(cfa, "std")
  lambda <- stdSolution$lambda 
  theta  <- stdSolution$theta
 
  avefunc <- \(lV) calcAVE(lV, theta.std = theta, lambda.std = lambda)
  relfunc <- \(lV) calcChronbach(R = R, x = indsLVs[[lV]])

  rel.std <- vapply(lVs, FUN.VALUE = numeric(1L), relfunc)
  AVE.std <- vapply(lVs, FUN.VALUE = numeric(1L), avefunc)
  res.std <- 1 - rel.std
  singleInds <- getSingleIndNames(lVs)

  rchar <- c(letters, 1:9)
  while(any(singleInds %in% colnames(data))) {
    # if matching names in data, add random characters, until there is no
    # match. Should in practive never be triggered
    singleInds <- paste0(singleInds, sample(rchar, length(singleInds)))
  }

  varNewInds <- structure(numeric(length(lVs)), names=singleInds)
  newData <- data

  for (lV in lVs) {
    indsLV     <- indsLVs[[lV]]
    newIndName <- singleInds[[lV]]
    X          <- as.matrix(data[indsLV])
    newInd     <- rowMeans(data[indsLV], na.rm = TRUE)

    newData[[newIndName]]    <- newInd
    varNewInds[[newIndName]] <- stats::var(newInd, na.rm = TRUE)
  }
  
  rel <- rel.std * varNewInds
  res <- res.std * varNewInds

  lhs <- unname(c(lVs, singleInds))
  rhs <- unname(c(singleInds, singleInds))
  op  <- c(rep("=~", k), rep("~~", k))
  mod <- unname(c(rep("1", k), res))

  newParTable <- rbind(
    parTableInner, data.frame(lhs = lhs, op = op, rhs = rhs, mod = mod)
  )

  newSyntax <- parTableToSyntax(newParTable)

  out <- list(syntax = newSyntax, data = newData, parTable = newParTable, 
              reliability = rel.std, AVE = AVE.std, lVs = lVs, single.items = singleInds)

  class(out) <- c("list", "modsem_relcorr")

  out
}


#' @export
print.modsem_relcorr <- function(x,...) {
  sep <- strrep(" ", 4)
  indent <- strrep(" ", 2)

  lVs <- x$lVs
  lVsf <- format(paste0(lVs, ":"), justify = "left")
  AVE_rel <- format(c(x$AVE, x$reliability), digits = 3, nsmall = 2, justify = "right")

  idxAVE <- seq_along(x$AVE)
  AVEs  <- AVE_rel[idxAVE]
  rels  <- AVE_rel[-idxAVE]

  cat("Average Variance Extracted:\n")
  for (i in seq_along(lVs)) {
    lV <- lVsf[[i]]
    AVE <- AVEs[[i]]
    cat(paste0(indent, lV, sep, AVE, "\n"))
  }

  cat(paste0("\nChronbach's Alpha:\n"))
  for (i in seq_along(lVs)) {
    lV <- lVsf[[i]]
    rel <- rels[[i]]
    cat(paste0(indent, lV, sep, rel, "\n"))
  }

  cat(paste0("\nGenerated Syntax:\n"))
  cat(padString(x$syntax, pad = indent))

  cat(paste0("Generated Items:\n"))
  datastr <- utils::capture.output(utils::str(x$data[x$single.items]))
  cat(padString(datastr, pad = indent))
}


getSingleIndNames <- function(lVs) {
  structure(paste0("composite_", lVs, "_"), names=lVs)
}


calcChronbach <- function(R, x) {
  Rx <- R[x, x] 
  k <- length(x)
  
  if (k == 1) return(1) # should be treated like an observed variable

  warnif(any(Rx < 0), "Some item covariances are negative! Consider ",
         "recoding your items!")
  
  Rx <- abs(Rx)
  varSum <- sum(diag(Rx))
  covSum <- sum(diag(Rx[lower.tri(Rx)]))

  (k / (k - 1)) * (1 - varSum / (varSum + 2 * covSum))
}


calcConstructRel <- function(lambda.std, lV) {
  lambda.std <- lambda.std[, lV, drop=FALSE]
  lambda.std <- lambda.std[lambda.std != 0]

  warnif(any(lambda.std < 0), "Some items have negative factor loadings!\n",
         "Consider recoding your items!")

  sum(lambda.std) ^ 2 / (sum(lambda.std) ^ 2 + sum(1 - lambda.std ^ 2))
}


calcAVE <- function(lambda.std, theta.std, lV) {
  inds <- rownames(lambda.std)

  lambda.std <- structure(lambda.std[, lV, drop=TRUE], names = inds)
  lambda.std <- lambda.std[lambda.std != 0]
  theta.std <- diag(theta.std)[names(lambda.std)]

  sum(lambda.std) ^ 2 / (sum(lambda.std) ^ 2 + sum(1 - theta.std ^ 2))
}
