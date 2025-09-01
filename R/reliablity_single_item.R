#' Reliability‑Corrected Single‑Item SEM
#'
#' Replace (some of) the first‑order latent variables in a lavaan measurement
#' model by single composite indicators whose error variances are fixed from
#' Cronbach's \eqn{\alpha}. The function returns a modified lavaan model
#' syntax together with an augmented data set that contains the newly created
#' composite variables, so that you can fit the full SEM in a single step.
#'
#' The resulting object can be fed directly into \code{modsem} or \code{lavaan::sem}
#' by supplying \code{syntax = $syntax} and \code{data = $data}.
#'
#' @param syntax A character string containing lavaan model syntax.  Must at
#'   least include the measurement relations (\code{=~}).
#'
#' @param data A \code{data.frame}, \code{tibble} or object coercible to a data frame that
#'   holds the raw observed indicators.
#'
#' @param choose \emph{Optional.} Character vector with the names of the latent
#'   variables you wish to replace by single indicators.  Defaults to \strong{all}
#'   first‑order latent variables in \code{syntax}.
#'
#' @param scale.corrected Should reliability-corrected items be scale-corrected? If \code{TRUE}
#'   reliability-corrected single items are corrected for differences in factor loadings between
#'   the items. Default is \code{TRUE}.
#'
#' @param warn.lav Should warnings from \code{lavaan::cfa} be displayed? If \code{FALSE}, they
#'   are suppressed.
#'
#' @return An object of S3 class \code{modsem_relcorr} (a named list) with elements:
#' \describe{
#'   \item{\code{syntax}}{Modified lavaan syntax string.}
#'   \item{\code{data}}{Data frame with additional composite indicator columns.}
#'   \item{\code{parTable}}{Parameter table corresponding to `syntax`.}
#'   \item{\code{reliability}}{Named numeric vector of reliabilities (one per latent
#'     variable).}
#'   \item{\code{AVE}}{Named numeric vector with Average Variance Extracted values.}
#'   \item{\code{lVs}}{Character vector of latent variables that were corrected.}
#'   \item{\code{single.items}}{Character vector with names for the generated reliability corrected items}
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
#' summary(est_lms)
#' }
#'
#' @export
relcorr_single_item <- function(syntax,
                                data,
                                choose = NULL,
                                scale.corrected = TRUE,
                                warn.lav = TRUE) {
  data <- as.data.frame(data)

  parTable      <- modsemify(syntax)
  intTerms      <- parTable$rhs[grepl(":", parTable$rhs)]

  parTableOuter <- parTable[parTable$op == "=~", ]
  parTableInner <- parTable[parTable$op != "=~", ]

  lVs            <- getLVs(parTable)
  allIndsLVs     <- getIndsLVs(parTable, lVs = lVs)
  allInds        <- unique(unlist(allIndsLVs))
  higherOrderLVs <- getHigherOrderLVs(parTable)
  firstOrderLVs  <- setdiff(lVs, higherOrderLVs)

  stopif(any(higherOrderLVs %in% allInds),
         "Third order (or higher) latent variables are not allowed (yet)!")

  choose <- if (is.null(choose)) lVs else choose
  ignore <- setdiff(lVs, choose)

  warnif(any(!choose %in% lVs), "Could not find latent variables:\n",
         paste(choose[!choose %in% lVs], collapse = ", "))
  choose <- choose[choose %in% lVs]

  ignoreRows <- parTableOuter$lhs %in% ignore

  parTableInner <- rbind(parTableInner, parTableOuter[ignoreRows, ])
  parTableOuter <- parTableOuter[!ignoreRows, ]

  stopif(!NROW(parTableOuter), "Cannot ignore all latent variables when creating ",
         "reliability corrected items!")

  warnif(any(parTableOuter$mod != ""),
         "Labels and modifiers in the measurement model used to create reliability ",
         "corrected single items will be ignored!\nThis may cause the model constraints ",
         "in the structural model to break!", immediate. = FALSE)

  lVs  <- setdiff(firstOrderLVs, ignore)
  indsLVs <- getIndsLVs(parTableOuter, lVs = lVs)
  inds <- unlist(indsLVs)
  R <- stats::cov(data[inds], use = "pairwise.complete.obs")
  k <- length(lVs)

  # Add residual Covariances
  rescovind <- parTableInner$lhs %in% inds &
               parTableInner$rhs %in% inds &
               parTableInner$op == "~~"
  parTableInner <- parTableInner[!rescovind, ]
  parTableOuter <- rbind(parTableOuter, parTableInner[rescovind, ])

  # Check if there are any indicators in the structural model
  indsInInner <- parTableInner$lhs %in% inds | parTableInner$rhs %in% inds
  if (any(indsInInner)) {
    warning2("Removing expressions containing indicators!\n",
             capturePrint(parTableInner[indsInInner, c("lhs", "op", "rhs")]),
             immediate. = FALSE)
    parTableInner <- parTableInner[!indsInInner, ]
  }

  wrap      <- if (warn.lav) \(x) x else suppressWarnings
  cfaSyntax <- parTableToSyntax(parTableOuter)
  cfa       <- wrap(lavaan::cfa(cfaSyntax, data = data, se = "none"))

  cov.lv      <- lavaan::lavInspect(cfa, "cov.lv")
  stdSolution <- lavaan::lavInspect(cfa, "std")
  lambda <- stdSolution$lambda
  theta  <- stdSolution$theta

  avefunc <- \(lV) calcAVE(lV, theta.std = theta, lambda.std = lambda)
  relfunc <- \(lV) calcChronbach(R = R, x = indsLVs[[lV]])

  rel.std <- vapply(lVs, FUN.VALUE = numeric(1L), relfunc)
  AVE.std <- vapply(lVs, FUN.VALUE = numeric(1L), avefunc)
  res.std <- 1 - rel.std
  singleInds <- getSingleIndNames(lVs)

  rchar <- c(letters, seq_len(9))
  while(any(singleInds %in% colnames(data))) {
    # if matching names in data, add random characters, until there is no
    # match. Should in practive never be triggered
    singleInds <- paste0(singleInds, sample(rchar, length(singleInds)))
  }

  varNewInds <- structure(numeric(length(lVs)), names=singleInds)
  newData <- data

  if (scale.corrected) {
    res <- stats::setNames(numeric(length(lVs)), nm = lVs)

    for (lV in lVs) {
      ITEM <- getScaleCorrectedItem(parTable = parTable, lV = lV,
                                    data = data, cfa = cfa)

      newIndName <- singleInds[[lV]]
      newInd     <- ITEM$item
      residual   <- ITEM$residual

      newData[[newIndName]] <- newInd
      res[[lV]] <- residual
    }

  } else {
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
  }

  if (length(higherOrderLVs)) {
    for (lV in intersect(lVs, allInds)) {
      # If lV is an indicator for a higher order variable
      # we add a correction to the residual, to account for
      # measurement error on both levels
      res[[lV]] <- res[[lV]] + getHigherOrderResidual(lV = lV, cfa = cfa)
    }
  }

  # first order lVs
  newParTable <- parTableInner
  firstOrder  <- !lVs %in% allInds
  secondOrder <- lVs %in% allInds

  if (any(firstOrder)) {
    lVs1        <- lVs[firstOrder]
    singleInds1 <- singleInds[firstOrder]
    k1          <- length(lVs1)
    res1        <- res[firstOrder]

    lhs <- unname(c(lVs1, singleInds1))
    rhs <- unname(c(singleInds1, singleInds1))
    op  <- c(rep("=~", k1), rep("~~", k1))
    mod <- unname(c(rep("1", k1), res1))

    newParTable <- rbind(
      newParTable, data.frame(lhs = lhs, op = op, rhs = rhs, mod = mod)
    )
  }

  if (any(secondOrder)) {
    message(
      "Correcting first order items/lVs for higher order measurement model!\n",
      "Note that the first order LVs no longer can be used directly in the\n",
      "structural model!"
    )

    lVs2        <- lVs[secondOrder]
    singleInds2 <- singleInds[secondOrder]
    k2          <- length(lVs2)
    res2        <- res[secondOrder]

    for (hlV in higherOrderLVs) {
      lVs2h <- intersect(lVs2, allIndsLVs[[hlV]])
      singleInds2h <- singleInds2[lVs2h]
      res2h <- res2[lVs2h]
      k2h   <- length(singleInds2h)

      lhs <- unname(c(rep(hlV, k2h), singleInds2h))
      rhs <- unname(c(singleInds2h, singleInds2h))
      op  <- c(rep("=~", k2h), rep("~~", k2h))
      mod <- unname(c(rep("", k2h), res2h))

      newParTable <- rbind(
        newParTable, data.frame(lhs = lhs, op = op, rhs = rhs, mod = mod)
      )
    }
  }

  newSyntax <- parTableToSyntax(newParTable)

  out <- list(syntax = newSyntax, data = newData,
              parTable = newParTable, reliability = rel.std,
              AVE = AVE.std, lVs = lVs, single.items = singleInds,
              intTerms = intTerms, residuals = res, res.std = res.std,
              cov.lv = cov.lv)

  class(out) <- c("list", "modsem_relcorr")

  out
}


getScaleCorrectedItem <- function(parTable, lV, data, cfa) {
  measr  <- parTable[parTable$lhs == lV & parTable$op == "=~", ]
  inds   <- unique(measr$rhs)

  matrices <- lavaan::lavInspect(cfa, what = "coef")

  theta  <- matrices$theta[inds, inds]
  lambda <- as.vector(matrices$lambda[inds, lV])
  I      <- rep(1, length(inds))

  X    <- as.matrix(data[, inds])
  item <- rowSums(X) / sum(lambda)

  residual <- (t(I) %*% theta %*% I) / sum(lambda)^2 # variance of (e1 + e2 + ... + en)

  list(X = X, item = item, residual = residual)
}


getHigherOrderResidual <- function(lV, cfa) {
  parTable <- lavaan::parameterEstimates(cfa)
  inds     <- getInds(parTable)

  if (!lV %in% inds) return(0)

  residual <- parTable[parTable$lhs == lV &
                       parTable$op == "~~" &
                       parTable$rhs == lV, "est"]

  ifelse(!length(residual), yes = 0, no = residual[[1L]])
}


simulateCrosssimulateCrossResCovRCS <- function(corrected,
                                                elemsInIntTerms,
                                                mc.reps = 3e4,
                                                parTable,
                                                include.normal.inds = FALSE) {
  EMPTY <- list(syntax = "", rows = NULL)

  if (!length(elemsInIntTerms))
    return(EMPTY)

  intTerms        <- names(elemsInIntTerms)
  single.items    <- corrected$single.items
  residuals       <- corrected$residuals
  constructs      <- names(single.items)

  intTermIsRCS <- vapply(elemsInIntTerms, FUN.VALUE = logical(1L),
                         FUN = \(elems) all(elems %in% constructs))
  elemsInIntTerms <- elemsInIntTerms[intTermIsRCS]

  if (!length(elemsInIntTerms))
    return(EMPTY)

  Phi   <- corrected$cov.lv
  Theta <- diag(residuals)
  nm    <- names(residuals)
  Phi   <- Phi[nm, nm]
  dimnames(Theta) <- list(nm, nm)

  XI  <- mvtnorm::rmvnorm(mc.reps, sigma = Phi)
  EPS <- mvtnorm::rmvnorm(mc.reps, sigma = Theta)
  Y   <- XI + EPS

  CONSTRUCTS <- as.data.frame(XI)
  ITEMS      <- as.data.frame(Y)

  colnames(CONSTRUCTS) <- nm
  colnames(ITEMS)      <- single.items[nm]

  if (include.normal.inds) {
    RESIDUALS  <- as.data.frame(EPS)
    colnames(RESIDUALS)  <- single.items[nm]
  } else RESIDUALS <- list()

  getres <- function(y, x) {
    missing <- is.na(y) | is.na(x)

    x.c <- x[!missing]
    y.c <- y[!missing]

    epsilon <- stats::residuals(stats::lm(x.c ~ y.c))

    out <- rep(NA, length(epsilon))
    out[!missing] <- epsilon

    out
  }

  for (intTerm in intTerms) {
    elems  <- elemsInIntTerms[[intTerm]]
    items  <- single.items[elems]
    nameXZ <- paste0(items, collapse = "")

    xz.construct <- apply(CONSTRUCTS[elems], MARGIN = 1, FUN = prod)
    xz.item      <- apply(ITEMS[items],  MARGIN = 1, FUN = prod)
    xz.residual  <- getres(y = xz.item, x = xz.construct)

    CONSTRUCTS[[intTerm]] <- xz.construct
    RESIDUALS[[nameXZ]]   <- xz.residual
    ITEMS[[intTerm]]      <- xz.item
  }

  Epsilon <- stats::cov(as.data.frame(RESIDUALS), use = "na.or.complete")
  newRows <- NULL

  for (i in seq_len(NROW(Epsilon))) {
    for (j in seq_len(i)) {
      mod <- Epsilon[i, j]
      lhs <- rownames(Epsilon)[[i]]
      rhs <- colnames(Epsilon)[[j]]

      matching <- parTable[((parTable$lhs == rhs & parTable$rhs == lhs) |
                            (parTable$lhs == rhs & parTable$rhs == lhs)) &
                           parTable$op == "~~", ]
      if (NROW(matching)) next

      newRows <- rbind(
        newRows,
        createParTableRow(c(lhs, rhs), op = "~~", mod = mod)
      )
    }
  }

  syntaxAdditions <- parTableToSyntax(newRows)

  list(syntax = syntaxAdditions, rows = newRows)
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
  covSum <- sum(Rx[lower.tri(Rx)])

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
