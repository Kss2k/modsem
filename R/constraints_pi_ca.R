# Functitions for specifying constraints in the constrained approach


labelFactorLoadings <- function(parTable) {
  # Firstly Label factor loadings in Partable
  loadingLabels <- apply(parTable[parTable$op == "=~", c("rhs", "lhs")],
                         MARGIN = 1,
                         FUN = function(x) createLabelLambda(x[1], x[2]))
  parTable[parTable$op == "=~", "mod"] <- loadingLabels
  parTable
}


specifyFactorLoadingsSingle <- function(parTable, relDf) {
  latentProdName <- stringr::str_c(rownames(relDf), collapse = "")

  for (indProd in colnames(relDf)) {
    indProdLabel     <- createLabelLambda(indProd, latentProdName)
    indsInProdLabels <- createLabelLambda(relDf[[indProd]], rownames(relDf))
    vecLhsRhs        <- c(indProdLabel, stringr::str_c(indsInProdLabels,
                                                       collapse = " * "))

    newRow   <- createParTableRow(vecLhsRhs, op = "==")
    parTable <- rbind(parTable, newRow)
  }

  parTable
}


specifyFactorLoadings <- function(parTable, relDfs) {
  parTable <- specifyFactorLoadingsSingle(parTable, relDfs[[1]])
  if (length(relDfs) <= 1) return(parTable)
  specifyFactorLoadings(parTable, relDfs[-1])
}


addVariances <- function(pt) {
  # Add variance-labels if missing for LVs
  LVs <- getLVs(pt)
  if (length(LVs) == 0) return(pt)

  specifiedLVs <- pt[pt$lhs %in% LVs &
                     pt$op == "~~" &
                     pt$rhs %in% LVs &
                     pt$lhs == pt$rhs, "lhs"] |> unique()
  toBeSpecifiedLvs <- LVs[!LVs %in% specifiedLVs]

  newRows <- lapply(toBeSpecifiedLvs, FUN = function(x)
                    createParTableRow(c(x, x), op = "~~")) |>
    purrr::list_rbind()
  pt <- rbind(pt, newRows)

  # Add variance for observed variables if missing
  observed <- getIndicators(pt, observed = TRUE)

  specifiedOvs <- pt[pt$lhs %in% observed &
                     pt$op == "~~" &
                     pt$rhs %in% observed &
                     pt$lhs == pt$rhs, "lhs"] |> unique()
  toBeSpecifiedOvs <- observed[!observed %in% specifiedOvs]

  newRows <- lapply(toBeSpecifiedOvs,
                    FUN = function(x) createParTableRow(c(x, x), op = "~~")) |>
    purrr::list_rbind()
  rbind(pt, newRows)
}


addCovariances <- function(pt) {
  # Add covariances for exogenous variables if missing
  pt       <- stripColonsParTable(pt)
  latents  <- getLVs(pt)
  if (length(latents) == 0) return(pt)

  higherOrder <- getHigherOrderLVs(pt)
  firstOrder  <- latents[!latents %in% higherOrder]

  if (length(higherOrder)) { # then we should not touch the first order covariances
    combosHO <- getUniqueCombos(higherOrder)
    ptHO     <- pt[pt$lhs %in% higherOrder & pt$rhs %in% higherOrder, ]
    combosHO$connected <- !is.na(apply(combosHO, MARGIN = 1, function(xy)
                                       trace_path(ptHO, xy[[1]], xy[[2]])))
    toBeSpecified <- combosHO[!combosHO$connected, c("V1", "V2")]
    if (nrow(toBeSpecified)) {
      newRows <- apply(toBeSpecified[c("V1", "V2")],
                       MARGIN = 1,
                       FUN = function(x) createParTableRow(x, op = "~~")) |>
        purrr::list_rbind()
      pt <- rbind(pt, newRows)
    }
  } else {
    combosFO <- getUniqueCombos(firstOrder)
    ptFO     <- pt[pt$lhs %in% firstOrder & pt$rhs %in% firstOrder, ]
    combosFO$connected <- !is.na(apply(combosFO, MARGIN = 1, function(xy)
                                       trace_path(ptFO, xy[[1]], xy[[2]])))
    toBeSpecified <- combosFO[!combosFO$connected, c("V1", "V2")]
    if (nrow(toBeSpecified)) {
      newRows <- apply(toBeSpecified[c("V1", "V2")],
                       MARGIN = 1,
                       FUN = function(x) createParTableRow(x, op = "~~")) |>
        purrr::list_rbind()
      pt <- rbind(pt, newRows)
    }
  }

  pt
}


labelParameters <- function(pt) {
  latents <- getLVs(pt)
  endogenous <- latents[latents %in% pt[pt$op == "~", "lhs"]]
  exogenous <- latents[!latents %in% endogenous]
  observed <- getOVs(pt)

  # Gamma
  pt[pt$op == "~", "mod"] <-
    apply(pt[pt$op == "~", c("rhs", "lhs")],
          MARGIN = 1, FUN = function(x)
            createLabelGamma(x[[1]], x[[2]]))

  # Variances of exogenous
  pt[pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% exogenous, "mod"] <-
    vapply(pt[pt$op == "~~" & pt$lhs == pt$rhs &
              pt$lhs %in% exogenous, "lhs"],
           FUN.VALUE = vector("character", length = 1L),
           FUN = createLabelVar)

  # Variances of endogenous
  pt[pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% endogenous, "mod"] <-
    vapply(pt[pt$op == "~~" & pt$lhs == pt$rhs &
              pt$lhs %in% endogenous, "lhs"],
           FUN.VALUE = vector("character", length = 1L),
           FUN = createLabelZeta)

  # Variance of Observed Variables
  pt[pt$op == "~~" & pt$rhs %in% observed & pt$lhs == pt$rhs, "mod"] <-
    vapply(pt[pt$op == "~~" & pt$rhs %in% observed &
           pt$lhs == pt$rhs, "rhs"],
           FUN.VALUE = vector("character", length = 1L),
           FUN = createLabelVar)

  # Covariances
  pt[pt$op == "~~" & pt$lhs != pt$rhs, "mod"] <-
    apply(pt[pt$op == "~~" & pt$lhs != pt$rhs, c("lhs", "rhs")],
          MARGIN = 1, FUN = function(x)
            createLabelCov(x[[1]], x[[2]]))

  pt
}


specifyVarCovSingle <- function(parTable, relDf) {
  # This function specifies variances for latents, indicators,
  # and indicator products. It will also specifies covariances for latent
  # products, and elements int those products.
  stopif(nrow(relDf) > 2, "Constraints for products with more than two ",
         " elements are not supported for this method")

  # General info
  elemsInProdTerm <- rownames(relDf)
  latentProd <- stringr::str_c(rownames(relDf), collapse = "")

  # Variance of latent product
  labelLatentProd <- createLabelVar(latentProd)
  labelsElemsInProd <- vapply(elemsInProdTerm,
                              FUN.VALUE = vector("character", length = 1L),
                              FUN = function(x) trace_path(parTable, x, x))

  labelCovElems <- trace_path(parTable, elemsInProdTerm[[1]],
                              elemsInProdTerm[[2]]) |> paste0(" ^ 2")

  lhs <- labelLatentProd
  rhs <- paste(stringr::str_c(labelsElemsInProd, collapse = " * "),
               labelCovElems, sep = " + ")
  varLatentProd <- createParTableRow(c(lhs, rhs), op = "==")

  # covariances between elems and latents
  labelsCovElemProd <-
    vapply(elemsInProdTerm, FUN.VALUE = vector("character", length = 1L),
           FUN = function(elem) createLabelCov(elem, latentProd)) # wrap in anonymous fun
                                                             # scope latentProd
  
  labelsCovElemProd <- labelsCovElemProd[labelsCovElemProd %in% parTable$mod]

  if (length(labelsCovElemProd)) { # should not be added in higher order models
    covsElemsProd <- lapply(labelsCovElemProd, FUN = function(x)
                              createParTableRow(c(x, "0"), op = "==")) |>
      purrr::list_rbind()
  } else covsElemsProd <- NULL

  # Variances of product indicators
  constrained.varProdInds <- vector("list", length = ncol(relDf))

  for (indProd in colnames(relDf)) {
    labelVarIndProd <- createLabelVar(indProd)

    labelsFactorLoadings <- vector("character", length = nrow(relDf))
    labelsVarLatents     <- vector("character", length = nrow(relDf))
    labelsVarInds        <- vector("character", length = nrow(relDf))

    for (latent in seq_len(nrow(relDf))) {
      labelsFactorLoadings[[latent]] <-
        createLabelLambdaSquared(relDf[latent, indProd],
                                 rownames(relDf)[[latent]])
      labelsVarLatents[[latent]] <-
        trace_path(parTable, rownames(relDf)[[latent]],
                   rownames(relDf)[[latent]])
      labelsVarInds[[latent]] <- createLabelVar(relDf[latent, indProd])
    }

    lhs  <- labelVarIndProd
    rhs1 <- paste(labelsFactorLoadings[[1]], labelsVarLatents[[1]],
                  labelsVarInds[[2]], sep = " * ")
    rhs2 <- paste(labelsFactorLoadings[[2]], labelsVarLatents[[2]],
                  labelsVarInds[[1]], sep = " * ")
    rhs3 <- paste(labelsVarInds[[1]], labelsVarInds[[2]], sep = " * ")
    rhs  <- paste(rhs1, rhs2, rhs3, sep = " + ")

    constrained.varProdInds[[indProd]] <- createParTableRow(c(lhs, rhs), op = "==")
  }

  constrained.varProdInds <- purrr::list_rbind(constrained.varProdInds)
  rbindParTable(parTable, rbind(varLatentProd,
                                covsElemsProd,
                                constrained.varProdInds))
}


specifyVarCov <- function(parTable, relDfs) {
  parTable <- specifyVarCovSingle(parTable, relDfs[[1]])
  if (length(relDfs) <= 1) return(parTable)
  specifyVarCov(parTable, relDfs[-1])
}
