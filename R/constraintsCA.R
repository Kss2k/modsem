labelFactorLoadings <- function(parTable) {
  # Firstly Label factor loadings in Partable
  loadingLabels <- apply(parTable[parTable$op == "=~", c("rhs", "lhs")],
                         MARGIN = 1,
                         FUN = function(x) createLabelLambda(x[1], x[2])
  )
  parTable[parTable$op == "=~", "mod"] <- loadingLabels
  parTable
}



specifyFactorLoadingsSingle <- function(parTable, relDf) {
  latentProdName <- stringr::str_c(rownames(relDf), collapse = "")
  for (indProd in colnames(relDf)) {
    indProdLabel <- createLabelLambda(indProd, latentProdName)
    indsInProdLabels <- createLabelLambda(
      relDf[[indProd]],
      rownames(relDf)
    )
    vecLhsRhs <- c(indProdLabel, stringr::str_c(indsInProdLabels,
      collapse = " * "
    ))
    newRow <- createParTableRow(vecLhsRhs, op = "==")
    parTable <- rbind(parTable, newRow)
  }
  parTable
}



specifyFactorLoadings <- function(parTable, relDfs) {
  parTable <- specifyFactorLoadingsSingle(parTable, relDfs[[1]])

  if (length(relDfs) <= 1) {
    return(parTable)
  } else {
    return(specifyFactorLoadings(parTable, relDfs[-1]))
  }
}


addVariances <- function(pt) {
  # Add variances if missing for latents
  latents <- pt[pt$op == "=~", "lhs"] |> unique()
  if (length(latents) == 0) return(pt)
  alreadySpecifed <- pt[pt$lhs %in% latents &
    pt$op == "~~" &
    pt$rhs %in% latents &
    pt$lhs == pt$rhs, "lhs"] |> unique()
  toBeSpecified <- latents[!latents %in% alreadySpecifed]
  newRows <- lapply(toBeSpecified,
                    FUN = function(x) createParTableRow(c(x, x), op = "~~")) |>
    purrr::list_rbind()
  pt <- rbind(pt, newRows)
  # Add variance for observed variables if missing 
  observed <- pt[pt$op == "=~", "rhs"] |> unique()
  alreadySpecifed <- pt[pt$lhs %in% observed &
    pt$op == "~~" &
    pt$rhs %in% observed &
    pt$lhs == pt$rhs, "lhs"] |> unique()
  toBeSpecified <- observed[!observed %in% alreadySpecifed]
  newRows <- lapply(toBeSpecified,
                    FUN = function(x) createParTableRow(c(x, x), op = "~~")) |>
    purrr::list_rbind()
  rbind(pt, newRows)
}


addCovariances <- function(pt) {
  # Add covariances for exogenous variables if missing
  latents <- pt[pt$op == "=~", "lhs"] |> unique()
  if (length(latents) == 0) return(pt)
  combos <- getUniqueCombinations(latents)
  combos$connected <- !is.na(apply(combos, MARGIN = 1, function(xy) {
    tracePath(pt, xy[1], xy[2])
  }))
  toBeSpecified <- combos[!combos$connected, c("V1", "V2")]
  newRows <- apply(toBeSpecified[c("V1", "V2")],
                   MARGIN = 1,
                   FUN = function(x) createParTableRow(x, op = "~~")) |>
    purrr::list_rbind()
  rbind(pt, newRows)
}

labelParameters <- function(pt) {
  latents <- pt[pt$op == "=~", "lhs"] |> unique()
  endogenous <- latents[latents %in% pt[pt$op == "~", "lhs"]]
  exogenous <- latents[!latents %in% endogenous]
  observed <- pt[pt$op == "=~", "rhs"] |> unique()

  # Gamma
  pt[pt$op == "~" & pt$rhs != "1", "mod"] <-
    apply(pt[pt$op == "~" & pt$rhs != "1", c("rhs", "lhs")],
      MARGIN = 1,
      FUN = function(x) {
        stringr::str_c(c("Gamma", x), collapse = "_")
      }
    )
  # Variances of exogenous
  pt[pt$op == "~~" & pt$lhs == pt$rhs &
    pt$lhs %in% exogenous, "mod"] <-
    vapply(pt[pt$op == "~~" & pt$lhs == pt$rhs &
          pt$lhs %in% exogenous, "lhs"],
          FUN.VALUE = vector("character", length = 1L),
          FUN = function(x) paste("Var", x, sep = "_")
    )
  # Variances of endogenous
  pt[pt$op == "~~" & pt$lhs == pt$rhs &
    pt$lhs %in% endogenous, "mod"] <-
    vapply(pt[pt$op == "~~" & pt$lhs == pt$rhs &
          pt$lhs %in% endogenous, "lhs"],
          FUN.VALUE = vector("character", length = 1L),
          FUN = function(x) paste("Zeta", x, sep = "_")
    )

  # Variance of Observed Variables  
  pt[pt$op == "~~" & pt$rhs %in% observed &
     pt$lhs == pt$rhs, "mod"] <-
    vapply(pt[pt$op == "~~" & pt$rhs %in% observed &
           pt$lhs == pt$rhs, "rhs"],
          FUN.VALUE = vector("character", length = 1L),
          FUN = function(x) paste("Var", x, sep = "_"))

  # Covariances
  pt[pt$op == "~~" & pt$lhs != pt$rhs, "mod"] <-
    apply(pt[pt$op == "~~" & pt$lhs != pt$rhs, c("lhs", "rhs")],
      MARGIN = 1,
      FUN = function(x) {
        stringr::str_c(c("Cov", x), collapse = "_")
      }
    )
  pt
}


specifyVarCovSingle <- function(parTable, relDf) {
  # This funciton will specify variances for latents, indicators,
  # and indicator products. It will also specify covariances for latent
  # products, and elements int those products.
  if (nrow(relDf) > 2) {
    stop(
      "Constraints for products with more than two elements are not supported yet:",
      capturePrint(print(relDf))
    )
  }
  # General info
  elemsInProdTerm <- rownames(relDf)
  latentProd <- stringr::str_c(rownames(relDf), collapse = "")

  # Variance of latent product
  labelLatentProd <- createLabelVar(latentProd)
  
  labelsElemsInProd <- vapply(elemsInProdTerm, 
                               FUN.VALUE = vector("character", length = 1L),
                               FUN = function(x) tracePath(parTable, x, x))

  labelCovElems <- tracePath(parTable,
                             elemsInProdTerm[[1]],
                             elemsInProdTerm[[2]]) |>
    paste0(" ^ 2")

  lhs <- labelLatentProd
  rhs <- paste(stringr::str_c(labelsElemsInProd, collapse = " * "),
    labelCovElems,
    sep = " + "
  )
  varLatentProd <- createParTableRow(c(lhs, rhs), op = "==")

  # covariances between elems and latents
  labelsCovElemProd <-
    vapply(elemsInProdTerm,
      FUN.VALUE = vector("character", length = 1L),
      FUN = function(elem, latentProd) {
        createLabelCov(elem, latentProd)
      },
      latentProd = latentProd
    )
  covsElemsProd <- lapply(labelsCovElemProd,
    FUN = function(x) {
      createParTableRow(c(x, "0"), op = "==")
    }
  ) |>
    purrr::list_rbind()

  # Variances of product indicators
  constrainedVarProdInds <- vector("list", length = ncol(relDf))

  for (indProd in colnames(relDf)) {
    labelVarIndProd <- createLabelVar(indProd)

    labelsFactorLoadings <- vector("character", length = nrow(relDf))
    labelsVarLatents <- vector("character", length = nrow(relDf))
    labelsVarInds <- vector("character", length = nrow(relDf))

    for (latent in 1:nrow(relDf)) {
      labelsFactorLoadings[[latent]] <- createLabelLambdaSquared(
        relDf[latent, indProd],
        rownames(relDf)[[latent]]
      )
      labelsVarLatents[[latent]] <- tracePath(parTable,
                                              rownames(relDf)[[latent]],
                                              rownames(relDf)[[latent]])
      labelsVarInds[[latent]] <- createLabelVar(relDf[latent, indProd])
    }
    lhs <- labelVarIndProd
    rhs1 <- paste(labelsFactorLoadings[[1]],
      labelsVarLatents[[1]],
      labelsVarInds[[2]],
      sep = " * "
    )
    rhs2 <- paste(labelsFactorLoadings[[2]],
      labelsVarLatents[[2]],
      labelsVarInds[[1]],
      sep = " * "
    )
    rhs3 <- paste(labelsVarInds[[1]], labelsVarInds[[2]], sep = " * ")

    rhs <- paste(rhs1, rhs2, rhs3, sep = " + ")

    constrainedVarProdInds[[indProd]] <- createParTableRow(c(lhs, rhs), op = "==")
  }
  constrainedVarProdInds <- purrr::list_rbind(constrainedVarProdInds)
  rbindParTable(parTable, rbind(
    varLatentProd,
    covsElemsProd,
    constrainedVarProdInds
  ))
}



specifyVarCov <- function(parTable, relDfs) {
  parTable <- specifyVarCovSingle(parTable, relDfs[[1]])

  if (length(relDfs) <= 1) {
    return(parTable)
  } else {
    return(specifyVarCov(parTable, relDfs[-1]))
  }
}

# DEPRECATEDlabelVarCovSingle <- function(parTable, relDf, isCovConflictCorrected = FALSE) {
#   # Some general Information
#   elemsInProdTerm <- rownames(relDf)
#   latentProd <- stringr::str_c(rownames(relDf), collapse = "")
#   prodInds <- colnames(relDf)
#   indsInProds <- unique(unlist(relDf))
#
#   # Adding variances to parTable -----------------------------------------------
#   if (!isCovConflictCorrected &&
#       NROW(parTable[parTable$lhs %in% elemsInProdTerm &
#                     parTable$op == "~", ]) > 0) {
#     warning("Warning adding a variance label to an endogenous variable \n",
#             "This will likely be a mistake since ̈́'eta ~~ eta' gives the ",
#             "variance of the disturbace term of the eta variable, ",
#             "not its variance. \n",
#             "Use another method, or alter the parTable ",
#             "using \nremoveFromParTable = ... \naddToParTable = ...\n")
#   }
#   variancesToBeSpecified <- c(latentProd,
#                               elemsInProdTerm,
#                               prodInds,
#                               indsInProds)
#   newRows <- lapply(variancesToBeSpecified,
#                     FUN = function(xi)
#                       createParTableRow(c(xi, xi), op = "~~")) |>
#     purrr::list_rbind()
#
#
#   newRows[["mod"]] <- vapply(newRows[["lhs"]],
#                           FUN.VALUE = vector("character", length = 1L),
#                           FUN = function(x)
#                             stringr::str_c(c("Var", x),
#                                            collapse = "_")
#                           )
#
#   # Covariances latents --------------------------------------------------------
#   covsElemsInProd <- apply(getUniqueCombinations(elemsInProdTerm),
#         MARGIN = 1,
#         FUN = createParTableRow,
#         op = "~~") |>
#     purrr::list_rbind()
#
#   covsElemToProd <-
#     lapply(elemsInProdTerm,
#            FUN = function(elem, latentProd)
#              createParTableRow(c(elem, latentProd), op = "~~"),
#            latentProd = latentProd) |>
#     purrr::list_rbind()
#
#   covs <- rbind(covsElemsInProd, covsElemToProd)
#   # Label the covs
#   covs[["mod"]] <- apply(covs[c("lhs", "rhs")],
#                        MARGIN = 1,
#                        FUN = function(x)
#                          stringr::str_c(c("Cov", x), collapse = "_"))
#   newRows <- rbind(newRows, covs)
#   rbindParTable(parTable, newRows)
# }



# DEPRECATEDlabelVarCov <- function(parTable, relDfs, isCovConflictCorrected = FALSE) {
#  parTable <- labelVarCovSingle(parTable, relDfs[[1]], isCovConflictCorrected)
#
#  if (length(relDfs) <= 1) {
#    return(parTable)
#  } else {
#    return(labelVarCov(parTable, relDfs[-1], isCovConflictCorrected))
#  }
# }


# DEPRECATEDspecifyVarCovSingle <- function(parTable, relDf) {
#   # This funciton will specify variances for latents, indicators,
#     # and indicator products. It will also specify covariances for latent
#     # products, and elements int those products.
#   if (nrow(relDf) > 2) {
#     stop("Constraints for products with more than two elements are not supported yet:",
#          capturePrint(print(relDf)))
#   }
#   # General info
#   elemsInProdTerm <- rownames(relDf)
#   latentProd <- stringr::str_c(rownames(relDf), collapse = "")

#   # Variance of latent product
#   labelLatentProd <- createLabelVar(latentProd)
#   labelsElemsInProd <- createLabelVar(elemsInProdTerm)
#
#   labelCovElems <- createLabelCov(elemsInProdTerm[[1]],
#                                   elemsInProdTerm[[2]]) |>
#     paste0(" ^ 2")
#   lhs <- labelLatentProd
#   rhs <- paste(stringr::str_c(labelsElemsInProd, collapse = " * "),
#                labelCovElems, sep = " + ")
#   varLatentProd <- createParTableRow(c(lhs, rhs), op = "==")

#   # covariances between elems and latents
#   labelsCovElemProd <-
#     vapply(elemsInProdTerm,
#            FUN.VALUE = vector("character", length = 1L),
#            FUN = function(elem, latentProd)
#              createLabelCov(elem, latentProd),
#            latentProd = latentProd)
#   covsElemsProd <- lapply(labelsCovElemProd,
#                           FUN = function(x)
#                             createParTableRow(c(x, "0"), op = "==")) |>
#     purrr::list_rbind()

#   # Variances of product indicators
#   constrainedVarProdInds <- vector("list", length = ncol(relDf))

#   for (indProd in colnames(relDf)) {
#     labelVarIndProd <- createLabelVar(indProd)

#     labelsFactorLoadings <- vector("character", length = nrow(relDf))
#     labelsVarLatents <- vector("character", length = nrow(relDf))
#     labelsVarInds  <- vector("character", length = nrow(relDf))

#     for (latent in 1:nrow(relDf)) {
#       labelsFactorLoadings[[latent]] <- createLabelLambdaSquared(relDf[latent, indProd],
#                                                           rownames(relDf)[[latent]])
#       labelsVarLatents[[latent]] <- createLabelVar(rownames(relDf)[[latent]])
#       labelsVarInds[[latent]] <- createLabelVar(relDf[latent, indProd])
#     }
#     lhs <- labelVarIndProd
#     rhs1 <- paste(labelsFactorLoadings[[1]],
#                   labelsVarLatents[[1]],
#                   labelsVarInds[[2]], sep = " * ")
#     rhs2 <- paste(labelsFactorLoadings[[2]],
#                   labelsVarLatents[[2]],
#                   labelsVarInds[[1]], sep = " * ")
#     rhs3 <- paste(labelsVarInds[[1]], labelsVarInds[[2]], sep = " * ")

#     rhs <- paste(rhs1, rhs2, rhs3, sep = " + ")

#     constrainedVarProdInds[[indProd]] <- createParTableRow(c(lhs, rhs), op = "==")
#   }
#   constrainedVarProdInds <- purrr::list_rbind(constrainedVarProdInds)

#   rbindParTable(parTable, rbind(varLatentProd,
#                                 covsElemsProd,
#                                 constrainedVarProdInds))
# }


# specifyVarCov <- function(parTable, relDfs) {
#   parTable <- specifyVarCovSingle(parTable, relDfs[[1]])
#
#   if (length(relDfs) <= 1) {
#     return(parTable)
#   } else {
#     return(specifyVarCov(parTable, relDfs[-1]))
#   }
# }
