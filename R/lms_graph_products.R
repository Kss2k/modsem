# Product terms for the recursive (graph) LMS backend.
#
# The legacy backend stores interactions in `omegaXiXi` / `omegaEtaXi`, which
# are stacked (numEtas * numXis) blocks laid out for Kronecker products. That
# layout can only express two-way terms, and it structurally forbids products
# between endogenous variables.
#
# The graph backend instead describes products with two matrices:
#
#   productDesign  P by k, integer exponents. Entry (p, v) is how many times
#                  latent variable v appears in product p, so `X:X` is a 2 and
#                  `X:Z:W` is three 1s. Purely structural -- never a free
#                  parameter, and deliberately not registered in DA_BLOCKS.
#
#   omega          numEtas by P coefficients. `NA` is free, a number is fixed,
#                  and 0 means the product does not enter that equation. This
#                  is the same NA/number/label convention every other parameter
#                  block uses, so it needs no special handling downstream.
#
# Column order of `productDesign` matches the graph state vector: xis first,
# then etas.


# Products are unordered, so `X:Z` and `Z:X` must map to one column. Sorting
# the factors by their position in the latent vector gives a canonical name;
# repeated factors survive the sort, which is what makes `X:X` an exponent.
lmsGraphProductName <- function(factors, latent) {
  index <- match(factors, latent)
  mod_stopif(anyNA(index), sprintf(
    paste0("The LMS graph backend only supports products of latent variables. ",
           "`%s` is not a latent variable in this model."),
    paste(unique(factors[is.na(index)]), collapse = "`, `")))
  paste(latent[sort(index)], collapse = ":")
}


lmsGraphProductFactors <- function(product) stringr::str_split_1(product, ":")


lmsGraphProductMatrices <- function(parTable, xis, etas) {
  latent <- c(xis, etas)
  intTerms <- getIntTermRows(parTable)
  intTerms <- intTerms[intTerms$lhs %in% etas, , drop = FALSE]

  if (!NROW(intTerms)) return(list(
    productDesign = matrix(0L, 0L, length(latent),
                           dimnames = list(NULL, latent)),
    omega = matrix(numeric(), length(etas), 0L, dimnames = list(etas, NULL)),
    omegaLabel = matrix(character(), length(etas), 0L,
                        dimnames = list(etas, NULL))
  ))

  factors <- lapply(intTerms$rhs, lmsGraphProductFactors)
  canonical <- vapply(factors, lmsGraphProductName, character(1L),
                      latent = latent)
  products <- unique(canonical)

  design <- matrix(0L, length(products), length(latent),
                   dimnames = list(products, latent))
  for (p in seq_along(products)) {
    counts <- table(lmsGraphProductFactors(products[[p]]))
    design[p, names(counts)] <- as.integer(counts)
  }

  omega <- matrix(0, length(etas), length(products),
                  dimnames = list(etas, products))
  label <- matrix("", length(etas), length(products),
                  dimnames = list(etas, products))

  for (i in seq_len(NROW(intTerms))) {
    eta <- intTerms$lhs[[i]]
    product <- canonical[[i]]
    modifier <- intTerms$mod[[i]]
    # A label leaves the numeric entry at 0 and records the label instead,
    # matching how `getFreeOrConstIntTerms()` treats the legacy omega blocks.
    if (canBeNumeric(modifier, includeNA = TRUE))
      omega[eta, product] <- as.numeric(modifier)
    else label[eta, product] <- modifier
  }

  lmsGraphCheckProductOrder(design, omega, label, xis, etas)
  list(productDesign = design, omega = omega, omegaLabel = label)
}


# Latent variables are evaluated in structural order, so a product may only
# reference etas that are already available when its equation is evaluated. The
# Kronecker layout made this impossible to violate; the product representation
# does not, so it has to be checked. `omega` rows are etas in topological
# order, so a product entering equation j may only contain etas before j.
lmsGraphCheckProductOrder <- function(design, omega, omegaLabel, xis, etas) {
  if (!NCOL(omega) || !length(etas) || !NROW(design))
    return(invisible(TRUE))

  etaColumns <- length(xis) + seq_along(etas)
  entered <- is.na(omega) | omega != 0 | omegaLabel != ""

  for (j in seq_along(etas)) for (p in which(entered[j, ])) {
    used <- which(design[p, etaColumns, drop = TRUE] > 0L)
    if (!length(used)) next
    mod_stopif(any(used >= j), sprintf(
      paste0("`%s` uses the product `%s`, which contains an endogenous ",
             "variable that is not evaluated before it. Products may only ",
             "reference endogenous variables preceding the equation they ",
             "enter."),
      etas[[j]], rownames(design)[[p]]))
  }
  invisible(TRUE)
}


# Swap a DA model over to the product representation. The Kronecker blocks are
# replaced rather than supplemented: they are left zero-dimensional, exactly as
# `lambdaY` and `thetaEpsilon` are for the LMS approach, so every loop- and
# vector-based routine skips them without needing to know which backend is in
# use. Mirrors `lmsGraphPrepareOrdered()`.
lmsGraphPrepareProducts <- function(model) {
  old.theta <- model$theta

  for (g in seq_along(model$models)) {
    submodel <- model$models[[g]]
    products <- lmsGraphProductMatrices(
      model$parTable, submodel$info$xis, submodel$info$etas
    )

    submodel$matrices$productDesign      <- products$productDesign
    submodel$matrices$omega              <- products$omega
    submodel$labelMatrices$omega         <- products$omegaLabel
    submodel$matrices$omegaXiXi          <- matrix(numeric(), 0L, 0L)
    submodel$matrices$omegaEtaXi         <- matrix(numeric(), 0L, 0L)
    submodel$labelMatrices$omegaXiXi     <- matrix(character(), 0L, 0L)
    submodel$labelMatrices$omegaEtaXi    <- matrix(character(), 0L, 0L)

    model$models[[g]] <- submodel
  }

  params <- createTheta(model, parTable.in = model$parTable)
  common <- intersect(names(old.theta), names(params$theta))
  params$theta[common] <- old.theta[common]
  model$params[names(params)] <- params
  model$theta <- params$theta
  model$params$bounds <- getParamBounds(model)
  model$params$gradientStruct <- getGradientStruct(model, model$theta,
                                                   method = "lms")
  model
}
