prepDataModsemDA <- function(data, allIndsXis, allIndsEtas, missing = "listwise",
                             cluster = NULL, sampling.weights = NULL, compress.data = FALSE) {

  if (is.null(data) || !NROW(data))
    return(list(data.full = NULL, n = 0, n.obs = 0, k = 0, p = 0,
                cluster = NULL, idx = NULL))

  if (!is.null(cluster)) {
    mod_stopif(length(cluster) > 1L, "`cluster` must be a single variable!")

    CLUSTER <- as.factor(data[, cluster])

  } else CLUSTER <- NULL

  if (!is.null(sampling.weights)) {
    mod_stopif(length(sampling.weights) > 1L, "`sampling.weights` must be a single variable!")
    WEIGHTS <- data[, sampling.weights]

  } else WEIGHTS <- NULL

  sortData(data, allIndsXis, allIndsEtas) |>
    castDataNumericMatrix() |>
    handleMissingData(missing = missing, CLUSTER = CLUSTER, WEIGHTS = WEIGHTS) |>
    compressData(compress = compress.data) |>
    patternizeMissingDataFIML()
}


indexData <- function(data) {
  attr(data, "idx")
}


# `weights` are stored in pattern-concatenated order (the stacked rows of
# `data.split`), whereas `data.full` keeps the input row order. Anything
# combining the two -- sample moments, N -- needs the weights mapped back.
# With a single missing-data pattern the permutation is the identity.
rowWeightsDA <- function(data) {
  w <- data$weights
  if (is.null(w)) return(NULL)

  perm <- unlist(data$rowidx)
  if (is.null(perm) || length(perm) != length(w)) return(w)

  out <- numeric(length(w))
  out[perm] <- w
  out
}


# Weighted analogue of `cov(X, use = "pairwise.complete.obs")`. Reduces to it
# exactly when every weight is 1: pairwise means, denominator sum(w) - 1.
weightedCovPairwise <- function(X, w) {
  p  <- NCOL(X)
  S  <- matrix(NA_real_, p, p, dimnames = list(colnames(X), colnames(X)))
  ok <- !is.na(X)

  for (i in seq_len(p)) for (j in i:p) {
    use <- ok[, i] & ok[, j]
    W   <- sum(w[use])
    if (W <= 1) next

    mi <- sum(w[use] * X[use, i]) / W
    mj <- sum(w[use] * X[use, j]) / W

    S[i, j] <- S[j, i] <-
      sum(w[use] * (X[use, i] - mi) * (X[use, j] - mj)) / (W - 1)
  }

  S
}


# Map from compressed rows back to input rows: `out[i]` is the compressed row
# that input row `i` was collapsed into. NULL when the data was not compressed,
# which is how callers detect that no expansion is needed. Like the weights,
# `index` is stored in pattern-concatenated order and must be permuted back.
decompressRowsDA <- function(data) {
  index <- data$index
  if (!is.list(index)) return(NULL)

  perm <- unlist(data$rowidx)
  if (!is.null(perm) && length(perm) == length(index)) {
    ordered <- vector("list", length(index))
    ordered[perm] <- index
    index <- ordered
  }

  out <- integer(sum(lengths(index)))
  for (j in seq_along(index)) out[index[[j]]] <- j
  out
}


weightedColMeans <- function(X, w) {
  vapply(seq_len(NCOL(X)), FUN.VALUE = numeric(1L), FUN = function(j) {
    use <- !is.na(X[, j])
    sum(w[use] * X[use, j]) / sum(w[use])
  })
}

sortData <- function(data, allIndsXis, allIndsEtas) {
  inds    <- unique(c(allIndsXis, allIndsEtas))
  ovs     <- colnames(data)
  missing <- inds[!inds %in% ovs]

  mod_stopif(!all(inds %in% ovs), paste0("Missing observed variables in data:\n  ",
         paste(missing, collapse = ", ")))

  data[unique(c(allIndsXis, allIndsEtas))]
}


anyAllNA <- function(data) {
  any(vapply(data, FUN.VALUE = logical(1L), function(x) all(is.na(x))))
}


castDataNumericMatrix <- function(data) {
  force(data) # evaluate to check for errors
  data <- tryCatch({
    numericData <- lapplyDf(data, FUN = as.numeric)
  },
  warning = function(w) {
    mod_msg_warn(paste0("Warning in converting data to numeric: \n", w))
    numericData <- suppressWarnings(lapplyDf(data, FUN = as.numeric))
    mod_stopif(anyAllNA(numericData), "Unable to convert data to type numeric")
    numeric
  },
  error = function(e) {
    mod_msg_stop("Unable to convert data to type numeric")
  })
  as.matrix(data)
}


# A double formatted with "%.17g" round-trips exactly, so rows are merged only
# when they are bit-for-bit identical. Plain `as.character()` would truncate to
# 15 significant digits and could silently merge distinct continuous rows.
compressKey <- function(x) {
  if (is.numeric(x)) sprintf("%.17g", x) else as.character(x)
}


# Collapse duplicate rows into a single row carrying their summed weight.
#
# The likelihood only ever sees the data through row patterns: `m` identical
# rows contribute `m` copies of the same term, so replacing them with one row
# of weight `m` is exact rather than an approximation. On categorical
# indicators the row space is finite and the saving is large (Mplus compresses
# for the same reason, which is also what makes the two directly comparable).
compressData <- function(data, compress = TRUE) {
  if (!compress || is.null(data) || NROW(data) < 2L)
    return(data)

  CLUSTER <- attr(data, "cluster")
  WEIGHTS <- attr(data, "weights")
  INDEX   <- attr(data, "index")

  # Rows may only be merged when they also agree on the cluster they belong
  # to -- collapsing across clusters would destroy the nesting structure the
  # cluster-robust machinery relies on.
  keyed <- as.data.frame(data)
  if (!is.null(CLUSTER)) keyed[["..cluster.."]] <- CLUSTER

  key   <- do.call(paste, c(lapply(keyed, FUN = compressKey), list(sep = "\r")))
  first <- !duplicated(key)
  u     <- sum(first)

  if (u >= NROW(data)) return(data)

  mod_msg_note("Compressing data from", NROW(data), "->", u, "rows..")

  if (is.null(WEIGHTS)) WEIGHTS <- rep(1, NROW(data))
  if (is.null(INDEX))   INDEX   <- seq_len(NROW(data))

  # `match()` against the first occurrences numbers the groups 1..u in order of
  # first appearance, which is exactly the row order of `udata`. Both
  # `rowsum(reorder = TRUE)` and `split()` order their output by that group
  # number, so no further alignment is needed.
  id    <- match(key, key[first])
  udata <- data[first, , drop = FALSE]

  # NOTE: attributes must be set after subsetting, as subsetting a matrix
  # drops any custom attributes.
  attr(udata, "weights") <- as.vector(rowsum(WEIGHTS, group = id, reorder = TRUE))
  attr(udata, "cluster") <- CLUSTER[first]

  # `index` is a plain vector on uncompressed data, but a list of the original
  # row numbers here. Whenever `index` is a list we know the data has been
  # compressed, which is what functions like `modsem_inspect()` and
  # `modsem_predict()` need in order to decompress back to the input rows.
  attr(udata, "index") <- unname(split(INDEX, id))

  # The number of rows actually estimated on is no longer the sample size, so
  # carry the latter separately for anything reporting N (AIC/BIC, nobs(), ...)
  attr(udata, "n.obs") <- attr(data, "n.obs") %||% NROW(data)

  udata
}


patternizeMissingDataFIML <- function(data) {
  # if we are not using fiml, the missing data should already have been removed...
  CLUSTER <- attr(data, "cluster")
  WEIGHTS <- attr(data, "weights")
  INDEX   <- attr(data, "index")
  N.OBS   <- attr(data, "n.obs")

  Y   <- as.matrix(data)
  obs <- !is.na(Y)

  rowMissingAll <- apply(obs, MARGIN = 1, FUN = \(x) !any(x))
  colMissingAll <- apply(obs, MARGIN = 2, FUN = \(x) !any(x))
  mod_stopif(any(colMissingAll),
         paste0("Please remove variables with only missing values:\n  ",
         paste0(colnames(obs)[colMissingAll], collapse = ", ")))

  patterns <- unique(obs, MARING = 2)

  if (any(rowMissingAll)) { # remove patterns where all are missing
    # This shouldn't really happen, as it should be handled already in
    # `handleMissingData()`. Regardless, we should handle it if it ever happens
    kept <- data[!rowMissingAll, , drop = FALSE]

    # subsetting a matrix drops custom attributes, so re-attach them
    attr(kept, "cluster") <- CLUSTER[!rowMissingAll]
    attr(kept, "weights") <- WEIGHTS[!rowMissingAll]
    attr(kept, "index")   <- INDEX[!rowMissingAll]
    attr(kept, "n.obs")   <- N.OBS

    return(patternizeMissingDataFIML(kept))
  }

  ids <- seq_len(NROW(patterns))
  n   <- NROW(Y)
  k   <- NCOL(Y)

  rowidx <- vector("list", length = NROW(ids))
  colidx <- vector("list", length = NROW(ids))
  data.split  <- vector("list", length = NROW(ids))
  n.pattern  <- numeric(NROW(ids))
  d.pattern  <- numeric(NROW(ids))

  for (id in ids) {
    mask  <- matrix(patterns[id, ], nrow = n, ncol = k, byrow = TRUE)
    match <- apply(obs == mask, MARGIN = 1, FUN = all)
    ridx  <- which(match)
    cidx  <- which(patterns[id, ])

    rowidx[[id]] <- ridx
    colidx[[id]] <- cidx
    data.split[[id]] <- Y[ridx, cidx, drop = FALSE]
    n.pattern[[id]] <- sum(match)
    d.pattern[[id]] <- length(cidx)
  }

  # Everything downstream (posterior matrices, densities, score contributions)
  # is in pattern-concatenated order (i.e., the stacked rows of `data.split`),
  # so the row-level vectors must be permuted accordingly. With a single
  # pattern, `perm` is the identity.
  permutation <- unlist(rowidx)

  list(
    ids        = ids,
    rowidx     = rowidx,
    colidx     = colidx,
    colidx0    = lapply(colidx, FUN = \(idx) idx - 1),
    patterns   = patterns,
    data.split = data.split,
    n.pattern  = n.pattern,
    d.pattern  = d.pattern,
    n          = NROW(data),

    # `n` counts the rows actually estimated on; after `compressData()` those
    # are distinct row patterns rather than observations, so anything reporting
    # a sample size (nobs(), AIC/BIC, ...) must use `n.obs` instead.
    n.obs      = N.OBS %||% NROW(data),
    k          = NCOL(data),
    p          = length(ids),
    data.full  = data,
    is.fiml    = length(ids) > 1L,
    cluster    = CLUSTER[permutation],
    weights    = WEIGHTS[permutation],
    index      = INDEX[permutation]
  )
}


handleMissingData <- function(data, missing = "listwise", CLUSTER = NULL, WEIGHTS = NULL) {
  missing       <- tolower(missing)
  completeCases <- stats::complete.cases(data)
  anyMissing    <- any(!completeCases)
  allMissing    <- all(!completeCases)
  INDEX         <- seq_len(NROW(data))

  if (!anyMissing){
    attr(data, "cluster") <- CLUSTER
    attr(data, "weights") <- WEIGHTS
    attr(data, "index")   <- INDEX
    return(data)

  } else if (allMissing) {
    missingAllCol <- apply(data, MARGIN = 2, FUN = \(x) all(is.na(x)))
    colsMissing   <- colnames(data)[missingAllCol]

    mod_msg_stop(paste0("Please remove variables with only missing values:\n  ",
          paste0(colsMissing, collapse = ", ")))
  }

  if (missing %in% c("listwise", "casewise", "complete")) {
    mod_msg_warn(paste0("Removing missing values list-wise!\n",
             "Consider using `missing=\"fiml\"`, `missing=\"impute\"`, ",
             "or the `modsem_mimpute()` function!\n"))

    out <- data[completeCases, ]
    attr(out, "cluster") <- CLUSTER[completeCases]
    attr(out, "weights") <- WEIGHTS[completeCases]
    attr(out, "index")   <- INDEX[completeCases]

    return(out)

  } else if (missing == "impute") {
    mod_msg_note("Imputing missing values. Consider using the `modsem_mimpute()` function!")

    imp  <- Amelia::amelia(data, m = 1, p2s = 0)

    imp1 <- as.matrix(as.data.frame(imp$imputations[[1]]))
    attr(imp1, "cluster") <- CLUSTER
    attr(imp1, "weights") <- WEIGHTS
    attr(imp1, "index")   <- INDEX

    return(imp1)

  } else if (missing %in% c("fiml", "ml", "direct")) {
    rowMissingAll <- apply(data, MARGIN = 1, FUN = \(x) all(is.na(x)))
    keep          <- !rowMissingAll

    data <- data[keep, , drop = FALSE] # we already know that !all(rowMissingAll)

    # NOTE: attributes must be set after subsetting, as subsetting a
    # matrix drops any custom attributes!
    attr(data, "cluster") <- CLUSTER[keep]
    attr(data, "weights") <- WEIGHTS[keep]
    attr(data, "index")   <- INDEX[keep]

    return(data)

  } else {
    mod_msg_stop(sprintf("Unrecognized value for `missing`: `%s`", missing))
  }
}
