prepDataModsemDA <- function(data, allIndsXis, allIndsEtas, missing = "listwise",
                             cluster = NULL, sampling.weights = NULL) {

  if (is.null(data) || !NROW(data))
    return(list(data.full = NULL, n = 0, k = 0, p = 0, cluster = NULL, idx = NULL))

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
    patternizeMissingDataFIML()
}


indexData <- function(data) {
  attr(data, "idx") 
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


patternizeMissingDataFIML <- function(data) {
  # if we are not using fiml, the missing data should already have been removed...
  CLUSTER <- attr(data, "cluster")
  WEIGHTS <- attr(data, "weights")
  INDEX   <- attr(data, "index")

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
    return(patternizeMissingDataFIML(data[!rowMissingAll, , drop = FALSE]))
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
    k          = NCOL(data),
    p          = length(ids),
    data.full  = data,
    is.fiml    = length(ids) > 1L,
    cluster    = CLUSTER,
    weights    = WEIGHTS,
    index      = INDEX
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
    attr(data, "index")  <- INDEX[completeCases]

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
    attr(data, "cluster") <- CLUSTER
    attr(data, "weights") <- WEIGHTS
    attr(data, "index")   <- INDEX

    rowMissingAll <- apply(data, MARGIN = 1, FUN = \(x) all(is.na(x)))
    data          <- data[!rowMissingAll, , drop = FALSE] # we've already know that
                                                          # all(rowMissingAll) != TRUE
    return(data)

  } else {
    mod_msg_stop(sprintf("Unrecognized value for `missing`: `%s`", missing))
  }
}
