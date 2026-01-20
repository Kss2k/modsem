modsemMatrix <- function(mat, symmetric = isSymmetric(mat), is.public = FALSE) {
  if (is.null(mat)) return(mat)

  if (is.public) {
    rn <- rownames(mat)
    cn <- colnames(mat)

    if (!is.null(rn)) {
      isTempRn <- startsWith(rn, prefix = TEMP_OV_PREFIX)
      rnClean <- stringr::str_remove_all(rn, pattern = TEMP_OV_PREFIX)
      isDupTempRn <- rnClean[isTempRn] %in% rnClean[!isTempRn]

      keepRn <- rep(TRUE, length(rn))
      keepRn[isTempRn][isDupTempRn] <- FALSE
      mat <- mat[keepRn, , drop = FALSE]

      rownames(mat) <- stringr::str_replace_all(
        string = rnClean[keepRn],
        pattern = OP_REPLACEMENTS_INV
      )
    }

    if (!is.null(cn)) {
      isTempCn <- startsWith(cn, prefix = TEMP_OV_PREFIX)
      cnClean <- stringr::str_remove_all(cn, pattern = TEMP_OV_PREFIX)
      isDupTempCn <- cnClean[isTempCn] %in% cnClean[!isTempCn]

      keepCn <- rep(TRUE, length(cn))
      keepCn[isTempCn][isDupTempCn] <- FALSE
      mat <- mat[ , keepCn, drop = FALSE]

      colnames(mat) <- stringr::str_replace_all(
        string = cnClean[keepCn],
        pattern = OP_REPLACEMENTS_INV
      )
    }
  }

  class(mat) <- unique(c("ModsemMatrix", class(mat)))

  if (symmetric)
    class(mat) <- unique(c("ModsemSymmetricMatrix", class(mat)))

  mat
}


#' @export
print.ModsemSymmetricMatrix <- function(x, digits = 3, sep = " ", ...) {
  if (nrow(x) == ncol(x)) x[upper.tri(x)] <- NA

  print.ModsemMatrix(x, digits = digits, sep = sep, ...)
}


#' @export
print.ModsemMatrix <- function(x, digits = 3, shift = 0L, ...) {
  y <- matrix(formatNumeric(x, digits = digits), nrow = nrow(x),
              ncol = ncol(x), dimnames = dimnames(x))

  y <- unclass(y)
  # Remove NAs
  y[is.na(x)] <- ""

  if (!is.null(colnames(x))) {
    colnames(y) <- abbreviate(colnames(x), minlength = digits + 3L) }
  if (shift > 0L) {
    empty.string <- rep(strrep(x = " ", times = shift), times = nrow(x))

    if (!is.null(rownames(x))) rownames(y) <- paste0(empty.string, rownames(x))
    else                       rownames(y) <- empty.string
  }

  print(y, ..., quote = FALSE, right = TRUE)

  invisible(x)
}


modsemParTable <- function(parTable, is.public = FALSE) {
  if (is.null(parTable)) return(parTable)

  if (is.public) {
    paramCols <- c("lhs", "op", "rhs", "label")
    paramCols <- intersect(paramCols, colnames(parTable))

    for (col in paramCols) {
      parTable[[col]] <- stringr::str_replace_all(
        string = parTable[[col]],
        pattern = OP_REPLACEMENTS_INV
      )
    }
  }

  class(parTable) <- unique(c("ModsemParTable", class(parTable)))
  parTable
}


modsemVector <- function(vec, is.public = FALSE) {
  if (is.null(vec)) return(vec)

  if (is.public) {
    nm <- names(vec)
    isTemp <- startsWith(nm, prefix = TEMP_OV_PREFIX)

    if (any(isTemp)) {
      clean <- stringr::str_remove_all(nm, pattern = TEMP_OV_PREFIX)

      isDupTemp <- clean[isTemp] %in% clean[!isTemp]

      keep <- rep(TRUE, length(nm))
      keep[isTemp][isDupTemp] <- FALSE

      vec <- stats::setNames(vec[keep], nm = clean[keep])
    }

    names(vec) <- stringr::str_replace_all(
      string = names(vec),
      pattern = OP_REPLACEMENTS_INV
    )
  }

  class(vec) <- unique(c("ModsemVector", class(vec)))

  vec
}


#' @export
print.ModsemVector <- function(x, digits = 3L, ...) {
  y <- formatNumeric(x, digits = digits)
  print(y, quote = FALSE)
}


#' @export
print.ModsemParTable <- function(x, nd = 3L, ...) {
  row.names <- rownames(x)
  y <- lapply(x, \(x) if (is.numeric(x)) round(x, nd) else x) |>
    as.data.frame()

  rownames(y) <- row.names

  print(y, ...)

  invisible(x)
}
