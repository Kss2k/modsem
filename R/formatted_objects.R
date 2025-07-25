modsemMatrix <- function(mat, symmetric = isSymmetric(mat)) {
  if (is.null(mat)) return(mat)

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


modsemParTable <- function(parTable) {
  if (is.null(parTable)) return(parTable)

  class(parTable) <- unique(c("ModsemParTable", class(parTable)))
  parTable
}


modsemVector <- function(vec) {
  if (is.null(vec)) return(vec)
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
