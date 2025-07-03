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
print.ModsemMatrix <- function(x, digits = 3, sep = " ", ...) {
  body <- matrix(formatNumeric(x, digits = digits), 
                 nrow = nrow(x), ncol = ncol(x))
  # Remove NAs
  body[is.na(x)] <- ""

  header <- colnames(x)
  index  <- format(c("", rownames(x)), justify = "right")

  X <- rbind(header, body)
  X <- matrix(apply(X, MARGIN = 2, FUN = format, justify = "right"), 
              nrow = nrow(X), ncol = ncol(X))
  
  if (MODSEM_COLORS$active) {
    neg <- !is.na(x) & x < 0
    pos <- !is.na(x) & !neg

    f.pos <- MODSEM_COLORS$f.numeric.positive
    f.neg <- MODSEM_COLORS$f.numeric.negative

    X[-1, ][pos] <- as.character(f.pos(X[-1, ][pos]))
    X[-1, ][neg] <- as.character(f.neg(X[-1, ][neg]))
  }

  width.out <- options("width")$width
  if (is.null(width.out) || is.na(width.out)) width.out <- 80L

  buffer   <- index
  i        <- 1L
  k        <- ncol(X)
  iter     <- 0L # If something goes wrong, we use iter to stop us from ending up in an endless loop
  max.iter <- ncol(X) # Worst case scenario, we print one column per outer while-loop
  width.used <- 0L

  while (i <= k &&iter <= max.iter) {
    iter <- iter + 1L
    j    <- 0L

    while (i <= k && width.used < width.out) {
      sep_i       <- if (max(nchar(buffer)) > 0) " " else ""
      buffer.next <- paste0(buffer, sep_i, X[, i])
      width.next  <- max(nchar(buffer.next)) # should all be same, but take max to be certain...

      if (j <= 0L || width.next < width.out) {
        buffer     <- buffer.next
        width.used <- width.next
      } else break

      i <- i + 1L
      j <- j + 1L
    }

    for (line in buffer) cat(line, "\n")

    buffer <- index
  }
}


modsemVector <- function(vec) {
  if (is.null(vec)) return(vec)
  class(vec) <- unique(c("ModsemVector", class(vec)))

  vec
}


#' @export
print.ModsemVector <- function(x, ...) {
  y <- matrix(x, nrow = 1, ncol = length(x),
              dimnames = list("", names(x)))
  
  print.ModsemMatrix(y, ...)
}


modsemParTable <- function(parTable) {
  if (is.null(parTable)) return(parTable)

  class(parTable) <- unique(c("ModsemParTable", class(parTable)))
  parTable
}


#' @export
print.ModsemParTable <- function(x, nd = 3L, ...) {
    row.names <- rownames(x)
    y <- lapply(x, \(x) if (is.numeric(x)) round(x, nd) else x) |>
      as.data.frame()

    rownames(y) <- row.names

    colorize(print(y, ...))

    invisible(x)
}
