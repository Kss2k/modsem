MODSEM_COLORS <- rlang::env(
  numeric.positive = "orange1",
  numeric.negative = "tomato",
  active = FALSE,
  f.numeric.positive = \(x) x,
  f.numeric.negative = \(x) x,
  available = grDevices::colors()
)


NUMERIC.P <- "(?<![A-Za-z0-9._-])([0-9]+(?:[.-][0-9]+)*(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)(?![A-Za-z])"
NUMERIC.N <- "(?<![A-Za-z0-9._-])(-[0-9]+(?:[.-][0-9]+)*(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)(?![A-Za-z])"


#' Set Color Scheme for Numeric Highlighting in \code{\link{modsem}}
#'
#' Sets the ANSI color styles for positive and negative numeric values to be used 
#' in text output functions in \code{\link{modsem}}. Also activates or deactivates color usage.
#'
#' @param numeric.positive A color name (as accepted by the \code{cli} package or base \code{R}) 
#'   to be used for positive numeric values. Default is \code{"orange1"}.
#' @param numeric.negative A color name for negative numeric values. Default is \code{"orange1"}.
#' @param active Logical; if \code{TRUE}, color highlighting will be activated. If \code{FALSE}, 
#'   any color-related behavior is disabled, regardless of color validity.
#'
#' @details 
#' This function sets both the string names of the color (\code{numeric.positive}, \code{numeric.negative})
#' and their corresponding compiled ANSI styles via \code{cli::make_ansi_style()}, storing them in
#' the `MODSEM_COLORS` environment.
#'
#' If either color is not a valid \code{R} color (i.e., not in \code{grDevices::colors()}), or if \code{active} is 
#' \code{FALSE}, then highlighting is disabled and the function returns \code{FALSE}.
#'
#' @return Logical: \code{TRUE} if colors were successfully set and activated, \code{FALSE} if deactivated.
#'
#' @examples
#' m1 <- "
#' # Outer Model
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#' # Inner Model
#'   Y ~ X + Z + X:Z
#' "
#'
#' est <- modsem(m1, data = oneInt)
#' 
#' # Summary with default colors
#' print(summary(est))
#'
#' # Change colors
#' set_modsem_colors(numeric.positive = "green", numeric.negative = "red", active = TRUE)
#' print(summary(est))
#'
#' # Disable colors
#' set_modsem_colors(active = FALSE)
#' print(summary(est))
set_modsem_colors <- function(numeric.positive = "orange1",
                              numeric.negative = "orange1",
                              active = TRUE) {
  MODSEM_COLORS$available <- grDevices::colors()

  colors <- c(numeric.positive, numeric.negative)
  
  if (!active || any(!colors %in% MODSEM_COLORS$available)) {
    MODSEM_COLORS$active <- FALSE
    return(FALSE)
  }

  MODSEM_COLORS$active <- TRUE
  MODSEM_COLORS$numeric.positive <- numeric.positive
  MODSEM_COLORS$numeric.negative <- numeric.negative
  MODSEM_COLORS$f.numeric.positive <- cli::make_ansi_style(numeric.positive)
  MODSEM_COLORS$f.numeric.negative <- cli::make_ansi_style(numeric.negative)

  TRUE
}
set_modsem_colors <- function(numeric.positive = "orange1",
                              numeric.negative = "orange1",
                              active = TRUE) {
  colors <- c(numeric.positive, numeric.negative)
  if (!active || any(!colors %in% MODSEM_COLORS$available)) {
    MODSEM_COLORS$active <- FALSE
    return(FALSE)
  }

  MODSEM_COLORS$active <- TRUE
  MODSEM_COLORS$numeric.positive <- numeric.positive
  MODSEM_COLORS$numeric.negative <- numeric.negative
  MODSEM_COLORS$f.numeric.positive <- cli::make_ansi_style(numeric.positive)
  MODSEM_COLORS$f.numeric.negative <- cli::make_ansi_style(numeric.negative)

  TRUE
}


getColorizedASCII <- function(color) {
  as.character(cli::make_ansi_style(color)("\\0"))
}


colorize <- function(expr, 
                     numeric.positive = MODSEM_COLORS$numeric.positive, 
                     numeric.negative = MODSEM_COLORS$numeric.negative,
                     append = "\n") {
  if (!MODSEM_COLORS$active) return(expr)

  output <- stringr::str_c(utils::capture.output(expr), collapse = "\n")

  rep.numeric.p <- getColorizedASCII(numeric.positive)
  rep.numeric.n <- getColorizedASCII(numeric.negative)

  mapping <- c(rep.numeric.p, rep.numeric.n)
  names(mapping) <- c(NUMERIC.P, NUMERIC.N)

  colorized <- tryCatch(
    stringr::str_replace_all(string  = output, pattern = mapping),
    error = function(e) {
      warning2("Colorization of output failed! Message:", e,
               immediate. = FALSE)
      output
    }
  )

  cat(colorized, append)
}


colorFormatNum <- function(x, f = "%f") {
  str <- sprintf(f, x)
  if (sign(x) == 1) 
    MODSEM_COLORS$f.numeric.positive(str)
  else 
    MODSEM_COLORS$f.numeric.negative(str)
}
