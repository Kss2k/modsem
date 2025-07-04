MODSEM_COLORS <- rlang::env(
  ## user-configurable colour names
  positive = NA_character_,
  negative = NA_character_,
  true     = NA_character_,
  false    = NA_character_,
  nan      = NA_character_,
  na       = NA_character_,
  inf      = NA_character_,
  string   = NA_character_,

  ## run-time flags & helpers
  active   = FALSE,
  available = NULL,

  ## pre-compiled styling functions (identity placeholders)
  f.positive = \(x) x,
  f.negative = \(x) x,
  f.true     = \(x) x,
  f.false    = \(x) x,
  f.nan      = \(x) x,
  f.na       = \(x) x,
  f.inf      = \(x) x,
  f.string   = \(x) x
)

PATTERN.pos    <- "(?<![-+A-Za-z0-9._~])([0-9]+(?:[.-][0-9]+)*(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)(?![A-Za-z])"
PATTERN.neg    <- "(?<![-+A-Za-z0-9._~])(-[0-9]+(?:[.-][0-9]+)*(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)(?![A-Za-z])"
PATTERN.TRUE   <- "(?<![A-Za-z0-9])TRUE(?![A-Za-z0-9])"
PATTERN.FALSE  <- "(?<![A-Za-z0-9])FALSE(?![A-Za-z0-9])"
PATTERN.NaN    <- "(?<![A-Za-z0-9])NaN(?![A-Za-z0-9])"
PATTERN.NA     <- "(?<![A-Za-z0-9])NA(?![A-Za-z0-9])"
PATTERN.Inf    <- "(?<![A-Za-z0-9])-?Inf(?![A-Za-z0-9])"
PATTERN.String <- "(\"(?:[^\"\\\\]|\\\\.)*\"|'(?:[^'\\\\]|\\\\.)*')"


setAvailableColors <- function() {
  MODSEM_COLORS$available <- grDevices::colors()
}


styleOrAsIs <- function(col) {
  if (is.na(col)) return(\(x) x)
  cli::make_ansi_style(col)
}


getColorizedASCII <- function(col) {
  if (is.na(col)) return("\\0")            # no colour → leave as-is
  as.character(cli::make_ansi_style(col)("\\0"))
}


resetModsemColors <- function() {
  MODSEM_COLORS$active <- FALSE
  fields <- c("positive","negative","true","false","nan","na","inf","string")
  for (f in fields) MODSEM_COLORS[[f]] <- NA_character_

  funs <- grep("^f\\.", names(MODSEM_COLORS), value = TRUE)
  for (f in funs)  MODSEM_COLORS[[f]] <- \(x) x

  setAvailableColors()
}


#' Define or disable the colour theme used by \code{modsem}
#'
#' All arguments are optional; omitted ones fall back to the defaults below.
#' Pass \code{active = FALSE} to turn highlighting off (and reset the palette).
#'
#' @param expr Expression or object with output which should be colorized.
#' @param positive Colour of positive numbers.
#' @param negative Colour of negative numbers.
#' @param true Colour of \code{TRUE}.
#' @param false Colour of \code{FALSE}.
#' @param nan Colour of \code{NaN}.
#' @param na Colour of \code{NA}.
#' @param inf Colour of \code{-Inf} and \code{Inf}.
#' @param string Colour of quoted strings.
#' @param append String appended after the coloured output (default `\\n`).
#'
#' @examples
#' 
#' set_modsem_colors(positive = "green", 
#'                   negative = "red",
#'                   true = "green", 
#'                   false = "red", 
#'                   na = "purple")
#' 
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
#' colorize_output(summary(est))
#' colorize_output(print(est))
#' 
#' \dontrun{
#' est_lms <- modsem(m1, data = oneInt, method = "lms")
#'
#' colorize_output(summary(est_lms))
#' 
#' colorize_output(modsem_inspect(est_lms))
#' }
#'
#' @return `TRUE` if colours are active afterwards, otherwise `FALSE`.
#' @export
set_modsem_colors <- function(positive = "green3",
                              negative = positive,
                              true     = "green3",
                              false    = "tomato",
                              nan      = "tomato",
                              na       = "yellow3",
                              inf      = "magenta",
                              string   = "cyan",
                              active   = TRUE) {
  if (is.null(MODSEM_COLORS$available)) setAvailableColors()

  cols <- c(positive, negative, true, false, nan, na, inf, string)
  if (!active || any(!cols %in% MODSEM_COLORS$available)) {
    reetModsemColors()
    return(FALSE)
  }

  MODSEM_COLORS$positive <- positive
  MODSEM_COLORS$negative <- negative
  MODSEM_COLORS$true     <- true
  MODSEM_COLORS$false    <- false
  MODSEM_COLORS$nan      <- nan
  MODSEM_COLORS$na       <- na
  MODSEM_COLORS$inf      <- inf
  MODSEM_COLORS$string   <- string

  MODSEM_COLORS$f.positive <- styleOrAsIs(positive)
  MODSEM_COLORS$f.negative <- styleOrAsIs(negative)
  MODSEM_COLORS$f.true     <- styleOrAsIs(true)
  MODSEM_COLORS$f.false    <- styleOrAsIs(false)
  MODSEM_COLORS$f.nan      <- styleOrAsIs(nan)
  MODSEM_COLORS$f.na       <- styleOrAsIs(na)
  MODSEM_COLORS$f.inf      <- styleOrAsIs(inf)
  MODSEM_COLORS$f.string   <- styleOrAsIs(string)

  MODSEM_COLORS$active <- TRUE

  MODSEM_COLORS$active
}

#' Capture, colourise, and emit console text
#' @param expr Expression or object with output which should be colorized.
#' @param positive Colour of positive numbers.
#' @param negative Colour of negative numbers.
#' @param true Colour of \code{TRUE}.
#' @param false Colour of \code{FALSE}.
#' @param nan Colour of \code{NaN}.
#' @param na Colour of \code{NA}.
#' @param inf Colour of \code{-Inf} and \code{Inf}.
#' @param string Colour of quoted strings.
#' @param append String appended after the coloured output (default `\\n`).
#' @return Invisibly returns the *plain* captured text.
#'
#' @examples
#' 
#' set_modsem_colors(positive = "green", 
#'                   negative = "red",
#'                   true = "green", 
#'                   false = "red", 
#'                   na = "purple")
#' 
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
#' colorize_output(summary(est))
#' colorize_output(print(est))
#' 
#' \dontrun{
#' est_lms <- modsem(m1, data = oneInt, method = "lms")
#'
#' colorize_output(summary(est_lms))
#' 
#' colorize_output(modsem_inspect(est_lms))
#' }
#' @export
colorize_output <- function(expr,
                            positive = MODSEM_COLORS$positive,
                            negative = MODSEM_COLORS$negative,
                            true     = MODSEM_COLORS$true,
                            false    = MODSEM_COLORS$false,
                            nan      = MODSEM_COLORS$nan,
                            na       = MODSEM_COLORS$na,
                            inf      = MODSEM_COLORS$inf,
                            string   = MODSEM_COLORS$string,
                            append = "\n") {
  out <- stringr::str_c(utils::capture.output(expr), collapse = "\n")

  if (!MODSEM_COLORS$active) {
    cat(out, append)
    return(invisible())
  }

  ## build mapping only for active colours --------------------------------
  mapping <- c()
  add <- function(pat, col) {
    if (!is.na(col)) mapping <<- c(mapping, setNames(getColorizedASCII(col), pat))
  }

  add(PATTERN.pos,    positive)
  add(PATTERN.neg,    negative)
  add(PATTERN.TRUE,   true)
  add(PATTERN.FALSE,  false)
  add(PATTERN.NaN,    nan)
  add(PATTERN.NA,     na)
  add(PATTERN.Inf,    inf)
  add(PATTERN.String, string)

  coloured <- tryCatch(
    stringr::str_replace_all(out, mapping),
    error = function(e) {
      warning("Colourization failed: ", conditionMessage(e), call. = FALSE)
      out
    }
  )

  cat(coloured, append)
  invisible(out)
}
