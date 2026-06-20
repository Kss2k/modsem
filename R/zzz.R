PKG_INFO <- rlang::env(version = NULL)


UNICODE_MSG_STRINGS <- list(
  updateStatusLog0 = "\rIter=%d Mode=%s LogLik=%.2f \u0394LL=%.2g rel\u0394LL=%.2g",
  incrementIterations0 = "\rEval=%d LogLik=%.2f \u0394LL=%.2g rel\u0394LL=%.2g"
)


ASCII_MSG_STRINGS <- list(
  updateStatusLog0 = "\rIter=%d Mode=%s LogLik=%.2f dLL=%.2g reldLL=%.2g",
  incrementIterations0 = "\rEval=%d LogLik=%.2f dLL=%.2g reldLL=%.2g"
)


MSG_STRINGS <- rlang::env(
  unicode = FALSE,
  strings = ASCII_MSG_STRINGS
)


# Box-drawing characters for table output (see simple_slopes.R)
UNICODE_BOX_CHARS <- list(
  V_LINE    = "\u2502", # vertical line
  H_LINE    = "\u2500", # horizontal line
  D_CROSS   = "\u252c", # down cross
  LU_CORNER = "\u250c", # left upper corner
  LL_CORNER = "\u2514", # left lower corner
  RU_CORNER = "\u2510", # right upper corner
  RL_CORNER = "\u2518", # right lower corner
  H_DLINE   = "\u2550", # horizontal double line
  F_DCROSS  = "\u256a", # full double cross
  R_DCROSS  = "\u2561", # right double cross
  L_DCROSS  = "\u255e", # left double cross
  U_CROSS   = "\u2534"  # up cross
)


ASCII_BOX_CHARS <- list(
  V_LINE    = "|",
  H_LINE    = "-",
  D_CROSS   = "+",
  LU_CORNER = "+",
  LL_CORNER = "+",
  RU_CORNER = "+",
  RL_CORNER = "+",
  H_DLINE   = "=",
  F_DCROSS  = "+",
  R_DCROSS  = "+",
  L_DCROSS  = "+",
  U_CROSS   = "+"
)


BOX_CHARS <- rlang::env(chars = ASCII_BOX_CHARS)


getPackageVersion <- function(pkgname) {
  tryCatch({
    c(read.dcf(
      file = system.file("DESCRIPTION", package = pkgname),
      fields = "Version"
    ))
  }, error = function(e) {
    mod_msg_warn("Failed to get package version")
    "??" # replace this with a hard-coded value?
  })
}


.onLoad <- function(libname, pkgname) {
  PKG_INFO$version <- getPackageVersion(pkgname)

  # Are unicode characters supported?
  utf8 <- tryCatch(isTRUE(l10n_info()[["UTF-8"]]), error = \(e) FALSE)
  MSG_STRINGS$unicode <- utf8

  if (utf8) MSG_STRINGS$strings <- UNICODE_MSG_STRINGS
  else      MSG_STRINGS$strings <- ASCII_MSG_STRINGS

  if (utf8) BOX_CHARS$chars <- UNICODE_BOX_CHARS
  else      BOX_CHARS$chars <- ASCII_BOX_CHARS

  resetModsemColors()
}


.onAttach <- function(libname, pkgname) {
  version <- getPackageVersion(pkgname)
  message <- sprintf("This is %s (%s). Please report any bugs!", pkgname, version)

  openMP_Enabled <- checkOpenMP_Cpp()

  if (!openMP_Enabled) {
    message <- paste0(message, "\n",
                     "OpenMP is not available! Multi-threading will not work properly!")
  }

  packageStartupMessage(message)
}
