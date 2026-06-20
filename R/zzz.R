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
