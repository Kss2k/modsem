PKG_INFO <- rlang::env(version = NULL)


getPackageVersion <- function(pkgname) {
  tryCatch({
    read.dcf(file = system.file("DESCRIPTION", package = pkgname),
             fields = "Version")
  }, error = function(e) {
    warning2("Failed to get package version")
    "??" # replace this with a hard-coded value?
  })
}


.onLoad <- function(libname, pkgname) {
  PKG_INFO$version <- getPackageVersion(pkgname)
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
