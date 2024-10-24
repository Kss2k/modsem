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
}
