getMethodSettingsPI <- function(method, args) {
  defaultResCov <- "simple"
  settings <- list(
    rca = list(
        center.before = FALSE,
        center.after = FALSE,
        residuals.prods = TRUE,
        residual.cov.syntax = TRUE,
        constrained.prod.mean = FALSE,
        constrained.loadings = FALSE,
        constrained.var = FALSE,
        res.cov.method = defaultResCov,
        res.cov.across = TRUE,
        match = FALSE,
        match.recycle = FALSE),
    uca  = list(
        center.before = TRUE,
        center.after = FALSE,
        residuals.prods = FALSE,
        residual.cov.syntax = TRUE,
        constrained.prod.mean = TRUE,
        constrained.loadings = FALSE,
        constrained.var = FALSE,
        res.cov.method = defaultResCov,
        res.cov.across = TRUE,
        match = FALSE,
        match.recycle = FALSE),
    pind  = list(
        center.before = FALSE,
        center.after = FALSE,
        residuals.prods = FALSE,
        residual.cov.syntax = FALSE,
        constrained.prod.mean = FALSE,
        constrained.loadings = FALSE,
        constrained.var = FALSE,
        res.cov.method = defaultResCov,
        res.cov.across = TRUE,
        match = FALSE,
        match.recycle = FALSE),
    dblcent  = list(
        center.before = TRUE,
        center.after = TRUE,
        residuals.prods = FALSE,
        residual.cov.syntax = TRUE,
        constrained.prod.mean = FALSE,
        constrained.loadings = FALSE,
        constrained.var = FALSE,
        res.cov.method = defaultResCov,
        res.cov.across = TRUE,
        match = FALSE,
        match.recycle = FALSE),
    ca = list(
        center.before = TRUE,
        center.after = FALSE,
        residuals.prods = FALSE,
        residual.cov.syntax = TRUE,
        constrained.prod.mean = TRUE,
        constrained.loadings = TRUE,
        constrained.var = TRUE,
        res.cov.method = "ca",
        res.cov.across = FALSE,
        match = TRUE,
        match.recycle = TRUE)
    )

  settingNames <- unique(unlist(lapply(settings, FUN = names)))
  args <- args[settingNames]

  if (is.null(args)) return(settings[method])

  isMissing <- vapply(args, FUN.VALUE = logical(1L), FUN = is.null)
  missingArgs <- settingNames[isMissing]
  stopif(!method %in% names(settings), "Unrecognized method")
  
  c(settings[[method]][missingArgs], args[!isMissing])
}


#' default arguments for product indicator approaches
#'
#' @param method which method to get the settings for
#' @return list
#' @export
#' @description
#' This function returns the default settings for the product indicator approaches
#' @examples
#' library(modsem)
#' default_settings_pi()
default_settings_pi <- function(method = c("rca", "uca", "pind", "dblcent", "ca")) {
  getMethodSettingsPI(method = method, args = NULL)
}
