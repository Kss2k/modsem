getMethodSettingsPI <- function(method, args) {
    settingNames <- c("center.before", "center.after",
                    "residuals.prods", "residual.cov.syntax",
                    "constrained.prod.mean", "constrained.loadings",
                    "constrained.var", "constrained.res.cov.method")
    args <- args[settingNames]
    isMissing <- vapply(args, FUN.VALUE = logical(1L), FUN = is.null)
    missingArgs <- settingNames[isMissing]
    defaultResCov <- "simple" # could be defaultResCov both there is not really 
                              # support for it in the literature
    settings <- list(
        rca = list(
            center.before = FALSE,
            center.after = FALSE,
            residuals.prods = TRUE,
            residual.cov.syntax = TRUE,
            constrained.prod.mean = FALSE,
            constrained.loadings = FALSE,
            constrained.var = FALSE,
            constrained.res.cov.method =  defaultResCov),
        uca  = list(
            center.before = TRUE,
            center.after = FALSE,
            residuals.prods = FALSE,
            residual.cov.syntax = TRUE,
            constrained.prod.mean = TRUE,
            constrained.loadings = FALSE,
            constrained.var = FALSE,
            constrained.res.cov.method =  defaultResCov),
        pind  = list(
            center.before = FALSE,
            center.after = FALSE,
            residuals.prods = FALSE,
            residual.cov.syntax = FALSE,
            constrained.prod.mean = FALSE,
            constrained.loadings = FALSE,
            constrained.var = FALSE,
            constrained.res.cov.method =  defaultResCov),
        dblcent  = list(
            center.before = TRUE,
            center.after = TRUE,
            residuals.prods = FALSE,
            residual.cov.syntax = TRUE,
            constrained.prod.mean = FALSE,
            constrained.loadings = FALSE,
            constrained.var = FALSE,
            constrained.res.cov.method =  defaultResCov),
        ca = list(
            center.before = TRUE,
            center.after = FALSE,
            residuals.prods = FALSE,
            residual.cov.syntax = TRUE,
            constrained.prod.mean = TRUE,
            constrained.loadings = TRUE,
            constrained.var = TRUE,
            constrained.res.cov.method =  "ca")
        )

    if  (!method %in% names(settings)) {
        stop2("Unrecognized method")
    }
    c(settings[[method]][missingArgs], args[!isMissing])
}

