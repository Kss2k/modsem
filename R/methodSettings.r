
getMethodSettings <- function(method, args) {
    settingNames <- c("centerBefore", "centerAfter",
                    "residualsProds", "residualCovSyntax",
                    "constrainedProdMean", "constrainedLoadings",
                    "constrainedVar", "constrainedResCovMethod")
    args <- args[settingNames]
    isMissing <- vapply(args, FUN.VALUE = logical(1L), FUN = is.null)
    missingArgs <- settingNames[isMissing]

    settings <- list(
        rca = list(
            centerBefore = FALSE,
            centerAfter = FALSE,
            residualsProds = TRUE,
            residualCovSyntax = TRUE,
            constrainedProdMean = FALSE,
            constrainedLoadings = FALSE,
            constrainedVar = FALSE,
            constrainedResCovMethod =  "equality"),
        uca  = list(
            centerBefore = TRUE,
            centerAfter = FALSE,
            residualsProds = FALSE,
            residualCovSyntax = TRUE,
            constrainedProdMean = TRUE,
            constrainedLoadings = FALSE,
            constrainedVar = FALSE,
            constrainedResCovMethod =  "equality"),
        pind  = list(
            centerBefore = FALSE,
            centerAfter = FALSE,
            residualsProds = FALSE,
            residualCovSyntax = FALSE,
            constrainedProdMean = FALSE,
            constrainedLoadings = FALSE,
            constrainedVar = FALSE,
            constrainedResCovMethod =  "equality"),
        dblcent  = list(
            centerBefore = TRUE,
            centerAfter = TRUE,
            residualsProds = FALSE,
            residualCovSyntax = TRUE,
            constrainedProdMean = FALSE,
            constrainedLoadings = FALSE,
            constrainedVar = FALSE,
            constrainedResCovMethod =  "equality"),
        ca = list(
            centerBefore = TRUE,
            centerAfter = FALSE,
            residualsProds = FALSE,
            residualCovSyntax = TRUE,
            constrainedProdMean = TRUE,
            constrainedLoadings = TRUE,
            constrainedVar = TRUE,
            constrainedResCovMethod =  "ca")
        )

    if  (!method %in% names(settings)) {
        stop("Unrecognized method")
    }
    c(settings[[method]][missingArgs], args[!isMissing])
}

