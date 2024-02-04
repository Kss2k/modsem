# ------------------------------------------------------------------------------
# Function for running the same model with multiple methods 
allMethods <- c("rca", "uca", "ca", "dblcent", "mplus", "pind")
allNativeMethods <- allMethods[allMethods != "mplus"]
fastMethods <- c("rca", "uca", "dblcent", "pind")

runMultipleMethods <- function(modelSyntax, 
                               data, 
                               methods = allNativeMethods,
                               ...) {
  estimates <- structure(vector("list", length = length(methods)),
                         names = methods)
  for (method in methods) {
    estimates[[method]] <- tryCatch(
      modsem(modelSyntax, data, method, ...),
      warning = function(w) {
        warning("Warning in ", method, "\n")
        capturePrint(w)
        modsem(modelSyntax, data, method, ...)
      },
      error = function(e) {
        warning("Error in ", method, "\n")
        capturePrint(e)
        NA
      }
    )
  }
  estimates
}
