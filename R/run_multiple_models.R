# ------------------------------------------------------------------------------
# Function for running the same model with multiple methods 
allMethods <- c("rca", "uca", "ca", "dblcent", "mplus", "pind")
allNativeMethods <- allMethods[allMethods != "mplus"]
fastMethods <- c("rca", "uca", "dblcent", "pind")

runMultipleMethods <- function(model.syntax, 
                               data, 
                               methods = allNativeMethods,
                               ...) {
  estimates <- structure(vector("list", length = length(methods)),
                         names = methods)
  for (method in methods) {
    estimates[[method]] <- tryCatch(
      modsem(model.syntax, data, method, ...),
      warning = function(w) {
        warning2("Warning in ", method, "\n", capturePrint(w), "\n")
        modsem(model.syntax, data, method, ...)
      },
      error = function(e) {
        warning2("Error in ", method, "\n", capturePrint(e))
        NA
      }
    )
  }
  estimates
}
