getFromBack <- function(x, i = 0) {
  x[(length(x) - i):length(x)]
}

#' @importFrom R6 R6Class
#' @importFrom dplyr if_else
PathTracer <- R6::R6Class("PathTracer", public = list(
  pt = NULL, 
  paths = NULL,
  x = NULL,
  y = NULL,
  initialize = function(pt, x, y, addCovPt = TRUE) {
    # Remove any potential ':' from the model
    pt <- lapplyDf(pt, stringr::str_remove_all, pattern = ":")
    structuralVars <- pt[pt$op == "~" & pt$rhs != "1", c("lhs", "rhs")] |>
      unlist() |> unique() 
    pt <- pt[pt$lhs %in% structuralVars & pt$rhs %in% structuralVars, ]
    pt$mod[pt$mod == ""] <- 
      apply(pt[pt$mod == "", c("lhs", "op", "rhs")],
            MARGIN = 1,
            stringr::str_c, 
            collapse = "")
    if (addCovPt) {
      reversedCovPaths <- pt[pt$lhs != pt$rhs & 
                             pt$op == "~~", ] 
      colnames(reversedCovPaths) <- c("rhs", "op", "lhs", "mod")
      pt <- rbind(pt, reversedCovPaths)
    }
    self$pt <- pt
    self$x <- x 
    self$y <- y
    self$paths <- list()
  },
  trace = function(var = self$x, currentPath = c(), covCount = 0, 
                   from = "lhs") {
    if (covCount > 1) return(FALSE)
    if (var == self$y && from == "rhs") {
      self$paths <- c(self$paths, character(1))
      self$paths[[length(self$paths)]] <- currentPath
      return(TRUE)
    }
    branches <- self$pt[self$pt[[from]] == var, ]
    if (NROW(branches) == 0) return(FALSE)
    
    rest <- logical(length(nrow(branches)))
    for (i in seq_len(nrow(branches))) {
      branch <- branches[i, ]
      nextFrom <- from
      nextCovCount <- covCount
      if (from == "lhs") nextVar <- branch$rhs else nextVar <- branch$lhs
      # reverse direction of travel if covariance (but not before selecting next variable)
      if (branch$op == "~~") {
        nextCovCount <- covCount + 1
        switch(from, lhs = nextFrom <- "rhs", rhs = nextFrom <- "lhs")
      }
      rest <- c(rest, self$trace(nextVar, currentPath = c(currentPath, branch$mod),
                                 covCount = nextCovCount, from = nextFrom))
      
    } 
    any(rest)
  },
  clean = function() {
    squaredExprs <- purrr::map_chr(self$paths, function(x) {
             reps <- sort(x) |> table() |> as.data.frame()
             dplyr::if_else(reps[[2]] > 1, 
                     true = paste(reps[[1]], 
                                  reps[[2]], 
                                  sep = " ^ "),
                     false = reps[[1]]) |>
                       stringr::str_c(collapse = " * ")
          })
    reps <- squaredExprs |> table() |> as.data.frame()
    dplyr::if_else(reps[[2]] > 1, 
            true = paste(reps[[2]], 
                         reps[[1]], 
                         sep = " * "),
            false = reps[[1]]) |>
      stringr::str_c(collapse = " + ")
  },
  generateSyntax = function(parenthesis = TRUE, ...) {
    self$paths <- list()
    solvedSelf <- self$trace(...) 
    
    if (!solvedSelf) return(NA)
    if (parenthesis) return(paste0("(", self$clean(), ")"))
    self$clean()
  }
))


#' Estimate formulas for (co-)variance paths using Wright's path tracing rules
#'
#' @param pt A data frame with columns lhs, op, rhs, and mod, from modsemify(syntax)
#' @param x source variable
#' @param y destination variable
#' @param parenthesis if TRUE, the output will be enclosed in parenthesis
#' @param ... additional arguments passed to tracePath
#'
#' @return A string with the estimated path (simplified if possible)
#' @export 
#' @description
#' This function estimates the path from x to y using the path tracing rules, 
#' note that it only works with structural parameters, so "=~" are ignored. If you
#' you want to use the measurement model, it should work if you replace it "=~" with
#' "~" in the mod column of pt.
#' @examples
#' library(modsem)
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 +x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#'
#'   # Inner model
#'   Y ~ X + Z + X:Z
#''
#' pt <- modsemify(m1)
#' tracePath(pt, "Y", "Y") # variance of Y
tracePath <- function(pt, x, y, parenthesis = TRUE, ...) {
  pathTracer <- PathTracer$new(pt, x, y)
  out <- pathTracer$generateSyntax(parenthesis = parenthesis, ...)
  rm(pathTracer)
  out
}
