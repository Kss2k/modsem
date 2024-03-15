getFromBack <- function(x, i = 0) {
  x[(length(x) - i):length(x)]
}

PathTracer <- R6::R6Class("PathTracer", public = list(
  pt = NULL,
  paths = NULL,
  x = NULL,
  y = NULL,
  nParams = NULL,
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
    self$nParams <- 0
    
  },
  

  trace = function(var = self$x, currentPath = c(), covCount = 0, 
                   lastDirection = "lhs->rhs") {
    if (!is.null(currentPath) && covCount > 1) return(FALSE)
    if (var == self$y && lastDirection == "rhs->lhs") {
      self$nParams <- self$nParams + length(currentPath)
      self$paths <- c(self$paths, character(1))
      self$paths[[length(self$paths)]] <- currentPath
      return(TRUE)
    }
    branches <- self$pt[self$pt$lhs == var | self$pt$rhs == var, ]
    
    if (lastDirection == "lhs->rhs") {
      branches <- branches[branches$lhs == var, ]#& 
                           #branches$lhs != branches$lhs, ]
      nextDirection <- "lhs->rhs"
    } else if (lastDirection == "rhs->lhs") {
      branches <- branches[branches$rhs == var, ]
      nextDirection <- "rhs->lhs"
    }
    if (NROW(branches) == 0) return(FALSE)
    rest <- logical(length(nrow(branches)))
    for (i in seq_len(nrow(branches))) {
      branch <- branches[i, ]
      if (nextDirection == "lhs->rhs") nextVar <- branch$rhs
      else nextVar <- branch$lhs
      # reverse direction of travel if covariance (but not before selecting next variable)
      if (branch$op == "~~") {
        covCount <- covCount + 1
        if (lastDirection == "lhs->rhs") nextDirection <- "rhs->lhs"
        else nextDirection <- "lhs->rhs"
      }
      rest <- c(rest, self$trace(nextVar, 
                                 currentPath = c(currentPath, branch$mod),
                                 covCount = covCount,
                                 lastDirection = nextDirection))
      if (branch$op == "~~") {
        # undo reversal and count (not pretty, but it works)
        covCount <- covCount - 1
        if (lastDirection == "lhs->rhs") nextDirection <- "lhs->rhs"
        else nextDirection <- "rhs->lhs"
      }
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
  generateSyntax = function(enclose = TRUE, ...) {
    self$paths <- list()
    self$nParams <- 0
    solvedSelf <- self$trace(...) 
    
    if (!solvedSelf) return(NA)
    out <- self$clean()
    if (enclose) out <- paste("(", out, ")", sep = "")
    out
  }
))


#' Estimate formulas for (co-)variance paths using Wright's path tracing rules
#'
#' @param pt A data frame with columns lhs, op, rhs, and mod, from modsemify(syntax)
#' @param x source variable
#' @param y destination variable
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
tracePath <- function(pt, x, y, ...) {
  pathTracer <- PathTracer$new(pt, x, y)
  out <- pathTracer$generateSyntax(...)
  rm(pathTracer)
  out
}
