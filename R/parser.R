evalToken <- function(token, lhs, rhs) {
  UseMethod("evalToken")
}


#' @export
evalToken.LavOperator <- function(token, lhs, rhs) {
  if (is.LavToken(rhs)) rhs <- list(rhs)
  if (is.LavToken(lhs)) lhs <- list(lhs)

  stopif(!is.atomic(lhs) && is.LavOperator(lhs$op), "Unexpected operator ", 
         highlightErrorToken(lhs$op))

  stopif(!is.atomic(rhs) && is.LavOperator(rhs$op), "Unexpected operator ", 
         highlightErrorToken(rhs$op))

  list(lhs = lhs, op = token, rhs = rhs)
}


#' @export
evalToken.LavToken <- function(token, lhs, rhs) {
  token
}


#' @export
evalToken.LavName <- function(token, lhs, rhs) {
  stopif(is.LavToken(rhs), "Unexpected token ", 
         highlightErrorToken(rhs))
  stopif(is.LavToken(lhs), "Unexpected token ",
         highlightErrorToken(lhs))
  token
}


#' @export
evalToken.LavNumeric <- function(token, lhs, rhs) {
  stopif(is.LavToken(rhs), "Unexpected token ", 
         highlightErrorToken(rhs))
  stopif(is.LavToken(lhs), "Unexpected token ",
         highlightErrorToken(lhs))
  token
}


#' @export
evalToken.LavAdd <- function(token, lhs, rhs) {
  if (is.LavToken(rhs)) {
    rhs <- list(rhs)
  }
  if (is.LavToken(lhs)) {
    lhs <- list(lhs)
  }
  c(lhs, rhs)
}


#' @export
evalToken.LavModify <- function(token, lhs, rhs) {
  structure(rhs,
            modifier = lhs)
}


#' @export
evalToken.LavBlank <- function(token, lhs, rhs) {
  NULL
}


#' @export
evalToken.LavInteraction <- function(token, lhs, rhs) {
  if (!"LavName" %in% class(lhs) || !"LavName" %in% class(rhs)) {
    stop2("Interactions are reserved for objects ", highlightErrorToken(token))
  }
  out <- paste0(lhs, token, rhs)
  attributes(out) <- attributes(lhs)
  out
}


#' @export
evalToken.LavComment <- function(token, lhs, rhs) {
  NULL
}


#' @export
evalToken.LavFunction <- function(token, lhs, rhs) {
  functionCall <- paste0(token, stringr::str_c(unlist(rhs), collapse = ","), ")")
  out <- eval(rlang::parse_expr(functionCall), envir = modEnv)
  attributes(out) <- attributes(token)
  out
}


#' @export
evalToken.LeftBracket <- function(token, lhs, rhs) {
  rhs
}


#' @export
evalToken.RightBracket <- function(token, lhs, rhs) {
  lhs
}


#' @export
evalToken.LavSeperator <- function(token, lhs, rhs) {
  if ("LavToken" %in% class(rhs)) {
    rhs <- list(rhs)
  }
  if ("LavToken" %in% class(lhs)) {
    lhs <- list(lhs)
  }
  c(lhs, rhs)
}


evalTokens <- function(syntaxTree) {
  if (is.null(syntaxTree) || length(syntaxTree) == 0) {
    return(NULL)
  }

  lhs <- evalTokens(syntaxTree$lhs)
  rhs <- evalTokens(syntaxTree$rhs)

  return(evalToken(syntaxTree$node,
                   lhs = lhs,
                   rhs = rhs))
}


parseSyntaxTrees <- function(syntaxTrees) {
  lapply(syntaxTrees, evalTokens)
}


createParTableBranch <- function(syntaxTree) {
  if (is.null(syntaxTree)) {
    return(NULL)
  }
  rhs <- vector("character", length(syntaxTree[["rhs"]]))
  mod <- rhs

  for (i in seq_along(syntaxTree[["rhs"]])) {
    rhs[[i]] <- getTokenString(syntaxTree[["rhs"]][[i]])
    modifier <- getTokenString(attr(syntaxTree[["rhs"]][[i]], "modifier"))
    if (!is.null(modifier)) {
      mod[[i]] <- modifier
    }
  }
  lhs <- vapply(syntaxTree[["lhs"]], FUN.VALUE = character(1L),
                FUN = getTokenString)
  lhs <- vapply(lhs, FUN.VALUE = character(length(rhs)),
                FUN = function(x, len)
                  rep(x, len),
                len = length(rhs)) |> as.vector()
  op <- rep(getTokenString(syntaxTree$op), length(rhs))
  parTable <- data.frame(lhs = lhs, op = op, rhs = rhs, mod = mod)

  # post-processing
  parTable[parTable$op == ":=", "mod"] <- parTable[parTable$op == ":=", "lhs"]
  intercepts <- parTable$op == "~" & parTable$rhs == "1"
  parTable[intercepts, "op"]  <- "~1"
  parTable[intercepts, "rhs"] <- ""

  parTable
}


#' Generate parameter table for \code{lavaan} syntax
#'
#' @param syntax model syntax
#'
#' @return \code{data.frame} with columns \code{lhs, op, rhs, mod}
#' @export modsemify
#'
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
#' modsemify(m1)
modsemify <- function(syntax) {
  stopif(!is.character(syntax) && length(syntax) > 1,
         "Syntax is not a string og length 1")
  syntaxTrees <- createSyntaxTreesSyntax(syntax)
  parsedTrees <- parseSyntaxTrees(syntaxTrees)
  purrr::list_rbind(lapply(parsedTrees,
                           FUN = createParTableBranch))
}


parTableToSyntax <- function(parTable, removeColon = FALSE) {
  intercepts <- parTable$op == "~1"
  parTable[intercepts, "rhs"] <- "1"
  parTable[intercepts, "op"]  <- "~"

  out <- ""
  if (removeColon) {
    parTable$lhs <- stringr::str_remove_all(parTable$lhs, ":")
    parTable$rhs <- stringr::str_remove_all(parTable$rhs, ":")
    parTable$mod <- stringr::str_remove_all(parTable$mod, ":")
  }
  for (i in 1:nrow(parTable)) {
    if (parTable[["mod"]][[i]] != "" && parTable[["op"]][[i]] != ":=") {
      modifier <- paste0(parTable[["mod"]][[i]], "*")
    } else {
      modifier <- ""
    }
    line <- paste0(parTable[["lhs"]][[i]], " ",
                   parTable[["op"]][[i]], " ",
                   modifier,
                   parTable[["rhs"]][[i]], "\n")
    out <- paste0(out, line)
  }

  out
}


mergeTokens <- function(x, y) {
  stopif(!"LavName" %in% class(x) || !"LavName" %in% class(x),
         "Interactions are reserved for objects ",
         highlightErrorToken(x))

  out <- paste0(x, y)
  attributes(out) <- attributes(x)
  out

}
