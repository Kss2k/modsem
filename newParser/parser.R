


evalToken <- function(token, lhs, rhs) {
  UseMethod("evalToken")
}



evalToken.LavOperator <- function(token, lhs, rhs) {

  if ("LavToken" %in% class(rhs)) {
    rhs <- list(rhs)
  }
  if (length(lhs) > 0 && !"LavToken" %in% class(lhs)) {
    stop("Incorrect number of elements on the lhs of operator",
         highlightErrorToken(token))
  } else if (!is.atomic(lhs)) {
    if (is.LavOperator(lhs$op)) {
      stop("Unexpected operator ", highlightErrorToken(lhs$op))
    }
  } else if (!is.atomic(rhs)) {
    if (is.LavOperator(rhs$op)) {
      stop("Unexpected operator ", highlightErrorToken(rhs$op))
    }
  }
  list(lhs = lhs, op = token, rhs = rhs)
}



evalToken.LavToken <- function(token, lhs, rhs) {
  token
}



evalToken.LavAdd <- function(token, lhs, rhs) {
  if ("LavToken" %in% class(rhs)) {
    rhs <- list(rhs)
  }
  if ("LavToken" %in% class(lhs)) {
    lhs <- list(lhs)
  }
  c(lhs, rhs)
}



evalToken.LavModify <- function(token, lhs, rhs) {
  structure(rhs,
            modifier = lhs)
}



evalToken.LavBlank <- function(token, lhs, rhs) {
  NULL
}

evalToken.LavInteraction <- function(token, lhs, rhs) {
  mergeTokens(lhs, rhs)
}

evalToken.LavComment <- function(token, lhs, rhs) {
  NULL
}

evalToken.LavFunction <- function(token, lhs, rhs) {
  browser()
  updateVariablesEnvir()
  functionCall <- paste0(token, stringr::str_c(unlist(rhs), collapse = ","), ")")
  functionCall
  eval(rlang::parse_expr(functionCall), envir = modVarsEnv)
}

evalToken.LeftBracket <- function(token, lhs, rhs) {
  # if (!is.null(rhs$rightBracket) || !rhs$rightBracket) {
  #   stop("Unmatched bracket ", highlightErrorToken(token))
  # }
  rhs
}



evalToken.RightBracket <- function(token, lhs, rhs) {
  #list(body = lhs, rightBracket = TRUE)
  lhs
}



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
  lapply(syntaxTrees,
         evalTokens)
}



createParTableBranch <- function(syntaxTree) {
  if (is.null(syntaxTree)) {
    return(NULL)
  }
  # Run a function to check wheter there is nesting in Rhs, which throws an error
  #maxDepth(syntaxTree[["rhs"]])
  rhs <- vector("character", length(syntaxTree[["rhs"]]))
  mod <- rhs

  for (i in seq_along(syntaxTree[["rhs"]])) {
    rhs[[i]] <- getTokenString(syntaxTree[["rhs"]][[i]])
    modifier <- getTokenString(attr(syntaxTree[["rhs"]][[i]], "modifier"))
    if (!is.null(modifier)) {
      mod[[i]] <- modifier
    }
  }

  lhs <- rep(getTokenString(syntaxTree$lhs), length(rhs))
  op <- rep(getTokenString(syntaxTree$op), length(rhs))
  data.frame(lhs = lhs, op = op, rhs = rhs, mod = mod)
}




modsemify <- function(syntax) {
  if (!is.character(syntax) && length(syntax) > 1) {
    stop("Syntax is not a string og length 1")
  }
  syntaxTrees <- createSyntaxTreesSyntax(syntax)
  parsedTrees <- parseSyntaxTrees(syntaxTrees)
  purrr::list_rbind(lapply(parsedTrees,
                           FUN = createParTableBranch))
}



parTableToSyntax <- function(parTable, removeColon = FALSE) {
  out <- ''
  for (i in 1:nrow(parTable)) {
    if (parTable[["mod"]][i] != "") {
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
  if (removeColon == TRUE) {
    out <- stringr::str_remove_all(out, ":")
  }
  out
}




mergeTokens <- function(x, y) {
  if (!"LavName" %in% class(x) || !"LavName" %in% class(x)) {
    stop("Interactions are reserved for objects ", highlightErrorToken(x))
  }


  out <- paste0(x, y)
  attributes(out) <- attributes(token)
  out

}
