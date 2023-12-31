
expressionType <- function(op, lhs, rhs) {
if (is.null(op) & is.null(rhs)) {
    if (is.null(lhs)) {
      stop("Missing lhs")
    }

  }
}



createExprNode <- function(node, lhs = NULL, rhs = NULL,
                           priority = 1) {
  structure(list(node = node,
                 lhs = lhs,
                 rhs = rhs),
            class = "exprNode",
            priority = priority)
}



getMinTokenPriority <- function(listTokens, min = NA) {
 # browser()
  if (is.null(listTokens) || length(listTokens) == 0) {
    if (is.na(min)) {
      stop("Unable to find minimum priority for tokens")
    }
    return(min)
  }

  tokenPriority <- getTokenPriority(listTokens[[1]])
  if (is.na(min) || (tokenPriority < min)) {
    min <- tokenPriority
  }
  getMinTokenPriority(listTokens[-1], min)
}



chooseToken <- function(listTokens, i = 1, chosenTokenIdx = NULL,
                        leftClosures = list()) {

  if (is.null(listTokens) || i > length(listTokens)) {
    if (is.null(chosenTokenIdx)) {
      if (length(leftClosures) > 0) {
        stop("Unmatched left bracket", last(leftClosures))
      }
    }
    return(chosenTokenIdx)
  }
  token <- listTokens[[i]]

  if (is.LeftClosure(token)) {
    leftClosures <- appendToList(leftClosures, token)

  } else if (is.RightClosure(token)) {
    if (length(leftClosures) == 0) {
      stop("Unmatched right bracket", highlightErrorToken(token))
    }
    leftClosures <- leftClosures[-1]
    if (length(leftClosures) == 0) {
      return(i)
    }
  }

  if (length(leftClosures) == 0) {
    if (is.null(chosenTokenIdx)) {
      chosenTokenIdx <- i
    }
    if (getTokenPriority(token) < getTokenPriority(listTokens[[chosenTokenIdx]])) {
      chosenTokenIdx <- i
    }
  }
  chooseToken(listTokens, i + 1, chosenTokenIdx, leftClosures)
}


createSyntaxTreeLine <- function(listTokens, i = 1) {
  minPriority <- getMinTokenPriority(listTokens)
  if (i > length(listTokens) ||
        length(listTokens) == 1 &
          getTokenPriority(listTokens[[1]]) < 0) {
    return(NULL)
  }

  if (getTokenPriority(listTokens[[i]]) == minPriority) {
    if (i != 1) {
      lhs <- createSyntaxTreeLine(listTokens[1:(i-1)], i = 1)
    } else lhs <- NULL

    if (i != length(listTokens)) {
      rhs <- createSyntaxTreeLine(listTokens[(i + 1):length(listTokens)], i = 1)
    } else rhs <- NULL

    return(
      createExprNode(listTokens[[i]], lhs = lhs, rhs = rhs,
                   priority = getTokenPriority(listTokens[[i]]))
      )
  }
  createSyntaxTreeLine(listTokens, i + 1)
}



#
# createSyntaxTreeLine2 <- function(listTokens, i = 0) {
#   if (i > length(listTokens) ||
#       length(listTokens) == 1 &
#       getTokenPriority(listTokens[[1]]) < 0) {
#     return(NULL)
#   }
#
#
#
#   if (is.LeftClosure(listTokens[[i]])) {
#     leftClosure <- list(listTokens[[i]])
#   } else leftClosure <- list()
#   i <- chooseToken(listTokens, leftClosures = leftClosures)
#   if (i != 1) {
#     lhs <- createSyntaxTreeLine2(listTokens[1:(i-1)], i = 0)
#   } else lhs <- NULL
#
#   if (i != length(listTokens)) {
#     rhs <- createSyntaxTreeLine2(listTokens[(i + 1):length(listTokens)], i = 0)
#   } else rhs <- NULL
#
#   return(
#     createExprNode(listTokens[[i]], lhs = lhs, rhs = rhs,
#                    priority = getTokenPriority(listTokens[[i]]))
#   )
#
# }

createSyntaxTreesSyntax <- function(syntax) {
  tokensLines <- tokenizeSyntax(syntax)
  lapply(tokensLines, createSyntaxTreeLine)
}



last <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  x[[length(x)]]
}
