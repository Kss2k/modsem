
getLines <- function(syntax) {
  # split syntax into individual lines and return as vector
  unlist(strsplit(syntax, "\n|;"))
}
get


getCharsLine <- function(line, i = 1) {
  if (i > nchar(line)) {
    return(NULL)
  }
  rest <- getCharsLine(line, i + 1)
  c(substr(line, i, i),rest)
}



tokenizeSyntax <- function(syntax) {
  lines <- lapply(getLines(syntax),
                FUN = getCharsLine)
  # remove NULL's
  isNull <- vapply(lines,
                   FUN = is.null,
                   FUN.VALUE = vector("logical", length = 1L))

  tokensLines <- lapply(lines[!isNull],
         getTokensLine)
  # remove lines which are list(0) since they only
  isListZeroLen <- vapply(tokensLines,
                          FUN = function(list) length(list) <= 0 ||
                            areAllLavBlank(list),
                          FUN.VALUE = vector("logical", length = 1L))
  tokensLines[!isListZeroLen]

}



# I want this to be a recursive function, which builds tokens as it goes,
  # and stops builing each token when certain inputs are given.
getTokensLine <- function(line, i = 1,
                          token = structure("", class = "LavBlank"),
                          vecTokens = list()) {

  if (!is.character(line)) {
    stop("Expected string in parseLine()")
  }
  if (grepl("\\s+", token)) {
    token <- structure("",
                       class = "LavBlank")
  }

  if (nchar(token) == 1) {
    if (grepl("[[:alpha:]_.]", token)) {
      class(token) <- "LavName"

    } else if (grepl("[\\=\\~\\*\\+\\(\\)\\<\\>\\-\\,\\:\\^]", token)) {
      class(token) <- "LavOperator"

    } else if (grepl("[[:alnum:]]", token)) {
      class(token) <- "LavNumeric"

    } else if (grepl('\\"' , token)) {
      class(token) <- "LavString"

    }else {
      stop("Unrecognized class of token", token)
    }

  }

  # Break and return if you reach the end of the line, or a comment
  if (i > length(line) || line[[i]] == "#") {
    vecTokens[[length(vecTokens) + 1]] <- token
    return(vecTokens)
  }

  # Whatever case, token should be finished if you encounter \\s+
  if (grepl("\\s+", line[[i]])) {
    vecTokens[[length(vecTokens) + 1]] <- token
    out <- getTokensLine(line, i + 1, token = line[[i]], vecTokens)

  } else if (!grepl("\\s+", line[[i]])) {

    if (fitsToken(token, nextChar = line[[i]])) {

      updatedToken <- structure(paste0(token, line[[i]]),
                                class = class(token))


      out <- getTokensLine(line, i + 1, updatedToken, vecTokens)

    } else if (!fitsToken(token, nextChar = line[[i]])) {
      vecTokens[[length(vecTokens) + 1]] <- token
      out <- getTokensLine(line, i + 1, token = line[[i]], vecTokens)
    }

  }
  # remove LavBlank from tokens before returning
  isNotLavBlank <- vapply(out,
         FUN = function(x) class(x) != "LavBlank",
         FUN.VALUE = vector("logical", length = 1L))
  out[isNotLavBlank]
}



fitsToken <- function(token, ...) {
  UseMethod("fitsToken")
}



fitsToken.LavName <- function(token, nextChar) {
  if (length(nextChar) != 1) {
    stop("Wrong length of nextChar", nextChar)
  }
  # if object name ends with ( it is a function,
    # and next char belongs to a new object
  if (grepl("\\($", token)) {
    return(FALSE)
  }
  grepl("[[:alpha:][:digit:]_.\\(]", nextChar)[[1]]
}



fitsToken.LavString <- function(token, nextChar) {
  if (length(nextChar) != 1) {
    stop("Wrong length of nextChar", nextChar)
  }
  # if object name ends with ( it is a function,
  # and next char belongs to a new object
  if (grepl('\\"$', token)) {
    return(FALSE)
  }
  grepl("[[:graph:][:space:]]", nextChar)[[1]]
}



fitsToken.LavOperator <- function(token, nextChar) {
  if (length(nextChar) != 1) {
    stop("Wrong length of nextChar", nextChar)
  }
  completeToken <- paste0(token, nextChar)
  switch(completeToken,
         "=~" = TRUE,
         "~~" = TRUE,
         "<-" = TRUE,
         "->" = TRUE,
         "==" = TRUE,
         "!=" = TRUE,
         FALSE)
  #grepl("[\\=\\~\\*\\+\\(\\)\\<\\>\\-]", nextChar)[[1]]
}



fitsToken.LavBlank <- function(token, nextChar) {
  TRUE
}



fitsToken.LavNumeric <- function(token, nextChar) {
  if (length(nextChar) != 1) {
    stop("Wrong length of nextChar", nextChar)
  }
  grepl("[[:digit:].]", nextChar)
}


getClassLavOp <- function(op) {
  if (length(op) != 1 || !is.character(op)) {
    stop("incorrect input in getClassLavOp(): ")
  }
  # does not support all types of operators yet
  switch (op,
    "=~" = "LavMeasure",
    "~"  = "LavPredict",
    "~~" = "LavCovar",
    "+"  = "LavAdd",
    "*"  = "LavModify",
    "<"  = "LavLessLeft",
    ">"  = "LavLessRight",
    "==" = "LavEquals",
    ":"  = "LavInteraction",
    "("  = "LeftBracket",
    ")"  = "RightBracket",
    "="  = "LavAssign",
    stop("Unrecognized operator")
  )
}



isSpecificationOperator <- function(op) {
  grepl("LavMeasure|LavPredict|LavCovar|LavAssign", getClassLavOp(op))
}



isEqualityOperator <- function(op) {
  grepl("LavLessLeft|LavLessRight|LavEquals", getClassLavOp(op))
}



isModificationOperator <- function(op) {
  grepl("LavModify|LavAdd|LeftBracket|RightBracket", getClassLavOp(op))
}


isInteractionOperator <- function(op) {
    grepl("LavInteraction", getClassLavOp(op))
}


isAssignmentOperator <- function(op) {
  grepl("LavAssign", getClassLavOp(op))
}


getExpressionType <- function(op) {
  if (isSpecificationOperator(op)) {
    "specification"
  } else if (isEqualityOperator(op)) {
    "equality"
  } else if (isInteractionOperator(op)) {
    "empty"
  } else {
    stop("Unrecognized expression type based on first operator in line, operator: ", op)
  }
}



getClassLavName <- function(name) {
  if (length(name) != 1 || !is.character(name)) {
    stop("incorrect input in getClassLavName()")
  }
  if(grepl("\\($", name)) {

    return("LavFunction")
  }
  "LavObject"
}



doesOperatorFitExprType <- function(op, expressionType) {
  if (is.null(expressionType)) {
    expressionType <- "empty"
  }
  switch (expressionType,
    "empty" = isSpecificationOperator(op) |
              isEqualityOperator(op) |
              isInteractionOperator(op),
    "specification" = isModificationOperator(op) |
                      isInteractionOperator(op),
    "equality" = !isSpecificationOperator(op),
    stop("unrecognizxed expression type"))
}


evalTokens <- function(listTokens,
                       lhs = NULL,
                       op = NULL,
                       rhs = NULL,
                       expressionType = "empty") {

  if (1 > length(listTokens)) {
    return(NULL)

  } else if (length(listTokens) == 1) {
    return(listTokens)
  }

  if (is.null(lhs)) {
    if (class(listTokens[[1]]) == "LavName") {
      # figure out wheter it is an object or a function call
      className <- getClassLavName(listTokens[[1]])

      if (className == "LavObject") {
        lhs <- listTokens[[1]]
        restParsed <- evalTokens(listTokens[-1, drop = FALSE], #i + 1,
                                 lhs = lhs,
                                 op = op,
                                 rhs = rhs,
                                 expressionType = expressionType)
        return(restParsed)
      } else if (className == "LavFunction") {

        seqWithFunctionCall <- 1:findFunctionEnd(listTokens)

        outputName <- evalLavFunction(listTokens[seqWithFunctionCall])
        # replace sequence with output name, and recursing on rest
        if (length(seqWithFunctionCall) < length(listTokens)) {
          restParsed <- evalTokens(listTokens[-seqWithFunctionCall],
                                   lhs = outputName,
                                   op = op,
                                   rhs = rhs,
                                   expressionType = expressionType)
          return(restParsed)
        }
        return(outputName)
      }


    } else if (class(listTokens[[1]]) == "LavNumeric") {
      lhs <- listTokens[[1]]
      restParsed <- evalTokens(listTokens[-1, drop = FALSE],
                               lhs = lhs,
                               op = op,
                               rhs = rhs,
                               expressionType = expressionType)
      return(restParsed)

    } else {
      stop("Expected name at the start of: ",listTokens)
    }

  } else if (is.null(op)) {
    if (class(listTokens[[1]]) == "LavOperator") {
      op <- listTokens[[1]]
      class(op) <- getClassLavOp(op)
      restExpression <- listTokens[-1, drop = FALSE]

      if (!doesOperatorFitExprType(op, expressionType)) {
        stop("Unexpected operator at: ", paste(lhs, listTokens))

      } else if (length(listTokens) == 1) {
        return(listTokens)
      }
      if (expressionType == "empty") {
        expressionType <- getExpressionType(op)
      }

      return(evalOp(op, lhs = lhs, rhs = restExpression, expressionType = expressionType))

    } else {
      stop("Expected operator after object name: ", listTokens[[1]])
    }
  }
}



evalOp <- function(op, lhs, rhs, expressionType) {
  UseMethod("evalOp")
}



evalOp.LavMeasure <- function(op, lhs, rhs, expressionType = "specification") {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs, expressionType = expressionType)

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType)
}



evalOp.LavPredict <- function(op, lhs, rhs, expressionType = "specification") {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs, expressionType = expressionType)

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType )
}



evalOp.LavCovar <- function(op, lhs, rhs, expressionType = "specification") {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs, expressionType = expressionType)

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType)
}



evalOp.LavEquals <- function(op, lhs, rhs, expressionType = "equality") {
  # op is just here to fetch the method
  rest <- stringr::str_c(unlist(rhs), collapse = "")

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType)
}



evalOp.LavLessRight <- function(op, lhs, rhs, expressionType = "equality") {
  # op is just here to fetch the method
  rest <- stringr::str_c(unlist(rhs), collapse = "")

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType)
}



evalOp.LavLessLeft <- function(op, lhs, rhs, expressionType = "equality") {
  # op is just here to fetch the method
  rest <- stringr::str_c(unlist(rhs), collapse = "")

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType)
}



evalOp.LavAdd <- function(op, lhs, rhs, expressionType = "Specification") {
  # op is just here to fetch method
  if (length(rhs) == 1) {
    return(list(lhs, rhs[[1]]))
  } else if (length(rhs) <= 0 ||
             !(class(rhs[[1]]) %in% c("LavName", "LavNumeric"))) {
    stop("expected name after +")
  }
  rest <- evalTokens(listTokens = rhs, expressionType = expressionType)

  c(list(lhs), rest)
}



evalOp.LavModify <- function(op, lhs, rhs, expressionType = "empty") {
  rest <- evalTokens(rhs, expressionType = expressionType)
  attr(rest[[1]], "LavMod") <- lhs
  rest
}



evalOp.LavInteraction <- function(op, lhs, rhs, expressionType = "empty") {
  rest <- evalTokens(listTokens = rhs, expressionType = expressionType)
  # combine lhs and rest into one, inheriting attributes from lhs
  combinedLhsRhs <- paste0(lhs, ":", rest[[1]])

  if (expressionType == "empty") {
    return(c(lhs = combinedLhsRhs, rest[-1, drop = FALSE]))
  }

  if (length(rest) <= 1) {
    return(list(combinedLhsRhs))
  }

  c(combinedLhsRhs, rest[-1, drop = FALSE])
}



createSyntaxTree <- function(syntax) {
  tokensLines <- tokenizeSyntax(syntax)
  lapply(tokensLines,
         FUN = evalTokens)
}



createParTableBranch <- function(syntaxBranch) {
  # Run a function to check wheter there is nesting in Rhs, which throws an error
  #maxDepth(syntaxBranch[["rhs"]])
  rhs <- vector("character", length(syntaxBranch[["rhs"]]))
  mod <- rhs

  for (i in seq_along(syntaxBranch[["rhs"]])) {
    rhs[[i]] <- syntaxBranch[["rhs"]][[i]]
    LavModify <- attr(syntaxBranch[["rhs"]][[i]], "LavMod")
    if (!is.null(LavModify)) {
      mod[[i]] <- LavModify
    }
  }

  lhs <- rep(syntaxBranch$lhs, length(rhs))
  op <- rep(syntaxBranch$op, length(rhs))
  data.frame(lhs = lhs, op = op, rhs = rhs, mod = mod)
}



modsemify <- function(syntax) {
  syntaxTree <- createSyntaxTree(syntax)
  purrr::list_rbind(lapply(syntaxTree,
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




areAllLavBlank <- function(list) {
  isLavBlank <- vapply(list, FUN = function(x) class(x) != "LavBlank",
         FUN.VALUE = vector("logical", length = 1L))
  !as.logical(sum(isLavBlank))
}
