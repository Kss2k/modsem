
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
  modsemParseEnv$lines <- lines[!isNull]
    tokensLines <- lapply(
      lines[!isNull],
      function(line) {
        modsemParseEnv$currentLine <- modsemParseEnv$currentLine + 1
        out <- getTokensLine(line)
        modsemParseEnv$tokensInLine <- list()
        out
        })
  # remove lines which are list(0) since they only
  isListZeroLen <- vapply(tokensLines,
                          FUN = function(list) length(list) <= 0 ||
                            areAllLavBlank(list),
                          FUN.VALUE = vector("logical", length = 1L))
  modsemParseEnv$lines <- tokensLines[!isListZeroLen]
  tokensLines[!isListZeroLen]
}



# I want this to be a recursive function, which builds tokens as it goes,
  # and stops builing each token when certain inputs are given.
getTokensLine <- function(line, i = 1,
                          token = structure("", class = "LavBlank"),
                          vecTokens = list()) {
  modsemParseEnv$tokensInLine <- c(modsemParseEnv$tokensInLine, token)

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

    } else {
      parseError("Unrecognized class of token: ")
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
    parseError("Unrecognized operator: ")
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
    parseError("Unrecognized expression type based on first operator in line")
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
    parseError("Unrecognized expression type"))
}


evalTokens <- function(listTokens,
                       lhs = NULL,
                       op = NULL,
                       rhs = NULL,
                       expressionType = "empty") {
  modsemParseEnv$tokensInLine <- c(modsemParseEnv$tokensInLine, listTokens[[1]])

  if (1 > length(listTokens)) {
    return(NULL)

  } else if (length(listTokens) == 1) {
    if (class(listTokens[[1]]) %in% c("LavName", "LavNumeric", "LavString")) {
      if (is.null(op)) {
        parseError("Expected operator")
      }
      return(listTokens)
    } else {
      parseError("Line should not end in operator")
    }
  }

  if (is.null(lhs)) {
    if (class(listTokens[[1]]) == "LavName") {
      # figure out wheter it is an object or a function call
      className <- getClassLavName(listTokens[[1]])

      if (className == "LavObject") {
        lhs <- listTokens[[1]]
        restParsed <- evalTokens(listTokens[-1, drop = FALSE], #i + 1,
                                 lhs = lhs,
                                 op = NULL,
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
                                   op = NULL,
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
      parseError("Expected name at the start of line ")
    }

  } else if (is.null(op)) {
    if (class(listTokens[[1]]) == "LavOperator") {
      op <- listTokens[[1]]
      class(op) <- getClassLavOp(op)
      restExpression <- listTokens[-1, drop = FALSE]

      if (!doesOperatorFitExprType(op, expressionType)) {
        # stop("Unexpected operator at: ", paste(lhs, listTokens))
        parseError("Unexpected operator at")

      } else if (length(listTokens) == 1) {
        return(listTokens)
      }
      if (expressionType == "empty") {
        expressionType <- getExpressionType(op)
      }
      return(evalOp(op, lhs = lhs, rhs = restExpression,
                    expressionType = expressionType))
      # return(
      #   tryCatch(evalOp(op, lhs = lhs, rhs = restExpression,
      #                   expressionType = expressionType),
      #            error = function(e) {
      #               parseError("Error in evaluating operator")
      #            }))

    } else {
      parseError("Invalid modifier symbol (should be '*' or '?') at line ")
      stop("ModSEM ERROR: \n",
           "invalid modifier symbol (should be '*' or '?') at line ",
           getLineNumber(), " pos ", getTokenPosition(), " token ",
           listTokens[[1]], "\n",
           highlightError())
    }
  }
}



evalOp <- function(op, lhs, rhs, expressionType) {
  UseMethod("evalOp")
}



evalOp.LavMeasure <- function(op, lhs, rhs, expressionType = "specification") {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs, op = "=~", expressionType = expressionType)

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType)
}



evalOp.LavPredict <- function(op, lhs, rhs, expressionType = "specification") {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs, op = "~", expressionType = expressionType)

  list(lhs = lhs, op = op, rhs = rest, expressionType = expressionType )
}



evalOp.LavCovar <- function(op, lhs, rhs, expressionType = "specification") {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs, op = "~~", expressionType = expressionType)

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
  if (length(rhs) <= 0 || !(class(rhs[[1]]) %in% c("LavName", "LavNumeric"))) {
    modsemParseEnv$tokensInLine <- c(modsemParseEnv$tokensInLine, rhs[[1]])
    parseError("Expected name or number after '+'")

  } else if (length(rhs) == 1) {
    return(list(lhs, rhs[[1]]))
  }
  rest <- evalTokens(listTokens = rhs,op = "+", expressionType = expressionType)

  c(list(lhs), rest)
}



evalOp.LavModify <- function(op, lhs, rhs, expressionType = "empty") {
  rest <- evalTokens(rhs, op = "*", expressionType = expressionType)
  attr(rest[[1]], "LavMod") <- lhs
  rest
}



evalOp.LavInteraction <- function(op, lhs, rhs, expressionType = "empty") {
  rest <- evalTokens(listTokens = rhs, op = ":", expressionType = expressionType)
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



evalOp.character <- function(op, lhs, rhs, expressionType) {
  parsError("Unexpected operator at line")
}



createSyntaxTree <- function(syntax) {
  tokensLines <- tokenizeSyntax(syntax)
  modsemParseEnv$currentLine <- 0
  lapply(tokensLines,
         function(tokensLine) {
           modsemParseEnv$currentLine <- modsemParseEnv$currentLine + 1
           modsemParseEnv$tokensInLine <- list()
           evalTokens(tokensLine)
         })
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
  resetModsemParseEnv()
  if (!is.character(syntax) && length(syntax) > 1) {
    stop("Syntax is not a string og length 1")
  }
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



getTokenPosition <- function() {
  modsemParseEnv$currentLine
}



getLineNumber <- function() {
  modsemParseEnv$currentLine
}


getCurrentToken <- function() {
  modsemParseEnv$tokensInLine[[length(modsemParseEnv$tokensInLine)]]
}



modsemParseEnv <- rlang::env(
  lines = NULL,
  currentLine = 0,
  tokensInLine = list()
)



resetModsemParseEnv <- function() {
  modsemParseEnv$lines <- NULL
  modsemParseEnv$currentLine <- 0
  modsemParseEnv$tokensInLine <- list()
}


highlightError <- function() {
  tokensInLine <- modsemParseEnv$tokensInLine
  tokensRegex <- vapply(tokensInLine,
                        FUN.VALUE = character(1L),
                        FUN = getRegexToken)
  lastElem <- length(tokensRegex)
  errorPattern <-
    paste0("(",stringr::str_c(tokensRegex[-lastElem],collapse = "\\s*"),
           "\\s*)(", tokensRegex[[lastElem]], ")")
  line <- modsemParseEnv$lines[[getLineNumber()]] |>
    stringr::str_c(collapse = "")
  regexMatch <- regexpr(errorPattern, line, perl = TRUE)
  location <- attr(regexMatch, "capture.start")[[2]]
  indent <- "      "
  message <- paste0(indent, line, "\n",
                    indent, stringr::str_c(rep(" ", location - 1),
                                           collapse = ""), "^")
  message
}



parseError <- function(message) {
  message <- paste0("ModSEM \n    ", message, " line ", getLineNumber(),
                    " token: '", getCurrentToken(), "'\n", highlightError())
  stop(message, call. = FALSE)
}



getRegexToken <- function(token) {
  out <- ""
  for (i in 1:nchar(token)) {
    if (grepl("[[:punct:]]", substr(token, i, i))) {
      out <- paste0(out, "\\", substr(token, i, i))
    } else {
      out <- paste0(out, substr(token, i, i))
    }
  }
  return(out)
}
