
modsemParseEnv <- rlang::env(
  syntaxLines = NULL
)


resetModsemParseEnv <- function() {
  modsemParseEnv$syntaxLines <- NULL
}



getCharsLine <- function(line, i = 1) {
  if (is.null(line) || nchar(line) == 0) {
    return(" ")
  } else if (i > nchar(line)) {
    return(NULL)
  }
  rest <- getCharsLine(line, i + 1)
  c(substr(line, i, i),rest)
}


getLines <- function(syntax) {
  # split syntax into individual lines and return as vector
  lines <- strsplit(syntax, "\n|;") |>
    unlist() |>
    as.list() |>
    lapply(getCharsLine)
  lines <- purrr::imap(lines,
                  function(x, lNum)
                    structure(x, lineNumber = lNum))
  lines
}

createTokensLine <- function(line, i = 1,
                             token = NULL, listTokens = list()) {
  # if (i == length(line)) browser()
  if (i > length(line)) {
    return(appendToList(listTokens, token))
  }

  if (length(listTokens) > 0 && is.EqualityOperator(last(listTokens))) {
    token <- buildMathExprToken(line[-(1:i)], pos = i)
    return(appendToList(listTokens, token))
  }


  if (is.null(token)) {
    token <- initializeToken(line[[i]], pos = i, line)
  } else {
    if (fitsToken(token, nextChar = line[[i]])) {
      token <- addCharToken(token, nextChar = line[[i]])

    } else {
      listTokens <- appendToList(listTokens, token)
      token <- initializeToken(line[[i]], pos = i, line)
    }
  }
  if ("LavComment" %in% class(token)) {
      return(listTokens)
    }
  return(createTokensLine(line, i + 1, token, listTokens = listTokens))
}

initializeToken <- function(char, pos, line) {
  if (grepl("#", char)) {
    type  <- "LavComment"
    priority <- 999
  } else if (grepl("\\s+", char)) {
    type <- "LavBlank"
    priority <- 999
  } else if (grepl("[[:alpha:]_.]", char)) {
    type <- "LavName"
    priority <- 10
  } else if (grepl("[\\(\\)]", char)) {
    type <- "LavClosure"
    priority <- 2
  } else if (grepl("[\\=\\~\\*\\+\\<\\>\\-\\,\\:\\^]", char)) {
    type <- "LavOperator"
    priority <- 0
  } else if (grepl("[[:alnum:]]", char)) {
    type <- "LavNumeric"
    priority <- 10
  } else if (grepl('\\"' , char)) {
    type <- "LavString"
    priority <- 10
  } else {
    stop("Unrecognized class of token in line ", attr(line, "lineNumber"),
         " pos ", pos, "\n",
         highlightError(line, pos = pos))
  }
  structure(char,
            pos = pos,
            lineNumber = attr(line, "lineNumber"),
            priority = priority,
            class = c(type, "LavToken"))
}



buildMathExprToken <- function(restLine, pos) {
  token <- stringr::str_c(restLine, collapse = "")
  structure(token,
            pos = pos,
            lineNumber = attr(restLine, "lineNumber"),
            priority = 10,
            class = c("LavMathExpr", "LavToken", "data.frame"))

}



fitsToken <- function(token, nextChar) {
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
  if (length(nextChar) != 1) {
    stop("Wrong length of nextChar", nextChar)
  }
  grepl("\\s+", nextChar)
}



is.EqualityOperator <- function(token) {
  switch(getTokenString(token),
         "==" = TRUE,
         "<" = TRUE,
         ">" = TRUE,
         FALSE)
  }



fitsToken.LavClosure <- function(token, nextChar) {
  FALSE
}



fitsToken.LavNumeric <- function(token, nextChar) {
  if (length(nextChar) != 1) {
    stop("Wrong length of nextChar", nextChar)
  }
  grepl("[[:digit:].]", nextChar)
}



addCharToken <- function(token, nextChar) {
  if (length(nextChar) != 1) {
    stop("Wrong length of nextChar ", nextChar)
  }

  out <- paste0(token, nextChar)
  attributes(out) <- attributes(token)
  out
}



fitsToken.LavComment <- function(token, nextChar, pos) {
  TRUE
}



assignSubClass <- function(token) {
  UseMethod("assignSubClass")
}



assignSubClass.LavOperator <- function(token) {
  switch (getTokenString(token),
          "=~" = {subClass <- "LavMeasure";     priority <- 0},
          "~"  = {subClass <- "LavPredict";     priority <- 0},
          "~~" = {subClass <- "LavCovar";       priority <- 0},
          "+"  = {subClass <- "LavAdd";         priority <- 1},
          "*"  = {subClass <- "LavModify";      priority <- 2},
          "<"  = {subClass <- "LavLessLeft";    priority <- 0},
          ">"  = {subClass <- "LavLessRight";   priority <- 0},
          "==" = {subClass <- "LavEquals";      priority <- 0},
          ":"  = {subClass <- "LavInteraction"; priority <- 2},
          "," =  {subClass <- "LavSeperator"; priority <- 0},
          stop("Unrecognized operator: ", highlightErrorToken(token))
  )
  structure(token,
            class = c(subClass, class(token)),
            priority = priority)
}



mergeTokensToString <- function(listTokens) {
  vapply(listTokens,
         FUN.VALUE = character(1L),
         FUN = getTokenString) |>
    stringr::str_c(collapse = "")
}



assignSubClass.LavClosure <- function(token) {
  switch (getTokenString(token),
          "("  = {subClass <- "LeftBracket";    priority <- 3},
          ")"  = {subClass <- "RightBracket";   priority <- 3},
          stop("Unrecognized operator: ", token)
  )
  structure(token,
            class = c(subClass, class(token)),
            priority = priority)
}



assignSubClass.LavName <- function(token) {
  if (grepl("\\($", getTokenString(token))) {
    subClass <- "LavFunction"
    priority <- 3
  } else {
    subClass <- "LavObject"
    priority <- 3
  }
  structure(token,
            class = c(subClass, class(token)),
            priority = priority)
}



assignSubClass.LavToken <- function(token) {
  token
}



appendToList <- function(list, elem) {
  list[[length(list) + 1]] <- elem
  list
}



prioritizeTokens <- function(listTokens, i = 1, brackets = list(),
                             nLeftBrackets = 0) {
  if (is.null(listTokens) || i > length(listTokens)) {
    if (nLeftBrackets != 0) {
        stop("Unmatched left bracket", highlightErrorToken(brackets[[1]]))
    }
    return(listTokens)
  }
  else if (is.RightClosure(listTokens[[i]])) {
    brackets <- brackets[-length(brackets)]
    nLeftBrackets <- nLeftBrackets - 1
    if (nLeftBrackets < 0) {
      stop("Unmatched right bracket ", highlightErrorToken(listTokens[[i]]))
    }
  }
  getTokenPriority(listTokens[[i]]) <-
    getTokenPriority(listTokens[[i]]) + nLeftBrackets*10

  if (is.LeftClosure(listTokens[[i]])) {
    brackets <- appendToList(brackets, listTokens[[i]])
    nLeftBrackets <- nLeftBrackets + 1
  }
  prioritizeTokens(listTokens, i + 1, brackets = brackets,
                   nLeftBrackets = nLeftBrackets)
}



removeLavBlankLine <- function(line, removeComments = TRUE) {
  if (is.null(line) || length(line) == 0) {
    return(line)
  }
  isBlankOrComment <- vapply(line,
                             FUN.VALUE = logical(1L),
                             FUN = function(token) is.LavBlankOrComment(token))
  line[!isBlankOrComment]
}



tokenizeSyntax <- function(syntax, optimize = TRUE) {
  resetModsemParseEnv()
  if (!is.character(syntax)) {
    stop("Syntax must be a string")
  }
  lines <- getLines(syntax)
  modsemParseEnv$syntaxLines <- lines

  tokenizedLines <- lines |>
    lapply(createTokensLine) |>
    lapply(removeLavBlankLine) |>
    lapply(FUN = function(tokens)
                   lapply(tokens, assignSubClass))
  modsemParseEnv$syntaxLines <- tokenizedLines
  tokenizedLines <- tokenizedLines |>
    lapply(prioritizeTokens)

  isEmpty <- vapply(tokenizedLines,
                    FUN.VALUE = logical(1L),
                    FUN = function(line) is.null(line) || length(line) == 0)

  tokenizedLines[!isEmpty]
}



highlightError <- function(line, pos) {
  line <- stringr::str_c(line, collapse = "")
  indent <- "      "
  message <- paste0(indent, line, "\n",
                    indent, stringr::str_c(rep(" ", pos - 1),
                                           collapse = ""), "^")
  message
}



highlightErrorToken <- function(token) {
  lineNumber <- attr(token, "lineNumber")
  line <- modsemParseEnv$syntaxLines[[lineNumber]]
  if (is.list(line)) {
    line <- vapply(line,
                   FUN.VALUE = character(1L),
                   FUN = getTokenString)
  }
  line <- stringr::str_c(line, collapse = "")
  pos <- getTokenPosition(token)
  indent <- "      "
  message <- paste0("\n", indent, line, "\n",
                    indent, stringr::str_c(rep(" ", pos - 1),
                                           collapse = ""), "^")
  message
}


getTokenString <- function(token) {
  token
}




getTokenPriority <- function(token) {
  attr(token, "priority")
}


`getTokenPriority<-` <- function(token, value) {
  attr(token, "priority") <- value
  token
}


getTokenPosition <- function(token) {
  attr(token, "pos")
}


is.LavOperator <- function(token) {
  "LavOperator" %in% class(token)
}



is.LavBlankOrComment <- function(token) {
  "LavBlank" %in% class(token) | "LavComment" %in% class(token)
}



is.LeftClosure <- function(token) {
  "LeftBracket" %in% class(token) | "LavFunction" %in% class(token)
}



is.RightClosure <- function(token) {
  "RightBracket" %in% class(token)
}



is.LavClosure <- function(token) {
  "LavClosure" %in% class(token)
}



is.LavOperator <- function(token) {
  "LavOperator" %in% class(token)
}



is.firstClassOperator <- function(token) {
  switch(getTokenString(token),
         "=~" = TRUE,
         "~"  = TRUE,
         "~~" = TRUE,
         "<"  = TRUE,
         ">"  = TRUE,
         "==" = TRUE,
         FALSE)
}


# LavToken        <- setClass("LavToken")
# setMethod("assignSubClass", "LavToken", assignSubClass.LavToken)
# LavName         <- setClass("LavName")
# setMethod("assignSubClass", "LavName", assignSubClass.LavName)
# LavObject       <- setClass("LavObject")
# LavFunciton     <- setClass("LavFunction")
# LavOperator     <- setClass("LavOperator")
# setMethod("assignSubClass", "LavOperator", assignSubClass.LavOperator)
# LavClosure      <- setClass("LavClosure")
# setMethod("assignSubClass", "LavClosure", assignSubClass.LavClosure)
# LeftBracket     <- setClass("LeftBracket")
# RightBracket    <- setClass("RightBracket")
# LavComment      <- setClass("LavComment")
# LavBlank        <- setClass("LavBlank")
# LavNumeric      <- setClass("LavNumeric")
# LavString       <- setClass("LavString")
# LavMeasure      <- setClass("LavMeasure")
# LavPredict      <- setClass("LavPredict")
# LavCovar        <- setClass("LavCovar")
# LavAdd          <- setClass("LavAdd")
# LavModify       <- setClass("LavModify")
# LavLessLeft     <- setClass("LavLessLeft")
# LavLessRight    <- setClass("LavLessRight")
# LavEquals       <- setClass("LavEquals")
# LavInteraction  <- setClass("LavInteraction")
# LavSeperator    <- setClass("LavSeperator")


