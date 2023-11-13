
getLines <- function(syntax) {
  # split syntax into individual lines and return as vector
  unlist(strsplit(syntax, "\n"))
}



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
  # remove lines which are list(0) since they only represent comments
  isListZeroLen <- vapply(tokensLines,
                          FUN = function(list) length(list) <= 0,
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

    } else if (grepl("[\\=\\~\\*\\+\\(\\)\\<\\>\\-\\,\\:]", token)) {
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
    stop("incorrect input in getClassLavOp()")
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
    ")"  = stop("Unmatched bracket for )"),
    stop("Unrecognized operator")
  )
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



evalTokens <- function(listTokens,
                       # i = 1,
                       lhs = NULL,
                       op = NULL,
                       rhs = NULL,
                       modifier = NULL) {

  if (1 > length(listTokens)) {
    return(NULL)
  } else if (length(listTokens) == 1) {
    return(listTokens[[1]])
  }

  if (is.null(lhs)) {
    if (class(listTokens[[1]]) == "LavName") {
      # figure out wheter it is an object or a function call
      className <- getClassLavName(listTokens[[1]])

      if (className == "LavObject") {
        lhs <- listTokens[[1]]
        restParsed <- evalTokens(listTokens[-1], #i + 1,
                                 lhs = lhs,
                                 op = op,
                                 rhs = rhs,
                                 modifier = modifier)
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
                                   modifier = modifier)
          return(restParsed)
        }
        return(outputName)
      }


    } else if (class(listTokens[[1]]) == "LavNumeric") {
      lhs <- listTokens[[1]]
      restParsed <- evalTokens(listTokens[-1],
                               lhs = lhs,
                               op = op,
                               rhs = rhs,
                               modifier = modifier)
      return(restParsed)

    } else {
      stop("Expected line to start with an object name or number: ",
           listTokens[[1]])
    }
  } else if (is.null(op)) {
    if (length(listTokens) == 1) {
      return(listTokens[[1]])

    } else if (class(listTokens[[1]]) == "LavOperator") {
      op <- listTokens[[1]]
      class(op) <- getClassLavOp(op)
      restExpression <- listTokens[-1]

      return(evalOp(op, lhs = lhs, rhs = restExpression))

    } else {
      stop("Expected operator after object name:")
    }
  }
}



evalOp <- function(op, lhs, rhs) {
  UseMethod("evalOp")
}



evalOp.LavMeasure <- function(op, lhs, rhs) {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs)

  list(lhs = lhs, op = op, rhs = rest )
}



evalOp.LavPredict <- function(op, lhs, rhs) {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs)

  list(lhs = lhs, op = op, rhs = rest )
}



evalOp.LavCovar <- function(op, lhs, rhs) {
  # op is just here to fetch the method
  rest <- evalTokens(listTokens = rhs)

  list(lhs = lhs, op = op, rhs = rest )
}



evalOp.LavEquals <- function(op, lhs, rhs) {
  # op is just here to fetch the method
  rest <- stringr::str_c(unlist(rhs), collapse = "")

  list(lhs = lhs, op = op, rhs = rest)
}



evalOp.LavLessRight <- function(op, lhs, rhs) {
  # op is just here to fetch the method
  rest <- stringr::str_c(unlist(rhs), collapse = "")

  list(lhs = lhs, op = op, rhs = rest)
}



evalOp.LavLessLeft <- function(op, lhs, rhs) {
  # op is just here to fetch the method
  rest <- stringr::str_c(unlist(rhs), collapse = "")

  list(lhs = lhs, op = op, rhs = rest)
}



evalOp.LavAdd <- function(op, lhs, rhs) {

  # op is just here to fetch method
  if (length(rhs) == 1) {
    return(list(lhs, rhs[[1]]))
  } else if (length(rhs) <= 0 ||
             !(class(rhs[[1]]) %in% c("LavName", "LavNumeric"))) {
    stop("expected name after +")
  }
  rest <- evalTokens(listTokens = rhs)

  c(list(lhs), rest)
}


evalOp.LavModify <- function(op, lhs, rhs) {
  attr(rhs[[1]], "LavMod") <- lhs

  evalTokens(rhs[-1], lhs = rhs[[1]])
}


evalOp.LavInteraction <- function(op, lhs, rhs) {

  rest <- evalTokens(listTokens = rhs)



  # combine lhs and rest into one, inheriting attributes from lhs
  combinedLhsRhs <- structure(paste0(lhs, ":", rest[[1]]),
                              class = class(lhs),
                              LavMod = attr(lhs, "LavMod"))
  if (length(rest) <= 1) {
    return(c(combinedLhsRhs))
  }

  c(combinedLhsRhs, rest[-1])
}



createSyntaxTree <- function(syntax) {
  tokensLines <- tokenizeSyntax(syntax)
  lapply(tokensLines,
         FUN = evalTokens)
}



createParTableBranch <- function(syntaxBranch) {

  # Run a function to check wheter there is nesting in Rhs, which throws an error
  maxDepth(syntaxBranch[["rhs"]])
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
