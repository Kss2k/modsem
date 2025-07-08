twostep_da <- function(syntax, data, method = "lms", ...) {
  data <- as.data.frame(data)
  parTable      <- modsemify(syntax)
  parTableOuter <- parTable[parTable$op == "=~", ]
  parTableInner <- parTable[parTable$op != "=~", ]

  etas <- getEtas(parTable, isLV = FALSE)
  syntaxCFA <- parTableToSyntax(parTableOuter)
  cfa       <- lavaan::cfa(syntaxCFA, data = data, meanstructure = TRUE)

  parTableCFA <- lavaan::parameterEstimates(cfa)

  # Remove (co-) variances for etas
  parTableCFA <- parTableCFA[!(parTableCFA$op == "~~" &
                               parTableCFA$lhs %in% etas |
                               parTableCFA$rhs %in% etas), ]

  parTableOuterFixed <- data.frame(
    lhs = parTableCFA$lhs,
    op = parTableCFA$op,
    rhs = parTableCFA$rhs,
    mod = as.character(parTableCFA$est)
  )

  if ("label" %in% names(parTableCFA)) {
    labelledParams <- parTableCFA[parTableCFA$label != "", ]
    parTableCustom <- data.frame(
      lhs = labelledParams$label,
      op = ":=",
      rhs = as.character(labelledParams$est),
      mod = ""
    )
  } else parTableCustom <- NULL
 
  etaIntercepts <- data.frame(
    lhs = etas,
    op = "~1",
    rhs = "",
    mod = ""
  )

  parTableTwoStep <- rbind(
    parTableCustom,
    parTableOuterFixed, 
    parTableInner,
    etaIntercepts
  )

  syntaxTwoStep <- parTableToSyntax(parTableTwoStep)

  twostep.spec <- list(
    syntax        = syntaxTwoStep,
    cfa           = cfa,
    parTable      = parTableTwoStep,
    parTableOuter = parTableOuterFixed,
    parTableInner = parTableInner
  )

  browser()
  fit <- modsem_da(
    model.syntax = syntaxTwoStep,
    data = data,
    method = method,
    ...
  )

  browser()
}
