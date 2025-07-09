twostep_da <- function(syntax, data, method = "lms", zero.tol = 1e-12,
                       fix.cov.xis = TRUE, ...) {
  data <- as.data.frame(data)
  parTable      <- modsemify(syntax)
  parTableOuter <- parTable[parTable$op == "=~", ]
  parTableInner <- parTable[parTable$op != "=~", ]

  etas <- getEtas(parTable, isLV = FALSE)
  lvs  <- getLVs(parTable)
  syntaxCFA <- parTableToSyntax(parTableOuter)
  cfa       <- lavaan::cfa(syntaxCFA, data = data, meanstructure = TRUE)

  parTableCFA <- lavaan::parameterEstimates(cfa)

  fixCovVars <- if (fix.cov.xis) etas else lvs
  # Remove (co-) variances for etas (and optionally xis)
  parTableCFA <- parTableCFA[!(parTableCFA$op == "~~" &
                               parTableCFA$lhs %in% fixCovVars |
                               parTableCFA$rhs %in% fixCovVars), ]

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

  fit <- modsem_da(
    model.syntax = syntaxTwoStep,
    data = data,
    method = method,
    ...
  )

  # prep partable from CFA
  parTableCFA <- rename(parTableCFA, se = "se.cfa")
  parTableCFA <- parTableCFA[c("lhs", "op", "rhs", "se.cfa")]
  parTableCFA$se.cfa[parTableCFA$se.cfa <= zero.tol] <- NA
  estParTable <- parameter_estimates(fit)

  # Add reversed order covariances, to ensure matching
  covrows <- parTableCFA$op == "~~" & parTableCFA$rhs != parTableCFA$lhs
  parTableCFACov     <- parTableCFA
  parTableCFACov$rhs <- parTableCFA$lhs
  parTableCFACov$lhs <- parTableCFA$rhs
  parTableCFA <- rbind(
    parTableCFA, parTableCFACov[covrows, ]
  )

  parTable.out <- leftJoin(
    x = estParTable,
    y = parTableCFA,
    by = c("lhs", "op", "rhs")
  )

  parTable.out$std.error <- ifelse(
    is.na(parTable.out$std.error) & !is.na(parTable.out$se.cfa),
    yes = parTable.out$se.cfa,
    no  = parTable.out$std.error
  )
  
  parTable.out$z.value  <- parTable.out$est / parTable.out$std.error
  parTable.out$p.value  <- 2 * stats::pnorm(-abs(parTable.out$z.value))
  parTable.out$ci.lower <- parTable.out$est - CI_WIDTH * parTable.out$std.error
  parTable.out$ci.upper <- parTable.out$est + CI_WIDTH * parTable.out$std.error

  fit$parTable <- parTable.out

  fit
}
