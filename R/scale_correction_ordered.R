correctScaleOrdered <- function(model.syntax, data = NULL, ordered = NULL, ...) {
  cols.ordered <- colnames(data) %in% ordered
  data[cols.ordered] <- lapply(data[cols.ordered], FUN = as.ordered)

  parTable  <- modsemify(model.syntax)
  outer     <- parTable[parTable$op == "=~", ]
  cfa.cat   <- lavaan::cfa(parTableToSyntax(outer), data = data, ...)
  estimates <- parameterEstimates(cfa.cat)

  thresholds <- estimates[estimates$op == "|", , drop = FALSE]

  if (!NROW(thresholds))
    return(data)

  rescaleVariable <- function(x, name, mean.structure = TRUE) {
    thresholds <- c(-Inf, sort(thresholds[thresholds$lhs == name, "est"]), Inf)
    levels <- levels(x)

    stopif(length(thresholds) != length(levels) + 1,
           "Unexpected number of thresholds and ordered levels!")

    y <- rep(NA, length(x))
    for (i in seq_along(levels)) {
      level <- levels(x)[[i]]

      t0 <- thresholds[i]
      t1 <- thresholds[i + 1]

      mu    <- (dnorm(t0) - dnorm(t1)) / (pnorm(t1) - pnorm(t0))
      z.mid <- qnorm((pnorm(t0) + pnorm(t1)) / 2)

      y[x == level] <- mu
    }

    if (mean.structure) offset <- mean(as.integer(as.factor(x)), na.rm = TRUE)
    else                offset <- 0

    y + offset
  }

  ordered.vars <- unique(thresholds$lhs)

  data.y <- data
  for (var in ordered.vars)
    data.y[[var]] <- rescaleVariable(data[[var]], name = var)

  data.y
}
