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

      y[x == level] <- z.mid
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


modsemOrderedScaleCorrection <- function(model.syntax,
                                         data,
                                         method = "lms",
                                         ordered = NULL,
                                         calc.se = TRUE,
                                         R = 5,
                                         N = 1e5, # 100,000
                                         verbose = interactive(),
                                         ...) {
  if (is.null(verbose))
    verbose <- TRUE # default

  cols <- colnames(data)
  cols.ordered <- cols[cols %in% ordered | sapply(data, is.ordered)]

  MU <- stats::setNames(vector("list", length(cols.ordered)), nm = cols.ordered)

  for (col in cols.ordered) {
    data[[col]] <- as.integer(as.ordered(data[[col]]))
    MU[[col]]   <- vector("list", length(unique(data[[col]])))
  }

  rescaleOrderedVariable <- function(name, data, sim.ov, k = 200, eps = 1e-3) {
    n <- NROW(data)
    N <- NROW(sim.ov)

    x <- as.integer(as.ordered(data[[name]]))
    y <- sim.ov[, name]
    y <- (y - mean(y)) / sd(y) # standardize
    t <- c(-Inf, seq(min(y) + eps, max(y) - eps, length.out = k), Inf)
    z <- cut(y, breaks = t, ordered_result = TRUE)

    exp.densities <- table(z) / N
    exp.cdf <- cumsum(exp.densities)

    obs.densities <- table(x) / n
    obs.cdf <- cumsum(obs.densities)

    z              <- as.integer(z) # more convenient in the loop
    names(exp.cdf) <- seq_len(length(exp.cdf))

    out <- rep(NA, n)
    for (i in seq_along(obs.cdf)) {
      quantile  <- obs.cdf[i]
      match.idx <- exp.cdf <= quantile
      match.cdf <- exp.cdf[match.idx]
      exp.cdf   <- exp.cdf[!match.idx]

      match.values <- as.integer(names(match.cdf))
      lower <- match.values[1L]
      upper <- last(match.values)

      mu_i <- mean(y[z >= lower & z <= upper])
      MU[[name]][[i]] <<- c(MU[[name]][[i]], mu_i)
      mu <- mean(MU[[name]][[i]])

      out[!is.na(x) & x == i] <- mu 
    }

    out
  }

  rescaleOrderedData <- function(data, sim.ov) {
    data.y <- data

    for (col in cols.ordered)
      data.y[[col]] <- rescaleOrderedVariable(name = col, data = data, 
                                              sim.ov = sim.ov)

    data.y
  }


  ERROR <- \(e) {warning2(e, immediate. = FALSE); NULL}

  data.x <- data
  data.y <- data
  for (r in seq_len(R)) {
    printedLines <- utils::capture.output(split = TRUE, {
      if (verbose) printf("Bootstrapping scale-correction %d/%d...\n", r, R)

      fit.naive <- modsem(
        model.syntax = model.syntax,
        data         = data.y,
        method       = method,
        calc.se      = FALSE,
        verbose      = verbose,
        ...
      )

      sim <- simulateDataParTable(
        parTable = parameter_estimates(fit.naive),
        N        = N
      )

      data.y <- rescaleOrderedData(data = data.x, sim.ov = sim$oV)
    })

    nprinted <- length(printedLines)
    eraseConsoleLines(nprinted)
  }

  modsem(
    model.syntax = model.syntax, 
    method       = method,
    data         = data.y, 
    calc.se      = calc.se,
    verbose      = verbose,
    ...
  )
}
