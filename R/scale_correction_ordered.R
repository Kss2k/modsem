modsemOrderedScaleCorrection <- function(model.syntax,
                                         data,
                                         method = "lms",
                                         ordered = NULL,
                                         calc.se = TRUE,
                                         R = 5,
                                         N = 1e5, # 100,000
                                         verbose = interactive(),
                                         ...) {
  message("Scale correcting ordinal variables. ",
          "This is an experimental feature!\n",
          "See `help(modsem_da)` for more information.")

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
    y <- (y - mean(y)) / stats::sd(y) # standardize
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
      mu_k <- c(MU[[name]][[i]], mu_i)

      MU[[name]][[i]] <<- mu_k

      k   <- length(mu_k)
      q20 <- floor(k / 5)
      mu  <- mean(mu_k[seq_len(k) > q20]) # drop first 20 percent

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
