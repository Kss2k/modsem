#' @export
modsem_gsem <- function(model.syntax = NULL,
                        data = NULL,
                        group = NULL,
                        method = "gsem",
                        verbose = NULL,
                        optimize = NULL,
                        nodes = NULL,
                        missing = NULL,
                        convergence.abs = NULL,
                        convergence.rel = NULL,
                        optimizer = NULL,
                        center.data = NULL,
                        standardize.data = NULL,
                        standardize.out = NULL,
                        standardize = NULL,
                        mean.observed = NULL,
                        cov.syntax = NULL,
                        double = NULL,
                        calc.se = NULL,
                        FIM = NULL,
                        EFIM.S = NULL,
                        OFIM.hessian = NULL,
                        EFIM.parametric = NULL,
                        robust.se = NULL,
                        R.max = NULL,
                        max.iter = NULL,
                        max.step = NULL,
                        start = NULL,
                        epsilon = NULL,
                        quad.range = NULL,
                        adaptive.quad = NULL,
                        adaptive.frequency = NULL,
                        adaptive.quad.tol = NULL,
                        n.threads = NULL,
                        algorithm = NULL,
                        em.control = NULL,
                        ordered = NULL,
                        ordered.probit.correction = FALSE,
                        cluster = NULL,
                        cr1s = FALSE,
                        rcs = FALSE,
                        rcs.choose = NULL,
                        rcs.scale.corrected = TRUE,
                        orthogonal.x = NULL,
                        orthogonal.y = NULL,
                        auto.fix.first = NULL,
                        auto.fix.single = NULL,
                        auto.split.syntax = NULL,
                        ...) {
  if (is.null(model.syntax)) {
    stop2("No model.syntax provided")
  } else if (!is.character(model.syntax)) {
    stop2("The provided model syntax is not a string!")
  } else if (length(model.syntax) > 1) {
    stop2("The provided model syntax is not of length 1")
  }

  if (is.null(data)) {
    stop2("No data provided")
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  args <-
    getMethodSettingsDA(method,
      args =
        list(
          verbose            = verbose,
          optimize           = optimize,
          nodes              = nodes,
          convergence.abs    = convergence.abs,
          convergence.rel    = convergence.rel,
          optimizer          = optimizer,
          center.data        = center.data,
          standardize.data   = standardize.data,
          standardize.out    = standardize.out,
          standardize        = standardize,
          mean.observed      = mean.observed,
          double             = double,
          calc.se            = calc.se,
          FIM                = FIM,
          EFIM.S             = EFIM.S,
          OFIM.hessian       = OFIM.hessian,
          EFIM.parametric    = EFIM.parametric,
          robust.se          = robust.se,
          R.max              = R.max,
          max.iter           = max.iter,
          max.step           = max.step,
          epsilon            = epsilon,
          quad.range         = quad.range,
          adaptive.quad      = adaptive.quad,
          adaptive.frequency = adaptive.frequency,
          adaptive.quad.tol  = adaptive.quad.tol,
          n.threads          = n.threads,
          algorithm          = algorithm,
          em.control         = em.control,
          missing            = missing,
          orthogonal.x       = orthogonal.x,
          orthogonal.y       = orthogonal.y,
          auto.fix.first     = auto.fix.first,
          auto.fix.single    = auto.fix.single,
          auto.split.syntax  = auto.split.syntax,
          cr1s               = cr1s,
          group              = group
        )
    )

  cont.cols <- setdiff(colnames(data), c(cluster, group))

  if (args$center.data)
    data[cont.cols] <- lapply(data[cont.cols], FUN = centerIfNumeric, scaleFactor = FALSE)

  if (args$standardize.data)
    data[cont.cols] <- lapply(data[cont.cols], FUN = scaleIfNumeric, scaleFactor = FALSE)

  group.info <- getGroupInfo(
    model.syntax       = model.syntax,
    cov.syntax         = NULL,
    data               = data,
    group              = group,
    auto.split.syntax  = args$auto.split.syntax
  )

  model <- specifyModelGsem(
    group.info         = group.info,
    method             = method,
    m                  = args$nodes,
    mean.observed      = args$mean.observed,
    double             = args$double,
    quad.range         = args$quad.range,
    adaptive.quad      = args$adaptive.quad,
    adaptive.frequency = args$adaptive.frequency,
    missing            = args$missing,
    orthogonal.x       = args$orthogonal.x,
    orthogonal.y       = args$orthogonal.y,
    auto.fix.first     = args$auto.fix.first,
    auto.fix.single    = args$auto.fix.single,
    cluster            = cluster
  )

  if (args$optimize) {
    model <- tryCatch({
      .optimize <- purrr::quietly(optimizeStartingParamsGsem)
      #.optimize <- optimizeStartingParamsDA
      result    <- .optimize(model, args = args, group = group, engine = "sam")
      warnings  <- result$warnings

      if (length(warnings)) {
        fwarnings <- paste0(
          paste0(seq_along(warnings), ". ", warnings),
          collapse = "\n"
        )

        warning2("warning when optimizing starting parameters:\n", fwarnings)
      }

      result$result

    }, error = function(e) {
      warning2("unable to optimize starting parameters:\n", e)
      model
    })
  }

  if (!is.null(start)) {
    checkStartingParams(start, model = model) # throws error if somethings wrong
    model$theta <- start
  }

  # We want to limit the number of threads available to OpenBLAS.
  # Depending on the OpenBLAS version, it might not be compatible with
  # OpenMP. If `n.blas > 1L` you might end up getting this message:
  #> OpenBLAS Warning : Detect OpenMP Loop and this application may hang.
  #>                    Please rebuild the library with USE_OPENMP=1 option.
  # We don't want to restrict OpenBLAS in any other setttings in other settings,
  # e.g., lavaan::sem, so we reset after the model has been estimated.
  setThreads(n = args$n.threads, n.blas = 1L)
  on.exit(resetThreads()) # clean up at end of function

  est <- tryCatch(emGsem(model,
      verbose           = args$verbose,
      convergence.abs   = args$convergence.abs,
      convergence.rel   = args$convergence.rel,
      calc.se           = args$calc.se,
      FIM               = args$FIM,
      EFIM.S            = args$EFIM.S,
      OFIM.hessian      = args$OFIM.hessian,
      EFIM.parametric   = args$EFIM.parametric,
      robust.se         = args$robust.se,
      max.iter          = args$max.iter,
      max.step          = args$max.step,
      epsilon           = args$epsilon,
      optimizer         = args$optimizer,
      R.max             = args$R.max,
      em.control        = args$em.control,
      algorithm         = args$algorithm,
      adaptive.quad     = args$adaptive.quad,
      quad.range        = args$quad.range,
      adaptive.quad.tol = args$adaptive.quad.tol,
      nodes             = args$nodes,
      cr1s              = args$cr1s,
      ...
  ),
  error = function(e) {
    if (args$verbose) cat("\n")
    message <- paste0("modsem [%s]: Model estimation failed!\n",
                      "Message: %s")
    stop2(sprintf(message, method, e$message))
  })

  # Finalize the model object
  # Expected means and covariances
  est$expected.matrices <- tryCatch(
    calcExpectedMatricesDA(
      parTable = est$parTable,
      xis  = getXisModelDA(model$models[[1L]]), # taking both the main model and cov model into account
      etas = getEtasModelDA(model$models[[1L]])  # taking both the main model and cov model into account
    ),
    error = function(e) {
      warning2("Failed to calculate expected matrices: ", e$message)
      NULL
    })

  # Arguments
  est$args <- args
  class(est) <- c("modsem_da", "modsem")

  # Check the results
  postCheckModel(est)

  # Return
  if (args$standardize.out) standardize_model(est) else est
}
