modsem_bayes <- function(model.syntax = NULL,
                         data = NULL,
                         chains = 2,
                         iter = 2000,
                         warmup = iter / 2,
                         ordered = NULL,
                         ordered.link = "logit",
                         parameterization = "centered",
                         ...) {
  promptInstallRstan()

  if (!requireNamespace("rstan", quietly = TRUE))
     mod_msg_stop("The 'rstan' package is required to use the `modsem_stan()` function!")

  if (rcs) { # use reliability-correct single items?
    corrected <- relcorr_single_item(
      syntax          = model.syntax,
      data            = data,
      choose          = rcs.choose,
      scale.corrected = rcs.scale.corrected,
      warn.lav        = FALSE
    )

    model.syntax <- corrected$syntax
    data         <- corrected$data
  }

  stopif(is.null(model.syntax) && rcs,
         "`model.syntax` argument is needed when `rcs=TRUE`!")

  if (is.null(compiled.model) || rcs) {
    stopif(is.null(model.syntax),
           "One of `model.syntax` or `compiled.model` has to be provided!")
    # pass ordered through so codegen knows which indicators are ordinal
    compiled.model <- compile_stan_model(model.syntax,
                                         ordered = ordered,
                                         ordered.link = ordered.link,
                                         parameterization = parameterization)
  
  } else {
    # normalize ordered for downstream data
    # prep even when compiled.model is provided
    if (is.null(ordered)) ordered <- character(0)
  }

  lVs     <- compiled.model$info$lVs
  indsLVs <- compiled.model$info$indsLVs
  inds    <- unique(unlist(indsLVs))
  etas    <- compiled.model$info$etas
  deps    <- c(inds, etas)

  # IMPORTANT: pass ordered to the data builder so it supplies INDICATORS_<ind> and K_<ind>
  stan_data <- getStanData(
    compiled.model = compiled.model,
    data = data,
    ordered = ordered
  )

  message("Sampling Stan model...")
  fit <- rstan::sampling(
    object  = compiled.model$stan_model,
    data    = stan_data,
    chains  = chains,
    iter    = iter,
    warmup  = warmup,
    pars    = compiled.model$info$exclude.pars,
    include = FALSE,
    ...
  )

  diagnostics <- rstan::summary(fit)$summary
  fit
}


promptInstallRstan <- function() {
  if (!interactive() || requireNamespace("rstan", quietly = TRUE))
    return(NULL)

  printf("The `rstan` package is needed to use `modsem_stan()`\n")
  printf("Do you want to install it? (y/n)\n")

  choice <- readLines(n = 1L)

  if (!nchar(choice))
    return(invisible(NULL))

  if (tolower(substr(choice, 1L, 1L)) == "y")
    install.packages("rstan")
}
