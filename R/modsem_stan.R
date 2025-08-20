#' Interaction between latent variables using Markow-Chain Monte-Carlo sampling via \code{STAN}
#'
#' @param model.syntax \code{lavaan} syntax. Overwritten by \code{compiled_model}.
#'
#' @param data A dataframe with observed variables used in the model.
#'
#' @param compiled_model Compiled model from \code{\link{compile_stan_model}}. Saves time if the
#'   same \code{model.syntax} has to be reused multiple times.
#'
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4.
#'
#' @param iter A positive integer specifying the number of iterations for
#' each chain (including warmup). The default is 2000.
#'
#' @param warmup A positive integer specifying the number of warmup (aka
#' burnin) iterations per chain. If step-size adaptation is on
#' (which it is by default), this also controls the number of
#' iterations for which adaptation is run (and hence these
#' warmup samples should not be used for inference). The number
#' of warmup iterations should be smaller than \code{iter} and the
#' default is \code{iter/2}.
#'
#' @param ... Arguments passed to \code{stan::sampling}.
#' @export
modsem_stan <- function(model.syntax = NULL,
                        data = NULL,
                        compiled_model = NULL,
                        chains = 2,
                        iter = 2000,
                        warmup = iter / 2,
                        ordered = NULL,
                        ...) {
  if (is.null(compiled_model)) {
    stopif(is.null(model.syntax),
           "One of `model.syntax` or `compiled_model` has to be provided!")
    # pass ordered through so codegen knows which indicators are ordinal
    compiled_model <- compile_stan_model(model.syntax, ordered = ordered)
  } else {
    # normalize ordered for downstream data prep even when compiled_model is provided
    if (is.null(ordered)) ordered <- character(0)
  }

  lVs     <- compiled_model$info$lVs
  indsLVs <- compiled_model$info$indsLVs
  inds    <- unique(unlist(indsLVs))
  etas    <- compiled_model$info$etas
  deps    <- c(inds, etas)

  # IMPORTANT: pass ordered to the data builder so it supplies INDICATORS_<ind> and K_<ind>
  stan_data <- getStanData(compiled_model = compiled_model, data = data, ordered = ordered)

  message("Sampling Stan model...")
  fit <- rstan::sampling(object  = compiled_model$stan_model,
                         data    = stan_data,
                         chains  = chains,
                         iter    = iter,
                         warmup  = warmup,
                         pars    = compiled_model$info$exclude.pars,
                         include = FALSE,
                         ...)

  diagnostics <- rstan::summary(fit)$summary

  samples     <- as.matrix(fit)
  # DO NOT drop indexed params globally; we need cutpoints like x2__CUTPOINTS[1].
  # Vectors of length N (e.g., latent time-series) are already excluded via `exclude.pars`.

  # Map Stan parameter labels -> lavaan-like labels
  cleanPars <- \(pars) stringr::str_replace_all(pars, STAN_OPERATOR_LABELS)

  # We will build lhs/op/rhs in two passes:
  #  (A) cutpoints:  <ind>__CUTPOINTS[i]  ->  lhs=<ind>, op="|", rhs=paste0("t", i)
  #  (B) everything else: use existing operator mapping

  par_names_raw <- colnames(samples)

  is_cut <- grepl("__CUTPOINTS\\[[0-9]+\\]$", par_names_raw)
  cut_raw <- par_names_raw[is_cut]

  # Extract lhs (indicator) and threshold index
  cut_lhs <- sub("__CUTPOINTS\\[[0-9]+\\]$", "", cut_raw)
  cut_k   <- as.integer(sub("^.*__CUTPOINTS\\[([0-9]+)\\]$", "\\1", cut_raw))
  cut_rhs <- paste0("t", cut_k)
  cut_op  <- rep("|", length(cut_raw))

  # Non-cutpoint parameters
  noncut_raw <- par_names_raw[!is_cut]
  noncut_clean <- cleanPars(noncut_raw)

  OP <- "~~|=~|~1|~"  # lavaan operators we already produce
  noncut_op <- stringr::str_extract(noncut_clean, pattern = OP)
  noncut_op[is.na(noncut_op)] <- ":="  # fall-back for any derived labels

  lr <- stringr::str_split_fixed(noncut_clean, pattern = OP, n = 2)
  noncut_lhs <- lr[, 1]
  noncut_rhs <- lr[, 2]

  # Square residual SDs (variances) for deps (as before)
  isSD <- noncut_lhs == noncut_rhs & noncut_op == "~~" & (noncut_lhs %in% deps | noncut_rhs %in% deps)
  if (any(isSD)) samples[, match(noncut_raw[isSD], par_names_raw)] <- samples[, match(noncut_raw[isSD], par_names_raw)]^2

  # Combine back the label pieces for all params (cutpoints first or interleaved – order does not matter)
  all_lhs <- c(cut_lhs, noncut_lhs)
  all_op  <- c(cut_op,  noncut_op)
  all_rhs <- c(cut_rhs, noncut_rhs)

  # Reorder samples consistently with the combined vectors
  samples <- samples[, c(cut_raw, noncut_raw), drop = FALSE]

  # Summaries
  coefs <- apply(samples, MARGIN = 2, FUN = mean)
  vcov  <- cov(samples)
  diagnostics <- diagnostics[colnames(samples), , drop = FALSE]
  rhat  <- diagnostics[, "Rhat"]
  neff  <- diagnostics[, "n_eff"]

  # Build parTable (lavaan-like), including thresholds
  pars_clean_for_table <- cleanPars(colnames(samples))  # human-friendly labels where relevant

  se <- sqrt(diag(vcov))
  parTable <- data.frame(
    lhs = all_lhs, op = all_op, rhs = all_rhs,
    est = coefs, std.error = se,
    z.value = coefs / se, p.value = 2 * stats::pnorm(-abs(coefs / se)),
    ci.lower = coefs - CI_WIDTH * se,
    ci.upper = coefs + CI_WIDTH * se,
    R.hat = rhat, n.eff = neff,
    row.names = NULL
  )

  # Add fixed first-loading rows (as before)
  parTable <- rbind(
    data.frame(lhs = lVs, op  = "=~",
               rhs = vapply(lVs, FUN = \(lV) indsLVs[[lV]][[1L]], FUN.VALUE = character(1)),
               est = 1, std.error = NA, z.value = NA, p.value = NA,
               ci.lower = NA, ci.upper = NA, R.hat = NA, n.eff = NA),
    parTable
  )

  # Nice ordering: put thresholds with other measurement items if you like; otherwise keep default
  parTable <- parTable[order(parTable$op), ]
  rownames(parTable) <- NULL
  parTable <- modsemParTable(parTable[parTable$op != ":=", , drop = FALSE])

  out <- list(
    fit = fit,
    parTable = parTable,
    coefs = coefs,
    vcov = vcov,
    samples = samples
  )
  class(out) <- "modsem_stan"
  out
}


#' @export
summary.modsem_stan <- function(object, ...) {
  parTable <- object$parTable
  parTable$n.eff <- as.character(round(parTable$n.eff)) # print as integer, not float
  summarize_partable(parTable)
}


#' @export
print.modsem_stan <- function(x, ...) {
  print(x$parTable)
}


#' @export
parameter_estimates.modsem_stan <- function(object, ...) {
  object$parTable
}


#' @export
standardized_estimates.modsem_stan <- function(object,
                                               monte.carlo = FALSE,
                                               mc.reps = 10000,
                                               tolerance.zero = 1e-10,
                                               ...) {
  stdSolution <- standardizedSolutionCOEFS(
    object,
    monte.carlo = monte.carlo,
    mc.reps = mc.reps,
    tolerance.zero = tolerance.zero,
    ...
  )

  stdSolution$parTable
}


#' @export
#' @importFrom stats vcov
vcov.modsem_stan <- function(object, ...) {
  object$vcov
}


#' @export
#' @importFrom stats coef
coef.modsem_stan <- function(object, ...) {
  object$coefs
}
