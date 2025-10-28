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
#' @param ordered Ordered (i.e., ordinal) variables.
#'
#' @param ordered.link Link function to be used for ordered variables.
#'
#' @param parameterization What parameterization to use? Default is \code{"centered"}.
#'
#' @param rcs Should latent variable indicators be replaced with reliability-corrected
#'   single item indicators instead? See \code{\link{relcorr_single_item}}.
#'
#' @param rcs.choose Which latent variables should get their indicators replaced with
#'   reliability-corrected single items? It is passed to \code{\link{relcorr_single_item}}
#'   as the \code{choose} argument.
#'
#' @param rcs.scale.corrected Should reliability-corrected items be scale-corrected? If \code{TRUE}
#'   reliability-corrected single items are corrected for differences in factor loadings between
#'   the items. Default is \code{TRUE}.
#'
#' @param ... Arguments passed to \code{stan::sampling}.
#'
#' @export
modsem_stan <- function(model.syntax = NULL,
                        data = NULL,
                        compiled_model = NULL,
                        chains = 2,
                        iter = 2000,
                        warmup = iter / 2,
                        ordered = NULL,
                        ordered.link = "logit",
                        parameterization = "centered",
                        rcs = FALSE,
                        rcs.choose = NULL,
                        rcs.scale.corrected = TRUE,
                        ...) {
  if (!requireNamespace("rstan", quietly = TRUE))
     stop2("The 'rstan' package is required to use the `modsem_stan()` function!")

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

  if (is.null(compiled_model) || rcs) {
    stopif(is.null(model.syntax),
           "One of `model.syntax` or `compiled_model` has to be provided!")
    # pass ordered through so codegen knows which indicators are ordinal
    compiled_model <- compile_stan_model(model.syntax,
                                         ordered = ordered,
                                         ordered.link = ordered.link,
                                         parameterization = parameterization)
  
  } else {
    # normalize ordered for downstream data
    # prep even when compiled_model is provided
    if (is.null(ordered)) ordered <- character(0)
  }

  lVs     <- compiled_model$info$lVs
  indsLVs <- compiled_model$info$indsLVs
  inds    <- unique(unlist(indsLVs))
  etas    <- compiled_model$info$etas
  deps    <- c(inds, etas)

  # IMPORTANT: pass ordered to the data builder so it supplies INDICATORS_<ind> and K_<ind>
  stan_data <- getStanData(compiled_model = compiled_model,
                           data = data,
                           ordered = ordered)

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
  cutRaw <- par_names_raw[is_cut]

  # Extract lhs (indicator) and threshold index
  cutLhs <- sub("__CUTPOINTS\\[[0-9]+\\]$", "", cutRaw)
  cut_k   <- as.integer(sub("^.*__CUTPOINTS\\[([0-9]+)\\]$", "\\1", cutRaw))
  cutRhs <- if (length(cut_k)) paste0("t", cut_k) else NULL
  cutOp  <- rep("|", length(cutRaw))

  # Non-cutpoint parameters
  noncutRaw <- par_names_raw[!is_cut]
  noncutClean <- cleanPars(noncutRaw)

  OP <- "~~~|~~|=~|~1|~"  # lavaan operators we already produce
  noncutOp <- stringr::str_extract(noncutClean, pattern = OP)
  noncutOp[is.na(noncutOp)] <- ":="

  lr <- stringr::str_split_fixed(noncutClean, pattern = OP, n = 2)
  noncutLhs <- lr[, 1]
  noncutRhs <- lr[, 2]

  # Square residual SDs (variances) for deps (as before)
  # isSD <- noncutLhs == noncutRhs & noncutOp == "~~~"
  # if (any(isSD)) samples[, match(noncutRaw[isSD], par_names_raw)] <- samples[, match(noncutRaw[isSD], par_names_raw)]^2

  # Combine back the label pieces for all params (cutpoints first or interleaved – order does not matter)
  allLhs <- c(cutLhs, noncutLhs)
  allOp  <- c(cutOp,  noncutOp)
  allRhs <- c(cutRhs, noncutRhs)
  allOp[allOp == "~~~"] <- "~~"

  # Reorder samples consistently with the combined vectors
  samples <- samples[, c(cutRaw, noncutRaw), drop = FALSE]

  # Remove :=
  keep   <- allOp != ":="
  lhs <- allLhs[keep]
  op  <- allOp[keep]
  rhs <- allRhs[keep]

  samples         <- samples[, keep, drop = FALSE]
  namesSamplesRaw <- colnames(samples)
  diagnostics     <- diagnostics[namesSamplesRaw, , drop = FALSE]

  # Detect duplicated residual covariance columns (order-insensitive)
  sym_mask <- op == "~~" & lhs != rhs
  dedup_key <- paste(
    op,
    ifelse(sym_mask, pmin(lhs, rhs), lhs),
    ifelse(sym_mask, pmax(lhs, rhs), rhs),
    sep = "::"
  )

  dup_idx <- duplicated(dedup_key)
  if (any(dup_idx)) {
    samples     <- samples[, !dup_idx, drop = FALSE]
    lhs         <- lhs[!dup_idx]
    op          <- op[!dup_idx]
    rhs         <- rhs[!dup_idx]
    diagnostics <- diagnostics[!dup_idx, , drop = FALSE]
    dedup_key   <- dedup_key[!dup_idx]
  }

  # Align residual covariance orientation with the compiled parTable, if possible
  sym_mask <- op == "~~" & lhs != rhs
  if (any(sym_mask)) {
    par_input <- compiled_model$info$parTable
    input_sym <- par_input$op == "~~" & par_input$lhs != par_input$rhs
    if (any(input_sym)) {
      input_key <- paste(
        par_input$op[input_sym],
        pmin(par_input$lhs[input_sym], par_input$rhs[input_sym]),
        pmax(par_input$lhs[input_sym], par_input$rhs[input_sym]),
        sep = "::"
      )

      current_key <- dedup_key[sym_mask]
      match_idx   <- match(current_key, input_key)
      has_match   <- !is.na(match_idx)

      if (any(has_match)) {
        idx <- which(sym_mask)[has_match]
        lhs[idx] <- par_input$lhs[input_sym][match_idx[has_match]]
        rhs[idx] <- par_input$rhs[input_sym][match_idx[has_match]]
      }

      if (any(!has_match)) {
        idx      <- which(sym_mask)[!has_match]
        lhs_new  <- pmin(lhs[idx], rhs[idx])
        rhs_new  <- pmax(lhs[idx], rhs[idx])
        lhs[idx] <- lhs_new
        rhs[idx] <- rhs_new
      }
    } else {
      idx      <- which(sym_mask)
      lhs_new  <- pmin(lhs[idx], rhs[idx])
      rhs_new  <- pmax(lhs[idx], rhs[idx])
      lhs[idx] <- lhs_new
      rhs[idx] <- rhs_new
    }
  }

  colnames(samples) <- paste0(lhs, op, rhs)
  rownames(diagnostics) <- colnames(samples)

  # Summaries
  coefs <- apply(samples, MARGIN = 2, FUN = mean)
  vcov  <- stats::cov(samples)
  rhat  <- diagnostics[, "Rhat"]
  neff  <- diagnostics[, "n_eff"]

  # Build parTable (lavaan-like), including thresholds
  pars_clean_for_table <- cleanPars(colnames(samples))  # human-friendly labels where relevant

  se <- sqrt(diag(vcov))
  ci.lower <- apply(samples, MARGIN = 2, FUN = quantile, probs = 0.05)
  ci.upper <- apply(samples, MARGIN = 2, FUN = quantile, probs = 0.95)
  # non-symmetric P-value (Mplus note: https://www.statmodel.com/download/FAQ-Bootstrap%20-%20Pvalue.pdf)
  B <- NROW(samples)
  M <- apply(samples, MARGIN = 2, FUN = \(x) sum(x < 0))
  p.value <- 2 * apply(cbind(M/B, 1 - M/B), MARGIN = 1, FUN = min) # Two-sided (non-symmetric)

  # handle NaNs and zero SEs
  se.zero <- se <= .Machine$double.eps
  se  [se.zero]                <- NA
  rhat[se.zero | is.nan(rhat)] <- NA
  neff[se.zero | is.nan(neff)] <- NA

  parTable <- data.frame(
    lhs = lhs, op = op, rhs = rhs,
    est = coefs, std.error = se,
    z.value = coefs / se,
    p.value = p.value,
    ci.lower = ci.lower,
    ci.upper = ci.upper,
    R.hat = rhat, n.eff = neff,
    row.names = NULL
  )

  parTable <- tryCatch(sortParTableStan(parTable, compiled_model$info$parTable),
                       error = \(e) parTable)
  parTable <- modsemParTable(parTable)

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
summary.modsem_stan <- function(object, standardized = FALSE, ...) {
  parTable <- parameter_estimates(object)

  if (standardized) {
    parTable <- addTransformedEstimatesPT(
      parTable    = parTable,
      values.from = "est",
      values.to   = "Std.all",
      merge.by    = c("lhs", "op", "rhs"),
      FUN         = standardized_estimates,
      object      = object,
      ...
    )
  }

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


sortParTableStan <- function(parTable, parTable.input) {
  etas <- getSortedEtas(parTable.input, isLV = TRUE)
  xis  <- getXis(parTable.input, checkAny = FALSE, etas = etas, isLV = TRUE)

  indsXis  <- unlist(getIndsLVs(parTable.input, lVs = xis))
  indsEtas <- unlist(getIndsLVs(parTable.input, lVs = etas))

  opOrder <- c("=~", "~", "~1", "~~", "|", ":=")
  varOrder <- unique(c(indsXis, indsEtas, xis, etas))

  getScore <- function(x, order.by) {
    order.by <- unique(c(order.by, x)) # ensure that all of x is in order.by
    mapping  <- stats::setNames(seq_along(order.by), nm = order.by)
    score    <- mapping[x]

    if (length(score) != length(x)) {
      warning2("Sorting of parameter estimates failed!\n",
               immediate. = FALSE)

      return(seq_along(x))
    }

    score
  }

  scoreLhs <- getScore(x = parTable$lhs, order.by = varOrder)
  scoreOp  <- getScore(x = parTable$op,  order.by = opOrder)
  scoreRhs <- getScore(x = parTable$rhs, order.by = varOrder)

  out <- parTable[order(scoreOp, scoreLhs, scoreRhs), , drop = FALSE]
  rownames(out) <- NULL

  out
}
