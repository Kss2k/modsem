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
modsem_stan <- function(model.syntax = NULL,
                        data = NULL,
                        compiled_model = NULL,
                        chains = 2, 
                        iter = 2000, 
                        warmup = iter / 2, 
                        ...) {
  if (is.null(compiled_model)) {
    stopif(is.null(model.syntax), 
           "One of `model.syntax` or `compiled_model` has to be provided!") 

    compiled_model <- compile_stan_model(model.syntax)
  }

  lVs     <- compiled_model$info$lVs
  indsLVs <- compiled_model$info$indsLVs
  inds    <- unique(unlist(indsLVs))
  etas    <- compiled_model$info$etas
  deps    <- c(inds, etas)

  stan_data <- get_stan_data(compiled_model = compiled_model, data = data)

  message("Sampling Stan model...")
  fit <- rstan::sampling(object  = compiled_model$stan_model,
                         data    = stan_data,
                         chains  = chains,
                         iter    = iter,
                         warmup  = warmup,
                         pars    = compiled_model$info$exclude.pars,
                         include = FALSE,
                         ...)

  samples <- as.matrix(fit)
  samples <- samples[, !grepl("\\[[0-9]+\\]", colnames(samples))]
  
  pars <- stringr::str_replace_all(colnames(samples), STAN_OPERATOR_LABELS)
  colnames(samples) <- pars
 
  OP <- "~~|=~|~1|~"
  op <- stringr::str_extract(pars, pattern = OP)
  op[is.na(op)] <- ":="
  lr <- stringr::str_split_fixed(pars, pattern = OP, n = 2)

  lhs <- lr[, 1]
  rhs <- lr[, 2]

  isSD <- lhs == rhs & op == "~~" & (lhs %in% deps | rhs %in% deps)
  samples[, isSD] <- samples[, isSD]^2

  coefs <- apply(samples, MARGIN = 2, FUN = mean)
  vcov  <- cov(samples)

  se <- sqrt(diag(vcov))
  parTable <- data.frame(
    lhs = lhs, op = op, rhs = rhs, est = coefs, std.error = se,
    z.value = coefs / se, p.value = 2 * stats::pnorm(-abs(coefs / se)),
    ci.lower = coefs - CI_WIDTH * se, ci.upper = coefs + CI_WIDTH * se
  )

  parTable <- rbind(
      data.frame(lhs = lVs, op  = "=~",
                 rhs = sapply(lVs, FUN = \(lV) indsLVs[[lV]][[1L]]),
                 est = 1, std.error = NA, z.value = NA, p.value = NA, 
                 ci.lower = NA, ci.upper = NA),
      parTable
  )

  parTable <- parTable[order(parTable$op), ]
  rownames(parTable) <- NULL
  parTable <- modsemParTable(parTable[parTable$op != ":=", , drop = FALSE])

  out <- list(fit = fit,
              parTable = parTable,
              coefs = coefs,
              vcov = vcov,
              samples = samples)


  class(out) <- "modsem_stan"

  out
}


#' @export
summary.modsem_stan <- function(object, ...) {
  summarize_partable(object$parTable)
}


#' @export
print.modsem_stan <- function(x, ...) {
  print(x$parTable)
}


#' @export
parameter_estimates <- function(object, ...) {
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
