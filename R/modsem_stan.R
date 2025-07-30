modsem_stan <- function(model.syntax = NULL,
                        data = NULL,
                        compiled_model = NULL,
                        chains = 2, 
                        iter = 2000, 
                        warmup = 1000, ...) {
  if (is.null(compiled_model)) {
    stopif(is.null(model.syntax), 
           "One of `model.syntax` or `compiled_model` has to be provided!") 

    compiled_model <- compile_stan_model(m1)
  }

  stan_data <- get_stan_data(compiled_model = compiled_model, data = data)
 
  message("Sampling STAN model...")
  fit <- rstan::sampling(object = compiled_model$stan_model,
                         data   = stan_data,
                         chains = chains,
                         iter   = iter,
                         warmup = warmup,
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
  
  isSD <- lhs == rhs & op == "~~"
  samples[, isSD] <- samples[, isSD]^2

  coefs <- apply(samples, MARGIN = 2, FUN = mean)
  vcov  <- cov(samples)

  se <- sqrt(diag(vcov))
  parTable <- data.frame(
    lhs = lhs, op = op, rhs = rhs, est = coefs, std.error = se,
    z.value = coefs / se, p.value = 2 * stats::pnorm(-abs(coefs / se)),
    ci.lower = coefs - CI_WIDTH * se, ci.upper = coefs + CI_WIDTH * se
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
