#' @export
var_interactions.data.frame <- function(object, ignore.means = FALSE, 
                                        monte.carlo = FALSE, mc.reps = 1e6, ...) {

  # Preparation
  parTable <- removeInteractionVariances(fillColsParTable(object))

  ## interaction (and square) terms explicitly declared in the model
  intTerms <- unique(parTable$rhs[parTable$op == "~" & grepl(":", parTable$rhs)])
  if (length(intTerms) == 0L)
    return(modsemParTable(parTable))                 # nothing to do

  getIntTermLength <- \(xz) length(stringr::str_split(xz, pattern = ":")[[1L]])
  nway <- max(vapply(intTerms, FUN.VALUE = integer(1L), FUN = getIntTermLength))


  if (nway <= 2L && !monte.carlo) {
    varInteractionsAnalytic(
      parTable = parTable,
      intTerms = intTerms,
      ignore.means = ignore.means,
      ...
    )

  } else {
    varInteractionsMonteCarlo(
      parTable = parTable,
      intTerms = intTerms,
      ignore.means = ignore.means,
      mc.reps = mc.reps,
      ...
    )

  }
}


varInteractionsAnalytic <- function(parTable, intTerms, ignore.means = FALSE, ...) {
  # Gather first‑ and second‑order moments for the base variables --
  vars <- sort(unique(unlist(strsplit(intTerms, ":", fixed = TRUE))))

  if (ignore.means) mu <- stats::setNames(numeric(length(vars)), vars)
  else              mu <-  sapply(vars, getMean, parTable = parTable)

  getCov <- function(a, b) calcCovParTable(a, b, parTable)
  Sigma  <- outer(vars, vars, getCov)
  dimnames(Sigma) <- list(vars, vars)

  # Helper functions for the Gaussian fourth‑moment identity
  cov_PX <- function(i, j, k)
    mu[i] * Sigma[j, k] + mu[j] * Sigma[i, k]

  cov_PP <- function(i, j, k, l)
    Sigma[i, k] * Sigma[j, l] + Sigma[i, l] * Sigma[j, k] +
    mu[i] * mu[j] * Sigma[k, l] + mu[k] * mu[l] * Sigma[i, j] +
    mu[i] * mu[k] * Sigma[j, l] + mu[j] * mu[k] * Sigma[i, l] +
    mu[i] * mu[l] * Sigma[j, k] + mu[j] * mu[l] * Sigma[i, k]

  # Split product names once to avoid repetition
  prodPairs  <- lapply(intTerms, function(s) sort(strsplit(s, ":", fixed = TRUE)[[1L]]))
  names(prodPairs) <- intTerms                # keep original order/labels

  # Assemble rows to append
  rows <- vector("list", 0L)

  for (P in names(prodPairs)) {

    ij <- prodPairs[[P]];  i <- ij[1L];  j <- ij[2L]

    ## variance of the product term itself
    varP <- cov_PP(i, j, i, j)
    rows[[length(rows) + 1L]] <-
      data.frame(lhs = P, op = "~~", rhs = P, label = "",
                 est = varP, std.error = NA, z.value = NA,
                 p.value = NA, ci.lower = NA, ci.upper = NA)

    ## covariances with the *original* variables
    for (k in vars) {
      rows[[length(rows) + 1L]] <-
        data.frame(lhs = k, op = "~~", rhs = P, label = "",
                   est = cov_PX(i, j, k), std.error = NA, z.value = NA,
                   p.value = NA, ci.lower = NA, ci.upper = NA)
    }
  }

  # Covariances among the declared product terms
  if (length(prodPairs) > 1L) {
    comb <- utils::combn(names(prodPairs), 2L)
    for (cc in seq_len(ncol(comb))) {
      P1 <- comb[1L, cc];  P2 <- comb[2L, cc]
      ij <- prodPairs[[P1]];  kl <- prodPairs[[P2]]
      rows[[length(rows) + 1L]] <-
        data.frame(lhs = P1, op = "~~", rhs = P2, label = "",
                   est = cov_PP(ij[1L], ij[2L], kl[1L], kl[2L]),
                   std.error = NA, z.value = NA,
                   p.value = NA, ci.lower = NA, ci.upper = NA)
    }
  }

  # Bind & return
  parTable <- rbind(parTable, do.call(rbind, rows))
  modsemParTable(parTable)
}


varInteractionsMonteCarlo <- function(parTable, intTerms, ignore.means = FALSE, mc.reps = 1e6, ...) {
  # Gather first‑ and second‑order moments for the base variables --
  intVars  <- unique(unlist(strsplit(intTerms, ":", fixed = TRUE)))
  xis      <- unique(parTable[parTable$op == "~", "rhs"])
  vars     <- sort(unique(c(intVars, xis)))

  if (ignore.means) mu <- stats::setNames(numeric(length(vars)), vars)
  else              mu <-  sapply(vars, getMean, parTable = parTable)

  getCov <- function(a, b) calcCovParTable(a, b, parTable)
  Sigma  <- outer(vars, vars, getCov)
  dimnames(Sigma) <- list(vars, vars)

  XI <- as.data.frame(mvtnorm::rmvnorm(mc.reps, mean = mu, sigma = Sigma))

  for (intTerm in intTerms) {
    elems <- stringr::str_split(intTerm, pattern = ":")[[1L]]
    
    XI[[intTerm]] <- multiplyIndicatorsCpp(XI[elems])
  }

  Sigma.xz <- cov(XI)
  isXZ <- \(x) x %in% intTerms
  isXi <- \(x) x %in% xis
  isOK <- \(x, z)  (isXZ(x) && isXi(z)) || (isXZ(z) && isXi(x))

  for (i in seq_len(nrow(Sigma.xz))) for (j in seq_len(i)) {
    x   <- rownames(Sigma.xz)[[i]]
    z   <- colnames(Sigma.xz)[[j]]
    est <- Sigma.xz[x, z]
    
    if (!isOK(x, z))
      next

    newRow <- data.frame(lhs = x, op = "~~", rhs = z, label = "",
                         est = est, std.error = NA, z.value = NA,
                         p.value = NA, ci.lower = NA, ci.upper = NA)

    parTable <- rbind(parTable, newRow)
  }

  parTable
}


