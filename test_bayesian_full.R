devtools::load_all()
library(rstan)
rstan_options(auto_write = TRUE)      # cache compiled models
options(mc.cores = parallel::detectCores()) 


pattern_to_idx <- function(M, drop_first = FALSE) {
  ## M is a numeric pattern matrix with NA for *free* parameters.
  ## Returns list(row, col) of free positions.
  idx <- which(is.na(M), arr.ind = TRUE)
  if (drop_first)                 # remove the "first loading" rows if wanted
    idx <- idx[-1, , drop = FALSE]
  list(row = idx[, "row"], col = idx[, "col"])
}



specifyModelSTAN <- function(syntax = NULL,
                             data = NULL,
                             auto.fix.first = TRUE,
                             auto.fix.single = TRUE,
                             orthogonal.y = FALSE,   # <- currently ignored
                             orthogonal.x = FALSE,   # <- currently ignored
                             mean.observed = TRUE,
                             impute.na = TRUE) {     # <- was implicit before
  ## -----------------------------------------------------------------------
  ## 1. Parse lavaan‑like syntax and grab the usual pieces
  ## -----------------------------------------------------------------------
  if (!is.null(syntax)) parTable <- modsemify(syntax)
  stopif(is.null(parTable), "No parTable found")
  
  ## endogenous LVs (η)
  etas    <- rev(getSortedEtas(parTable,  isLV = TRUE, checkAny = TRUE))
  numEtas <- length(etas)
  
  ## exogenous LVs (ξ)
  xis     <- getXis(parTable, checkAny = TRUE)
  numXis  <- length(xis)
  
  ## interaction terms
  intTerms <- unique(parTable[grepl(":", parTable$rhs), "rhs"])
  numInts  <- length(intTerms)
  varsInts <- stringr::str_split(intTerms, ":", simplify = FALSE) |>
                stats::setNames(intTerms)
  
  ## bookkeeping
  lVs       <- c(xis, etas)
  numLVs    <- length(lVs)
  indsLVs   <- getIndsLVs(parTable, lVs)
  allInds   <- unique(unlist(indsLVs))
  numAllInds <- length(allInds)
  
  ## -----------------------------------------------------------------------
  ## 2. Clean / order raw data  (same helper you already use)
  ## -----------------------------------------------------------------------
  data <- cleanAndSortData(data,                       # user’s data.frame
                           unlist(indsLVs[xis]),       # indicator names for ξ
                           unlist(indsLVs[etas]),      # indicator names for η
                           impute.na = impute.na)
  
  ## -----------------------------------------------------------------------
  ## 3. Pattern matrices with NA for free parameters
  ##    (construct* helpers are the same as before)
  ## -----------------------------------------------------------------------
  LAMBDA <- constructLambda(lVs, indsLVs,
                            parTable = parTable,
                            auto.fix.first = auto.fix.first)$numeric
  
  GAMMA  <- constructGamma(etas, lVs, parTable = parTable)$numeric
  OMEGA  <- constructGamma(etas, intTerms, parTable = parTable)$numeric
  
  ## -----------------------------------------------------------------------
  ## 4. Translate those pattern matrices into *index vectors*
  ## -----------------------------------------------------------------------
  lambda_idx <- pattern_to_idx(LAMBDA, drop_first = TRUE)  # drop first loading
  gamma_idx  <- pattern_to_idx(GAMMA)
  omega_idx  <- pattern_to_idx(OMEGA)
  
  ## “first loading” row of each LV – used to fix λ = 1
  first_loading_row <- vapply(lVs, function(v) {
    indsLVs[[v]][1] |>                        # first indicator of that LV
      match(allInds)                          # row index in LAMBDA/Y
  }, FUN.VALUE = integer(1L))
  
  ## -----------------------------------------------------------------------
  ## 5. Interaction helpers:   prod_left / prod_right
  ## -----------------------------------------------------------------------
  lv_col <- setNames(seq_along(lVs), lVs)     # name → column index
  prod_left  <- integer(numInts)
  prod_right <- integer(numInts)
  for (q in seq_len(numInts)) {
    pair            <- varsInts[[q]]
    prod_left[q]    <- lv_col[pair[1]]
    prod_right[q]   <- lv_col[pair[2]]
  }
  
  ## -----------------------------------------------------------------------
  ## 6. Build the final list expected by sem_latent_interaction_fast.stan
  ## -----------------------------------------------------------------------
  list(
    ## sizes ---------------------------------------------------------------
    N        = nrow(data),
    K        = numAllInds,
    N_XIS    = numXis,
    N_ETAS   = numEtas,
    N_LVS    = numLVs,
    N_INT    = numInts,
    
    ## λ free positions ----------------------------------------------------
    N_FREE_LAMBDA = length(lambda_idx$row),
    lambda_row    = as.array(lambda_idx$row),
    lambda_col    = as.array(lambda_idx$col),
    
    ## Γ free positions ----------------------------------------------------
    N_FREE_GAMMA  = length(gamma_idx$row),
    gamma_row     = as.array(gamma_idx$row),
    gamma_col     = as.array(gamma_idx$col),
    
    ## Ω free positions ----------------------------------------------------
    N_FREE_OMEGA  = length(omega_idx$row),
    omega_row     = as.array(omega_idx$row),
    omega_col     = as.array(omega_idx$col),
    
    ## identification: first loading fixed to 1 ----------------------------
    first_loading_row = first_loading_row,
    
    ## latent‑interaction column pairs -------------------------------------
    prod_left  = as.array(prod_left),
    prod_right = as.array(prod_right),
    
    ## observed data -------------------------------------------------------
    Y = as.matrix(data[, allInds, drop = FALSE])
  )
}

model.syntax <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

stan_data <- specifyModelSTAN(model.syntax, data = oneInt)
stan_model <- stan_model("optimized_interaction_model.stan")

fit <- sampling(
  object = stan_model,
  data   = stan_data,
  chains = 2,
  iter   = 2000,
  warmup = 1000
)

summary(fit, c("gamma"))
stan_rhat(fit)
stan_trace(fit, "beta_X")
stan_trace(fit, "beta_XZ")
stan_trace(fit, "beta_Z")


# Three-way
model <- '
 X =~ x1 + x2 + x3
 Z =~ z1 + z2 + z3
 W =~ w1 + w2 + w3
 Y =~ y1 + y2 + y3

 Y ~ X + Z + W + X:Z + X:W + Z:W + X:Z:W
'

stan_data <- specifyModelSTAN(model, data = oneInt)

fit <- sampling(
  object = stan_model,
  data   = stan_data,
  chains = 2,
  iter   = 2000,
  warmup = 1000
)

summary(fit, c("gamma"))
