expandModsemParTable <- function(parTable,
                                 ordered = NULL,
                                 auto.fix.first = TRUE,
                                 orthogonal.x = FALSE,
                                 orthogonal.y = FALSE) {
  # Before moving the higher order measurment model to the structural model
  # we should check for factor loadings to be fixed
  if (auto.fix.first) {
    lvs <- getLVs(parTable)

    for (lv in lvs) {
      idx <- which(parTable$op == "=~" & parTable$lhs == lv)
      
      if (length(idx) && parTable$mod[idx[[1L]]] == "")
        parTable[idx[[1L]], "mod"] <- "1"
    }
  }

  parTable <- higherOrderMeasr2Struct(parTable)

  lvs  <- getLVs(parTable)
  inds <- getIndicators(parTable)
  ovs  <- getStructOVs(parTable)
  inds <- getIndicators(parTable)
  vars <- unique(c(lvs, ovs, inds))

  for (v in vars) {

    # Variance
    if (!any(parTable$lhs == v & parTable$op == "~~" & parTable$rhs == v)) {
      parTable <- rbind(
        parTable,
        data.frame(lhs = v, op = "~~", rhs = v, mod = "")
      )
    }
   
    # Intercept
    if (!any(parTable$lhs == v & parTable$op == "~1")) {
      mod <- if (v %in% lvs) "0" else ""
      parTable <- rbind(
        parTable,
        data.frame(lhs = v, op = "~~", rhs = v, mod = mod)
      )
    }

  }

  svars <- union(lvs, ovs)
  pure.xi <- isPureXi(svars, parTable = parTable)
  pure.eta <- isPureEta(svars, parTable = parTable)

  for (i in seq_along(svars)) for (j in seq_len(i - 1L)) {
    x <- svars[[i]]
    y <- svars[[j]]

    cond <- any(
      parTable$lhs == x & parTable$op == "~~" & parTable$rhs == y &
      parTable$lhs == y & parTable$op == "~~" & parTable$rhs == x
    )

    xx <- !orthogonal.x && pure.xi[[i]] && pure.xi[[j]]
    yy <- !orthogonal.y && pure.eta[[i]] && pure.eta[[j]]

    if (!cond && (xx || yy)) {
      parTable <- rbind(
        parTable,
        data.frame(lhs = x, op = "~~", rhs = y, mod = "")
      )
    }
  }

  # add priors
  parTable$prior <- NA # missing for now

  # finally remove duplicates
  removeParTableDuplicates(parTable)
}


isPureXi <- function(lvs, parTable) {
  etas <- unique(parTable[parTable$op == "~", "lhs"])
  !lvs %in% etas
}


isPureEta <- function(lvs, parTable) {
  preds <- unique(parTable[parTable$op == "~", "rhs"])
  etas <- unique(parTable[parTable$op == "~", "lhs"])
  lvs %in% etas & !lvs %in% preds
}


removeCovarianceDuplicates <- function(parTable) {
  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs
  isCov <- op == "~~" & lhs == rhs

  # order invariant labels
  sortedLabels <- paste0(pmax(lhs, rhs), "~~", pmin(lhs, rhs))

  parTable[!(duplicated(sortedLabels) & isCov), , drop = FALSE]
}


removeNaturalDuplicates <- function(parTable) {
  par <- paste0(parTable$lhs, parTable$op, parTable$rhs)
  parTable[!duplicated(par), , drop = FALSE]
}


removeParTableDuplicates <- function(parTable) {
  removeNaturalDuplicates(parTable) |>
    removeCovarianceDuplicates()
}


STAN_TEMPLATE <- "
data {
%s
}

parameters {
%s
}

transformed parameters {
%s
}

model {
%s
}
"


buildStanSyntaxFromParTable <- function(parTable) {

  inds    <- getIndicators(parTable)
  xis     <- getXis(parTable, isLV = FALSE)
  etas    <- getSortedEtas(parTable, isLV = FALSE)
  ovs     <- setdiff(getStructOVs(parTable), inds) # ovs only in the structural model
  lvs     <- getLVs(parTable)
  ordered <- unique(parTable[parTable$op == "|", "lhs"])
  cont    <- setdiff(inds, ordered)
  vars    <- unique(c(inds, xis, etas, ovs, lvs))

  notOk <- grepl("__", vars)
  mod_stopif(any(notOk),
    "`__` is a reserved pattern and can not be used in variable names!",
    "Please rename these variables:", paste0(vars[notOk], collapse = ", ")
  )

  ordOvs <- intersect(ovs, ordered)
  mod_stopif(length(ordOvs),
    "Observed variables which aren't indicators are not allowed to be ordered!",
    "Please redefine these as latent variables with a single indicator:",
    paste0(ordOvs, collapse = ", ")
  )

  IND_PREFIX <- "IND__"
  RES_PREFIX <- "RES__" # technically these are disturbances not residuals, but who cares?
  LV_PREFIX  <- "LV__"
  COV        <- "__COV__"
  MSR        <- "__MSR__"
  INTR       <- "__INTR__1"
  REG        <- "__REG__"
  MOD        <- "__MOD__"

  DATA <- c("int<lower=0> N;") # , "int<lower=0> N_CONT_INDS;")
  TRANSFORMED_DATA <- character(0L)
  PARAMETERS <- character(0L)
  TRANSFORMED_PARAMETERS <- character(0L)
  MODEL <- character(0L)
  GENERATED_QUANTITIES <- character(0L)

  # WE potentially need to pre-define labels/repeated parameters here
  isLab <- !canBeNumeric(parTable$mod) & parTable$mod != ""
  # parTable[isLab, "mod"] <- paste0(LAB, parTable[isLab, "mod"])
  labels <- unique(parTable[isLab, "mod"])

  for (lab in labels) {
    PARAMETERS <- c(PARAMETERS, sprintf("real %s;", lab))

    prior.idx <- which(parTable$lhs == lab & parTable$op == ":=") # use the := operator (for now)
    if (length(prior.idx)) {
      MODEL <- c(MODEL,
        sprintf("%s ~ %s", lab, parTable[prior.idx, "rhs"])
      )
    }
  }

  k <- length(lvs)
  PARAMETERS <- c(PARAMETERS,
    sprintf("vector[N] MAT__ZETA[%d];", k)
  )

  TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
    sprintf("matrix[%d,%d] MAT__PSI;", k, k),
    sprintf("vector[%d] VEC_ZETA_MU;", k)
  )

  for (i in seq_along(lvs)) {
    lv <- lvs[[i]]

    # Intercept
    intr <- paste0(lv, INTR)
    b0.idx <- which(parTable$op == "~1" & parTable$lhs == lv)

    if (length(b0.idx) != 1 || parTable[b0.idx, "mod"] == "") {
      PARAMETERS <- c(PARAMETERS, paste0("real ", intr, ";"))

    } else {
      mod <- parTable[b0.idx, "mod"]
      PARAMETERS <- c(TRANSFORMED_PARAMETERS,
        sprintf("real %s = %s;", intr, mod, ";")
      )
    }

    TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
      sprintf("VEC_ZETA_MU[%d] = %s;", i, intr)
    )
  }

  idxPsi <- which(
    parTable$lhs %in% lvs & parTable$op == "~~" & parTable$rhs %in% lvs
  )

  psiScalarTemplate <- "MAT__PSI[%d,%d] = %s;"
  for (idx in idxPsi) {
    lhs <- parTable[idx, "lhs"]
    rhs <- parTable[idx, "rhs"]
    mod <- parTable[idx, "mod"]
    par <- paste0(lhs, COV, rhs)
    i   <- which(lvs == lhs)
    j   <- which(lvs == rhs)

    if (mod == "") PARAMETERS  <- c(PARAMETERS, sprintf("real %s;", par))
    else           PARAMETERS  <- c(PARAMETERS, sprintf("real %s = %s;", par, par))

    TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
      sprintf(psiScalarTemplate, i, j, par),
      sprintf(psiScalarTemplate, j, i, par)
    )

  }
  
  MODEL <- c(MODEL,
    "MAT__ZETA ~ multi_normal(VEC_ZETA_MU, MAT__PSI);"
  )

  for (ov in ovs) {
    mod_msg_stop("Observed variables in the structural model are not allowed (yet)!")
  }

  for (xi in xis) {
    lxi <- paste0(LV_PREFIX, xi)
    i <- which(lvs == xi)
    TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
      sprintf("vector[N] %s = to_vector(MAT__ZETA[:,%d]);", lxi, i)
    )
  }

  for (eta in etas) {
    leta <- paste0(LV_PREFIX, eta)
    reta <- paste0(RES_PREFIX, eta)
    i    <- which(lvs == eta)

    TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
      sprintf("vector[N] %s = to_vector(MAT__ZETA[:,%d]);", reta, i)
    )

    eq <- reta
    idx <- which(parTable$lhs == eta & parTable$op == "~")

    for (j in idx) {
      pred <- parTable[j, "rhs"]
      mod  <- parTable[i, "mod"]
      reg  <- paste0(eta, REG, stringr::str_replace(pred, ":", MOD))

      if (mod == "") {
        PARAMETERS <- c(PARAMETERS, sprintf("real %s;", reg))
      } else {
        TRANSFORMED_PARAMETERS  <- c(TRANSFORMED_PARAMETERS,
          sprintf("real %s = %s;", reg, mod)
        )
      }

      pred <- stringr::str_split_1(pred, ":")
      prod <- paste0(paste0(LV_PREFIX, pred), collapse = ".*")

      eq <- sprintf("%s + %s", eq, paste0(reg, "*", prod))
    }

    TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
      sprintf("vector[N] %s = %s;", leta, eq)
    )
  }

  for (ind in inds) {
    # parameter names
    nm  <- paste0(IND_PREFIX, ind)
    dnm <- paste0(RES_PREFIX, nm)
    vnm <- paste0(nm, COV, nm)
    inm <- paste0(nm, INTR)

    DATA <- c(DATA,
      sprintf("vector[N] %s;", nm)
    )

    # Intercept
    b0.idx <- which(parTable$op == "~1" & parTable$lhs == ind)

    if (length(b0.idx) != 1 || parTable[b0.idx, "mod"] == "") {
      PARAMETERS  <- c(PARAMETERS, paste0("real ", inm, ";"))

    } else {
      mod <- parTable[b0.idx, "mod"]
      TRANSFORMED_PARAMETERS  <- c(TRANSFORMED_PARAMETERS,
        sprintf("real %s = %s;", inm, mod)
      )
    }
    
    # Residual variance 
    # This should probably be moved if we ever add residual covariances
    rv.idx <- which(
      parTable$lhs == ind & parTable$op == "~~" & parTable$rhs == ind
    )

    if (length(rv.idx) != 1 || parTable[rv.idx, "mod"] == "") {
      PARAMETERS <- c(PARAMETERS, paste0("real ", vnm, ";"))

    } else {
      mod <- parTable[b0.idx, "mod"]
      TRANSFORMED_PARAMETERS  <- c(TRANSFORMED_PARAMETERS,
        sprintf("real %s = %s;", vnm, mod)
      )
    }

    idx <- which(parTable$op == "=~" & parTable$rhs == ind)

    eq <- inm
    for (i in idx) {
      lv.i  <- parTable[i, "lhs"]
      mod.i <- parTable[i, "mod"]
      msr.i <- paste0(lv.i, MSR, ind)

      if (mod.i == "") {
        PARAMETERS <- c(PARAMETERS, sprintf("real %s;", msr.i))
      } else {
        TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
          sprintf("real %s = %s;", msr.i, mod.i)
        )
      }

      eq <- sprintf("%s + %s", eq, paste0(msr.i, "*", LV_PREFIX, lv.i))
    }
   
    # no residual covariances and continous indicators
    TRANSFORMED_PARAMETERS <- c(TRANSFORMED_PARAMETERS,
      sprintf("vector[N] %s = %s - (%s);", dnm, nm, eq)
    )

    MODEL <- c(MODEL,
      sprintf("%s ~ normal(0.0, sqrt(%s));", dnm, vnm)
    )
  }

  sprintf(STAN_TEMPLATE,
    paste0(DATA, collapse = "\n"),
    paste0(PARAMETERS, collapse = "\n"),
    paste0(TRANSFORMED_PARAMETERS, collapse = "\n"),
    paste0(MODEL, collapse = "\n")
  )
}
