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
  ovs     <- getStructOVs(parTable)
  lvs     <- getLVs(parTable)
  ordered <- unique(parTable[parTable$op == "|", "lhs"])

  INDICATOR_PREFIX <- "INDICATOR_"
  DISTURBANCE_PREFIX  <- "DISTURBANCE_"
  INTERCEPT_PREFIX <- "INTERCEPT_"
  VAR_PREFIX <- "VARIANCE_"
  LAMBDA_PREFIX <- "LAMBDA_"
  LV_PREFIX <- "LATENT__"
  COV  <- "__COVARIANCE__"
  MSR  <- "__MEASUREMENT__"
  INTR <- "__INTERCEPT__1"
  REG  <- "__REGRESSION__"
  LAB  <- "LABEL__"

  DATA <- "int<lower=0> N;"
  TRANSFORMED_DATA <- character(0L)
  PARAMETERS <- character(0L)
  TRANSFORMED_PARAMETERS <- character(0L)
  MODEL <- character(0L)
  GENERATED_QUANTITIES <- character(0L)

  # WE potentially need to pre-define labels/repeated parameters here
  isLab <- !canBeNumeric(parTable$mod)
  # parTable[isLab, "mod"] <- paste0(LAB, parTable[isLab, "mod"])
  labels <- unique(parTable[isLab, "mod"])

  for (lab in labels) {
    PARAMETERS <- c(PARAMETERS, sprintf("real %s;", lab))

    prior.idx <- which(parTable$lhs == lab & parTable$op == ":=")
    if (length(prior.idx)) {
      MODEL <- c(MODEL,
        sprintf("%s ~ %s", lab, parTable[prior.idx, "rhs"])
      )
    }
  }

  for (ind in inds) {
    # parameter names
    nm  <- paste0(INDICATOR_PREFIX, ind)
    dnm <- paste0(DISTURBANCE_PREFIX, nm)
    vnm <- paste0(nm, COV, dnm)
    inm <- paste0(nm, INTR)

    data.i <- paste0("array[N] real ", nm, ";")
    par.i <- character(0L)
    tpar.i <- character(0L)
    model.i <- character(0L)

    # Intercept
    b0.idx <- which(parTable$op == "~1" & parTable$lhs == ind)

    if (length(b0.idx) != 1 || parTable[b0.idx, "mod"] == "") {
      par.i  <- c(par.i, paste0("real ", inm, ";"))

    } else {
      mod.i <- parTable[b0.idx, "mod"]
      tpar.i  <- c(tpar.i, sprintf("real %s = %s;", inm, mod.i))
    }

    idx <- which(parTable$op == "=~" & parTable$rhs == ind)
    tpar.i <- paste0(
    )

    eq <- inm
    for (i in idx) {
      lv.i  <- parTable[i, "lhs"]
      mod.i <- parTable[i, "mod"]
      msr.i <- paste0(lv.i, MSR, ind)

      if (mod.i == "") par.i  <- c(par.i, sprintf("real %s;", msr.i))
      else tpar.i  <- c(tpar.i, sprintf("real %s = %s;", inm, mod.i))

      eq <- sprintf("%s + %s", eq, paste0(msr.i, "*", LV_PREFIX, lv.i))
    }
   
    # no residual covariances and continous indicators
    tpar.i <- c(tpar.i,
      sprintf("array[N] real %s = %s - (%s);", dnm, nm, eq)
    )

    model.i <- c(model.i,
      sprintf("%s ~ normal(0.0, sqrt(%s))", dnm, vnm)
    )

    DATA <- c(
      DATA, paste0(data.i, collapse = "\n")
    )

    PARAMETERS <- c(
      PARAMETERS, paste0(par.i, collapse = "\n")
    )

    TRANSFORMED_PARAMETERS <- c(
      TRANSFORMED_PARAMETERS, paste0(tpar.i, collapse = "\n")
    )

    MODEL <- c(
      MODEL, paste0(model.i, collapse = "\n")
    )
  }

  sprintf(STAN_TEMPLATE,
    paste0(DATA, collapse = "\n"),
    paste0(PARAMETERS, collapse = "\n"),
    paste0(TRANSFORMED_PARAMETERS, collapse = "\n"),
    paste0(MODEL, collapse = "\n")
  )
}
