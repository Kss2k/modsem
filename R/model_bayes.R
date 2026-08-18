expandModsemParTable <- function(parTable, ordered = NULL,
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



