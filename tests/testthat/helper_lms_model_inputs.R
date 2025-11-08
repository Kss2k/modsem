lms_model_inputs <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) return(cache)

    data("TPB", package = "modsem", envir = environment())

    fit <- suppressWarnings(
      suppressMessages(
        modsem(
          "
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC + PBC:SN

  BEH ~ 1
  ATT ~ 1
  INT ~ 1
  PBC ~ 1
  SN  ~ 1
        ",
          data = TPB,
          method = "lms",
          nodes = 8,
          optimize = FALSE,
          standardize.data = TRUE,
          mean.observed = FALSE,
          verbose = FALSE
        )
      )
    )

    submodel <- fit$model$models[[1L]]
    P_group  <- fit$estep$P_GROUPS[[1L]]
    dataR    <- submodel$data

    cache <<- list(
      modelR = submodel,
      P      = P_group,
      colidx = dataR$colidx0,
      n      = as.integer(dataR$n.pattern),
      d      = as.integer(dataR$d.pattern)
    )

    cache
  }
})
