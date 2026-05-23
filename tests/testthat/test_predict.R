devtools::load_all()

tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
'


testthat::expect_no_condition({
  TPB2 <- TPB
  TPB2[23:150, "att1"] <- NA

  fit_fiml <- modsem(tpb, TPB2, method = "lms", nodes = 32, missing = "fiml")

  ebm     <- modsem_predict(fit_fiml, method = "EBM")
  ml      <- modsem_predict(fit_fiml, method = "ML", type = "lv")
  ml_all  <- modsem_predict(fit_fiml, method = "ML", type = "all")

  # newdata: use first 20 rows (may contain missing)
  ebm_new <- modsem_predict(fit_fiml, newdata = TPB2[1:20, ], method = "EBM")
  ml_new  <- modsem_predict(fit_fiml, newdata = TPB2[1:20, ], method = "ML")
})


testthat::expect_no_condition({
  fit_qml <- modsem(tpb, TPB, method = "qml")

  ebm     <- modsem_predict(fit_qml, method = "EBM")
  ml      <- modsem_predict(fit_qml, method = "ML")
  ml_all  <- modsem_predict(fit_qml, method = "ML", type = "all")

  # newdata: use last 20 rows
  ebm_new <- modsem_predict(fit_qml, newdata = TPB[181:200, ], method = "EBM", type = "lv")
  ml_new  <- modsem_predict(fit_qml, newdata = TPB[181:200, ], method = "ML")
})


tpb.comp <- '
  ATT <~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH <~ b1 + b2

  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
'

testthat::expect_no_condition({
  fit_comp <- modsem(tpb.comp, TPB, method = "lms")

  ebm     <- modsem_predict(fit_comp, method = "EBM")
  ml      <- modsem_predict(fit_comp, method = "ML", type = "lv")

  ebm_new <- modsem_predict(fit_comp, newdata = TPB[1:20, ], method = "EBM", type = "all")
  ml_new  <- modsem_predict(fit_comp, newdata = TPB[1:20, ], method = "ML")
})
