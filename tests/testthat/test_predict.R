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

  ebm <- modsemPredictDA(fit_fiml, method = "EBM")
  ml  <- modsemPredictDA(fit_fiml, method = "ML")

  # newdata: use first 20 rows (may contain missing)
  ebm_new <- modsemPredictDA(fit_fiml, newdata = TPB2[1:20, ], method = "EBM")
  ml_new  <- modsemPredictDA(fit_fiml, newdata = TPB2[1:20, ], method = "ML")
})

testthat::expect_no_condition({
  fit_qml <- modsem(tpb, TPB, method = "qml")

  ebm <- modsemPredictDA(fit_qml, method = "EBM")
  ml  <- modsemPredictDA(fit_qml, method = "ML")

  # newdata: use last 20 rows
  ebm_new <- modsemPredictDA(fit_qml, newdata = TPB[181:200, ], method = "EBM")
  ml_new  <- modsemPredictDA(fit_qml, newdata = TPB[181:200, ], method = "ML")
})
