devtools::load_all()
m1 <- '
# Outer Model
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2 

# Inner model
  Y ~ X + Z + X:Z
'
run <- tryCatch({
    MplusAutomation::detectMplus()
    TRUE
  },
  error = function(e) FALSE
)
if (run) {
  mplus <- modsem(m1, oneInt, method = "mplus", estimator = "MLR")
  print(summary(mplus))
  plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = mplus)
}


tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC + INT:PBC
"

if (run) {
  mplus_tpb <- modsem(tpb, data = TPB, method = "mplus", rcs = TRUE)
  standardized_estimates(mplus_tpb, type = "stdyx")
}

tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  SUBJECTIVE_NORMS =~ sn1 + sn2 
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + PBC + INT:PBC:SUBJECTIVE_NORMS + INT:PBC + SUBJECTIVE_NORMS
"

if (run) {
  mplus_tpb_3way <- modsem(tpb, data = TPB[1:250, ], method = "mplus", rcs = TRUE)
  standardized_estimates(mplus_tpb_3way)
}
