devtools::load_all()
m1 <- '
# Outer Model
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2

# Inner model
  Y ~ b1 * X + b2 * Z + b3 * X:Z

  ccoef := (b1 + b2 * b3) / 2
  ccoef2 := b1 * b2

  b1 == b2
'

run <- tryCatch({
    MplusAutomation::detectMplus()
    TRUE
  },
  error = function(e) FALSE
)
if (run) {
  mplus <- modsem(m1, oneInt, method = "mplus", estimator = "MLR")
  standardized_estimates(mplus, type = "modsem")
  print(summary(mplus))
  print(vcov(mplus))
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
  print(standardized_estimates(mplus_tpb, type = "modsem"))
  print(standardized_estimates(mplus_tpb, type = "stdyx"))
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
  mplus_tpb_3way <- modsem(tpb, data = TPB[1:250, ], method = "mplus", rcs = TRUE,
                           integration = 8)
  standardized_estimates(mplus_tpb_3way)
}


# Test robust std.errors
mod <- '
# X1-3 are Level 1 variables
X1 =~ x1 
X2 =~ x2
X3 =~ x3

# W1-2 are Level 2 variables
W1 =~ w1
W2 =~ w2

fw =~ y1 + y2 + y3
fw ~ X1 + X2 + X3 + W1 + W2
'

# Standard errors corrected for clustering
fit.rc <- modsem(mod, lavaan::Demo.twolevel, rcs = TRUE, method = "mplus",
                 cluster = "cluster")
summary(fit.rc)
