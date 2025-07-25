devtools::load_all()

m1 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
# Inner Model
  Y ~ X + Z + X:Z
"

est <- modsem(m1, data = oneInt, method = "lms", mean.observed = FALSE)
est <- modsem(m1, data = oneInt, method = "qml", mean.observed = FALSE)
summary(est)


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

est <- modsem(tpb, data = TPB, method = "lms", mean.observed = FALSE,
              nodes = 32)
summary(est)


tpb_uk <- "
# Outer Model (Based on Hagger et al., 2007)
 ATT =~ att3 + att2 + att1 + att4
 SN =~ sn4 + sn2 + sn3 + sn1
 PBC =~ pbc2 + pbc1 + pbc3 + pbc4
 INT =~ int2 + int1 + int3 + int4
 BEH =~ beh3 + beh2 + beh1 + beh4

# Inner Model (Based on Steinmetz et al., 2011)
 # Causal Relationsships
 INT ~ ATT + SN + PBC
 BEH ~ INT + PBC
 BEH ~ INT:PBC
"

est <- modsem(tpb_uk, data = TPB_UK, method = "lms", mean.observed = FALSE,
              nodes = 32)
summary(est)
