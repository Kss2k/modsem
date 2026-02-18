devtools::load_all()

tpb <- "
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  INT ~ ATT + SN + PBC + ATT:PBC
  BEH ~ INT + PBC + ATT + INT:PBC + ATT:PBC
"

est_lms <- modsem(tpb, TPB, method = "lms", nodes = 32)
plot_jn(x = "ATT", z = "PBC", y = "BEH", model = est_lms, type = "indirect")
plot_jn(x = "ATT", z = "PBC", y = "BEH", model = est_lms, type = "total")

est_dblcent <- modsem(tpb, TPB, method = "dblcent")
p


model <- '
  X =~ att1 + att2 + att3 + att4 + att5
  W =~ pbc1 + pbc2 + pbc3
  M =~ int1 + int2 + int3
  Y =~ b1 + b2

  M ~ a*X
  Y ~ c1*X + c2*W + c3*X:W + b1*M + b2*M:W
'

fit <- modsem(model, data = TPB, method = "lms")
plot_jn(x = "ATT", z = "PBC", y = "BEH", model = est_dblcent, type = "indirect")
