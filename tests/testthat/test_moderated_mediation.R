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
