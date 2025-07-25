devtools::load_all()
tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Covariances
  ATT ~~ SN + PBC
  PBC ~~ SN
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * PBC
  BEH ~ INT + PBC
  BEH ~ d * INT:PBC

  a_b := a * b * d # check that labels arent overwritten without redefining
'

method <- c("ca", "rca", "uca", "dblcent")
ests <- lapply(method, function(m) modsem(tpb, data = TPB, method = m))
estsMatch <- lapply(method, function(m) modsem(tpb, data = TPB, method = m,
                                               match = TRUE))

print(summary(ests[[1]]))
