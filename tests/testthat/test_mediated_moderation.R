devtools::load_all()

tpb <- "
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
  INT ~ ATT + SN + PBC
  BEH ~ b * INT + PBC
  INT ~ a * PBC:SN
  BEH ~ c * PBC:SN

  IndirEffect := a * b
  TotalEffect := a * b + c
"

est_qml <- modsem(tpb, TPB, method = "qml")
standardized_estimates(est_qml)


est_lms <- modsem(tpb, TPB, method = "lms")
standardized_estimates(est_lms)
