devtools::load_all()


tpb <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3

  # Higher order constructs
  HIGH1 =~ PBC + ATT
  HIGH2 =~ INT + SN
  # Higher order interaction
  INTERACTION =~ PBC:INT + ATT:SN + PBC:SN + ATT:INT
  
  BEH =~ b1 + b2

  BEH ~ HIGH1 + HIGH2 + INTERACTION
'

est  <- modsem(tpb, TPB, method = "rca")
