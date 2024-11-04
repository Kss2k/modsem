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

  # Adding some constraints
  INTERACTION ~~ 0*HIGH1 + 0*HIGH2
'

est  <- modsem(tpb, data = TPB, method = "ca",
               suppress.warnings.match = TRUE)
summary(est)

testthat::expect_error(modsem(tpb, TPB, method = "lms"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")
testthat::expect_error(modsem(tpb, TPB, method = "qml"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")

tpb <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3

  # Higher order constructs
  HIGH1 =~ PBC + ATT
  HIGH2 =~ INT + SN
  # Higher order interaction
  
  BEH =~ b1 + b2

  BEH ~ HIGH1 + HIGH2 + HIGH1:HIGH2
'

testthat::expect_error(modsem(tpb, TPB, method = "rca"),
                       regexp = "The ':' operator is not allowed for higher order *")
