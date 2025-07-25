devtools::load_all()

tpb1 <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:pbc1
'

testthat::expect_error(modsem(tpb1, data = TPB, method = "lms"),
                       regexp = "Element in product term .*pbc1.*")

tpb2 <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC2
'

testthat::expect_error(modsem(tpb2, data = TPB, method = "lms"),
                       regexp = "Element in product term .*PBC2.*")

testthat::expect_error(modsem(tpb2, data = TPB, method = "dblcent"),
                       regexp = "Unknown element in product term: .*PBC2.*")
