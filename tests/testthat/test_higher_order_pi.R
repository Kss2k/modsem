devtools::load_all()

get_se <- function(x) sqrt(diag(vcov(x)))

# example 1
tpb <- '
  # First order constructs
  ATT =~ att1 + att2 + att3
  SN  =~ sn1 + sn2 + sn3
  PB =~ pb1 + pb2 + pb3
  PC =~ pc1 + pc2 + pc3
  BEH =~ b1 + b2

  # Higher order constructs
  INT =~ ATT + SN
  PBC =~ PC + PB

  # Higher order interaction
  INTxPBC =~ ATT:PC + ATT:PB + SN:PC + SN:PB

  # Structural model
  BEH ~ PBC + INT + INTxPBC
'

est_ca <- modsem(tpb, data = TPB_2SO, method = "ca")
testthat::expect_true(all(get_se(est_ca) < 0.15))

est_dblcent <- modsem(tpb, data = TPB_2SO, method = "dblcent")
testthat::expect_true(all(get_se(est_dblcent) < 0.15))

testthat::expect_error(modsem(tpb, TPB_2SO, method = "lms"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")
testthat::expect_error(modsem(tpb, TPB_2SO, method = "qml"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")

standardized_estimates(est_ca)
testthat::expect_error(
  standardized_estimates(est_ca, correction = TRUE, std.error = "delta"),
  regexp = "Correction of higher-order .*"
)

# example 2
tpb <- '
  # First order constructs
  ATT =~ att1 + att2 + att3
  SN  =~ sn1 + sn2 + sn3
  PBC =~ pbc1 + pbc2 + pbc3
  BEH =~ b1 + b2

  # Higher order constructs
  INT =~ ATT + PBC + SN

  # Higher order interaction
  INTxPBC =~ ATT:PBC + SN:PBC + PBC:PBC

  # Structural model
  BEH ~ PBC + INT + INTxPBC
'

est_ca <- modsem(tpb, data = TPB_1SO, method = "ca")
testthat::expect_true(all(get_se(est_ca) < 0.1))

est_dblcent <- modsem(tpb, data = TPB_1SO, method = "dblcent")
testthat::expect_true(all(get_se(est_dblcent) < 0.1))

testthat::expect_error(modsem(tpb, TPB_1SO, method = "lms"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")
testthat::expect_error(modsem(tpb, TPB_1SO, method = "qml"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")

# example 3
tpb <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  # Higher order constructs
  HIGH1 =~ PBC + ATT
  HIGH2 =~ INT + SN

  BEH ~ HIGH1 + HIGH2 + HIGH1:HIGH2
'


testthat::expect_error(modsem(tpb, TPB, method = "rca"),
                       regexp = "The ':' operator is not allowed for higher order *")
