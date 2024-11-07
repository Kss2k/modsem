devtools::load_all()

get_se <- function(x) sqrt(diag(vcov(x)))

# example 1
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
testthat::expect_true(all(get_se(est) < 0.6))

testthat::expect_error(modsem(tpb, TPB, method = "lms"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")
testthat::expect_error(modsem(tpb, TPB, method = "qml"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")

# example 2
tpb <- '
  # First order constructs
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  BEH =~ b1 + b2

  # Higher order constructs
  INT =~ PBC + ATT + SN

  # Higher order interaction
  INTxPBC =~ PBC:PBC + ATT:PBC + SN:PBC
  
  # Structural model
  BEH ~ PBC + INT + INTxPBC
'

est  <- modsem(tpb, data = TPB, method = "ca",
               suppress.warnings.match = TRUE)
summary(est)
testthat::expect_true(all(get_se(est) < 0.1))

est  <- modsem(tpb, data = TPB, method = "rca")
testthat::expect_true(all(get_se(est) < 0.1))

testthat::expect_error(modsem(tpb, TPB, method = "lms"),
                       regexp = "Higher-order latent variables are not supported in the lms and qml approaches.")
testthat::expect_error(modsem(tpb, TPB, method = "qml"),
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
