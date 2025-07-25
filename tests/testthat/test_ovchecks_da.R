devtools::load_all()


m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3

  # Inner model
  Y ~ X + Z + X:Z + x1
'

testthat::expect_error(modsem(m1, oneInt, method = "lms"),
                       regexp = "Observed variables are not allowed in .*")
testthat::expect_error(modsem(m1, oneInt, method = "qml"),
                       regexp = "Observed variables are not allowed in .*")


m2 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3

  # Inner model
  Y ~ X + Z + X:Z + jk
'

testthat::expect_error(modsem(m2, oneInt, method = "lms"),
                       regexp = "Observed variables are not allowed in .*")
testthat::expect_error(modsem(m2, oneInt, method = "qml"),
                       regexp = "Observed variables are not allowed in .*")


m3 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3 + jk

  # Inner model
  Y ~ X + Z + X:Z
'

testthat::expect_error(modsem(m3, oneInt, method = "lms"),
                       regexp = "Missing observed variables in data:.*jk.*")
testthat::expect_error(modsem(m3, oneInt, method = "qml"),
                       regexp = "Missing observed variables in data:.*jk.*")


m4 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3 + y1

  # Inner model
  Y ~ X + Z + X:Z
'

testthat::expect_error(modsem(m4, oneInt, method = "lms"),
                       regexp = "The same indicator cannot be used .* y1")
testthat::expect_error(modsem(m4, oneInt, method = "qml"),
                       regexp = "The same indicator cannot be used .* y1")


tpb_main <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT2 =~ a1 * att1 + a2 * att2 + att3 + att4 + att5
  SN =~ s1 * sn1 + sn2
  PBC =~ p1 * pbc1 + pbc2 + pbc3
  INT =~ i1 * int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + a * PBC
"

tpb_cov <- "
  INT ~ ATT + b * SN + c * PBC
"

expect_error(modsem(tpb_main, data = TPB, method = "lms",
                    calc.se=FALSE, cov.syntax = tpb_cov),
             regexp = "All latent variables in the cov-model must be an exogenous variable in the main model")
