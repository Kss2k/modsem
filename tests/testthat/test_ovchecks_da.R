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
