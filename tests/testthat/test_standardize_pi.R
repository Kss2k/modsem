devtools::load_all()

m1 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
# Inner Model
  Y ~ X + Z + a * X:Z

  c := a
"

est <- modsem(m1, data = oneInt)
parTable <- standardized_estimates(est)
parTable <- standardized_estimates(est, correction = TRUE)
varXZ <- parTable[parTable$lhs == "XZ" & parTable$rhs == "XZ", "est"]
testthat::expect_true(varXZ > 1.02)

a <- parTable[parTable$label == "a", "est"]
c <- parTable[parTable$label == "c", "est"]

testthat::expect_equal(a, c)

m2 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
# Inner Model
  Y ~ X + Z + X:Z + X:X
"

est <- modsem(m2, data = oneInt)
parTable <- standardized_estimates(est, correction = TRUE)
varXZ <- parTable[parTable$lhs == "XZ" & parTable$rhs == "XZ", "est"]
varXX <- parTable[parTable$lhs == "XX" & parTable$rhs == "XX", "est"]
testthat::expect_true(varXZ > 1.02)
testthat::expect_true(varXX > 1.7)
