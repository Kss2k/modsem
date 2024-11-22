devtools::load_all()
m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3

  # Inner model
  Y ~ X + Z + X:Z
'

oneInt2 <- oneInt
oneInt2[c(176, 176, 258, 1900),
        c(1, 2, 3, 7)] <- NA

# Double centering approach
testthat::expect_warning(modsem(m1, oneInt2),
                         regex = "Removing missing values case-wise")

est <- modsem(m1, oneInt2, na.rm=TRUE)
testthat::expect_true(!any(is.na(est$data)))

est <- modsem(m1, oneInt2, na.rm=FALSE)
testthat::expect_true(any(is.na(est$data)))

# Residual Centering Approach
testthat::expect_warning(modsem(m1, oneInt2, method = "rca"),
                         regex = "Removing missing values case-wise")

est <- modsem(m1, oneInt2, method = "rca", na.rm=TRUE)
testthat::expect_true(!any(is.na(est$data)))

est <- modsem(m1, oneInt2, method = "rca", na.rm=FALSE)
testthat::expect_true(any(is.na(est$data)))
