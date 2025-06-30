devtools::load_all()

m1 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner model
  Y ~ X + c * Z

  y1 ~~ 10000*y1
"

testthat::expect_warning(modsem(m1, oneInt, method = "lms", calc.se=FALSE), 
                         regexp = "Some estimated ov variances.*")

m2 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ 0.0001 * y1 + y2 + y3

# Inner model
  Y ~ X + c * Z

  Y ~~ 10000*Y
"

testthat::expect_warning(modsem(m2, oneInt, method = "lms", calc.se=FALSE), 
                         regexp = "Some estimated lv variances.*")

m3 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner model
  Y ~ X + c * Z
"

# we have to be a bit fancy in how we catch the errors here, since multiple
# are emitted...
quiet_modsem <- purrr::quietly(modsem)
results_lms <- quiet_modsem(m3, oneInt, method = "lms", auto.fix.first = FALSE)
pat1 <- "Standard errors for some .* could not be computed.*"
pat2 <- "The variance-covariance.*positive definite.*eigenvalue.*zero.*"
testthat::expect_true(any(grepl(pat1, results_lms$warnings)))
testthat::expect_true(any(grepl(pat2, results_lms$warnings)))

