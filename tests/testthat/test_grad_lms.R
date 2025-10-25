devtools::load_all()


m1 <- '
# Outer Model
  X =~ x1 + lx2 * x2 + lx3 * x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z

# Constraints
  x1 ~~ vx1*x1
  x2 ~~ vx2*x2
  x3 ~~ vx3*x3

  proj_vx1 := 0.8
  proj_vx2 := (lx2 ^ 2) * 0.8
  proj_vx3 := (lx3 ^ 2) * 0.8

  vx1 == 1 - proj_vx1
  vx2 == 1 - proj_vx2
  vx3 == 1 - proj_vx3
'

lms.no.opt <- modsem(m1, oneInt, method = "lms",
                     optimize = FALSE,
                     standardize.data = TRUE,
                     mean.observed = FALSE)
lms.opt <- modsem(m1, oneInt, method = "lms",
                  optimize = TRUE,
                  standardize.data = TRUE,
                  mean.observed = FALSE)

testthat::expect_equal(coef(lms.no.opt),
                       coef(lms.opt), tol = 2e-4)

m2 <- '
# Outer Model
  X =~ x1 + lx2 * x2 + lx3 * x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z

# Constraints
  x1 ~~ vx1*x1
  x2 ~~ vx2*x2
  x3 ~~ vx3*x3

  X ~ proj_vx1 * 1
  Z ~ proj_vx2 * 1
  Y ~ proj_vx3 * 1

  x1 ~ mu_x1 * 1
  z1 ~ mu_z1 * 1
  y1 ~ mu_y1 * 1

  proj_vx1 == 0.8
  proj_vx2 == (lx2 ^ 2) * 0.8
  proj_vx3 == (lx3 ^ 2) * 0.8

  mu_x1 == -proj_vx1
  mu_z1 == -proj_vx2
  mu_y1 == -proj_vx3

  vx1 == 1 - proj_vx1
  vx2 == 1 - proj_vx2
  vx3 == 1 - proj_vx3
'

lms.no.opt <- modsem(m2, oneInt, method = "lms",
                     optimize = FALSE,
                     standardize.data = TRUE)
lms.opt <- modsem(m2, oneInt, method = "lms",
                  optimize = TRUE,
                  standardize.data = TRUE)

testthat::expect_equal(coef(lms.no.opt),
                       coef(lms.opt), tol = 3e-4)


m3 <- '
# Outer Model
  X =~ x1 + lx2 * x2 + lx3 * x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z

# Constraints
  x1 ~~ a*x1
  x2 ~~ aa*x2
  x3 ~~ aaa*x3

  X ~ aaa.a   * 1
  Z ~ aaa.aa  * 1
  Y ~ aaa.aaa * 1

  x1 ~ aaaaaa.a   * 1
  z1 ~ aaaaaa.aa  * 1
  y1 ~ aaaaaa.aaa * 1

  aaa.a   == 0.8
  aaa.aa  == (lx2 ^ 2) * 0.8
  aaa.aaa == (lx3 ^ 2) * 0.8

  aaaaaa.a   == -aaa.a
  aaaaaa.aa  == -aaa.aa
  aaaaaa.aaa == -aaa.aaa

  a   == 1 - aaa.a
  aa  == 1 - aaa.aa
  aaa == 1 - aaa.aaa
'

lms.no.opt <- modsem(m3, oneInt, method = "lms",
                     optimize = FALSE,
                     standardize.data = TRUE)
lms.opt <- modsem(m3, oneInt, method = "lms",
                  optimize = TRUE,
                  standardize.data = TRUE)

testthat::expect_equal(coef(lms.no.opt),
                       coef(lms.opt), tol = 3e-4)
