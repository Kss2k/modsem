devtools::load_all()
set.seed(123)

m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + a * X:Z
'

est1 <- modsem(m1, data = oneInt, convergence = 1e-2, method = "qml",
               robust.se = TRUE)
print(summary(est1, scientific = TRUE))
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = est1)
plot_surface(x = "X", z = "Z", y = "Y", model = est1)


tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  # Covariances
  ATT ~~ SN + PBC
  PBC ~~ SN 
# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ b * INT + PBC 
  INT ~ ATT + SN + a * PBC
  BEH ~ PBC:INT
  VERY_LONG_LABEL_THAT_SHOULD_BE_SHORTENED := a * b
'

est2 <- modsem(tpb, data = TPB, method = "qml", 
               robust.se = TRUE,
               standardize = TRUE, convergence = 1e-2)
print(summary(est2, H0 = FALSE))
expect_warning(plot_jn(x = "INT", z = "PBC", y = "BEH", model = est2,
                       min_z = -1.5, max_z = -0.5),
               regex = "Truncating.*right.*left.*")
expect_warning(plot_jn(x = "INT", z = "PBC", y = "BEH", model = est2,
                       min_z = -1.5, max_z = 2.5) ,
               regex = "Truncating SD-range on the left!")
expect_warning(plot_jn(x = "INT", z = "PBC", y = "BEH", model = est2,
                       min_z = -2, max_z = 1) ,
               regex = "Truncating SD-range on the right!")

# check that standardizing twice does not change the estimates
testthat::expect_equal(standardized_estimates(est2), 
                       parameter_estimates(est2))

# check correct standardization
calcCovParTable("BEH", "BEH", parameter_estimates(est2))[[1]] |>
  testthat::expect_equal(1)

# check that estimates are standardized the same way 
# in vcov, coef and parameter_estimates
pt <- parameter_estimates(est2)
vcov_sd <- sqrt(vcov(est2)["BEH~PBC:INT", "BEH~PBC:INT"])
cond <- pt$lhs == "BEH" & pt$op == "~" & pt$rhs == "PBC:INT"
pt_sd <- pt[cond, "std.error"]
pt_est <- pt[cond, "est"]
coef_est <- coef(est2)[["BEH~PBC:INT"]]

expect_equal(vcov_sd, pt_sd)
expect_equal(pt_est, coef_est)

modsem_inspect(est2) 
coefficients(est2)

# Observed Variables
m3 <- '
X =~ x1
Y =~ y1
Z =~ z1
Y ~ X + Z + X:Z
'

est3 <- modsem(m3, data = oneInt, convergence = 1e-2, method = "qml",
               robust.se = TRUE)
print(summary(est3, scientific = TRUE))
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = est3)

tpb2 <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 #+ int2 + int3
  BEH =~ b1  + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + PBC 
  INT ~ ATT + SN + PBC
  BEH ~ PBC:INT
'

est4 <- modsem(tpb2, data = TPB, method = "qml", 
               robust.se = TRUE,
               standardize = TRUE, convergence = 1e-2)
print(summary(est4, H0 = FALSE))

testthat::expect_equal(standardized_estimates(est4), 
                       parameter_estimates(est4))

calcCovParTable("BEH", "BEH", parameter_estimates(est4))[[1]] |>
  testthat::expect_equal(1)

vcov(est4) 
modsem_inspect(est4) 
coef(est4)
coefficients(est4)
