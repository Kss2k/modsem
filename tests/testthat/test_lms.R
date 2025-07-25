devtools::load_all()
m1 <- "
# Outer Model
  X =~ x1
  Z =~ z1
  x1 ~~ 0.1 * x1
  Y =~ y1

# Inner model
  Y ~~ Y

  Y ~ a * X + c * Z + b * X:Z

  a == b ^ 2 + 0.01 * c

  a2 := sqrt(a - 0.0001)^2
  coef := sqrt(a * 2 + b^2)
"

est1 <- modsem(m1, oneInt, method = "lms", optimize = TRUE, verbose = TRUE,
               convergence.abs = 1e-2, R.max = 50000)
plot_interaction("X", "Z", "Y", "X:Z", -3:3, c(-0.5, 0.5), est1, standardized=TRUE)
print(summary(est1, adjusted.stat = TRUE))
plot_surface(x = "X", z = "Z", y = "Y", model = est1)
std_est1_mc <- standardize_model(est1, monte.carlo=TRUE)
std_est1_delta <- standardize_model(est1, monte.carlo=FALSE)
mdata <- modsem_predict(est1, center.data = FALSE)
ndata <- modsem_predict(est1, newdata = oneInt[1:100, ], center.data = FALSE)
testthat::expect_true(all(mdata[1:100, ] == ndata))

# test constraints and custom parameters
coefs <- coef(est1)
a <- coefs[["a"]]
b <- coefs[["b"]]
c <- coefs[["c"]]
testthat::expect_equal(a, b ^ 2 + 0.01 * c)

test_coefs <- function(est) {
  coefs <- coef(est)
  a <- coefs[["a"]]
  b <- coefs[["b"]]
  c <- coefs[["coef"]]
  testthat::expect_equal(c, sqrt(a * 2 + b^2))
}

test_coefs(est1)
test_coefs(std_est1_delta)
test_coefs(std_est1_mc)

m1 <- "
# Outer Model
  X =~ x1
  Z =~ z1
  x1 ~~ 0.1 * x1
  Y =~ y1

# Inner model
  Y ~ a * X + a * Z
  X ~~ varX * X
  Y ~~ Y
  Y ~ b * X:Z + 0.05 * X:X
  b == a * 1.2

  coef := sqrt(a * 2 + b^2)
"

testthat::expect_warning(
  modsem(m1, oneInt, method = "lms") ,
  regexp = "Variances and covariances .*"
)

testthat::expect_no_warning(
  modsem(m1, oneInt, method = "lms", cov.syntax = "")
)

# PROBLEM:
#   I have no clue why, but changing the ordering of how the interaction terms
#   are specified, ends up changing the number of iterations (and results ever
#   so slightly) -- even though the matrices are exactly the same. This can be
#   seen through the fact that the starting loglikelihoods are the same (if optimized)
#   indicating that the matrices are the same (i.e,. produce the same results, when
#   given the same values).
# ANSWER:
#   slightly different results from lavaan, giving slightly different
#   starting parameters
tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  LATENT_VAR_ATT =~ a1 * att1 + a2 * att2 + att3 + att4 + att5
  SN =~ s1 * sn1 + sn2
  PBC =~ p1 * pbc1 + pbc2 + pbc3
  INT =~ i1 * int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Covariances
  LATENT_VAR_ATT ~~ SN + PBC
  PBC ~~ SN
  # Causal Relationsships
  INT ~ gamma_int_att * LATENT_VAR_ATT + b * SN + b * PBC
  BEH ~ 0.2 * INT + a * PBC
  BEH ~ PBC:INT
  gamma_int_att == a
  p1 == 1 + 0.1
  a2 == 1
  s1 == 1
  i1 == 1
  my_custom_parameter := a * 2
"

covModel <- '
PBC ~ a * LATENT_VAR_ATT + SN
'

testthat::expect_warning({
  est2 <- modsem(tpb, TPB, method = "lms", verbose = TRUE, convergence.abs = 1,
                 cov.syntax = covModel, nodes = 16, robust.se = TRUE)
}, regexp = "It is recommended .* between endogenous variables .*")

ust_pt <- parameter_estimates(est2)
std_pt <- standardized_estimates(est2, monte.carlo=FALSE)
testthat::expect_true(all(is.na(ust_pt[ust_pt$label == "p1", "std.error"])))

# Check constraints
label <- "my_custom_parameter"
u_est <- ust_pt[ust_pt$label == label, "est"]
s_est <- std_pt[std_pt$label == label, "est"]
u_a   <- unique(ust_pt[ust_pt$label == "a", "est"])
s_a   <- unique(std_pt[std_pt$label == "a", "est"])
expect_true(round(u_est, 5) < round(s_est, 5))
testthat::expect_true(length(s_a) > 1)
testthat::expect_true(length(u_a) == 1)
testthat::expect_equal(u_est, 2 * u_a)
testthat::expect_equal(s_est, 2 * s_a[1])
testthat::expect_equal(
  unname(calcVarParTable(c("PBC", "INT", "BEH"), std_pt)), rep(1, 3)
)

plot_interaction(x = "INT", z = "PBC", y = "BEH", vals_z = c(-0.5, 0.5), model = est2)
print(summary(est2))
var_interactions(est2)
print(vcov(est2)[8:14, 8:14])
modsem_inspect(est2)
print(coef(est2))
coefficients(est2)
cat("Number of observations using 'nobs()':", nobs(est2), "\n")

tpb2 <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Covariances
  ATT ~~ SN + PBC
  PBC ~~ SN
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
'

testthat::expect_warning(modsem(tpb, TPB, method = "lms", convergence.abs = 1000,
                                nodes = 16, calc.se = FALSE),
                         regexp = "It is recommended .* between exogenous and endogenous .*")
modsem_predict(est2)


testthat::expect_true({
  parTable <- summary(est2, H0 = FALSE, standardized = TRUE,
                      intercepts = TRUE)$parTable
  all(is.na(parTable[parTable$op == "~1", "std.all"]))
})


testthat::expect_true({
  parTable <- summary(est2, H0 = FALSE, standardized = TRUE)$parTable
  all(is.na(parTable[parTable$op == "~1", "std.all"]))
})

tpb_uk <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att3 + att2 + att1 + att4
  SN =~ sn4 + sn2 + sn3 + sn1
  PBC =~ pbc2 + pbc1 + pbc3 + pbc4
  INT =~ int2 + int1 + int3 + int4
  BEH =~ beh3 + beh2 + beh1 + beh4

# Inner Model (Based on Steinmetz et al., 2011)
# Causal Relationships
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
'

testthat::expect_warning(
modsem(tpb_uk, TPB_UK, method = "lms", nodes=32, max.iter=10, calc.se=FALSE)
)


m1 <- "
# Outer Model
  X =~ x1
  Z =~ z1
  x1 ~~ 0.1 * x1
  Y =~ y1

# Inner model
  Y ~ X + Z
  Y ~ X:Z
"

testthat::expect_error(
  modsem(m1, oneInt, method = "lms", optimizer = "ssjj"),
  regexp = "*Model estimation failed!*"
)

m1 <- "
# Outer Model
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2

# Inner model
  Y ~ X + Z + X:Z + X:X
"

est3 <- modsem(m1, oneInt, method = "lms")
std <- standardized_estimates(est3)
testthat::expect_equal(unname(calcVarParTable("Y", std)), 1)
summarize_partable(std)


m1 <- "
# Outer Model
  X =~ x1 + x2
  Z =~ z1 + z2
  .Y =~ y1 + y2

# Inner model
  .Y ~ c * X + Z + X:Z

  X~~ a*X + b*Z

  a == 1
  b == .2
  c > 0.4
"

est_m1 <- modsem(m1, oneInt, method = "lms", cov.syntax = "")
coef_m1 <- coef(est_m1)
summary(est_m1)
print(modsem_inspect(est_m1, what = c("vcov.all", "theta", "lambda", "psi", "information", "coefficients.all")))

# test constraints
testthat::expect_equal(coef_m1[["a"]], 1)
testthat::expect_equal(coef_m1[["b"]], .2)
testthat::expect_equal(est_m1$model$info$bounds$lower[["c"]], 0.4)

m2 <- "
# Outer Model
  X =~ x1 + x2
  Z =~ z1 + z2
  .Y =~ y1 + y2

# Inner model
  .Y ~ c * X + Z + X:Z

  X~~ a*X + b*Z

  a == 1
  b == .2
  c > a * 2
"

testthat::expect_error(
  modsem(m2, oneInt, method = "lms", cov.syntax = ""),
  regexp = "Dynamic constraints .*"
)
