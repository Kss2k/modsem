devtools::load_all()
m1 <- "
# Outer Model
  X =~ x1
  Z =~ z1 
  x1 ~~ 0.1 * x1
  Y =~ y1

# Inner model
  Y ~ a * X + a * Z
  Y ~~ Y
  Y ~ b * X:Z + 0.05 * X:X
  b == a * 1.2

  coef := a * 2 + b^2
"

est1 <- modsem(m1, oneInt, method = "lms", optimize = TRUE, verbose = TRUE,
               convergence = 1e-2, R.max = 50000)
plot_interaction("X", "Z", "Y", "X:Z", -3:3, c(-0.5, 0.5), est1, standardized=TRUE)
print(summary(est1, adjusted.stat = TRUE))
plot_surface(x = "X", z = "Z", y = "Y", model = est1)
standardize_model(est1, monte.carlo=TRUE)
standardize_model(est1, monte.carlo=FALSE)

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
  est2 <- modsem(tpb, TPB, method = "lms", verbose = TRUE, convergence = 1, 
                 cov.syntax = covModel, nodes = 16, robust.se = TRUE)
}, regexp = "It is recommended .* between endogenous variables .*")

ust_pt <- parameter_estimates(est2)
std_pt <- standardized_estimates(est2, monte.carlo=FALSE)
testthat::expect_true(all(is.na(ust_pt[ust_pt$label == "p1", "std.error"])))

label <- "my_custom_parameter"
u_est <- ust_pt[ust_pt$label == label, "est"]
s_est <- std_pt[std_pt$label == label, "est"]
expect_true(round(u_est, 5) < round(s_est, 5))

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

testthat::expect_warning(modsem(tpb, TPB, method = "lms", convergence = 1000,
                                nodes = 16, calc.se = FALSE),
                         regexp = "It is recommended .* between exogenous and endogenous .*")


testthat::expect_true({
  parTable <- summary(est2, H0 = FALSE, standardized = TRUE,
                      intercepts = TRUE)$parTable
  all(parTable[parTable$op == "~1", "est"] == 0)
})


testthat::expect_true({
  parTable <- summary(est2, H0 = FALSE, standardized = TRUE)$parTable
  length(parTable[parTable$op == "~1", "est"]) == 0
})

tpb_uk <- '
# Outer Model (Based on Hagger et al., Citation2007)
  ATT =~ att3 + att2 + att1 + att4
  SN =~ sn4 + sn2 + sn3 + sn1
  PBC =~ pbc2 + pbc1 + pbc3 + pbc4
  INT =~ int2 + int1 + int3 + int4
  BEH =~ beh3 + beh2 + beh1 + beh4

# Inner Model (Based on Steinmetz et al., Citation2011)
# Causal Relationships
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
'

testthat::expect_warning(
modsem(tpb_uk, TPB_UK, method = "lms", nodes=32, max.iter=10, calc.se=FALSE)
)


testthat::expect_error(
  modsem(m1, oneInt, method = "lms", optimizer = "ssjj"),
  regexp = "*Model estimation failed!*"
)
