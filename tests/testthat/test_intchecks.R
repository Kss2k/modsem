devtools::load_all()

tpb1 <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:pbc1
'

# pbc1 is an indicator of PBC, so it becomes an endogenous variable in the
# interaction term, which requires `auto.split.syntax = TRUE`
testthat::expect_error(modsem(tpb1, data = TPB, method = "lms"),
                       regexp = "two endogenous variables")

# NOTE: the SAM engine can't handle the interaction term with the (split) pbc1,
# so we use the linear engine for the starting parameters
fit.tpb1 <- modsem(tpb1, data = TPB, method = "lms", auto.split.syntax = TRUE,
                   nodes = 32, optimize.engine = "linear")
pe.tpb1 <- parameter_estimates(fit.tpb1)
testthat::expect_true(is.finite(pe.tpb1[pe.tpb1$rhs == "INT:pbc1" & pe.tpb1$op == "~", "est"]))

tpb2 <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC2
'

testthat::expect_error(modsem(tpb2, data = TPB, method = "lms"),
                       regexp = "Missing observed variables in data: .*PBC2.*")

testthat::expect_error(modsem(tpb2, data = TPB, method = "dblcent"),
                       regexp = "Missing variables in data: .*PBC2.*")

tpb3 <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  INT:pbc1 ~ INT
'

testthat::expect_error(
  suppressWarnings(modsem(tpb3, data = TPB, method = "lms")),
  regexp = ".*Interaction terms.*dependent.*"
)

