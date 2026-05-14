devtools::load_all()


tpb.lin <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT <~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH <~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
'

testthat::expect_no_condition({
fit_lav <- lavaan::sem(tpb.lin, TPB, optim.gradient = NULL)
fit_mod <- modsem(tpb.lin, TPB, method = "lms")
})

coef_lav <- lavaan::coef(fit_lav)
coef_mod <- coef(fit_mod)[names(coef_lav)]

msr_lav <- lavaan::fitMeasures(fit_lav)
msr_mod <- fit_modsem_da(fit_mod)

testthat::expect_equal(c(coef_lav), c(coef_mod), tol = 1e-3)
testthat::expect_equal(msr_lav[["chisq"]], msr_mod[["chisq.value"]], tol = 1e-3)
testthat::expect_equal( # modsem does not fix composite variances like lavaan
  msr_mod[["chisq.df"]] + msr_mod[["npar"]] - msr_lav[["npar"]] - 15, # -15 for mean structure
  msr_lav[["df"]]
)

tpb.int <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT <~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH <~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC + INT:PBC
'

testthat::expect_no_condition(
  modsem(tpb.int, TPB, method = "lms", nodes = 32)
)

testthat::expect_error(
  modsem(tpb.int, TPB, method = "qml"),
  regexp = "Composite constructs .*"
)
