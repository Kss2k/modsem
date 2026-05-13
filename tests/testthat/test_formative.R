devtools::load_all()


tpb <- '
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


fit_lav <- lavaan::sem(tpb, TPB)
fit_mod <- modsem(tpb, TPB, method = "lms")
fit_modsem_da(fit_mod)

coef_lav <- lavaan::coef(fit_lav)
coef_mod <- coef(fit_mod)[names(coef_lav)]

testthat::expect_equal(c(coef_lav), c(coef_mod), tol = 1e-3)
