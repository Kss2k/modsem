devtools::load_all()
m1 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner model
  Y ~ X + c * Z
"

est_m1_lms <- modsem(m1, oneInt, method = "lms", calc.se=FALSE)
testthat::expect_true(est_m1_lms$iterations == 2L) # lowest possible number of iterations for lms
testthat::expect_warning(summary(est_m1_lms), "Comparative fit to H0 will not be calculated.")

est_m1_qml <- modsem(m1, oneInt, method = "qml", calc.se=FALSE)
testthat::expect_true(est_m1_qml$iterations == 1L) # lowest possible number of iterations for qml
testthat::expect_warning(summary(est_m1_qml), "Comparative fit to H0 will not be calculated.")


tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + a2 * att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEHAVIOUR_VARIABLE_LONG =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ a * ATT + a * SN + a * PBC
  BEHAVIOUR_VARIABLE_LONG ~ INT + b * PBC

  INT ~~ HELLO_THERE * INT
  BEHAVIOUR_VARIABLE_LONG ~~ HELLO_THERE2 * BEHAVIOUR_VARIABLE_LONG
"

est_tpb_lms <- modsem(tpb, TPB, method = "lms", calc.se=TRUE)
testthat::expect_true(est_tpb_lms$iterations == 2L)
testthat::expect_warning(summary(est_tpb_lms), "Comparative fit to H0 will not be calculated.")
print(summary(est_tpb_lms, H0 = FALSE, standardized = TRUE))

est_tpb_qml <- modsem(tpb, TPB, method = "qml", calc.se=FALSE)
testthat::expect_true(est_tpb_qml$iterations <= 1L)
testthat::expect_warning(summary(est_tpb_qml), "Comparative fit to H0 will not be calculated.")


tpb_main <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ a1 * att1 + a2 * att2 + att3 + att4 + att5
  SN =~ s1 * sn1 + sn2
  PBC =~ p1 * pbc1 + pbc2 + pbc3
  INT =~ i1 * int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + a * PBC
"

tpb_cov <- "
  INT ~ ATT + b * SN + c * PBC
"

est_split_lms <- modsem(tpb_main, data = TPB, method = "lms",
                           calc.se=FALSE, cov.syntax = tpb_cov)
testthat::expect_true(est_split_lms$iterations == 2)
testthat::expect_warning(summary(est_split_lms), "Comparative fit to H0 will not be calculated.")

est_split_qml <- modsem(tpb_main, data = TPB, method = "qml",
                           calc.se=FALSE, cov.syntax = tpb_cov)
testthat::expect_true(est_split_qml$iterations == 1)
testthat::expect_warning(summary(est_split_qml), "Comparative fit to H0 will not be calculated.")


dupsn1 <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ a1 * att1 + a2 * att2 + att3 + att4 + att5
  SN =~ s1 * sn1 + sn2
  PBC =~ p1 * pbc1 + pbc2 + pbc3 + sn1
  INT =~ i1 * int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + b * SN + c * PBC
  BEH ~ INT + a * PBC
"

est_dupsn1_lms <- modsem(dupsn1, TPB, method = "lms", calc.se=FALSE)
testthat::expect_true(est_dupsn1_lms$iterations == 2)
testthat::expect_warning(summary(est_tpb_lms), "Comparative fit to H0 will not be calculated.")

est_dupsn1_qml <- modsem(dupsn1, TPB, method = "qml", calc.se=FALSE)
testthat::expect_true(est_dupsn1_qml$iterations == 1)
testthat::expect_warning(summary(est_tpb_qml), "Comparative fit to H0 will not be calculated.")


tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  SN ~  ATT + PBC
  INT ~ ATT + PBC
  BEH ~ PBC
'
tpb_lms_2 <- modsem(tpb, TPB, method = "lms")
testthat::expect_true(tpb_lms_2$iterations == 2)
summary(tpb_lms_2, H0 = FALSE)


tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  att4 ~~ att5
  att4 ~~ int1

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + PBC + SN
  BEH ~ PBC
'

tpb_lms_3 <- modsem(tpb, TPB, method = "lms")
testthat::expect_true(tpb_lms_3$iterations == 2)
summary(tpb_lms_3, H0 = FALSE)


tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + a2 * att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC + INT:PBC
"

syntax_pi <- get_pi_syntax(tpb, method = "dblcent")
data_pi   <- get_pi_data(tpb, data = TPB, method = "dblcent")
est_lms_pi <- modsem(syntax_pi, data_pi, method = "lms", calc.se = FALSE)
est_qml_pi <- modsem(syntax_pi, data_pi, method = "qml", calc.se = FALSE)
testthat::expect_true(est_lms_pi$iterations == 2)
testthat::expect_true(est_qml_pi$iterations == 1)
