devtools::load_all()

m <- '
y1 ~ x1 + z1 + x1:z1
'

testthat::expect_no_error({
  fit_qml <- modsem(m, oneInt, method = "qml")
  fit_lms <- modsem(m, oneInt, method = "lms")
  modsem_predict(fit_lms)
  summary(fit_lms, standardized = TRUE)
  estimate_h0(fit_lms)
})

tpb <- '
int1 ~ att1 + sn1 + pbc1
b1 ~ int1 + pbc1 + int1:pbc1
'

testthat::expect_no_error({
  fit_qml <- modsem(tpb, TPB, method = "qml")
  fit_lms <- modsem(tpb, TPB, method = "lms")
  summary(fit_lms, standardized = TRUE)
})

m2 <- '
Z =~ z1 + z2
y1 ~ x1 + Z + x1:Z
'

testthat::expect_no_error({
  # test selection of integrating variable
  fit_lms <- modsem(m2, oneInt, method = "lms")
})
