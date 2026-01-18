devtools::load_all()

m <- '
y1 ~ x1 + z1 + x1:z1
'

fit_qml <- modsem(m, oneInt, method = "qml")
fit_lms <- modsem(m, oneInt, method = "lms")
summary(fit_lms, standardized = TRUE)


tpb <- '
int1 ~ att1 + sn1 + pbc1
b1 ~ int1 + pbc1 + int1:pbc1
'

fit_qml <- modsem(tpb, TPB, method = "qml")
fit_lms <- modsem(tpb, TPB, method = "lms")
summary(fit_lms, standardized = TRUE)
