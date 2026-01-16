devtools::load_all()

m <- '
y1 ~ x1 + z1 + x1:z1
'

fit_qml <- modsem(m, oneInt, method = "qml")
fit_lms <- modsem(m, oneInt, method = "lms")
