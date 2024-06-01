set.seed(123)
devtools::load_all()
m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + X:Z
'

est1 <- modsem(m1, data = oneInt, convergence = 1e-2, method = "qml")
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = est1)
