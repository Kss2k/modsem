devtools::load_all()
m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + X:Z
'

est1 <- modsem(m1, data = oneInt,
               method = "qml", maxIter = 10000)
