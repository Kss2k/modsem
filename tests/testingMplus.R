

m1 <- '
# Outer Model
X =~ x1
X =~x2 +x3
Z =~ z1 + z2 + z3
Y =~ y1 + y2 + y3


# Inner model
Y ~ X + Z
Y ~ X:Z
'

mplus <- modsem(m1, oneInt, method = "mplus")
