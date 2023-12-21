library(modsem)

m1 <- '
X =~ x1 + x2 + x3
Z =~ z1 + z2 + z3
Y =~ y1 + y2 + y3

'

m2 <- '
X =~ x1 + x2 + + x3
Z =~ z1 + 1 z2 + z3
Y =~ y1 + y2 + y3

'

m3 <-'
X =~ x1
X =~ x2 + x3
Y ~ X+Z
Y ~ X:Z + x1
Z =~ z1 + z2 + z3:v1
Y =~ y1 + y2 + y3
Y ~~ Y
Y >-= Y + 2 +3
'


tokenizeSyntax(m2)
tokenizeSyntax(m3)
modsemify(m3)
