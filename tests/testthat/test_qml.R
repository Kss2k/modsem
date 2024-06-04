set.seed(123)
devtools::load_all()

m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + X:Z
'

est1 <- modsem(m1, data = oneInt, convergence = 1e-2, method = "qml")
print(summary(est1))
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = est1)

tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + PBC 
  BEH ~ PBC:INT
'

covModel <- ' 
INT ~ ATT + SN + PBC
'

esttpb <- modsem(tpb, data = TPB, method = "qml", cov_syntax = covModel)
plot_interaction(x = "INT", z = "PBC", y = "BEH", 
                 xz = "PBC:INT", vals_z = c(-0.5, 0.5), 
                 model = esttpb)
summary(esttpb)


tpbNoInt <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + PBC 
  INT ~ ATT + SN + PBC
'

estTpbNoInt <- modsem(tpbNoInt, data = TPB, method = "qml")
print(summary(estTpbNoInt, H0 = FALSE))
