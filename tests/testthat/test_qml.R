set.seed(123)
devtools::load_all()
tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  LATT =~ att1 + att2 + att3 + att4 + att5
  LSN =~ sn1 + sn2
  LPBC =~ pbc1 + pbc2 + pbc3
  LINT =~ int1 + int2 + int3
  LBEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Covariances
  LATT ~~ cAsn * LSN + cApbc * LPBC
  LPBC ~~ cPbcSn * LSN 
  # Causal Relationsships
  LINT ~ gIa * LATT + gIsn * LSN + gIpbc * LPBC
  LBEH ~ LINT + LPBC 
  LBEH ~ LPBC:LATT
'
esttpb <- modsem(tpb, data = TPB, method = "qml")
summary(esttpb)

m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + X:Z
'

est1 <- modsem(m1, data = oneInt, convergence = 1e-2, method = "qml")
print(summary(est1))
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = est1)
