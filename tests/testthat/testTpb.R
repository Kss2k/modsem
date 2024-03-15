
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
  LBEH ~ LINT:LPBC  
'

method <- c("ca", "rca", "uca", "dblcent")
ests <- lapply(method, function(m) modsem(tpb, data = TPB, method = m))
estsMatch <- lapply(method, function(m) modsem(tpb, data = TPB, method = m,
                                               match = TRUE))
