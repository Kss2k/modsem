devtools::load_all()
m1 <- '
# Outer Model
  X =~ x1
  X =~ x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner model
  Y ~ X + Z
  Y ~ X:Z
'
# funnily enough, the starting parameters from the double centering approach
# give better loglikelihoods than the ones arrived at by the EM algorithm
# i.e., the loglikelihood decreases from the starting parameters
startTime1 <- Sys.time()
est1 <- modsem(m1, oneInt, method = "lms", 
               optimize = TRUE,
               convergence = 1e-2, sampleGrad = NULL, maxstep = 1)
duration1 <- Sys.time() - startTime1

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

startTime2 <- Sys.time()
est2 <- modsem(tpb, TPB, method = "lms", optimize = TRUE, 
               convergence = 1e-2, sampleGrad = NULL,
               maxstep = 1)
duration2 <- Sys.time() - startTime2


