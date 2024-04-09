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
# i.e., the loglikelyhodd decreases from the starting parameters
est1 <- modsem(m1, oneInt, method = "lms", 
               optimize = TRUE, verbose = FALSE,
               suppressWarnings = TRUE)
