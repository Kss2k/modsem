m1 <- '
# Outer Model
  X =~ x1
  X =~x2
  Z =~ z1 + z2
  Y =~ y1 + y2 

# Inner model
  Y ~ X + Z
  Y ~ X:Z
'

est1 <- modsem(m1, oneInt, method = "lms", optimize = TRUE)
