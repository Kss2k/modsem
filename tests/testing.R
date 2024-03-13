library(lavaan)
library(semTools)
library(tidyverse)
devtools::load_all()
# One Interaction --------------------------------------------------------------
  # Real Model
m1 <- '
# Outer Model
X =~ x1
X =~x2 +x3
Z =~ z1 + z2 + z3
Y =~ y1 + y2 + y3

# Inner model
Y ~ X + Z
Y ~ X:Z
X ~~ Z
X ~~ X
Z ~~ Z
Y ~~ Y
'

m2 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  X:Z =~ x1:z1 + x2:z2 + x3:z3 + x2:z3

  # Inner model
  Y ~ X + Z + X:Z
'
startVals <- c(
Lambda.x2  = 0.778,
Lambda.x3  = 0.888,
Lambda.x11 = 0.813,
Lambda.x12 = 0.891,
Lambda.y2  = 0.801,
Lambda.y3  = 0.905,
Gamma1 = 0.741,
Gamma2 = 0.473,
Theta.d1  = 0.157,
Theta.d8  = 0.158,
Theta.d15 = 0.170,
Theta.d22 = 0.174,
Theta.d29 = 0.159,
Theta.d36 = 0.154,
Theta.e1  = 0.178,
Theta.e5  = 0.159,
Theta.e9  = 0.156,
Psi = 3.302 ,
Phi1 = 1.057,
Phi2 = 0.238 ,
Phi4 = 1.013 ,
nu.x2 = 0,
nu.x3 = 0,
nu.x5 = 0, 
nu.x6 = 0,
nu.y2 = 0,
nu.y3 = 0, 
alpha = 0,
tau1 = 0.0,
tau2 = 0,
Omega3 = 0.633
)

eLms <- modsem(m1, oneInt, "lms", qml = FALSE, max.iter = 0, convergence = 0.5,
               optimize = TRUE)


scaledoneInt <- lapplyDf(oneInt, scale)
realModel <- lm(realY ~ realX*realZ, scaledoneInt)
summary(realModel)

  # RCA ------------------------------------------------------------------------
    # Using modsem()
m1Rca <- modsem(m1, oneInt, method = "rca", centerData = TRUE)
summary(m1Rca)
coefs <- m1Rca$coefParTable
coefs |>
  filter(!lhs %in% coefs$rhs[coefs$lhs == "XZ"] & !rhs %in% coefs$rhs[coefs$lhs == "XZ"] )

startVals2 <- 
    # Using semtools
oneIntST <- indProd(scaledoneInt,
                 var1 = c("x1", "x2", "x3"),
                 var2 = c("z1", "z2", "z3"),
                 match = FALSE,
                 residualC = TRUE)


m1semTools <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  X.Z =~ x1.z1 + x1.z2 + x1.z3 + x2.z1 + x2.z2 + x2.z3 + x3.z1 + x3.z2 + x3.z3
  # Inner model
  Y ~ X + Z + X.Z
  # Residual (Co)Variances: XZ
  x1.z1 ~~ x2.z1
  x1.z1 ~~ x3.z1
  x1.z1 ~~ x1.z2
  x1.z1 ~~ x1.z3
  x2.z1 ~~ x3.z1
  x2.z1 ~~ x2.z2
  x2.z1 ~~ x2.z3
  x3.z1 ~~ x3.z2
  x3.z1 ~~ x3.z3
  x1.z2 ~~ x2.z2
  x1.z2 ~~ x3.z2
  x1.z2 ~~ x1.z3
  x2.z2 ~~ x3.z2
  x2.z2 ~~ x2.z3
  x3.z2 ~~ x3.z3
  x1.z3 ~~ x2.z3
  x1.z3 ~~ x3.z3
  x2.z3 ~~ x3.z3
  x1.z1 ~~ 0*x2.z2
  x1.z1 ~~ 0*x3.z2
  x1.z1 ~~ 0*x2.z3
  x1.z1 ~~ 0*x3.z3
s x2.z1 ~~ 0*x1.z2
  x2.z1 ~~ 0*x3.z2
  x2.z1 ~~ 0*x1.z3
  x2.z1 ~~ 0*x3.z3
  x3.z1 ~~ 0*x1.z2
  x3.z1 ~~ 0*x2.z2
  x3.z1 ~~ 0*x1.z3
  x3.z1 ~~ 0*x2.z3
  x1.z2 ~~ 0*x2.z3
  x1.z2 ~~ 0*x3.z3
  x2.z2 ~~ 0*x1.z3
  x2.z2 ~~ 0*x3.z3
  x3.z2 ~~ 0*x1.z3
  x3.z2 ~~ 0*x2.z3



'

m1RcaST <- sem(m1semTools, oneIntST)
summary(m1RcaST)

  # UCA ------------------------------------------------------------------------
    # modsem
m1Uca <- modsem(m1, oneInt, method = "uca", standardizeData = TRUE)
summary(m1Uca)

    # semtools
oneIntSTuca <- indProd(scaledoneInt,
                    var1 = c("x1", "x2", "x3"),
                    var2 = c("z1", "z2", "z3"),
                    match = FALSE,
                    meanC = TRUE,
                    doubleMC = FALSE)


m1semTools <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  X.Z =~ x1.z1 + x1.z2 + x1.z3 + x2.z1 + x2.z2 + x2.z3 + x3.z1 + x3.z2 + x3.z3
  # Inner model
  Y ~ X + Z + X.Z
  # Residual (Co)Variances: XZ
  x1.z1 ~~ x2.z1
  x1.z1 ~~ x3.z1
  x1.z1 ~~ x1.z2
  x1.z1 ~~ x1.z3
  x2.z1 ~~ x3.z1
  x2.z1 ~~ x2.z2
  x2.z1 ~~ x2.z3
  x3.z1 ~~ x3.z2
  x3.z1 ~~ x3.z3
  x1.z2 ~~ x2.z2
  x1.z2 ~~ x3.z2
  x1.z2 ~~ x1.z3
  x2.z2 ~~ x3.z2
  x2.z2 ~~ x2.z3
  x3.z2 ~~ x3.z3
  x1.z3 ~~ x2.z3
  x1.z3 ~~ x3.z3
  x2.z3 ~~ x3.z3
  x1.z1 ~~ 0*x2.z2
  x1.z1 ~~ 0*x3.z2
  x1.z1 ~~ 0*x2.z3
  x1.z1 ~~ 0*x3.z3
  x2.z1 ~~ 0*x1.z2
  x2.z1 ~~ 0*x3.z2
  x2.z1 ~~ 0*x1.z3
  x2.z1 ~~ 0*x3.z3
  x3.z1 ~~ 0*x1.z2
  x3.z1 ~~ 0*x2.z2
  x3.z1 ~~ 0*x1.z3
  x3.z1 ~~ 0*x2.z3
  x1.z2 ~~ 0*x2.z3
  x1.z2 ~~ 0*x3.z3
  x2.z2 ~~ 0*x1.z3
  x2.z2 ~~ 0*x3.z3
  x3.z2 ~~ 0*x1.z3
  x3.z2 ~~ 0*x2.z3

  # Mean constraint
  X ~~ Cov_X_Z*Z
  X.Z ~ Cov_X_Z*1
'


m1semTools2 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  X.Z =~ x1.z1 + x1.z2 + x1.z3 + x2.z1 + x2.z2 + x2.z3 + x3.z1 + x3.z2 + x3.z3
  # Inner model
  Y ~ X + Z + X.Z
  # Residual (Co)Variances: XZ

  x1.z1 ~~ share_x1*x1.z2 + share_x1*x1.z3
  x1.z2 ~~ share_x1*x1.z3

  x2.z1 ~~ share_x2*x2.z2 + share_x2*x2.z3
  x2.z2 ~~ share_x2*x2.z3

  x3.z1 ~~ share_x3*x3.z2 + share_x3*x3.z3
  x3.z2 ~~ share_x3*x3.z3

  x1.z1 ~~ share_z1*x2.z1 + share_z1*x3.z1
  x2.z1 ~~ share_z1*x3.z1

  x1.z2 ~~ share_z2*x2.z2 + share_z2*x3.z2
  x2.z2 ~~ share_z2*x3.z2

  x1.z3 ~~ share_z3*x2.z3 + share_z3*x3.z3
  x2.z3 ~~ share_z3*x3.z3

  x1.z1 ~~ 0*x2.z2
  x1.z1 ~~ 0*x3.z2
  x1.z1 ~~ 0*x2.z3
  x1.z1 ~~ 0*x3.z3

  x2.z1 ~~ 0*x1.z2
  x2.z1 ~~ 0*x3.z2
  x2.z1 ~~ 0*x1.z3
  x2.z1 ~~ 0*x3.z3

  x3.z1 ~~ 0*x1.z2
  x3.z1 ~~ 0*x2.z2
  x3.z1 ~~ 0*x1.z3
  x3.z1 ~~ 0*x2.z3

  x1.z2 ~~ 0*x2.z3
  x1.z2 ~~ 0*x3.z3

  x2.z2 ~~ 0*x1.z3
  x2.z2 ~~ 0*x3.z3

  x3.z2 ~~ 0*x1.z3
  x3.z2 ~~ 0*x2.z3

  # Mean constraint
  X ~~ Cov_X_Z*Z
  X.Z ~ Cov_X_Z*1
'

m1UcaST <- sem(m1semTools2, oneIntSTuca)
summary(m1UcaST)

  # CA
m1Ca <- modsem(m1, oneInt, method = "ca", standardizeData = TRUE)
summary(m1Ca)

  #Mplus
m1Mplus <- modsem(m1, oneInt, method = "mplus", standardizeData = TRUE)
  # UCA
  # Product Indicator/ regression approach --------------------------------------
m1Reg <- modsem(m1, oneInt, method = "pind", standardizeData = TRUE)
summary(m1Reg)

  # double centering approach --------------------------------------------------
    # Using modsem
m1DblCent <- modsem(m1, oneInt, method = "dblcent", standardizeData = TRUE)
summary(m1DblCent)

    # Using semtools
### using semTools
oneIntSTdbl <- indProd(scaledoneInt,
                    var1 = c("x1", "x2", "x3"),
                    var2 = c("z1", "z2", "z3"),
                    match = FALSE,
                    meanC = TRUE,
                    doubleMC = TRUE)


m1semTools <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  X.Z =~ x1.z1 + x1.z2 + x1.z3 + x2.z1 + x2.z2 + x2.z3 + x3.z1 + x3.z2 + x3.z3
  # Inner model
  Y ~ X + Z + X.Z
  # Residual (Co)Variances: XZ
  x1.z1 ~~ x2.z1
  x1.z1 ~~ x3.z1
  x1.z1 ~~ x1.z2
  x1.z1 ~~ x1.z3
  x2.z1 ~~ x3.z1
  x2.z1 ~~ x2.z2
  x2.z1 ~~ x2.z3
  x3.z1 ~~ x3.z2
  x3.z1 ~~ x3.z3
  x1.z2 ~~ x2.z2
  x1.z2 ~~ x3.z2
  x1.z2 ~~ x1.z3
  x2.z2 ~~ x3.z2
  x2.z2 ~~ x2.z3
  x3.z2 ~~ x3.z3
  x1.z3 ~~ x2.z3
  x1.z3 ~~ x3.z3
  x2.z3 ~~ x3.z3
  x1.z1 ~~ 0*x2.z2
  x1.z1 ~~ 0*x3.z2
  x1.z1 ~~ 0*x2.z3
  x1.z1 ~~ 0*x3.z3
  x2.z1 ~~ 0*x1.z2
  x2.z1 ~~ 0*x3.z2
  x2.z1 ~~ 0*x1.z3
  x2.z1 ~~ 0*x3.z3
  x3.z1 ~~ 0*x1.z2
  x3.z1 ~~ 0*x2.z2
  x3.z1 ~~ 0*x1.z3
  x3.z1 ~~ 0*x2.z3
  x1.z2 ~~ 0*x2.z3
  x1.z2 ~~ 0*x3.z3
  x2.z2 ~~ 0*x1.z3
  x2.z2 ~~ 0*x3.z3
  x3.z2 ~~ 0*x1.z3
  x3.z2 ~~ 0*x2.z3
'

m1dblST <- sem(m1semTools, oneIntSTdbl)
summary(m1dblST)

  # Two interactions -----------------------------------------------------------

    # Real model
m2 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  G =~ g1 + g2 + g3
  H =~ h1 + h2 + h3
  # Inner model
  Y ~ X + Z + G + H + X:Z + G:H
'

scaledtwoInt <- lapplyDf(twoInt, scale)
realModel2 <- lm(realY ~ realX*realZ + realG*realH, scaledtwoInt)
summary(realModel2)

  # RCA ------------------------------------------------------------------------
# m2 <- '
#   # Outer Model
#   X =~ x1 + x2 +x3
#   Y =~ y1 + y2 + y3
#   Z =~ z1 + z2 + z3
#   G =~ g1 + g2 + g3
#   H =~ h1 + h2 + h3
#   X:Z =~ x1:z1
#   # Inner model
#   Y ~ X + Z + G + H + X:Z + G:H + g1:x1
# '

latentModelRca2 <- modsem(m2, twoInt, method = "rca")
summary(latentModelRca2)

  # UCA
latentModelUncent2 <- modsem(m2, twoInt, method = "uca", standardizeData = TRUE)
summary(latentModelUncent2$lavaan)

  # Constrained Approach
latentModelCA <- modsem(m2, twoInt, method = "ca", standardizeData = TRUE)

# Parceling -----------------------------------------------------------



# RCA ------------------------------------------------------------------------
m3 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  X:Z =~ x1:z1 + x2:z2 + x3:z3
  H =~ mean(h1, h2, h3)
  G =~ mean(g1, g2, g3)
  # Inner model
  Y ~ X + Z + X:Z + G + H + G:H
'


latentModelRcaParceling <- modsem(m3, twoInt, method = "pind")
summary(latentModelRcaParceling)

latentModelRcaParceling2 <- modsem(m3, twoInt, method = "rca")
summary(latentModelRcaParceling2)

latentModelRcaParceling <- modsem(m3, twoInt, method = "rca")
summary(latentModelRcaParceling)
# UCA
latentModelUncentParceling <- modsem(m3, twoInt, method = "uca")
summary(latentModelUncentParceling)
