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
Lambda.x2  = 0.778
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


# Testing "ca" 

caSyntax <- '
X =~ lambda_x1_X*x1
X =~ lambda_x2_X*x2
X =~ lambda_x3_X*x3
Z =~ lambda_z1_Z*z1
Z =~ lambda_z2_Z*z2
Z =~ lambda_z3_Z*z3
Y =~ lambda_y1_Y*y1
Y =~ lambda_y2_Y*y2
Y =~ lambda_y3_Y*y3
Y ~ Gamma_X_Y*X
Y ~ Gamma_Z_Y*Z
Y ~ Gamma_XZ_Y*XZ
XZ =~ lambda_x1z1_XZ*x1z1
XZ =~ lambda_x2z1_XZ*x2z1
XZ =~ lambda_x3z1_XZ*x3z1
XZ =~ lambda_x1z2_XZ*x1z2
XZ =~ lambda_x2z2_XZ*x2z2
XZ =~ lambda_x3z2_XZ*x3z2
XZ =~ lambda_x1z3_XZ*x1z3
XZ =~ lambda_x2z3_XZ*x2z3
XZ =~ lambda_x3z3_XZ*x3z3
X ~~ Var_X*X
Z ~~ Var_Z*Z
Y ~~ Zeta_Y*Y
XZ ~~ Var_XZ*XZ
x1 ~~ Var_x1*x1
x2 ~~ Var_x2*x2
x3 ~~ Var_x3*x3
z1 ~~ Var_z1*z1
z2 ~~ Var_z2*z2
z3 ~~ Var_z3*z3
y1 ~~ Var_y1*y1
y2 ~~ Var_y2*y2
y3 ~~ Var_y3*y3
x1z1 ~~ Var_x1z1*x1z1
x2z1 ~~ Var_x2z1*x2z1
x3z1 ~~ Var_x3z1*x3z1
x1z2 ~~ Var_x1z2*x1z2
x2z2 ~~ Var_x2z2*x2z2
x3z2 ~~ Var_x3z2*x3z2
x1z3 ~~ Var_x1z3*x1z3
x2z3 ~~ Var_x2z3*x2z3
x3z3 ~~ Var_x3z3*x3z3
X ~~ Cov_X_Z*Z
X ~~ Cov_X_XZ*XZ
Z ~~ Cov_Z_XZ*XZ
x2z1 ~~ Cov_x2z1_x1z1*x1z1
x3z1 ~~ Cov_x3z1_x1z1*x1z1
x3z1 ~~ Cov_x3z1_x2z1*x2z1
x1z2 ~~ Cov_x1z2_x1z1*x1z1
x1z2 ~~ 0*x2z1
x1z2 ~~ 0*x3z1
x2z2 ~~ 0*x1z1
x2z2 ~~ Cov_x2z2_x2z1*x2z1
x2z2 ~~ 0*x3z1
x2z2 ~~ Cov_x2z2_x1z2*x1z2
x3z2 ~~ 0*x1z1
x3z2 ~~ 0*x2z1
x3z2 ~~ Cov_x3z2_x3z1*x3z1
x3z2 ~~ Cov_x3z2_x1z2*x1z2
x3z2 ~~ Cov_x3z2_x2z2*x2z2
x1z3 ~~ Cov_x1z3_x1z1*x1z1
x1z3 ~~ 0*x2z1
x1z3 ~~ 0*x3z1
x1z3 ~~ Cov_x1z3_x1z2*x1z2
x1z3 ~~ 0*x2z2
x1z3 ~~ 0*x3z2
x2z3 ~~ 0*x1z1
x2z3 ~~ Cov_x2z3_x2z1*x2z1
x2z3 ~~ 0*x3z1
x2z3 ~~ 0*x1z2
x2z3 ~~ Cov_x2z3_x2z2*x2z2
x2z3 ~~ 0*x3z2
x2z3 ~~ Cov_x2z3_x1z3*x1z3
x3z3 ~~ 0*x1z1
x3z3 ~~ 0*x2z1
x3z3 ~~ Cov_x3z3_x3z1*x3z1
x3z3 ~~ 0*x1z2
x3z3 ~~ 0*x2z2
x3z3 ~~ Cov_x3z3_x3z2*x3z2
x3z3 ~~ Cov_x3z3_x1z3*x1z3
x3z3 ~~ Cov_x3z3_x2z3*x2z3
Cov_x2z1_x1z1 == lambda_x2_X * lambda_x1_X * (Var_X) * Var_z1
Cov_x3z1_x1z1 == lambda_x3_X * lambda_x1_X * (Var_X) * Var_z1
Cov_x3z1_x2z1 == lambda_x3_X * lambda_x2_X * (Var_X) * Var_z1
Cov_x1z2_x1z1 == lambda_z2_Z * lambda_z1_Z * (Var_Z) * Var_x1
Cov_x2z2_x2z1 == lambda_z2_Z * lambda_z1_Z * (Var_Z) * Var_x2
Cov_x2z2_x1z2 == lambda_x2_X * lambda_x1_X * (Var_X) * Var_z2
Cov_x3z2_x3z1 == lambda_z2_Z * lambda_z1_Z * (Var_Z) * Var_x3
Cov_x3z2_x1z2 == lambda_x3_X * lambda_x1_X * (Var_X) * Var_z2
Cov_x3z2_x2z2 == lambda_x3_X * lambda_x2_X * (Var_X) * Var_z2
Cov_x1z3_x1z1 == lambda_z3_Z * lambda_z1_Z * (Var_Z) * Var_x1
Cov_x1z3_x1z2 == lambda_z3_Z * lambda_z2_Z * (Var_Z) * Var_x1
Cov_x2z3_x2z1 == lambda_z3_Z * lambda_z1_Z * (Var_Z) * Var_x2
Cov_x2z3_x2z2 == lambda_z3_Z * lambda_z2_Z * (Var_Z) * Var_x2
Cov_x2z3_x1z3 == lambda_x2_X * lambda_x1_X * (Var_X) * Var_z3
Cov_x3z3_x3z1 == lambda_z3_Z * lambda_z1_Z * (Var_Z) * Var_x3
Cov_x3z3_x3z2 == lambda_z3_Z * lambda_z2_Z * (Var_Z) * Var_x3
Cov_x3z3_x1z3 == lambda_x3_X * lambda_x1_X * (Var_X) * Var_z3
Cov_x3z3_x2z3 == lambda_x3_X * lambda_x2_X * (Var_X) * Var_z3
Var_XZ == (Var_X) * (Var_Z) + ((Cov_X_Z) ^ 2)
Cov_X_XZ == 0
Cov_Z_XZ == 0
Var_x1z1 == (lambda_x1_X ^ 2) * (Var_X) * Var_z1 + (lambda_z1_Z ^ 2) * (Var_Z) * Var_x1 + Var_x1 * Var_z1
Var_x2z1 == (lambda_x2_X ^ 2) * (Var_X) * Var_z1 + (lambda_z1_Z ^ 2) * (Var_Z) * Var_x2 + Var_x2 * Var_z1
Var_x3z1 == (lambda_x3_X ^ 2) * (Var_X) * Var_z1 + (lambda_z1_Z ^ 2) * (Var_Z) * Var_x3 + Var_x3 * Var_z1
Var_x1z2 == (lambda_x1_X ^ 2) * (Var_X) * Var_z2 + (lambda_z2_Z ^ 2) * (Var_Z) * Var_x1 + Var_x1 * Var_z2
Var_x2z2 == (lambda_x2_X ^ 2) * (Var_X) * Var_z2 + (lambda_z2_Z ^ 2) * (Var_Z) * Var_x2 + Var_x2 * Var_z2
Var_x3z2 == (lambda_x3_X ^ 2) * (Var_X) * Var_z2 + (lambda_z2_Z ^ 2) * (Var_Z) * Var_x3 + Var_x3 * Var_z2
Var_x1z3 == (lambda_x1_X ^ 2) * (Var_X) * Var_z3 + (lambda_z3_Z ^ 2) * (Var_Z) * Var_x1 + Var_x1 * Var_z3
Var_x2z3 == (lambda_x2_X ^ 2) * (Var_X) * Var_z3 + (lambda_z3_Z ^ 2) * (Var_Z) * Var_x2 + Var_x2 * Var_z3
Var_x3z3 == (lambda_x3_X ^ 2) * (Var_X) * Var_z3 + (lambda_z3_Z ^ 2) * (Var_Z) * Var_x3 + Var_x3 * Var_z3
lambda_x1z1_XZ == lambda_x1_X * lambda_z1_Z
lambda_x2z1_XZ == lambda_x2_X * lambda_z1_Z
lambda_x3z1_XZ == lambda_x3_X * lambda_z1_Z
lambda_x1z2_XZ == lambda_x1_X * lambda_z2_Z
lambda_x2z2_XZ == lambda_x2_X * lambda_z2_Z
lambda_x3z2_XZ == lambda_x3_X * lambda_z2_Z
lambda_x1z3_XZ == lambda_x1_X * lambda_z3_Z
lambda_x2z3_XZ == lambda_x2_X * lambda_z3_Z
lambda_x3z3_XZ == lambda_x3_X * lambda_z3_Z
XZ ~ Mean_XZ*1
Mean_XZ == (Cov_X_Z)'
newDat <- modsem(m1, oneInt, method = "uca")$newData
caEst <- lavaan::sem(caSyntax, data = newDat)

