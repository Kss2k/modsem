setwd("C:/Users/slupp/OneDrive/Skrivebord/MasterOppgaveMehmet/data")
df1 <- readRDS("exampleData5correlatedPredictors.rds")

model <- '
X =~ x1 + x2 + x3
Z =~ z1 + z2 + z3
Y =~ y1 + y2 + y3

Y ~ X + Z + X:Z
'

# Real Model -------------------------------------------------------------------
centeredDf <- lapplyDf(df1, function(x) x - mean(x))
realModel <- lm(realY ~ realX*realZ, centeredDf)
summary(realModel)


# Recidual Centering -----------------------------------------------------------
rca <- modsem(model, df1, method = "rca")
summary(rca)


# Unconstrained approach -------------------------------------------------------
uca <- modsem(model, df1, method = "uca")
summary(uca)


# Basic Product indicator approach ---------------------------------------------
  # The regression approach struggles if the data isnt standardized/centered
pind <- modsem(model, df1, method = "pind")
summary(pind)

pind <- modsem(model, df1, method = "pind", centerData = TRUE)
summary(pind)

# LMS
lms <- modsem(model, centeredDf, method = "lms")
summary(lms)
