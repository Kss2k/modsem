library(lavaan)
library(modsem)
library(semTools)

setwd("C:/Users/slupp/OneDrive/Skrivebord/MasterOppgaveMehmet/data")
df1 <- readRDS("exampleData1.rds")
df2 <- readRDS("exampleData4.rds")
df3 <- readRDS("exampleData3.rds")

# Testing one interaction
m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3

  # Inner model
  Y ~ X + Z + X:Z
'

scaledDf1 <- lapplyDf(df1, scale)

realModel <- lm(realY ~ realX*realZ,scaledDf1)
summary(realModel)
m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  X:Z =~ mean(x1, x2, x3):mean(z1, z2, z3)
  # Inner model
  Y ~ X + Z + X:Z
'

m1Rca <- modsem(m1, df1, method = "rca",
                standardizeData = FALSE, isMeasurementSpecified = TRUE)
summary(m1Rca)


m1Uca <- modsem(m1, df1, method = "uca", standardizeData = TRUE)
summary(m1Uca)
