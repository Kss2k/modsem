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
m1 <-'
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  X:Z =~ mean(x1, x2, x3):mean(z1, z2, z3)
  # Inner model
  Y ~ X + Z + X:Z

'


test <- '
# Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + start(0.1)*z2 + z3

  # Inner model
  Y ~ X + Z
'



# Manual

df1$meanX <- rowMeans(df1[grepl("^x", colnames(df1))])
df1$meanY <- rowMeans(df1[grepl("^y", colnames(df1))])
df1$meanZ <- rowMeans(df1[grepl("^z", colnames(df1))])

  # Real result using lm()
realModel <- lm(meanY ~ meanX*meanZ,df1)
summary(realModel)
  # Using modsem
observedModel <- '
meanY ~ meanX + equal("meanY~meanZ")*meanZ + meanX:meanZ
'
observed <- modsem(observedModel, df1, method = "pind")
summary(observed)


# Using parceling in modsem
parcelingModel <- '
mean(y1, y2, y3) ~ start(0.1)*mean(x1, x2, x3) + equal("MEAN_x1_x2_x3")*mean(z1, z2, z3) + mean(x1, x2, x3):mean(z1, z2, z3)

'
parcM <- modsem(parcelingModel, df1, method = "pind")
summary(parcM)
