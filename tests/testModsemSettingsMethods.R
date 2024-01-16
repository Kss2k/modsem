#library(modsem)
devtools::load_all(".")
m1 <- '
# Outer Model
X =~ x1
X =~x2 +x3
Z =~ z1 + z2 + z3
Y =~ y1 + y2 + y3
X:Z =~ x1:z1 + z2:x2

# Inner model
Y ~ X + Z
Y ~ X:Z
'

caManual <- modsem(m1, oneInt, method = "pind", standardizeData = TRUE,
                   residualCovSyntax = TRUE, centerBefore = TRUE,
                   constrainedResCovMethod = "ca", constrainedVar = TRUE,
                   constrainedProdMean = TRUE, constrainedLoadings = TRUE)
summary(caManual)
caAuto <- modsem(m1, oneInt, method = "ca")
summary(caAuto)
caMod <- modsem(m1, simTest, method = "ca", constrainedResCovMethod = "simple",
                centerData = TRUE, centerAfter = TRUE, constrainedProdMean = FALSE)
summary(caMod)


 rcaManual <- modsem(m1, oneInt, method = "pind", residualsProds = TRUE)
summary(rcaManual)
rcaAuto <- modsem(m1, oneInt, method = "rca")
summary(rcaAuto)
rcaMod <- modsem(m1, oneInt, method = "rca", constrainedResCovMethod = "simple")
summary(rcaMod)


