devtools::load_all()

m1 <- '
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z + a * X:Z
  c := 10 * a
'

lms1 <- modsem(m1, oneInt, method = "lms", adaptive.quad=TRUE)

tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
'

lms2 <- modsem(tpb, TPB, method = "lms", nodes = 32, adaptive.quad=TRUE)
summary(lms2)


tpb_uk <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2 + sn3
  PBC =~ pbc1 + pbc2 + pbc3 + pbc4
  INT =~ int1 + int2 + int3 + int4
  BEH =~ beh1 + beh2 + beh3

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
"

lms3 <- modsem(tpb_uk, data = TPB_UK, "lms",
               nodes=32, FIM="observed", robust.se = FALSE,
               adaptive.quad=TRUE, algorithm ="EMA", optimize = TRUE,
               convergence.abs = 1e-8, convergence.rel = 1e-15)
print(summary(lms3))
#> Regressions:
#>                  Estimate  Std.Error  z.value  P(>|z|)
#>   INT ~         
#>     ATT             0.283      0.032    8.915    0.000
#>     SN             -0.023      0.021   -1.102    0.271
#>     PBC             3.069      0.315    9.751    0.000
#>   BEH ~         
#>     PBC             0.675      0.175    3.855    0.000
#>     INT             0.658      0.041   16.170    0.000
#>     INT:PBC         0.390      0.044    8.808    0.000
# Compared with Mplus
#> Regressions:
#>                  Estimate  Std.Error  z.value  P(>|z|)
#>   INT ~         
#>     ATT    (INT<)   0.283      0.032    8.844    0.000
#>     SN  (INT<-SN)  -0.023      0.021   -1.095    0.271
#>     PBC    (INT<)   3.065      0.238   12.878    0.000
#>   BEH ~         
#>     INT    (BEH<)   0.658      0.041   16.049    0.000
#>     PBC    (BEH<)   0.674      0.168    4.012    0.000
#>     INT:PB (BEH<)   0.389      0.036   10.806    0.000


nlsem <- '
ENJ =~ enjoy1 + enjoy2 + enjoy3 + enjoy4 + enjoy5
CAREER =~ career1 + career2 + career3 + career4
SC =~ academic1 + academic2 + academic3 + academic4 + academic5 + academic6
CAREER ~ ENJ + SC + ENJ:ENJ + SC:SC + ENJ:SC
'

# For such a small number of nodes it doesn't really matter if you use an
# adaptive quadrature, as all the nodes bring some value
# no warning for low number of nodes here, when using adaptive quadrature
lms4 <- modsem(nlsem, data = jordan, method = "lms",
               adaptive.quad=TRUE, OFIM.hessian = FALSE,
               nodes = 15, mean.observed = FALSE)
summary(lms4)
