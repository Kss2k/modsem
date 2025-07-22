devtools::load_all()

m1 <- '
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z + X:Z
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
  ATT =~ att3 + att2 + att1 + att4
  SN =~ sn4 + sn2 + sn3 + sn1
  PBC =~ pbc2 + pbc1 + pbc3 + pbc4
  INT =~ int2 + int1 + int3 + int4
  BEH =~ beh3 + beh2 + beh1 + beh4

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
"

lms3 <- modsem(tpb_uk, data = TPB_UK, "lms", 
               nodes=32, FIM="observed",
               adaptive.quad=TRUE, algorithm ="EMA")
summary(lms3)
#> Regressions:
#>                   Estimate  Std.Error  z.value  P(>|z|)
#>    INT ~         
#>      PBC             1.037      0.036    28.46    0.000
#>      ATT            -0.060      0.030    -2.04    0.041
#>      SN              0.051      0.033     1.55    0.120
#>    BEH ~         
#>      PBC             0.398      0.052     7.63    0.000
#>      INT             0.594      0.049    12.24    0.000
#>      PBC:INT         0.141      0.008    17.66    0.000

# Compared with Mplus
#> Regressions:
#>                    Estimate  Std.Error  z.value  Pr(>|z|)
#>   INT ~   
#>     ATT              -0.053      0.031    -1.71     0.089
#>     SN               -0.065      0.024    -2.71     0.008
#>     PBC               1.090      0.036    30.28     0.000
#>   BEH ~   
#>     PBC               0.405      0.052     7.79     0.000
#>     INT               0.588      0.048    12.25     0.000
#>     INT:PBC           0.141      0.008    17.62     0.000


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
