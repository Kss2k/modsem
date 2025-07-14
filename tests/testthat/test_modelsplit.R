devtools::load_all()

tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + PBC + SN
  BEH ~ INT + PBC 
  BEH ~ INT:PBC  
'

splitm <- function(s) splitParTable(modsemify(s))
s1 <- splitm(tpb)
testthat::expect_true(unique(s1$parTableCov$lhs) == "INT")

testthat::expect_warning(
  est_qml <- modsem(tpb, TPB, method = "qml", auto.split.syntax = FALSE)
)
summary(est_qml)


tpb2 <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  HAPPY =~ h1 + h2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + PBC + SN
  BEH ~ INT + PBC 
  BEH ~ INT:PBC
  HAPPY ~ BEH:INT
'

s2 <- splitm(tpb2)
testthat::expect_true(unique(s2$parTableCov$lhs) == "INT")

tpb3 <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  HAPPY =~ h1 + h2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + PBC + SN
  BEH ~ INT + PBC 
  INT ~ SN:PBC
  HAPPY ~ BEH:INT
'

s3 <- splitm(tpb3)
testthat::expect_true(is.null(s3$parTableCov))
