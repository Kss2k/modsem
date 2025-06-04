devtools::load_all()

m1 <- '
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z + X:Z
'

lms1 <- modsem(m1, oneInt, method = "lms")


tpb_uk <- "
# Outer Model (Based on Hagger et al., 2007)
 ATT =~ att3 + att2 + att1 + att4
 SN =~ sn4 + sn2 + sn3 + sn1
 PBC =~ pbc2 + pbc1 + pbc3 + pbc4
 INT =~ int2 + int1 + int3 + int4
 BEH =~ beh3 + beh2 + beh1 + beh4

# Inner Model (Based on Steinmetz et al., 2011)
 # Causal Relationsships
 INT ~ ATT + SN + PBC
 BEH ~ INT + PBC
 BEH ~ INT:PBC
"

devtools::load_all()
est <- modsem(tpb_uk, data = TPB_UK, "lms", nodes=64)

library(ggplot2)

quad_f <- model$quad
quad_a <- quad
quad_o <- order(quad_f$n)
quad_a$n <- quad_f$n[quad_o]
quad_a$w <- quad_f$w[quad_o]

m <- quad_f$m

d <- data.frame(
  x = c(quad_f$n, quad_a$n),
  w = c(quad_f$w, quad_a$w),
  p = c(cumsum(quad_f$w), cumsum(quad_a$w)),
  type = c(rep("fixed", m), rep("adaptive", m))
)


ggplot(d, aes(x = x, y = w, color = type)) +
  geom_point() +
  labs(x = "Nodes", y = "Weights", title = "Nodes and Weights for LMS") +
  theme_minimal()


ggplot(d, aes(x = x, y = p, color = type)) +
  geom_point(position = position_dodge(width=0.1)) +
  labs(x = "Nodes", y = "Weights", title = "Nodes and Weights for LMS") +
  theme_minimal()
# Standardized estimates
summary(lms1, standardized = TRUE)


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

lms2 <- modsem(tpb, TPB, method = "lms", nodes = 32)
summary(lms2)
