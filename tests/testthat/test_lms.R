devtools::load_all()
m1 <- "
# Outer Model
  X =~ x1
  Z =~ z1 
  x1 ~~ 0.1 * x1
  Y =~ y1

# Inner model
  Y ~ a * X + a * Z
  Y ~~ Y
  Y ~ b * X:Z + 0.05 * X:X
  b == a * 1.2
"
# funnily enough, the starting parameters from the double centering approach
# give better loglikelihoods than the ones arrived at by the EM algorithm
# i.e., the loglikelihood decreases from the starting parameters
startTime1 <- Sys.time()
est1 <- modsem(m1, oneInt, 
  method = "lms",
  optimize = TRUE, verbose = TRUE,
  convergence = 1e-2,
)
duration1 <- Sys.time() - startTime1
plot_interaction("X", "Z", "Y", "X:Z", -3:3, c(-0.5, 0.5), est1)
print(summary(est1))
# I have no clue why, but changing the ordering of how the interaction terms 
# are specified, ends up changing the number of iterations (and results ever 
# so slightly) -- even though the matrices are exactly the same. This can be 
# seen through the fact that the starting loglikelihoods are the same (if optimized) 
# indicating that the matrices are the same (i.e,. produce the same results, when
# given the same values). 
# Solution: slightly different results from lavaan, giving slightly different 
# starting parameters
tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Causal Relationsships
  INT ~ a * ATT + b * SN + b * PBC
  BEH ~ 0.2 * INT + a * PBC
  BEH ~ PBC:INT
"

covModel <- '
PBC ~ a * ATT + a * SN
'

startTime2 <- Sys.time()
est2 <- modsem(tpb, TPB, 
  method = "lms", optimize = TRUE, verbose = TRUE, 
  convergence = 1, sampleGrad = NULL, cov_syntax = covModel,
  nodes = 16
  # closer to mplus when tweaking the number of nodes and convergence criterion
  # nodes = 100, convergence = 1e-7 is very very close to mplus
)
duration2 <- Sys.time() - startTime2
plot_interaction(x = "INT", z = "PBC", y = "BEH", xz = "PBC:INT", vals_z = c(-0.5, 0.5), model = est2)
