library(mvtnorm)

tpbSyntax <- ' 
# Outer Model (Based on Hagger et al., 2007)
  Attitude =~ att1 + att2 + att3 + att4 + att5
  SubjectiveNorm =~ sn1 + sn2
  PerceivedBehaviouralControl =~ pbc1 + pbc2 + pbc3
  Intention =~ int1 + int2 + int3
  Behaviour =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Covariances
  Attitude ~~ cAsn * SubjectiveNorm + cApbc * PerceivedBehaviouralControl
  PerceivedBehaviouralControl ~~ cPbcSn * SubjectiveNorm 
  # Causal Relationsships
  Intention ~ gIa * Attitude + gIsn * SubjectiveNorm + gIpbc * PerceivedBehaviouralControl
  Behaviour ~ Intention + PerceivedBehaviouralControl 
  Behaviour ~ Intention:PerceivedBehaviouralControl  
'

combosParams <- expand.grid(
  repetition = 1,
  N = c(2000),
  # Disturbance coefficients
  sigAtt1 = 0.4,
  sigAtt2 = 0.4,
  sigAtt3 = 0.4,
  sigAtt4 = 0.4,
  sigAtt5 = 0.4,
  sigSn1  = 0.4,
  sigSn2  = 0.4,
  sigPbc1 = 0.4,
  sigPbc2 = 0.4,
  sigPbc3 = 0.4,
  sigInt1 = 0.4,
  sigInt2 = 0.4,
  sigInt3 = 0.4,
  sigB1   = 0.4,
  sigB2   = 0.4,
  # Intercepts Indicators 
  b0Att1 = 1,
  b0Att2 = 1,
  b0Att3 = 1,
  b0Att4 = 1,
  b0Att5 = 1,
  b0Sn1  = 1,
  b0Sn2  = 1,
  b0Pbc1 = 1,
  b0Pbc2 = 1,
  b0Pbc3 = 1,
  b0Int1 = 1,
  b0Int2 = 1,
  b0Int3 = 1,
  b0B1   = 1,
  b0B2   = 1,
  # Loadings
  b1Att1 = 1,
  b1Att2 = 0.9,
  b1Att3 = 0.8,
  b1Att4 = 0.7,
  b1Att5 = 0.9,
  b1Sn1  = 1,
  b1Sn2  = 0.9,
  b1Pbc1 = 1,
  b1Pbc2 = 0.9,
  b1Pbc3 = 0.8,
  b1Int1 = 1,
  b1Int2 = 0.9,
  b1Int3 = 0.8,
  b1B1   = 1,
  b1B2   = 0.9,
  # Covariances,
  covAttPbc = 0.6, 
  covAttSn  = 0.7,
  covPbcSn  = 0.7,
  # Innner Model
    # Intention ~ ...
  gammaIntAtt = 0.2,
  gammaIntSn  = 0.2, 
  gammaIntPbc = 0.2,
    # Behaviour ~ ...
  gammaBehInt = 0.2,
  gammaBehPbc = 0.2,
  gammaBehIntPbc = 0.2,
  # Zeta 
  zetaInt = 0.7,
  zetaBeh = 0.7
)


createDataParams <- function(prow) {
  Sigma <- diag(1, nrow = 3) 
  colnames(Sigma) <- rownames(Sigma) <- c("Att", "Pbc", "Sn") 
  Sigma[1, 2] <- Sigma[2, 1] <- prow$covAttPbc
  Sigma[1, 3] <- Sigma[3, 1] <- prow$covAttSn
  Sigma[2, 3] <- Sigma[3, 2] <- prow$covPbcSn
  distAttPbcSn <- rmvnorm(prow$N, c(0, 0, 0), Sigma) |> 
    as.data.frame()
  # Attitude 
  Attitude <- distAttPbcSn$V1 
  att1 <- prow$b0Att1 + prow$b1Att1 * Attitude + rnorm(prow$N, 0, prow$sigAtt1)
  att2 <- prow$b0Att2 + prow$b1Att2 * Attitude + rnorm(prow$N, 0, prow$sigAtt2)
  att3 <- prow$b0Att3 + prow$b1Att3 * Attitude + rnorm(prow$N, 0, prow$sigAtt3)
  att4 <- prow$b0Att5 + prow$b1Att4 * Attitude + rnorm(prow$N, 0, prow$sigAtt4)
  att5 <- prow$b0Att5 + prow$b1Att5 * Attitude + rnorm(prow$N, 0, prow$sigAtt5)

  # Subjective Norm 
  SubjectiveNorm <- distAttPbcSn$V2
  sn1 <- prow$b0Sn1 + prow$b1Sn1 * 
    SubjectiveNorm + 
    rnorm(prow$N, 0, prow$sigSn1)
  sn2 <- prow$b0Sn2 + prow$b1Sn2 * 
    SubjectiveNorm + 
    rnorm(prow$N, 0, prow$sigSn2)

  # Perceived Behavioural Control 
  PerceivedBehaviouralControl <- distAttPbcSn$V3
  pbc1 <- 
    prow$b0Pbc1 + 
    prow$b1Pbc1 * PerceivedBehaviouralControl + 
    rnorm(prow$N, 0, prow$sigPbc1)
  pbc2 <- 
    prow$b0Pbc2 + 
    prow$b1Pbc2 * PerceivedBehaviouralControl + 
    rnorm(prow$N, 0, prow$sigPbc2)
  pbc3 <- 
    prow$b0Pbc3 + 
    prow$b1Pbc3 * PerceivedBehaviouralControl + 
    rnorm(prow$N, 0, prow$sigPbc3)

  # Intention
  Intention <- 
    prow$gammaIntAtt * Attitude + 
    prow$gammaIntPbc * PerceivedBehaviouralControl +
    prow$gammaIntSn * SubjectiveNorm +
    rnorm(prow$N, 0, prow$zetaInt)
  int1 <- 
    prow$b0Int1 + 
    prow$b1Int1 * Intention + 
    rnorm(prow$N, 0, prow$sigInt1)
  int2 <- prow$b0Int2 + 
    prow$b1Int2 * Intention + 
    rnorm(prow$N, 0, prow$sigInt2)
  int3 <- 
    prow$b0Int3 + 
    prow$b1Int3 * Intention + 
    rnorm(prow$N, 0, prow$sigInt3)
  
  # Behaviour
  Behaviour <- 
    prow$gammaBehInt * Intention +
    prow$gammaBehPbc * PerceivedBehaviouralControl +
    prow$gammaBehIntPbc * Intention * PerceivedBehaviouralControl +
    rnorm(prow$N, 0, prow$zetaBeh)
  b1 <- 
    prow$b0B1 + 
    prow$b1B1 * Behaviour + 
    rnorm(prow$N, 0, prow$sigB1)
  b2 <- 
    prow$b0B2 + 
    prow$b1B2 * Behaviour + 
    rnorm(prow$N, 0, prow$sigB2)
  data.frame(att1, att2, att3, att4, att5, 
             sn1, sn2, 
             pbc1, pbc2, pbc3, 
             int1, int2, int3, 
             b1, b2)
}

theoryOfPlannedBehaviour <- createDataParams(combosParams[1, ])
save(theoryOfPlannedBehaviour, file = "../data/theoryOfPlannedBehaviour.rda")

