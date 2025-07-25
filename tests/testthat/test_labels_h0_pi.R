devtools::load_all()

# Test case for bug reported by Prof. Gordon Cheung (06.12.2025, CRAN 1.0.9)

# load data file (Example_C2.csv)
#> Data_C2 <- read.csv(file = 'Example_C2.csv')
#> S <- cov(Data_C2)
#> ltri <- S[lower.tri(S, diag = TRUE)]
#> ovs <- colnames(S)
# Specify the structural model with latent interaction (Example.C2.Model.3) #
ovs <- c("x1", "x2", "x3", "w1", "w2", "w3",
         "m1", "m2", "m3", "y1", "y2", "y3", "y4",
         "y5", "y6", "y7", "y8", "y9", "SAL",
         "sprom1", "sprom2", "sprom3", "perfev", "leader")
ltri <- c(3.91, 3.46, 3.25, -0.86, -0.63, -0.83, 0.84, 1.04, 0.95, -0.11, -0.19, -0.45,
          -0.16, -0.19, -0.34, -0.22, 0.06, -0.16, -0.5, 0.04, -0.47, -0.15, -0.04,
          -16.63, 3.79, 3.28, -0.81, -0.62, -0.77, 0.84, 0.98, 0.88, -0.15, -0.18, -0.32,
          -0.18, -0.19, -0.29, -0.28, 0.12, -0.2, -0.57, 0.09, -0.44, -0.26, -0.01,
          -18.97, 3.82, -0.96, -0.6, -0.85, 1.12, 1.16, 1.13, -0.18, -0.21, -0.37, -0.29,
          -0.27, -0.46, -0.32, 0, -0.16, -0.08, 0.13, -0.42, -0.3, -0.02, -9.41, 3.82,
          3.11, 3.28, -2.06, -2.05, -2, 0.41, 0.44, 0.5, 0.36, 0.34, 0.27, 0.08, 0.22,
          -0.13, 0.61, 0.21, -0.34, -0.09, 0.07, 6.34, 3.19, 2.98, -1.66, -1.63, -1.7,
          0.33, 0.34, 0.34, 0.19, 0.21, 0.15, 0.06, 0.13, -0.01, 0.49, 0.16, -0.48,
          -0.16, 0.04, 7.38, 3.32, -1.86, -1.78, -1.8, 0.25, 0.32, 0.43, 0.24, 0.25, 0.24,
          0.05, 0.16, -0.16, 0.56, 0.17, -0.45, -0.18, 0.05, 8.99, 3.61, 3.3, 2.81,
          -0.46, -0.46, -0.51, -0.61, -0.57, -0.79, -0.37, -0.39, -0.03, -0.4, -0.21, 0.36,
          0.24, -0.03, -3.67, 3.68, 2.71, -0.39, -0.47, -0.51, -0.53, -0.48, -0.76,
          -0.41, -0.43, -0.12, -0.76, -0.19, 0.13, 0.17, 0.01, 1.52, 3.19, -0.55, -0.56,
          -0.68, -0.62, -0.6, -0.71, -0.52, -0.54, -0.26, -0.67, -0.2, 0.21, 0.22,
          -0.09, -1.83, 1.52, 1.34, 1.24, 0.91, 1.09, 1.09, 0.82, 0.76, 0.55, 1.11, 0.69,
          0.6, 0.99, 0.11, 3.3, 1.46, 1.2, 0.88, 1.08, 1.1, 0.79, 0.78, 0.53, 1.04,
          0.66, 0.6, 1.03, 0.1, 0.6, 1.64, 0.98, 1.12, 1.19, 0.87, 0.79, 0.55, 1.03, 0.68,
          0.66, 0.79, 0.14, 4.43, 1.26, 1.12, 1.16, 0.89, 0.84, 0.48, 0.52, 0.47, 0.41,
          0.63, 0.13, 3.13, 1.42, 1.29, 0.9, 0.87, 0.54, 0.71, 0.64, 0.44, 0.85, 0.12,
          9.69, 1.76, 1.13, 1.04, 0.57, 0.77, 0.63, 0.58, 0.84, 0.12, 4.58, 1.58, 0.91,
          0.72, 0.79, 0.33, 0.49, 0.73, 0.14, 3.92, 1.34, 0.68, 0.72, 0.56, 0.33,
          0.65, 0.07, -0.43, 2.2, 1.21, 0.36, 0.05, 0.31, 0.1, 2.56, 12.29, 0.94, 0.43,
          1.02, 0.13, -12.93, 1.28, 0.81, 0.99, 0.06, 7.85, 3.39, 1.83, 0.07, 16.26, 2.51,
          0.1, 8.2, 0.13, 0.64, 5389.66)


set.seed(123)
S <- matrix(NA, nrow=length(ovs), ncol=length(ovs))
S[lower.tri(S, diag = TRUE)] <- ltri
S[upper.tri(S)] <- t(S)[upper.tri(S)]
Data_C2 <- as.data.frame(mvtnorm::rmvnorm(n=400, sigma=S))
colnames(Data_C2) <- ovs

Example.C2.Model.3 <- '
  # Measurement model #
  WFC =~ x1 + x2 + x3
  Burnout =~ m1 + m2 + m3
  Autonomy =~ w1 + w2 + w3
  Vigor =~ y1 + y2 + y3
  Dedi =~ y4 + y5 + y6
  Absorp =~ y7 + y8 + y9
  Promote =~ sprom1 + sprom2 + sprom3
  Salary =~ 1*SAL
  SAL ~~ 0*SAL
  Perform =~ 1* perfev
  perfev ~~ 0*perfev
  Engage =~ Vigor + Dedi + Absorp

  # a1 path #
  Burnout ~ a1*WFC

  # a2 path #
  Engage ~ a2*Burnout

  # a3 paths #
  Promote ~ a3a*Engage
  Perform ~ a3b*Engage
  Salary ~ a3c*Engage

  # z1 path #
  Burnout ~ z1*Autonomy:WFC

  # Main effect of moderator #
  Burnout ~ Autonomy

  # Direct Effects #
  Engage ~ WFC
  Promote ~ WFC
  Perform ~ WFC
  Salary ~ WFC

  # Variance of moderator #
  Autonomy ~~ varAuto*Autonomy

  L1stdz := sqrt(varAuto)
  L2a1 := a1
  L3a2 := a2
  L4a3a := a3a
  L5a3b := a3b
  L6a3c := a3c
  L7z1 := z1
'

est_Model.C3 <- modsem(Example.C2.Model.3, data = Data_C2, method="dblcent", estimator="MLR")
# method = dblcent, rca, uca, ca, pind
summary(est_Model.C3, fit.measure=TRUE, standardized=TRUE, rsq=TRUE)
