# ModSEM <img src="man/figures/ModSEM.png" alt="Logo" align = "right" height="139" class="logo">
This is a package which allows you to perform interactions between latent variables (i.e., moderation) in CB-SEM. See https://kss2k.github.io/intro_modsem/ for a tutorial.

# To Install 
```
# From CRAN 
install.packages("modsem")

# Latest version from Github
install.packages("devtools")
devtools::install_github("kss2k/modsem", build_vignettes = TRUE)
```

# Methods/Approaches

There are a number of approaches for estimating interaction effects in SEM. In `modsem()`, the `method = "method"` argument allows you to choose which to use.

- `"ca"` = constrained approach (Algina & Moulder, 2001)
    - Note that constraints can become quite complicated for complex models, 
      particularly when there is an interaction including enodgenous variables.
      The method can therefore be quite slow. 
- `"uca"` = unconstrained approach (Marsh, 2004)
- `"rca"` = residual centering approach (Little et al., 2006)
- `"dblcent"` = double centering approach (Marsh., 2013)
  - default 
- `"pind"` = basic product indicator approach (not recommended)
- `"lms"` = The latent moderated structural equations approach
- `"qml"` = The Quasi Maximum Likelihood estimation of latent moderated structural equations. 
- `"mplus"` 
  - estimates model through Mplus, if it is installed

# New Features (01.06.2024-04-06-2024)
- New function for plotting interaction effects (`plot_interaction()`)
- Interaction effects between two endogenous variables are now possible with the LMS 
  and QML approach, using the 'cov_syntax' argument, see the [vignette](file:///home/kjell/Dropbox/code/R/repos/modsem/docs/articles/interaction_two_etas.html) for more information.
  (only available on GitHub, not on CRAN)

# Examples 

## One interaction
```
library(modsem)
m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  
  # Inner model
  Y ~ X + Z + X:Z 
'

# Double centering approach
est1Dblcent <- modsem(m1, oneInt)
summary(est1Dblcent)

# Constrained approach
est1Ca <- modsem(m1, oneInt, method = "ca")
summary(est1Ca)

# QML approach 
est1Qml <- modsem(m1, oneInt, method = "qml")
summary(est1Qml) 

# LMS approach 
est1Lms <- modsem(m1, oneInt, method = "lms") 
summary(est1Lms)
```

## Theory Of Planned Behavior
```
tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Causal Relationsships
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ PBC:INT
"

# double centering approach
estTpbDblCent <- modsem(tpb, data = TPB, method = "dblcent")
summary(estTpbDblCent)

# Constrained approach using Wrigths path tracing rules for generating
# the appropriate constraints
estTpbCa <- modsem(tpb, data = TPB, method = "ca") 
summary(estTpbCa)

# LMS approach 
estTpbLms <- modsem(tpb, data = TPB, method = "lms")
summary(estTpbLms)
```
## Interactions between two observed variables
```
est2 <- modsem('y1 ~ x1 + z1 + x1:z1', data = oneInt, method = "pind")
summary(est2)

## Interaction between an obsereved and a latent variable 
m3 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  
  # Inner model
  Y ~ X + z1 + X:z1 
'

est3 <- modsem(m3, oneInt, method = "pind")
summary(est3)
```

## Multiple interaction terms
```
m4 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  G =~ g1 + g2 + g3
  H =~ h1 + h2 + h3
  
  # Inner model
  Y ~ X + Z + G + H + X:Z + G:H
'

# Using unconstrained approach
est4 <- modsem(m4, twoInt, method = "uca")
summary(est4)
```

## Interactionterms with more than two variables
```
m5 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  G =~ g1 + g2 + g3
  
  # Inner model
  Y ~ X + Z + G + X:Z:G
'

# Residual centering approach
est5 <- modsem(m5, tripleInt, standardizeData = TRUE, method = "rca")
summary(est5)
```
