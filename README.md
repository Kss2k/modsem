# ModSEM
This is a package which allows you to perform interactions between latent variables in CB-SEM. See https://bookdown.org/slupphaugkjell/quartomodsem/ for a tutorial. Note: Examples using `data = oneInt` won't give the same output for you, as the dataset has been replaced in the package 

# To Install 
```
install.packages("devtools")
devtools::install_github("kss2k/modsem")
```

# Methods/Approaches

There are a number of approaches for estimating interaction effects in SEM. In `modsem()`, the `method = "method"` argument allows you to choose which to use.

- "ca" = constrained approach (Algina & Moulder, 2001)
  - not recommended for cases where there is a main effect between variables in the interaction term, unless you know what you are doing.
  - e.g., 
  - use `removeFromParTable = "X ~~ Z"`
  - and `addToParTable = "new formula for covariance, with label Cov_X_Z"`
- "uca" = unconstrained approach (Marsh, 2004)
  - not recommended for cases where there is a main effect between variables in the interaction term, unless you know what you are doing.
  
- "rca" = residual centering approach (Little et al., 2006)
- "dblcent" = double centering approach (Marsh., 2013)
  - default
- "pind" = basic product indicator approach (not recommended)
- "lms" = The latent moderated structural equations approach through the nlsem package
  - do `qml = TRUE` for the quasi maximum likelihood version
  - do `optimize = TRUE` for faster convergence (experimental feature)
  - note: can only be done if you have a single endogenous (dependent) variable. 
- "mplus" 
  - estimates model through Mplus, if it is installed


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

est1 <- modsem(m1, oneInt)
summary(est1)
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

est4 <- modsem(m4, twoInt)
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

est5 <- modsem(m5, tripleInt, standardizeData = TRUE)
summary(est5)
```
