# ModSEM
latent variable interaction in CB-SEM. See https://bookdown.org/slupphaugkjell/an_introduction_to_latent_interaction_using_modsem/ for a tutorial.

# To Install 
```
install.packages("devtools")
devtools::install_github("kss2k/modsem")
```

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
