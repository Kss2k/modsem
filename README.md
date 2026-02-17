# `modsem` <img src="man/figures/modsem.png" alt="Logo" align = "right" height="139" class="logo">
<!-- badges: start -->
[![R-CMD-check](https://github.com/kss2k/modsem/actions/workflows/checks.yml/badge.svg)](https://github.com/kss2k/modsem/actions/workflows/checks.yml)
[![Tests](https://github.com/kss2k/modsem/actions/workflows/tests.yml/badge.svg)](https://github.com/kss2k/modsem/actions/workflows/tests.yml)
[![CRAN](https://www.r-pkg.org/badges/version/modsem)](https://cran.r-project.org/package=modsem)
[![PKGDOWN-Build](https://github.com/kss2k/modsem/actions/workflows/pkgdown.yml/badge.svg)](https://github.com/kss2k/modsem/actions/workflows/pkgdown.yml)
[![](https://cranlogs.r-pkg.org/badges/grand-total/modsem)](https://cran.r-project.org/package=modsem)
[![GitHub Clones](https://img.shields.io/badge/dynamic/json?color=success&label=Clone&query=count&url=https://gist.githubusercontent.com/Kss2k/131dd94e938508da4f039ec2a4ecb256/raw/clone.json&logo=github)](https://github.com/MShawon/github-clone-count-badge)

<!-- badges: end -->
`modsem` is an `R`-package for estimating interaction (i.e., moderation) effects between latent variables
in structural equation models (SEMs). See https://www.modsem.org for a tutorial.

# Installation
`modsem` is available on `CRAN` and `GitHub`, and can be installed as follows:

```R
# From CRAN 
install.packages("modsem")

# Latest version from GitHub
install.packages("devtools")
devtools::install_github("kss2k/modsem", build_vignettes = TRUE)
```
**Note**: The package needs to be compiled from source on `macOS` (if installing via `GitHub`) and `Linux`.
If you have issues installing the package on `macOS`, you might need to install the `gfortran` compiler.
A `C++` compiler is also required, but should be installed by default on most systems.
See the [R for macOs](https://cran.r-project.org/bin/macosx/tools/) page for more information.

If you're using `Windows`, consider installing [`OpenBLAS in R for Windows`](https://github.com/david-cortes/R-openblas-in-windows) 
for better perfmance. If you're using a `Linux` distribution, consider installing
the [`ropenblas` package](https://CRAN.R-project.org/package=ropenblas)

# Methods/Approaches

There are a number of approaches for estimating interaction effects in SEM. 
In `modsem()`, the `method = "method"` argument allows you to choose which to use.
Different approaches can be categorized into two groups: 
Product Indicator (PI) and Distribution Analytic (DA) approaches.

## Product Indicator (PI) Approaches:
- `"ca"` = constrained approach (Algina & Moulder, 2001)
    - Note that constraints can become quite complicated for complex models, 
      particularly when there is an interaction including enodgenous variables.
      The method can therefore be quite slow. 
- `"uca"` = unconstrained approach (Marsh, 2004)
- `"rca"` = residual centering approach (Little et al., 2006)
- `"dblcent"` = double centering approach (Marsh., 2013)
  - default 
- `"pind"` = basic product indicator approach (not recommended)

## Distribution Analytic (DA) Approaches
- `"lms"` = The Latent Moderated Structural equations (LMS) approach, see the [vignette](https://modsem.org/articles/lms_qml.html)
- `"qml"` = The Quasi Maximum Likelihood (QML) approach, see the [vignette](https://modsem.org/articles/lms_qml.html)
- `"mplus"` = `Mplus`
  - estimates model through `Mplus`, if it is installed

# Examples 

## Elementary Interaction Model (Kenny & Judd, 1984; Jaccard & Wan, 1995)
```R
library(modsem)

m1 <- '
  # Outer Model
  X =~ x1 + x2 + x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  
  # Inner model
  Y ~ X + Z + X:Z 
'

# Double centering approach
est1_dca <- modsem(m1, oneInt)
summary(est1_dca)

# Constrained approach
est1_ca <- modsem(m1, oneInt, method = "ca")
summary(est1_ca)

# QML approach 
est1_qml <- modsem(m1, oneInt, method = "qml")
summary(est1_qml, standardized = TRUE) 

# LMS approach 
est1_lms <- modsem(m1, oneInt, method = "lms") 
summary(est1_lms)
```

## Theory Of Planned Behavior
```R
tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ PBC:INT
"

# double centering approach
est_tpb_dca <- modsem(tpb, data = TPB, method = "dblcent")
summary(est_tpb_dca)

# Constrained approach using Wrigths path tracing rules for generating
# the appropriate constraints
est_tpb_ca <- modsem(tpb, data = TPB, method = "ca") 
summary(est_tpb_ca)

# LMS approach 
est_tpb_lms <- modsem(tpb, data = TPB, method = "lms")
summary(est_tpb_lms, standardized = TRUE) 

# QML approach 
est_tpb_qml <- modsem(tpb, data = TPB, method = "qml") 
summary(est_tpb_qml, standardized = TRUE)
```
## Interactions between two observed variables
```R
est2 <- modsem('y1 ~ x1 + z1 + x1:z1', data = oneInt, method = "dblcent")
summary(est2)
```

## Interaction between an obsereved and a latent variable 
```R
m3 <- '
  # Outer Model
  X =~ x1 + x2 + x3
  Y =~ y1 + y2 + y3
  
  # Inner model
  Y ~ X + z1 + X:z1
'

est3 <- modsem(m3, oneInt, method = "dblcent", 
               res.cov.method = "none") # res.cov.method = "simple" will lead
                                        # to an unidentifiable model. Instead we
                                        # constrain them to zero
summary(est3)
```
