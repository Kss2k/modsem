devtools::load_all()
set.seed(123)
oneIntMG <- oneInt
oneIntMG$group <- sample(c("G1", "G2"), nrow(oneInt), replace = TRUE)


m2 <- '
  X =~ x1 + lx2 * x2 + lx3 * x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ c(a1, a2) * X + Z + X:Z

  a1==a2
'

est <- modsem(m2, oneIntMG, method = "qml", group = "group", robust.se = TRUE)
summary(est)
#> Error in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, 
#>  : 
#>   invalid value 0 for 'digits' argument
#> In addition: Warning message:
#> Comparative fit to H0 will not be calculated. 


HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9
              visual ~ c(v_t1, v_t2) * textual + speed
              interaction := v_t2 - v_t1'

fit_lavaan <- lavaan::sem(HS.model,
                          data = lavaan::HolzingerSwineford1939,
                          group = "school")

# using modsem
fit_modsem <- modsem(HS.model,
                     data = lavaan::HolzingerSwineford1939,
                     group = "school", method = "lms", max.step = 1)


lavaanEst <- lavaan::parameterEstimates(fit_lavaan)
lavaanEst[is.na(lavaanEst)] <- -999
modsemEst <- lavaan::parameterEstimates(fit_modsem$lavaan)
modsemEst[is.na(modsemEst)] <- -999
testthat::expect_equal(lavaanEst, modsemEst)

set.seed(123)
oneIntMG <- oneInt
oneIntMG$group <- sample(c("G1", "G2"), nrow(oneInt), replace = TRUE)

m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

est.mg <- modsem(m1, oneIntMG, method = "lms", group = "group", robust.se = TRUE)
est.g1 <- modsem(m1, subset(oneIntMG, group == "G1"), method = "qml", robust.se = TRUE)
est.g2 <- modsem(m1, subset(oneIntMG, group == "G2"), method = "qml", robust.se = TRUE)

#> Iter=161 Mode=EM LogLik=-17477.99 ΔLL=3.2e-05 relΔLL=1.8e-09              
#> 
#> Warning: Model estimation failed, returning starting values!
#> Message: replacement has length zero
#> 
#> Error: modsem [lms]: Model estimation failed!
#> Message: replacement has length zero
standardized_estimates(est, correction = TRUE)
standardized_estimates(est, correction = TRUE, std.errors = "delta")
summary(est)

m2 <- '
  X =~ x1 + lx2 * x2 + lx3 * x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ c(a1, a2) * X + Z + X:Z

  a1==a2
'

est <- modsem(m2, oneIntMG, method = "lms", group = "group")
#> Error: Unknown labels in constraints: a1
standardized_estimates(est, correction = TRUE)
standardized_estimates(est, correction = TRUE, std.errors = "delta")
summary(est)
