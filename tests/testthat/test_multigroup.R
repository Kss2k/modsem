devtools::load_all()
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
                     group = "school")


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

est <- modsem(m1, oneIntMG, group = "group")
standardized_estimates(est, correction = TRUE)
standardized_estimates(est, correction = TRUE, std.errors = "delta")
summary(est)
