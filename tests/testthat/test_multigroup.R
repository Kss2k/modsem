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

est.mg <- modsem(m1, oneIntMG, method = "lms", group = "group", robust.se = TRUE)
est.mg.rcs <- modsem(m1, oneIntMG, method = "lms", group = "group", robust.se = TRUE, rcs = TRUE)
standardized_estimates(est.mg, correction = TRUE)
standardized_estimates(est.mg, correction = TRUE, std.errors = "delta")
summary(standardize_model(est.mg))

est.dca.mg <- modsem(m1, oneIntMG, method = "dblcent", group = "group")
est.dca.rcs.mg <- modsem(m1, oneIntMG, method = "dblcent", group = "group", rcs = TRUE)

m2 <- '
  X =~ x1 + lx2 * x2 + lx3 * x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ c(a1, a2) * X + Z + X:Z

  a1==a2
'

est <- modsem(m2, oneIntMG, method = "lms", group = "group")
standardized_estimates(est, correction = TRUE)
standardized_estimates(est, correction = TRUE, std.errors = "delta")
summary(est, standardized = TRUE, center = TRUE)
plot_jn(x = "X", z = "Z", y = "Y", model = est)
plot_jn(x = "X", z = "Z", y = "Y", model = est, plot.jn.points = FALSE)
plot_interaction(x = "X", z = "Z", y = "Y", model = est, vals_z = c(1, 0))
testthat::expect_warning(
  plot_surface(x = "X", z = "Z", y = "Y", model = est),
  regexp = "Plotting of surface.*"
)

oneIntMG2 <- oneIntMG
oneIntMG2$group[c(2, 388, 291, 1502)] <- NA

testthat::expect_error(
  modsem(m1, oneIntMG2, method = "lms", group = "group"),
  regexp = ".*group.* cannot contain missing values.*"
)
