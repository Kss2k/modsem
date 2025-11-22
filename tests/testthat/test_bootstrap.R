devtools::load_all()

m1 <- '
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2

  Y ~ X + Z + X:Z
'

fit_pi <- modsem(m1, oneInt)

bootstrap_modsem(fit_pi, FUN = coef, R = 10)

rsqr_diff <- function(est_h1) {
  est_h0 <- estimate_h0(est_h1, reduced = FALSE)
  r2.h1 <- modsem_inspect(est_h1, what = "r2")
  r2.h0 <- modsem_inspect(est_h0, what = "r2")

  r2.h1[names(r2.h0)] - r2.h0
}

bootstrap_modsem(fit_pi, FUN = rsqr_diff, R = 10)


fit_da <- modsem(m1, oneInt, method = "lms")
summary(fit_da)

bootstrap_modsem(fit_da, FUN = coef, R = 10L, type = "parametric")
bootstrap_modsem(fit_da, FUN = coef, R = 10L, type = "nonparametric")

tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC + INT:PBC
"

boot <- bootstrap_modsem(model = modsem,
                         model.syntax = tpb, data = TPB,
                         method = "dblcent", rcs = TRUE,
                         rcs.scale.corrected = TRUE,
                         rcs.mc.reps = 10000,
                         R = 10L,
                         FUN = "coef")
coef <- apply(boot, MARGIN = 2, FUN = mean, na.rm = TRUE)
se   <- apply(boot, MARGIN = 2, FUN = sd, na.rm = TRUE)

cat("Parameter Estimates:\n")
print(coef)

cat("Standard Errors: \n")
print(se)


for (method in c("lms", "dblcent")) {
  fit_rcs <- modsem(m1, oneInt, rcs = TRUE, method = method)

  testthat::expect_warning(
    bootstrap_modsem(fit_rcs, R = 10L),
    regexp = "bootstrapping a model .*rcs=TRUE.*naive"
  )
}
