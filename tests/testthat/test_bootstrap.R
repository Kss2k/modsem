devtools::load_all()

m1 <- '
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2

  Y ~ X + Z + X:Z
'

fit_pi <- modsem(m1, oneInt)

modsem_bootstrap(fit_pi, FUN = coef, R = 10)

rsqr_diff <- function(est_h1) {
  est_h0 <- estimate_h0(est_h1, reduced = FALSE)
  r2.h1 <- modsem_inspect(est_h1, what = "r2")
  r2.h0 <- modsem_inspect(est_h0, what = "r2")

  r2.h1[names(r2.h0)] - r2.h0
}

modsem_bootstrap(fit_pi, FUN = rsqr_diff, R = 10)


fit_da <- modsem(m1, oneInt, method = "lms")
summary(fit_da)

modsem_bootstrap(fit_da, FUN = coef, R = 10L)

