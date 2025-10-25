devtools::load_all()

m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3

  # Inner model
  Y ~ X + Z + X:Z
'

oneInt2 <- oneInt

set.seed(123)
k <- 200
I <- sample(nrow(oneInt2), k, replace = TRUE)
J <- sample(ncol(oneInt2), k, replace = TRUE)
for (k_i in seq_along(I)) oneInt2[I[k_i], J[k_i]] <- NA

# double centering approach
est <- modsem(m1, oneInt2)
testthat::expect_true(any(is.na(est$data)))

est <- modsem(m1, oneInt2, na.rm=TRUE)
testthat::expect_true(!any(is.na(est$data)))

est <- modsem(m1, oneInt2, na.rm=FALSE)
testthat::expect_true(any(is.na(est$data)))

# Residual Centering Approach
est <- modsem(m1, oneInt2, method = "rca")
testthat::expect_true(any(is.na(est$data)))

est <- modsem(m1, oneInt2, method = "rca", na.rm=TRUE)
testthat::expect_true(!any(is.na(est$data)))

est <- modsem(m1, oneInt2, method = "rca", na.rm=FALSE)
testthat::expect_true(any(is.na(est$data)))

# lms
testthat::expect_warning(
  modsem(m1, oneInt2, method = "lms", missing = "complete", convergence.abs = 1,
         calc.se = FALSE),
  regexp = "Removing.*Consider.*"
)


testthat::expect_message(
  modsem(m1, oneInt2, method = "lms", missing = "impute", convergence.abs = 1,
         calc.se = FALSE),
  regexp = "Imputing.*"
)

# qml
testthat::expect_message(
  modsem(m1, oneInt2, method = "qml", missing = "impute", convergence.rel =1e-1,
         calc.se = FALSE),
  regexp = "Imputing.*"
)


# test multiple imputation
mimp_lms <- modsem_mimpute(m1, oneInt2, method = "lms", m = 5, cov.syntax = "")
mimp_qml <- modsem_mimpute(m1, oneInt2, method = "qml", m = 5)

print(summary(mimp_lms))
summary(mimp_qml)

testthat::expect_false(
  mimp_lms$imputations$fitted[[1]]$logLik ==
    mimp_lms$logLik
)

sd.uncorrected <- sqrt(diag(vcov(mimp_lms$imputations$fitted[[1]])))
sd.corrected   <- sqrt(diag(vcov(mimp_lms)))
testthat::expect_true(all(sd.corrected > sd.uncorrected))

# test fiml
fiml_lms <- modsem(m1, oneInt2, method = "lms", missing = "fiml")
fiml_h0  <- estimate_h0(fiml_lms, calc.se = FALSE)
testthat::expect_true(fiml_h0$iterations == 2L)

testthat::expect_error(modsem(m1, oneInt2, method = "qml", missing = "fiml") ,
                       regexp = "Using FIML with QML is not available.*yet.*")

oneInt3 <- oneInt2
oneInt3$x1 <- NA
oneInt3$y2 <- NA

missing.options <- c("fiml", "ml", "direct", "casewise", "listwise",
                     "impute", "complete")
for (missing in missing.options) {
  testthat::expect_error(modsem(m1, oneInt3, method = "lms", missing = "fiml"),
                         regexp = "Please remove .*x1.*y2")
}
