devtools::load_all()

# --- Helpers ------------------------------------------------------------------

m_fully_comp <- '
  X <~ x1 + x2
  Z <~ z1 + z2
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

m_mixed <- '
  X =~ x1 + x2 + x3
  Z <~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

m_linear_comp <- '
  X <~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

m_latent_only <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

# --- parseLavaan: composite detection -----------------------------------------

spec_full  <- modsem:::parseLavaan(m_fully_comp, colnames(oneInt))
spec_mixed <- modsem:::parseLavaan(m_mixed,      colnames(oneInt))
spec_lin   <- modsem:::parseLavaan(m_linear_comp, colnames(oneInt))
spec_lv    <- modsem:::parseLavaan(m_latent_only, colnames(oneInt))

testthat::test_that("parseLavaan detects composites correctly", {
  testthat::expect_equal(spec_full$composites,  c("X", "Z"))
  testthat::expect_equal(spec_mixed$composites, "Z")
  testthat::expect_equal(spec_lin$composites,   "X")
  testthat::expect_length(spec_lv$composites,   0L)
})

testthat::test_that("isCompositeProd is TRUE only when ALL elements are composites", {
  testthat::expect_true(spec_full$isCompositeProd[["XZ"]])
  testthat::expect_false(spec_mixed$isCompositeProd[["XZ"]])
  testthat::expect_false(spec_lin$isCompositeProd[["XZ"]])
  testthat::expect_false(spec_lv$isCompositeProd[["XZ"]])
})

# --- Generated syntax: operator selection ------------------------------------

syntax_full  <- modsem_pi(m_fully_comp,  oneInt, run = FALSE)$syntax
syntax_mixed <- modsem_pi(m_mixed,       oneInt, run = FALSE)$syntax
syntax_lin   <- modsem_pi(m_linear_comp, oneInt, run = FALSE)$syntax

testthat::test_that("fully-composite product uses <~ in generated syntax", {
  testthat::expect_match(syntax_full, "XZ <~", fixed = TRUE)
  testthat::expect_no_match(syntax_full, "XZ =~", fixed = TRUE)
})

testthat::test_that("mixed (latent x composite) product uses =~ in generated syntax", {
  testthat::expect_match(syntax_mixed, "XZ =~", fixed = TRUE)
  testthat::expect_no_match(syntax_mixed, "XZ <~", fixed = TRUE)
})

testthat::test_that("linear composite (not in product) is passed through with <~", {
  testthat::expect_match(syntax_lin, "X <~", fixed = TRUE)
})

# --- Residual covariance toggle -----------------------------------------------

testthat::test_that("no residual covariance syntax for composite product by default", {
  testthat::expect_no_match(syntax_full, "~~", fixed = TRUE)
})

syntax_full_rcov <- modsem_pi(m_fully_comp, oneInt,
                              composite.int.res.cov = TRUE, run = FALSE)$syntax

testthat::test_that("composite.int.res.cov=TRUE adds residual covariances", {
  testthat::expect_match(syntax_full_rcov, "~~", fixed = TRUE)
})

testthat::test_that("mixed product still gets residual covariances by default", {
  testthat::expect_match(syntax_mixed, "~~", fixed = TRUE)
})

# --- method="ca" error --------------------------------------------------------

m_ca_comp_prod <- '
  X <~ x1 + x2
  Z <~ z1 + z2
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

m_ca_mixed_prod <- '
  X =~ x1 + x2 + x3
  Z <~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

# Composite is only used linearly; the interaction is purely between latent variables
m_ca_linear_only <- '
  ATT <~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC + INT:PBC
'

testthat::test_that('method="ca" errors on fully-composite product', {
  testthat::expect_error(
    modsem_pi(m_ca_comp_prod, oneInt, method = "ca", run = FALSE),
    regexp = 'method="ca".*composite'
  )
})

testthat::test_that('method="ca" errors on mixed (latent+composite) product', {
  testthat::expect_error(
    modsem_pi(m_ca_mixed_prod, oneInt, method = "ca", run = FALSE),
    regexp = 'method="ca".*composite'
  )
})

testthat::test_that('method="ca" allows linear composites with latent products', {
  testthat::expect_no_error(
    modsem_pi(m_ca_linear_only, TPB, method = "ca", run = FALSE)
  )
})

# --- Estimation: fully-composite interaction ----------------------------------

testthat::test_that("fully-composite interaction converges with dblcent", {
  fit <- suppressWarnings(modsem_pi(m_fully_comp, oneInt, method = "dblcent",
                                    optim.gradient = "numerical"))
  testthat::expect_true(lavaan::lavInspect(fit$lavaan, "converged"))
})

testthat::test_that("fully-composite interaction converges with rca", {
  fit <- suppressWarnings(modsem_pi(m_fully_comp, oneInt, method = "rca",
                                    optim.gradient = "numerical"))
  testthat::expect_true(lavaan::lavInspect(fit$lavaan, "converged"))
})

testthat::test_that("fully-composite interaction converges with uca", {
  fit <- suppressWarnings(modsem_pi(m_fully_comp, oneInt, method = "uca",
                                    optim.gradient = "numerical"))
  testthat::expect_true(lavaan::lavInspect(fit$lavaan, "converged"))
})

# --- Estimation: mixed interaction --------------------------------------------

testthat::test_that("mixed (latent x composite) interaction converges with dblcent", {
  fit <- suppressWarnings(modsem_pi(m_mixed, oneInt, method = "dblcent",
                                    optim.gradient = "numerical"))
  testthat::expect_true(lavaan::lavInspect(fit$lavaan, "converged"))
})

# --- Regression: pure latent models still work --------------------------------

testthat::test_that("pure latent models are unaffected by composite changes", {
  for (m in c("dblcent", "rca", "uca", "ca")) {
    fit <- suppressWarnings(modsem_pi(m_latent_only, oneInt, method = m))
    testthat::expect_true(
      lavaan::lavInspect(fit$lavaan, "converged"),
      info = paste("method =", m)
    )
  }
})
