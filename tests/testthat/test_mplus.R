devtools::load_all()
m1 <- '
# Outer Model
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2

# Inner model
  Y ~ b1 * X + b2 * Z + b3 * X:Z

  ccoef := (b1 + b2 * b3) / 2
  ccoef2 := b1 * b2

  b1 == b2
'

run <- tryCatch({
    MplusAutomation::detectMplus()
    TRUE
  },
  error = function(e) FALSE
)
if (run) {
  mplus <- modsem(m1, oneInt, method = "mplus", estimator = "MLR")
  standardized_estimates(mplus, type = "modsem")
  print(summary(mplus))
  print(vcov(mplus))

  # Standardized estimates are added as a separate column
  sum.std <- summary(mplus, standardized = TRUE)
  std     <- standardized_estimates(mplus)
  testthat::expect_equal(sum.std$format$extra.cols, "Std.all")
  testthat::expect_equal(sum.std$parTable$est, mplus$parTable$est)
  testthat::expect_equal(
    sum.std$parTable$Std.all[sum.std$parTable$rhs == "X:Z"],
    std$est[std$rhs == "X:Z" & std$op == "~"]
  )
  testthat::expect_output(print(sum.std), "Std.all")
  testthat::expect_output(print(summary(mplus, standardized = TRUE, type = "std")), "Std.lv")
  plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = mplus)
}


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

if (run) {
  mplus_tpb <- modsem(tpb, data = TPB, method = "mplus", rcs = TRUE)
  print(standardized_estimates(mplus_tpb, type = "modsem"))
  print(standardized_estimates(mplus_tpb, type = "stdyx"))
}

tpb <- "
# Outer Model (Based on Hagger et al., 2007)
  SUBJECTIVE_NORMS =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + PBC + INT:PBC:SUBJECTIVE_NORMS + INT:PBC + SUBJECTIVE_NORMS
"

if (run) {
  mplus_tpb_3way <- modsem(tpb, data = TPB[1:250, ], method = "mplus", rcs = TRUE,
                           integration = 8)
  standardized_estimates(mplus_tpb_3way)
}


# Test robust std.errors
mod <- '
# X1-3 are Level 1 variables
X1 =~ x1
X2 =~ x2
X3 =~ x3

# W1-2 are Level 2 variables
W1 =~ w1
W2 =~ w2

fw =~ y1 + y2 + y3
fw ~ X1 + X2 + X3 + W1 + W2
'

if (run) {
  # Standard errors corrected for clustering
  fit.rc <- modsem(mod, lavaan::Demo.twolevel, rcs = TRUE, method = "mplus",
                   cluster = "cluster", estimator = "MLR")
  summary(fit.rc)
}


# ──────────────────────────────────────────────────────────────────────────────
# Model specification (doesn't require Mplus)
# ──────────────────────────────────────────────────────────────────────────────
# Mplus names are case insensitive, and limited to 8 characters
names <- getMplusNames(c("x1", "X1", "a_very_long_name", "a_very_long_name2", "with"))
testthat::expect_equal(unname(names), c("x1", "X11", "a_very_l", "a_very_1", "with1"))
testthat::expect_true(!anyDuplicated(toupper(names)))

ml.syntax <- '
level: 1
  X =~ x1 + x2 + x3
  Y =~ y1 + y2 + y3
  Y ~ rv("s1")*X + x4 + x5 + rv("s2")*X:x4 + x4:x5
level: 2
  Yb =~ y1 + y2 + y3
  Yb + s1 + s2 ~ w1 + w2 + w1:w2
'

ml.data <- as.data.frame(matrix(rnorm(100 * 11), nrow = 100,
                                dimnames = list(NULL, c("x1", "x2", "x3", "x4", "x5",
                                                        "y1", "y2", "y3", "w1", "w2",
                                                        "cluster"))))

spec <- specifyMplusModel(modsemify(ml.syntax), data = ml.data, cluster = "cluster")
testthat::expect_true(spec$multilevel)
testthat::expect_true(spec$has.xwith)
testthat::expect_equal(spec$slopes$name, c("s1", "s2"))
testthat::expect_true(grepl("WITHIN = x1 x2 x3 x4 x5 x4x5;", spec$VARIABLE, fixed = TRUE))
testthat::expect_true(grepl("BETWEEN = w1 w2 w1w2;", spec$VARIABLE, fixed = TRUE))
testthat::expect_true(grepl("Xx4 | X XWITH x4;", spec$MODEL, fixed = TRUE))
testthat::expect_true(grepl("s1 | Y ON X;", spec$MODEL, fixed = TRUE))
testthat::expect_true(grepl("s2 | Y ON Xx4;", spec$MODEL, fixed = TRUE))
testthat::expect_true(grepl("Y ON x4x5;", spec$MODEL, fixed = TRUE))
testthat::expect_true(grepl("Yb ON w1w2;", spec$MODEL, fixed = TRUE))
testthat::expect_true(grepl("x4x5 = x4*x5;", spec$DEFINE, fixed = TRUE))
testthat::expect_true(grepl("w1w2 = w1*w2;", spec$DEFINE, fixed = TRUE))
testthat::expect_false(any(c("x4x5", "w1w2") %in% colnames(spec$data)))
testthat::expect_true(grepl("USEVARIABLES = .* x4x5 w1w2;", spec$VARIABLE))

analysis <- specifyMplusAnalysis(spec)
testthat::expect_equal(analysis$type, "twolevel random")
testthat::expect_equal(analysis$algorithm, "integration")
testthat::expect_false(grepl("STANDARDIZED", specifyMplusOutput(spec, analysis = analysis)))

# Bayesian estimation (no numerical integration), and additional ANALYSIS options
analysis.b <- specifyMplusAnalysis(spec, estimator = "bayes",
                                   options = list(biterations = "(5000)", processors = 4))
testthat::expect_null(analysis.b$algorithm)
testthat::expect_false(grepl("ALGORITHM|INTEGRATION", analysis.b$ANALYSIS))
testthat::expect_true(grepl("BITERATIONS = (5000);", analysis.b$ANALYSIS, fixed = TRUE))
testthat::expect_true(grepl("PROCESSORS = 4;", analysis.b$ANALYSIS, fixed = TRUE))
testthat::expect_false(grepl("PROCESSORS = 2;", analysis.b$ANALYSIS, fixed = TRUE))
testthat::expect_true(grepl("STANDARDIZED", specifyMplusOutput(spec, analysis = analysis.b)))
testthat::expect_error(specifyMplusAnalysis(spec, options = list(5000)), "named list")

# Higher order interactions are specified recursively
spec3 <- specifyMplusModel(modsemify("
  X =~ x1 + x2
  Z =~ x3 + x4
  W =~ y1 + y2
  y3 ~ X + Z + W + X:Z + X:Z:W
"), data = ml.data)
testthat::expect_false(spec3$multilevel)
testthat::expect_true(grepl("XZ | X XWITH Z;\nXZW | XZ XWITH W;", spec3$MODEL, fixed = TRUE))
testthat::expect_equal(specifyMplusAnalysis(spec3)$type, "random")

# Errors
testthat::expect_error(
  specifyMplusModel(modsemify(ml.syntax), data = ml.data),
  regexp = "cluster.*must be specified"
)
testthat::expect_error(
  specifyMplusModel(modsemify("level: 1\n Y =~ y1 + y2\nlevel: 2\n Yb =~ y1 + y2\n Yb ~ rv('s')*w1"),
                    data = ml.data, cluster = "cluster"),
  regexp = "Random slopes.*within level"
)
testthat::expect_error(
  specifyMplusModel(modsemify("Y =~ y1 + y2\n Y ~ rv('s')*x1"), data = ml.data),
  regexp = "Random slopes.*only available in multilevel"
)
testthat::expect_error(
  specifyMplusModel(modsemify("level: 1\n Y =~ y1 + y2\n Y ~ x9\nlevel: 2\n Yb =~ y1 + y2"),
                    data = ml.data, cluster = "cluster"),
  regexp = "Missing variables in data: x9"
)
testthat::expect_error(
  specifyMplusModel(modsemify("level: 1\n Y =~ y1 + y2\nlevel: 3\n Yb =~ y1 + y2"),
                    data = ml.data, cluster = "cluster"),
  regexp = "Unrecognized level"
)


# ──────────────────────────────────────────────────────────────────────────────
# Multilevel models
# ──────────────────────────────────────────────────────────────────────────────
compareLavaan <- function(fit, fit.lav, tol = 2e-3) {
  pt.mplus <- parameter_estimates(fit)
  pt.lav   <- lavaan::parameterEstimates(fit.lav)

  key   <- \(pt) paste(pt$lhs, pt$op, pt$rhs, pt$level)
  match <- match(key(pt.mplus), key(pt.lav))

  testthat::expect_false(anyNA(match))
  testthat::expect_lt(max(abs(pt.mplus$est - pt.lav$est[match])), tol)
}


two.level <- '
  level: 1
    fw =~ y1 + y2 + y3
    fw ~ x1 + x2 + x3
  level: 2
    fb =~ y1 + y2 + y3
    fb ~ w1 + w2
'

random.slopes <- '
  level: 1
    fw =~ y1 + y2 + y3
    fw ~ rv("s1")*x1 + rv("s2")*x2 + x3
  level: 2
    fb =~ y1 + y2 + y3
    fb ~ w1 + w2
    s1 + s2 ~ w1 + w2
'

if (run) {
  fit.ml <- modsem_mplus(two.level, data = lavaan::Demo.twolevel,
                         cluster = "cluster", cleanup = TRUE)
  compareLavaan(fit.ml, lavaan::sem(two.level, data = lavaan::Demo.twolevel,
                                    cluster = "cluster"))
  testthat::expect_equal(NROW(vcov(fit.ml)), length(coef(fit.ml)))
  testthat::expect_true("fb=~y2.l2" %in% names(coef(fit.ml)))
  testthat::expect_no_error(standardized_estimates(fit.ml))
  testthat::expect_no_error(capture.output(summary(fit.ml)))
  testthat::expect_error(standardized_estimates(fit.ml, type = "modsem"),
                         regexp = "not available for\\s+multilevel")

  fit.rs <- modsem_mplus(random.slopes, data = lavaan::Demo.twolevel,
                         cluster = "cluster", cleanup = TRUE)
  fit.rs.lav <- tryCatch( # random slopes are only available in newer versions of lavaan
    lavaan::sem(random.slopes, data = lavaan::Demo.twolevel, cluster = "cluster"),
    error = \(e) NULL
  )

  if (!is.null(fit.rs.lav)) compareLavaan(fit.rs, fit.rs.lav)
  testthat::expect_equal(NROW(vcov(fit.rs)), length(coef(fit.rs)))
  testthat::expect_true(all(c("s1~1.l2", "s2~~s2.l2") %in% names(coef(fit.rs))))
  testthat::expect_no_error(capture.output(summary(fit.rs)))
}


# Latent interaction at the within level, with a random intercept and a random slope
set.seed(2934)
simMultilevelInteraction <- function(K = 100, nk = 20) {
  N  <- K * nk
  cl <- rep(seq_len(K), each = nk)
  w  <- rnorm(K)
  u  <- rnorm(K, sd = 0.5)                  # random intercept of Y
  s  <- 0.3 + 0.2 * w + rnorm(K, sd = 0.2)  # random slope of Y on X

  X <- rnorm(N)
  Z <- 0.3 * X + rnorm(N, sd = sqrt(1 - 0.09))
  Y <- u[cl] + s[cl] * X + 0.3 * Z + 0.25 * X * Z + 0.3 * w[cl] + rnorm(N, sd = 0.6)

  ind <- \(f, l) sapply(l, \(li) li * f + rnorm(N, sd = 0.5))
  data <- data.frame(ind(X, c(1, .8, .9)), ind(Z, c(1, .7, .8)), ind(Y, c(1, .9, .8)), w[cl], cl)
  names(data) <- c("x1", "x2", "x3", "z1", "z2", "z3", "y1", "y2", "y3", "w", "cluster")
  data
}

ml.int.data <- simMultilevelInteraction()

ml.int <- '
  level: 1
    X =~ x1 + x2 + x3
    Z =~ z1 + z2 + z3
    Y =~ y1 + y2 + y3
    Y ~ rv("s")*X + Z + b*X:Z
  level: 2
    Yb =~ y1 + y2 + y3
    Yb + s ~ w
    b2 := 2 * b
'

if (run) {
  fit.int <- modsem_mplus(ml.int, data = ml.int.data, cluster = "cluster",
                          integration = 8, cleanup = TRUE)
  pt <- parameter_estimates(fit.int)
  est <- \(lhs, op, rhs, level) pt$est[pt$lhs == lhs & pt$op == op & pt$rhs == rhs & pt$level == level]

  testthat::expect_equal(est("Y", "~", "X:Z", 1), 0.25, tolerance = 0.1)
  testthat::expect_equal(est("Y", "~", "Z", 1), 0.3, tolerance = 0.1)
  testthat::expect_equal(est("s", "~1", "", 2), 0.3, tolerance = 0.1)
  testthat::expect_equal(est("s", "~", "w", 2), 0.2, tolerance = 0.1)
  testthat::expect_equal(est("b2", ":=", "2 * b", 0), 2 * est("Y", "~", "X:Z", 1))

  testthat::expect_equal(NROW(vcov(fit.int)), length(coef(fit.int)))
  testthat::expect_false(anyNA(vcov(fit.int)))
  testthat::expect_no_error(capture.output(summary(fit.int)))
  testthat::expect_no_error(var_interactions(fit.int))
}


# Observed and between level interactions (computed as product terms using DEFINE)
ml.obs.int <- '
  level: 1
    fw =~ y1 + y2 + y3
    fw ~ x1 + x2 + x3 + x2:x3
  level: 2
    fb =~ y1 + y2 + y3
    fb ~ w1 + w2 + w1:w2
'

if (run) {
  fit.obs <- modsem_mplus(ml.obs.int, data = lavaan::Demo.twolevel,
                          cluster = "cluster", cleanup = TRUE)
  pt <- parameter_estimates(fit.obs)
  testthat::expect_true(all(c("x2:x3", "w1:w2") %in% pt$rhs[pt$op == "~"]))
  testthat::expect_equal(pt$level[pt$rhs == "x2:x3"], 1L)
  testthat::expect_equal(pt$level[pt$rhs == "w1:w2"], 2L)
}


# Bayesian estimation
if (run) {
  m.bayes <- '
    X =~ x1 + x2 + x3
    Z =~ z1 + z2 + z3
    Y =~ y1 + y2 + y3
    Y ~ X + Z + X:Z
  '

  fit.bayes <- suppressWarnings(
    modsem_mplus(m.bayes, data = oneInt, estimator = "bayes", cleanup = TRUE,
                 analysis.options = list(bseed = 123, biterations = "(5000)"))
  )

  pt <- parameter_estimates(fit.bayes)
  b  <- pt[pt$op == "~", ]
  testthat::expect_equal(b$est, c(0.67, 0.57, 0.72), tolerance = 0.05)
  testthat::expect_true(all(is.na(b$z.value)))
  testthat::expect_true(all(b$ci.lower < b$est & b$est < b$ci.upper))
  testthat::expect_true(all(is.na(pt$std.error[pt$op == "=~" & pt$rhs %in% c("x1", "z1", "y1")])))
  testthat::expect_equal(sort(unname(sqrt(diag(vcov(fit.bayes))))),
                         sort(pt$std.error[!is.na(pt$std.error)]), tolerance = 1e-2)
  testthat::expect_no_error(capture.output(summary(fit.bayes, standardized = TRUE)))
}
