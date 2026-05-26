devtools::load_all()
set.seed(123)

m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + a * X:Z

X ~~ varX * X # check that no warning is thrown
              # (one should be thrown for LMS but not QML)
'

est1 <- modsem(m1, data = oneInt, convergence.rel = 1e-2, method = "qml",
               robust.se = TRUE)
print(summary(est1, scientific = TRUE))
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = est1)
plot_surface(x = "X", z = "Z", y = "Y", model = est1)

testthat::test_that("QML honors finite-difference gradient fallback", {
  model <- est1$start.model
  if (is.null(model)) model <- est1$model

  model$params$gradientStruct$useFDGradient <- TRUE
  model$params$gradientStruct$locations <- NULL
  model$params$gradientStruct$Jacobian <- NULL

  grad <- gradientLogLikQml(theta = model$theta, model = model, sum = TRUE)
  scores <- gradientLogLikQml(theta = model$theta, model = model, sum = FALSE)

  testthat::expect_length(grad, length(model$theta))
  testthat::expect_equal(ncol(scores), length(model$theta))
})

testthat::test_that("QML FD mapping includes composite variance blocks", {
  model <- est1$start.model
  if (is.null(model)) model <- est1$model

  submodel <- fillModel(theta = model$theta, model = model,
                        method = "qml")$models[[1L]]
  submodel$matrices$W <- matrix(1, 1, 1)
  submodel$matrices$T <- matrix(1, 1, 1)

  grad <- gradLogLikQmlFDCpp(
    submodel = submodel,
    block = c(14L, 15L),
    row = c(0L, 0L),
    col = c(0L, 0L),
    symmetric = c(0L, 1L)
  )

  testthat::expect_true(all(is.finite(grad)))
})


tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  # Covariances
  ATT ~~ SN + PBC
  PBC ~~ SN
# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ b * INT + PBC
  INT ~ ATT + SN + a * PBC
  BEH ~ PBC:INT

  # INT ~~ BEH
  VERY_LONG_LABEL_THAT_SHOULD_BE_SHORTENED := a * b
'

testthat::expect_no_condition(
est2 <- modsem(tpb, data = TPB, method = "qml",
               robust.se = TRUE,
               standardize = TRUE, convergence.rel = 1e-2)
)

print(summary(est2, H0 = FALSE))
expect_warning(plot_jn(x = "INT", z = "PBC", y = "BEH", model = est2,
                       min_z = -1.5, max_z = -0.5),
               regex = "Truncating.*right.*left.*")
expect_warning(plot_jn(x = "INT", z = "PBC", y = "BEH", model = est2,
                       min_z = -1.5, max_z = 2.5) ,
               regex = "Truncating SD-range on the left!")
expect_warning(plot_jn(x = "INT", z = "PBC", y = "BEH", model = est2,
                       min_z = -2, max_z = 1) ,
               regex = "Truncating SD-range on the right!")

# check that standardizing twice does not change the estimates
deterministicCols <- c("lhs", "op", "rhs", "est") # std.error based cols are subject to rng
testthat::expect_equal(standardized_estimates(est2)[deterministicCols],
                       parameter_estimates(est2)[deterministicCols], tol = 1e-3)

# check correct standardization
calcCovParTable("BEH", "BEH", parameter_estimates(est2))[[1]] |>
  testthat::expect_equal(1)

# check that estimates are standardized the same way
# in vcov, coef and parameter_estimates
pt <- parameter_estimates(est2)
vcov_sd <- sqrt(vcov(est2)["BEH~PBC:INT", "BEH~PBC:INT"])
cond <- pt$lhs == "BEH" & pt$op == "~" & pt$rhs == "PBC:INT"
pt_sd <- pt[cond, "std.error"]
pt_est <- pt[cond, "est"]
coef_est <- coef(est2)[["BEH~PBC:INT"]]

expect_equal(vcov_sd, pt_sd)
expect_equal(pt_est, coef_est)

modsem_inspect(est2)
coefficients(est2)

# Observed Variables
m3 <- '
X =~ x1
Y =~ y1
Z =~ z1
Y ~ X + Z + X:Z
'

est3 <- modsem(m3, data = oneInt, convergence.rel = 1e-2, method = "qml",
               robust.se = TRUE)
print(summary(est3, scientific = TRUE))
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = est3,
                 standardized = TRUE)

tpb2 <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 #+ int2 + int3
  BEH =~ b1  + b2

# Inner Model (Based on Steinmetz et al., 2011)
  BEH ~ INT + PBC
  INT ~ ATT + SN + PBC
  BEH ~ PBC:INT
'

est4 <- modsem(tpb2, data = TPB, method = "qml",
               robust.se = TRUE,
               standardize = TRUE, convergence.rel = 1e-2)
print(summary(est4, H0 = FALSE))

testthat::expect_equal(standardized_estimates(est4)[deterministicCols],
                       parameter_estimates(est4)[deterministicCols])

calcCovParTable("BEH", "BEH", parameter_estimates(est4))[[1]] |>
  testthat::expect_equal(1)

vcov(est4)
modsem_inspect(est4)
coef(est4)
coefficients(est4)

testthat::test_that("simple QML analytical gradient matches C++ FD for supported blocks", {
  syntax <- '
  X =~ x1 + x2 + x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  Y ~ X + Z + X:Z
  '

  group.info <- parseModelArgumentsByGroupDA(
    model.syntax = syntax,
    cov.syntax = NULL,
    method = "qml",
    data = oneInt,
    group = NULL,
    auto.split.syntax = FALSE
  )

  model <- specifyModelDA(
    group.info = group.info,
    method = "qml",
    m = 16,
    mean.observed = TRUE,
    double = FALSE,
    quad.range = Inf,
    adaptive.quad = FALSE,
    adaptive.frequency = 3,
    missing = "complete",
    orthogonal.x = FALSE,
    orthogonal.y = FALSE,
    auto.fix.first = TRUE,
    auto.fix.single = TRUE,
    fix.composite.var = TRUE,
    cluster = NULL,
    sampling.weights = NULL
  )

  model.filled <- fillModel(model = model, theta = model$theta, method = "qml")
  locations <- subset(
    model$params$gradientStruct$locations,
    group == 1 & block %in% c(0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 16)
  )

  grad.analytical <- analyticalGradQmlCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric
  )

  grad.fd <- gradLogLikQmlFDCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric,
    eps = 1e-7,
    ncores = 1L
  )

  testthat::expect_true(all(is.finite(grad.analytical)))
  testthat::expect_equal(grad.analytical, grad.fd, tolerance = 3e-3)

  tpb.linear <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  '

  group.info <- parseModelArgumentsByGroupDA(
    model.syntax = tpb.linear,
    cov.syntax = NULL,
    method = "qml",
    data = TPB,
    group = NULL,
    auto.split.syntax = FALSE
  )

  model <- specifyModelDA(
    group.info = group.info,
    method = "qml",
    m = 16,
    mean.observed = TRUE,
    double = FALSE,
    quad.range = Inf,
    adaptive.quad = FALSE,
    adaptive.frequency = 3,
    missing = "complete",
    orthogonal.x = FALSE,
    orthogonal.y = FALSE,
    auto.fix.first = TRUE,
    auto.fix.single = TRUE,
    fix.composite.var = TRUE,
    cluster = NULL,
    sampling.weights = NULL
  )

  model.filled <- fillModel(model = model, theta = model$theta, method = "qml")
  locations <- subset(
    model$params$gradientStruct$locations,
    group == 1 & block == 11
  )

  grad.analytical <- analyticalGradQmlCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric
  )

  grad.fd <- gradLogLikQmlFDCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric,
    eps = 1e-6,
    ncores = 1L
  )

  testthat::expect_true(all(is.finite(grad.analytical)))
  testthat::expect_equal(grad.analytical, grad.fd, tolerance = 2e-3)

  tpb.interaction <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
  BEH ~ ATT:PBC
  ATT ~ 1
  PBC ~ 1
  '

  group.info <- parseModelArgumentsByGroupDA(
    model.syntax = tpb.interaction,
    cov.syntax = NULL,
    method = "qml",
    data = TPB,
    group = NULL,
    auto.split.syntax = FALSE
  )

  model <- specifyModelDA(
    group.info = group.info,
    method = "qml",
    m = 16,
    mean.observed = TRUE,
    double = FALSE,
    quad.range = Inf,
    adaptive.quad = FALSE,
    adaptive.frequency = 3,
    missing = "complete",
    orthogonal.x = FALSE,
    orthogonal.y = FALSE,
    auto.fix.first = TRUE,
    auto.fix.single = TRUE,
    fix.composite.var = TRUE,
    cluster = NULL,
    sampling.weights = NULL
  )

  model.filled <- fillModel(model = model, theta = model$theta, method = "qml")
  locations <- subset(
    model$params$gradientStruct$locations,
    group == 1 & block %in% c(2, 9, 11, 12, 13)
  )

  grad.analytical <- analyticalGradQmlCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric
  )

  grad.fd <- gradLogLikQmlFDCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric,
    eps = 1e-6,
    ncores = 1L
  )

  testthat::expect_true(model.filled$models[[1]]$info$kOmegaEta > 0L)
  testthat::expect_equal(sort(unique(locations$block)),
                         c(2, 9, 11, 12, 13))
  testthat::expect_true(all(is.finite(grad.analytical)))
  testthat::expect_equal(grad.analytical, grad.fd, tolerance = 2e-3)
})

testthat::test_that("split QML uses hybrid gradient through covModel Jacobian", {
  tpb <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
  '

  group.info <- parseModelArgumentsByGroupDA(
    model.syntax = tpb,
    cov.syntax = NULL,
    method = "qml",
    data = TPB,
    group = NULL,
    auto.split.syntax = TRUE
  )

  model <- specifyModelDA(
    group.info = group.info,
    method = "qml",
    m = 16,
    mean.observed = TRUE,
    double = FALSE,
    quad.range = Inf,
    adaptive.quad = FALSE,
    adaptive.frequency = 3,
    missing = "complete",
    orthogonal.x = FALSE,
    orthogonal.y = FALSE,
    auto.fix.first = TRUE,
    auto.fix.single = TRUE,
    fix.composite.var = TRUE,
    cluster = NULL,
    sampling.weights = NULL
  )

  grad.hybrid <- gradientLogLikQml(
    theta = model$theta,
    model = model,
    sign = 1,
    epsilon = 1e-6
  )

  objective <- function(theta) {
    model.filled <- fillModel(model = model, theta = theta, method = "qml")
    logLikQmlCpp(model.filled$models[[1L]], ncores = 1L)
  }

  f0 <- objective(model$theta)
  grad.fd <- numeric(length(model$theta))
  indices <- c(
    model$params$SELECT_THETA_COV[[1L]],
    model$params$SELECT_THETA_MAIN[[1L]]
  )

  for (i in indices) {
    theta.i <- model$theta
    theta.i[[i]] <- theta.i[[i]] + 1e-6
    grad.fd[[i]] <- (objective(theta.i) - f0) / 1e-6
  }

  testthat::expect_true(isTRUE(model$params$gradientStruct$hasCovModel))
  testthat::expect_true(all(is.finite(grad.hybrid)))
  testthat::expect_equal(grad.hybrid[indices], grad.fd[indices],
                         tolerance = 3e-3)
})

testthat::test_that("split QML uses hybrid observation scores through covModel Jacobian", {
  tpb <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
  '

  group.info <- parseModelArgumentsByGroupDA(
    model.syntax = tpb,
    cov.syntax = NULL,
    method = "qml",
    data = TPB,
    group = NULL,
    auto.split.syntax = TRUE
  )

  model <- specifyModelDA(
    group.info = group.info,
    method = "qml",
    m = 16,
    mean.observed = TRUE,
    double = FALSE,
    quad.range = Inf,
    adaptive.quad = FALSE,
    adaptive.frequency = 3,
    missing = "complete",
    orthogonal.x = FALSE,
    orthogonal.y = FALSE,
    auto.fix.first = TRUE,
    auto.fix.single = TRUE,
    fix.composite.var = TRUE,
    cluster = NULL,
    sampling.weights = NULL
  )

  scores.hybrid <- gradientLogLikQml_i(
    theta = model$theta,
    model = model,
    sign = 1,
    epsilon = 1e-6
  )

  obs_objective <- function(theta) {
    model.filled <- fillModel(model = model, theta = theta, method = "qml")
    logLikQmlGroup(model.filled$models[[1L]], sign = 1, sum = FALSE)
  }

  f0 <- obs_objective(model$theta)
  indices <- c(
    model$params$SELECT_THETA_COV[[1L]],
    model$params$SELECT_THETA_MAIN[[1L]]
  )
  indices <- indices[seq_len(min(12L, length(indices)))]

  scores.fd <- matrix(0, nrow = NROW(scores.hybrid), ncol = length(indices))
  colnames(scores.fd) <- names(model$theta)[indices]

  for (j in seq_along(indices)) {
    theta.j <- model$theta
    theta.j[[indices[[j]]]] <- theta.j[[indices[[j]]]] + 1e-6
    scores.fd[, j] <- (obs_objective(theta.j) - f0) / 1e-6
  }

  testthat::expect_true(isTRUE(model$params$gradientStruct$hasCovModel))
  testthat::expect_true(all(is.finite(scores.hybrid[, indices])))
  testthat::expect_equal(scores.hybrid[, indices], scores.fd,
                         tolerance = 1e-4)
})

testthat::test_that("unsplit QML uses analytical non-constant-Binv observation scores", {
  tpb <- '
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
  BEH ~ ATT:PBC
  ATT ~ 1
  PBC ~ 1
  '

  group.info <- parseModelArgumentsByGroupDA(
    model.syntax = tpb,
    cov.syntax = NULL,
    method = "qml",
    data = TPB,
    group = NULL,
    auto.split.syntax = FALSE
  )

  model <- specifyModelDA(
    group.info = group.info,
    method = "qml",
    m = 16,
    mean.observed = TRUE,
    double = FALSE,
    quad.range = Inf,
    adaptive.quad = FALSE,
    adaptive.frequency = 3,
    missing = "complete",
    orthogonal.x = FALSE,
    orthogonal.y = FALSE,
    auto.fix.first = TRUE,
    auto.fix.single = TRUE,
    fix.composite.var = TRUE,
    cluster = NULL,
    sampling.weights = NULL
  )

  model.filled <- fillModel(model = model, theta = model$theta, method = "qml")
  locations <- subset(
    model$params$gradientStruct$locations,
    group == 1 & block %in% c(2, 9, 11, 12, 13)
  )

  scores.analytical <- analyticalObsGradQmlCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric
  )

  scores.fd <- gradObsLogLikQmlFDCpp(
    submodel = model.filled$models[[1]],
    block = locations$block,
    row = locations$row,
    col = locations$col,
    symmetric = locations$symmetric,
    eps = 1e-6,
    ncores = 1L
  )

  testthat::expect_true(model.filled$models[[1]]$info$kOmegaEta > 0L)
  testthat::expect_equal(sort(unique(locations$block)),
                         c(2, 9, 11, 12, 13))
  testthat::expect_true(all(is.finite(scores.analytical)))
  testthat::expect_equal(scores.analytical, scores.fd, tolerance = 1e-4)
})
