gaussianDensity <- function(nodes, ...) {
  # Three "observations" whose posteriors sit at different places, so a
  # quasi-adaptive rule has something to adapt to.
  centres <- c(-1, 0, 1)
  t(vapply(centres, function(mu)
    exp(-0.5 * rowSums((nodes - mu)^2)), numeric(NROW(nodes))))
}


testthat::test_that("fixed product rules integrate polynomials of the normal", {
  for (integration in c("gh", "rect")) {
    spec <- quadSpec(integration = integration, adaptive = "none",
                     nodes = 25L, k = 2L)
    rule <- buildQuadRule(spec)

    testthat::expect_equal(NROW(rule$n), 625L)
    testthat::expect_equal(NCOL(rule$n), 2L)
    testthat::expect_true(rule$common)
    testthat::expect_false(rule$packed)

    # A rule over the standard normal must reproduce its first two moments.
    testthat::expect_equal(sum(rule$w), 1, tolerance = 1e-10)
    testthat::expect_equal(as.numeric(colSums(rule$n * rule$w)), c(0, 0),
                           tolerance = 1e-10)
    testthat::expect_equal(as.numeric(colSums(rule$n^2 * rule$w)), c(1, 1),
                           tolerance = if (integration == "gh") 1e-10 else 1e-3)
    testthat::expect_equal(sum(rule$n[, 1L] * rule$n[, 2L] * rule$w), 0,
                           tolerance = 1e-10)
  }
})


testthat::test_that("the rectangular rule spans the requested range", {
  rule <- buildQuadRule(quadSpec(integration = "rect", adaptive = "none",
                                 nodes = 15L, k = 1L, rect.range = 5))
  testthat::expect_equal(range(rule$n[, 1L]), c(-5, 5))
  testthat::expect_equal(length(unique(round(diff(sort(rule$n[, 1L])), 10))), 1L)
})


testthat::test_that("quasi-adaptive rules stay near the requested node budget", {
  gh <- buildQuadRule(
    quadSpec(integration = "gh", adaptive = "quasi", nodes = 12L, k = 2L),
    density = gaussianDensity
  )
  # The pruner targets `nodes` per dimension; it is allowed to miss by a
  # couple of nodes per axis but must not blow past the fixed-rule budget.
  testthat::expect_lte(NROW(gh$n), 16L^2L)
  testthat::expect_gte(NROW(gh$n), 8L^2L)
  testthat::expect_true(gh$common)

  rect <- buildQuadRule(
    quadSpec(integration = "rect", adaptive = "quasi", nodes = 12L, k = 2L),
    density = gaussianDensity
  )
  # The rectangular rule adapts its range and then drops corner nodes -- those
  # extreme in every dimension at once, which carry almost none of the integral
  # -- so it lands at or below the requested budget rather than exactly on it.
  testthat::expect_lte(NROW(rect$n), 144L)
  testthat::expect_gte(NROW(rect$n), 64L)
  testthat::expect_lte(max(abs(rect$n)), 5)
  testthat::expect_true(rect$common)
})


testthat::test_that("product rules drop corner nodes", {
  # A corner of a Cartesian grid is extreme in every dimension at once, so its
  # weight is the product of the extreme one-dimensional weights and it carries
  # essentially none of the integral. Both product rules must shed them.
  for (integration in c("gh", "rect")) {
    spec <- quadSpec(integration = integration, adaptive = "quasi",
                     nodes = 10L, k = 3L)
    rule <- buildQuadRule(spec, density = gaussianDensity)
    testthat::expect_lt(NROW(rule$n), 10L^3L)

    # An adapted rule covers a trimmed region, so its weights sum to the mass
    # it spans rather than to 1. Renormalising them would turn the estimate
    # into a conditional expectation -- see the closed-form test below.
    testthat::expect_lte(sum(rule$w), 1 + 1e-8)
    testthat::expect_gt(sum(rule$w), 0.5)
  }
})


testthat::test_that("the rectangular quasi rule shrinks onto concentrated mass", {
  # A posterior far tighter than the +-5 default leaves most of the fixed grid
  # empty, which is exactly the case range adaptation exists to fix.
  concentrated <- function(nodes, ...)
    t(vapply(c(-0.3, 0, 0.3), function(mu)
      exp(-0.5 * rowSums((nodes - mu)^2) / 0.04), numeric(NROW(nodes))))

  spec <- quadSpec(integration = "rect", adaptive = "quasi", nodes = 12L,
                   k = 1L, adaptive.quad.tol = 1e-6)
  rule <- buildQuadRule(spec, density = concentrated)

  testthat::expect_equal(NROW(rule$n), 12L)
  testthat::expect_lt(max(abs(rule$n)), 5)
  testthat::expect_lte(sum(rule$w), 1 + 1e-8)
  # Nodes stay equispaced within the trimmed range.
  testthat::expect_equal(length(unique(round(diff(sort(rule$n[, 1L])), 8))), 1L)
})


testthat::test_that("a quasi-adaptive rule needs a density function", {
  testthat::expect_error(
    buildQuadRule(quadSpec(integration = "gh", adaptive = "quasi", k = 2L)),
    "needs a density function"
  )
})


testthat::test_that("Monte-Carlo draws are reproducible and leave the RNG alone", {
  spec <- quadSpec(integration = "mc", adaptive = "none", nodes = 500L, k = 3L)

  # The same spec must always give the same draws, so the EM likelihood is a
  # deterministic function of the parameters.
  first <- buildQuadRule(spec)
  second <- buildQuadRule(spec)
  testthat::expect_identical(first$n, second$n)
  testthat::expect_equal(NROW(first$n), 500L)
  testthat::expect_equal(NCOL(first$n), 3L)
  testthat::expect_equal(sum(first$w), 1)

  # Drawing the rule must not disturb the caller's random stream.
  set.seed(99)
  expected <- stats::runif(2L)
  set.seed(99)
  observed <- c(stats::runif(1L), { buildQuadRule(spec); stats::runif(1L) })
  testthat::expect_equal(observed, expected)

  spec$seed <- 4321L
  testthat::expect_false(identical(buildQuadRule(spec)$n, first$n))
})


testthat::test_that("Monte-Carlo node counts are totals, not per dimension", {
  testthat::expect_equal(quadTotalNodes(
    quadSpec(integration = "mc", adaptive = "none", nodes = 5000L, k = 4L)), 5000)
  testthat::expect_equal(quadTotalNodes(
    quadSpec(integration = "gh", adaptive = "none", nodes = 15L, k = 3L)), 15^3)
})


testthat::test_that("an oversized product rule errors and names the alternatives", {
  spec <- quadSpec(integration = "gh", adaptive = "none", nodes = 15L, k = 8L)
  testthat::expect_error(buildQuadRule(spec), "exceeds the limit")
  testthat::expect_error(buildQuadRule(spec), "integration = \"mc\"")
})


testthat::test_that("the availability matrix rejects unsupported combinations", {
  full.gh <- quadSpec(integration = "gh", adaptive = "full", k = 2L)
  testthat::expect_error(checkQuadAvailable(full.gh, "legacy"),
                         "lms.backend = \"graph\"")
  testthat::expect_silent(checkQuadAvailable(full.gh, "graph"))

  mc.legacy <- quadSpec(integration = "mc", adaptive = "none", k = 2L)
  testthat::expect_error(checkQuadAvailable(mc.legacy, "legacy"),
                         "not available for the legacy LMS backend")

  mc.quasi <- quadSpec(integration = "mc", adaptive = "quasi", k = 2L)
  testthat::expect_error(checkQuadAvailable(mc.quasi, "graph"),
                         "adaptive = \"none\" or \"full\"")

  for (backend in c("legacy", "graph")) for (integration in c("gh", "rect"))
    testthat::expect_silent(checkQuadAvailable(
      quadSpec(integration = integration, adaptive = "quasi", k = 2L), backend))
})


testthat::test_that("node defaults follow integration, adaptive and backend", {
  nodesFor <- function(...) getMethodSettingsDA(
    "lms", args = c(list(...), list(nodes = NULL)))$nodes

  # Shared rules: the legacy backend integrates only the nonlinear dimensions,
  # so it can afford -- and needs -- more points than the graph backend.
  for (integration in c("gh", "rect")) for (adaptive in c("none", "quasi")) {
    testthat::expect_equal(nodesFor(lms.backend = "legacy",
      integration = integration, adaptive = adaptive), 24)
    testthat::expect_equal(nodesFor(lms.backend = "graph",
      integration = integration, adaptive = adaptive), 15)
  }
  # Per-observation and Monte-Carlo counts do not depend on the backend.
  for (backend in c("legacy", "graph")) {
    testthat::expect_equal(nodesFor(lms.backend = backend,
      integration = "gh", adaptive = "full"), 5)
    testthat::expect_equal(nodesFor(lms.backend = backend,
      integration = "mc", adaptive = "none"), 500)
    testthat::expect_equal(nodesFor(lms.backend = backend,
      integration = "mc", adaptive = "full"), 250)
  }

  # An explicit value always wins.
  testthat::expect_equal(
    getMethodSettingsDA("lms", args = list(integration = "mc",
                                           adaptive = "full",
                                           nodes = 42))$nodes, 42)
})


testthat::test_that("the model stores a quadrature spec rather than a node grid", {
  spec <- quadSpec(integration = "gh", adaptive = "quasi", nodes = 15L, k = 3L)
  testthat::expect_null(spec$n)
  testthat::expect_null(spec$w)
  testthat::expect_true(isAdaptiveQuad(spec))
  testthat::expect_false(isAdaptiveQuad(setAdaptiveQuad(spec, "none")))
})


testthat::test_that("rules reproduce closed-form integrals of the normal", {
  # Agreement with another approximation cannot catch a bias that affects both.
  # These integrands have exact values, so they can.
  #
  #   exp(a'z)                      E[f] = exp(a'a / 2)      smooth, fat tails
  #   exp(-|z - mu|^2 / (2 s2))     E[f] = (s2/(1+s2))^(k/2) *
  #                                        exp(-|mu|^2 / (2(1+s2)))
  #
  # The second is shaped like a real posterior: narrow and off-centre. It is the
  # case that exposed the weight-renormalisation bias in the rectangular
  # quasi-adaptive rule, where adaptation made the rule WORSE than not adapting.
  k <- 2L
  exponential <- list(
    f = function(z, a = rep(1, k)) exp(as.vector(z %*% a)),
    true = exp(sum(rep(1, k)^2) / 2)
  )
  bump <- local({
    mu <- rep(1, k); s2 <- 0.25
    list(f = function(z) exp(-0.5 * rowSums(sweep(z, 2L, mu)^2) / s2),
         true = (s2 / (1 + s2))^(k / 2) * exp(-sum(mu^2) / (2 * (1 + s2))))
  })

  for (case in list(exponential, bump)) {
    density <- function(nodes, ...) matrix(case$f(nodes), nrow = 1L)
    estimate <- function(rule) sum(case$f(rule$n) * rule$w)

    for (integration in c("gh", "rect")) {
      fixed <- buildQuadRule(quadSpec(integration = integration,
                                      adaptive = "none", nodes = 15L, k = k))
      quasi <- buildQuadRule(quadSpec(integration = integration,
                                      adaptive = "quasi", nodes = 15L, k = k),
                             density = density)
      error <- function(rule) abs(estimate(rule) - case$true) / abs(case$true)

      # Every rule must at least be in the right ballpark.
      testthat::expect_lt(error(fixed), 1e-1)
      testthat::expect_lt(error(quasi), 1e-1)
      # Adapting must never make a rule worse than not adapting. This is the
      # signature the renormalisation bug produced.
      testthat::expect_lte(error(quasi), max(error(fixed), 1e-3))
    }
  }
})
