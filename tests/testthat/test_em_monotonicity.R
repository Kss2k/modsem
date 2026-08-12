# EM increases the observed log-likelihood at every iteration. That is the
# defining property of the algorithm, not a quality-of-fit nicety: a decrease
# means the M-step lowered Q, and every convergence guarantee is gone.
#
# This exists because it did happen. The ECM M-step split (#20) put the
# structural parameters in a block optimised through a surrogate built from
# E-step sufficient statistics, on the reading that they were separable from
# the measurement kernel once the latent nodes were fixed. They are not: the
# graph backend's nodes are (xi, zeta), so eta is RECONSTRUCTED from them
# through alpha/Gamma/Omega, and every structural parameter moves the indicator
# means. Perturbing structural parameters alone changed the measurement
# objective in 39 of 39 probes, and the joint objective's entire response to
# such a move lived in the measurement half (d.joint + d.measurement = 0 to
# 7e-11). The surrogate missed that by a median of 205%.
#
# The symptom was a plain-EM run decreasing the observed log-likelihood in 134
# of 299 iterations (#27). Nothing caught it for the whole life of the split,
# which is why this file exists.
#
# Keep this cheap enough to run every time. It is a guard on the M-step, not a
# fit test, so a small model and few iterations are the point.

skip_if(Sys.getenv("MODSEM_SKIP_SLOW_TESTS") == "true")

monotonicityViolations <- function(fit, tolerance = 1e-8) {
  history <- fit$iteration.history
  expect_true(!is.null(history) && NROW(history) > 2L)
  steps <- diff(history$loglik)
  steps[steps < -tolerance]
}

# Observed-log-likelihood monotonicity is the SYMPTOM, and it is too weak to
# catch this on its own: the default path kept a monotone log-likelihood while
# lowering Q in 28 of 139 M-steps, because Delta-logLik >= Delta-Q only removes
# the guarantee when Delta-Q is negative rather than forcing a decrease. So
# assert the generalised-EM property directly -- the M-step must not decrease
# the complete-data objective it was handed.
#
# Q is evaluated with ONE function (backend$complete) at the M-step's own input
# and output against the M-step's own posterior, so the two numbers are
# comparable. Values from different iterations are NOT: Q(theta | theta_t) uses
# each iteration's own weights.
# Both candidates are instrumented so this keeps testing the M-step that is
# actually wired up, rather than one the backend may no longer call. The depth
# guard matters because the ECM step delegates to the joint one when the block
# split is degenerate, which would otherwise record the same step twice.
mstepObjectiveGains <- function(expr, steps = c("mstepLms", "mstepLmsGraphEcm")) {
  namespace <- asNamespace("modsem")
  gains <- numeric(0)
  depth <- 0L

  originals <- stats::setNames(lapply(steps, get, envir = namespace), steps)
  on.exit({
    for (name in steps) {
      assign(name, originals[[name]], envir = namespace)
      lockBinding(name, namespace)
    }
  }, add = TRUE)

  for (name in steps) {
    original <- originals[[name]]
    recording <- local({
      inner <- original
      function(theta, model, P, max.step, ...) {
        depth <<- depth + 1L
        outermost <- depth == 1L
        on.exit(depth <<- depth - 1L, add = TRUE)
        out <- inner(theta, model, P, max.step, ...)
        if (outermost) {
          backend <- getLmsBackend("graph")
          Q <- function(x) backend$complete(theta = x, model = model, P = P,
                                            sign = -1, epsilon = 1e-6)
          # sign = -1, so this is a minimisation: positive gain = improved
          gains <<- c(gains, Q(theta) - Q(out$par))
        }
        out
      }
    })
    unlockBinding(name, namespace)
    assign(name, recording, envir = namespace)
  }

  force(expr)
  gains
}

set.seed(123)

test_that("plain EM never decreases the observed log-likelihood (continuous)", {
  fit <- modsem(
    'X =~ x1 + x2
     Z =~ z1 + z2
     Y =~ y1 + y2
     Y ~ X + Z + X:Z',
    oneInt[c("x1", "x2", "z1", "z2", "y1", "y2")],
    method = "lms", lms.backend = "graph",
    adaptive = "full", integration = "gh", nodes = 8L,
    algorithm = "EM", max.iter = 40L, calc.se = FALSE, verbose = FALSE
  )

  violations <- monotonicityViolations(fit)
  expect_equal(length(violations), 0L,
               info = sprintf("%d decreasing EM steps, worst %.6g",
                              length(violations),
                              if (length(violations)) min(violations) else 0))
})


test_that("plain EM never decreases the observed log-likelihood (ordered)", {
  # Ordered indicators are where this was found: the threshold parameters make
  # the measurement block much stiffer, so a structural step that lowers Q
  # shows up immediately rather than being absorbed.
  rthreshold <- function(k, offset = stats::runif(1, -1, 1), sigma = 0.35) {
    t <- seq_len(k) - mean(seq_len(k)) + offset
    c(-Inf, t + stats::runif(k, -sigma, sigma), Inf)
  }
  vars <- c("x1", "x2", "z1", "z2", "y1", "y2")
  data <- oneInt[vars]
  for (v in vars) {
    z <- (data[[v]] - mean(data[[v]])) / stats::sd(data[[v]])
    data[[v]] <- cut(z, breaks = rthreshold(4L), ordered_result = TRUE)
  }

  fit <- modsem(
    'X =~ x1 + x2
     Z =~ z1 + z2
     Y =~ y1 + y2
     Y ~ X + Z + X:Z',
    data, method = "lms", lms.backend = "graph", ordered = vars,
    adaptive = "full", integration = "gh", nodes = 8L,
    algorithm = "EM", max.iter = 40L, calc.se = FALSE, verbose = FALSE
  )

  violations <- monotonicityViolations(fit)
  expect_equal(length(violations), 0L,
               info = sprintf("%d decreasing EM steps, worst %.6g",
                              length(violations),
                              if (length(violations)) min(violations) else 0))
})


test_that("the graph M-step never decreases the complete-data objective", {
  rthreshold <- function(k, offset = stats::runif(1, -1, 1), sigma = 0.35) {
    t <- seq_len(k) - mean(seq_len(k)) + offset
    c(-Inf, t + stats::runif(k, -sigma, sigma), Inf)
  }
  vars <- c("x1", "x2", "z1", "z2", "y1", "y2")
  data <- oneInt[vars]
  for (v in vars) {
    z <- (data[[v]] - mean(data[[v]])) / stats::sd(data[[v]])
    data[[v]] <- cut(z, breaks = rthreshold(4L), ordered_result = TRUE)
  }

  gains <- mstepObjectiveGains(modsem(
    'X =~ x1 + x2
     Z =~ z1 + z2
     Y =~ y1 + y2
     Y ~ X + Z + X:Z',
    data, method = "lms", lms.backend = "graph", ordered = vars,
    adaptive = "full", integration = "gh", nodes = 8L,
    algorithm = "EM", max.iter = 40L, calc.se = FALSE, verbose = FALSE
  ))

  expect_true(length(gains) > 5L)
  decreases <- gains[gains < -1e-9]
  expect_equal(length(decreases), 0L,
               info = sprintf("M-step lowered Q in %d of %d steps, worst %.6g",
                              length(decreases), length(gains),
                              if (length(decreases)) min(decreases) else 0))
})
