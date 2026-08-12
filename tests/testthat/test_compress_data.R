devtools::load_all()

# `compress.data` collapses duplicate rows into one weighted row. The claim it
# rests on is that this is exact rather than an approximation: identical rows
# contribute identical terms to the likelihood. These tests pin down both the
# bookkeeping (weights, index, cluster, n.obs) and that claim end to end.

withAttrs <- function(m, cluster = NULL, weights = NULL) {
  attr(m, "cluster") <- cluster
  attr(m, "weights") <- weights
  attr(m, "index")   <- seq_len(NROW(m))
  m
}


test_that("compressData collapses duplicates and conserves weight", {
  m <- withAttrs(matrix(c(1, 2, 1, 3, 2, 1, 2, 4), ncol = 2,
                        dimnames = list(NULL, c("a", "b"))))
  # rows: (1,2) (2,1) (1,2) (3,4)  -> three distinct rows
  cm <- suppressMessages(compressData(m, compress = TRUE))

  expect_equal(NROW(cm), 3L)
  expect_equal(attr(cm, "weights"), c(2, 1, 1))
  expect_equal(attr(cm, "n.obs"), 4L)

  # the index must partition the original rows, and the listed rows must be
  # exactly the ones that were merged
  index <- attr(cm, "index")
  expect_true(is.list(index))
  expect_equal(sort(unlist(index)), 1:4)
  expect_equal(lengths(index), attr(cm, "weights"))
  expect_equal(index[[1L]], c(1L, 3L))

  # expanding by the weights must reproduce the input rows exactly
  rebuilt <- cm[rep(seq_len(NROW(cm)), attr(cm, "weights")), , drop = FALSE]
  expect_equal(unname(m[unlist(index), ]), unname(rebuilt))
})


test_that("compressData is a no-op when nothing is duplicated", {
  m <- withAttrs(matrix(1:8, ncol = 2))
  expect_identical(compressData(m, compress = TRUE), m)
  expect_identical(compressData(m, compress = FALSE), m)
})


test_that("compressData carries existing sampling weights through", {
  m <- withAttrs(matrix(c(1, 1, 2, 5, 5, 6), ncol = 2),
                 weights = c(0.5, 1.5, 2))
  cm <- suppressMessages(compressData(m, compress = TRUE))

  expect_equal(NROW(cm), 2L)
  expect_equal(attr(cm, "weights"), c(2, 2)) # 0.5 + 1.5, and 2
})


test_that("compressData never merges rows from different clusters", {
  # identical indicator rows, but two clusters: collapsing them would destroy
  # the nesting the cluster-robust variance estimator depends on
  m <- withAttrs(matrix(rep(c(1, 2), each = 4), ncol = 2),
                 cluster = factor(c("a", "a", "b", "b")))
  cm <- suppressMessages(compressData(m, compress = TRUE))

  expect_equal(NROW(cm), 2L)
  expect_equal(attr(cm, "weights"), c(2, 2))
  expect_equal(as.character(attr(cm, "cluster")), c("a", "b"))
})


test_that("compressData does not merge rows that merely look equal", {
  # differ in the 16th significant digit: as.character() would print these
  # identically, so a naive key would silently merge two distinct rows
  x <- c(1, 1 + 2 * .Machine$double.eps, 1)
  m <- withAttrs(cbind(x, 0))
  cm <- suppressMessages(compressData(m, compress = TRUE))

  expect_equal(NROW(cm), 2L)
  expect_equal(attr(cm, "weights"), c(2, 1))
})


test_that("compress.data defaults to TRUE only for the lms graph backend", {
  graph  <- getMethodSettingsDA("lms", args = list(lms.backend = "graph"))
  legacy <- getMethodSettingsDA("lms", args = list(lms.backend = "legacy"))
  qml    <- getMethodSettingsDA("qml", args = list())

  expect_true(graph$compress.data)
  expect_false(legacy$compress.data)
  expect_false(qml$compress.data)

  # an explicit value always wins over the backend default
  expect_false(
    getMethodSettingsDA("lms", args = list(lms.backend = "graph",
                                           compress.data = FALSE))$compress.data
  )
})


test_that("compressing duplicate rows leaves the fit unchanged", {
  set.seed(1234)
  data <- oneInt[sample(NROW(oneInt), 300L), ]
  data <- rbind(data, data[1:80, ]) # 80 exact duplicates

  syntax <- "X =~ x1 + x2 + x3
             Z =~ z1 + z2 + z3
             Y =~ y1 + y2 + y3
             Y ~ X + Z + X:Z"

  fit <- function(compress) suppressWarnings(suppressMessages(modsem(
    syntax, data, method = "lms", lms.backend = "graph", adaptive = "full",
    nodes = 5L, calc.se = TRUE, verbose = FALSE, compress.data = compress)))

  uncompressed <- fit(FALSE)
  compressed   <- fit(TRUE)

  # the two runs sum the same terms in a different order, so they agree to
  # floating-point reassociation rather than bit-for-bit
  expect_equal(as.numeric(compressed$logLik), as.numeric(uncompressed$logLik),
               tolerance = 1e-6)

  coef.u <- parameter_estimates(uncompressed)
  coef.c <- parameter_estimates(compressed)
  expect_equal(coef.c$est, coef.u$est, tolerance = 1e-4)
  expect_equal(coef.c$std.error, coef.u$std.error, tolerance = 1e-4)

  # estimation used 300 rows, but 380 observations were supplied
  expect_equal(NROW(compressed$model$models[[1L]]$data$data.full), 300L)
  expect_equal(nobs(compressed), 380)
  expect_equal(nobs(uncompressed), 380)

  # ... and the sample size that reaches the fit measures is the observation
  # count, not the number of distinct rows
  fit.u <- modsem:::fit_modsem_da(uncompressed, chisq = FALSE)
  fit.c <- modsem:::fit_modsem_da(compressed, chisq = FALSE)
  expect_equal(fit.c$BIC, fit.u$BIC, tolerance = 1e-4)
  expect_equal(fit.c$aBIC, fit.u$aBIC, tolerance = 1e-4)
})


test_that("predictions are expanded back to one row per observation", {
  set.seed(1234)
  data <- oneInt[sample(NROW(oneInt), 300L), ]
  data <- rbind(data, data[1:80, ])

  syntax <- "X =~ x1 + x2 + x3
             Z =~ z1 + z2 + z3
             Y =~ y1 + y2 + y3
             Y ~ X + Z + X:Z"

  fit <- function(compress) suppressWarnings(suppressMessages(modsem(
    syntax, data, method = "lms", nodes = 8L, calc.se = FALSE,
    verbose = FALSE, compress.data = compress)))

  pred.u <- modsem_predict(fit(FALSE))
  pred.c <- modsem_predict(fit(TRUE))

  expect_equal(NROW(pred.c), NROW(data))
  expect_equal(unname(as.matrix(pred.c)), unname(as.matrix(pred.u)),
               tolerance = 1e-6)

  # the duplicated block must come back with the predictions of the rows it
  # was collapsed into
  expect_equal(unname(as.matrix(pred.c)[1:80, ]),
               unname(as.matrix(pred.c)[301:380, ]))
})


test_that("ordered data compresses hard and gives the same fit", {
  # the case this feature is aimed at: a finite row space, so most rows repeat
  rthreshold <- function(k, offset = stats::runif(1, -1, 1), sigma = 0.35) {
    t <- seq_len(k) - mean(seq_len(k)) + offset
    c(-Inf, t + stats::runif(k, -sigma, sigma), Inf)
  }
  cut_data <- function(d, k = 5, choose) {
    z <- function(x) (x - mean(x)) / stats::sd(x)
    for (v in choose) d[[v]] <- cut(z(d[[v]]), breaks = rthreshold(k),
                                    ordered_result = TRUE)
    d
  }

  vars <- colnames(oneInt)[!grepl("3", colnames(oneInt))]
  set.seed(2837290)
  data <- cut_data(oneInt, choose = vars)[vars]

  syntax <- "X =~ x1 + x2
             Z =~ z1 + z2
             Y =~ y1 + y2
             Y ~ X + Z + X:Z"

  fit <- function(compress) suppressWarnings(suppressMessages(modsem(
    syntax, data, method = "lms", lms.backend = "graph", ordered = vars,
    adaptive = "full", nodes = 10L, algorithm = "EM", max.iter = 5L,
    calc.se = FALSE, verbose = FALSE, compress.data = compress)))

  uncompressed <- fit(FALSE)
  compressed   <- fit(TRUE)

  expect_equal(NROW(compressed$model$models[[1L]]$data$data.full), 733L)
  expect_equal(nobs(compressed), NROW(data))
  expect_equal(as.numeric(compressed$logLik), as.numeric(uncompressed$logLik),
               tolerance = 1e-6)
  expect_equal(parameter_estimates(compressed)$est,
               parameter_estimates(uncompressed)$est, tolerance = 1e-3)
})
