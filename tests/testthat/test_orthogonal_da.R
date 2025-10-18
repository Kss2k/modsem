devtools::load_all()

set.seed(123)
n <- 250
X <- rnorm(n, mean = 0, sd = sqrt(1.2))
Z <- rnorm(n, mean = 0, sd = sqrt(1.4))

Y1 <- 0.7 * X + 0.3 * Z + 0.5 * X * Z + rnorm(n, mean = 0, sd = sqrt(0.4))
Y2 <- 0.3 * X + 0.7 * Z + 0.1 * X * Z + rnorm(n, mean = 0, sd = sqrt(0.6))

x1 <-       X + 1.2 + rnorm(n, mean = 0, sd = sqrt(0.2))
x2 <- 0.8 * X + 0.9 + rnorm(n, mean = 0, sd = sqrt(0.2))

z1 <-       Z + 1.2 + rnorm(n, mean = 0, sd = sqrt(0.2))
z2 <- 0.8 * Z + 0.9 + rnorm(n, mean = 0, sd = sqrt(0.2))

y1 <-       Y1 + 1.2 + rnorm(n, mean = 0, sd = sqrt(0.2))
y2 <- 0.8 * Y1 + 0.9 + rnorm(n, mean = 0, sd = sqrt(0.2))

y3 <-       Y2 + 1.2 + rnorm(n, mean = 0, sd = sqrt(0.2))
y4 <- 0.8 * Y2 + 0.9 + rnorm(n, mean = 0, sd = sqrt(0.2))

data <- data.frame(x1, x2, z1, z2, y1, y2, y3, y4)

m_linear <- '
 X =~ x1 + x2
 Z =~ z1 + z2
 Y1 =~ y1 + y2
 Y2 =~ y3 + y4

 Y1 ~ X + Z
 Y2 ~ a * X + b * Z

 a < .33
 b > .721
 a > 0
 b < 1
'

lavaan::summary(lavaan::sem(m_linear, data))
testthat::expect_no_warning({
  lms_linear_obl <- modsem(m_linear, data = data, method = "lms")
  lms_linear_ort <- modsem(m_linear, data = data, method = "lms",
                           orthogonal.x = TRUE, orthogonal.y = TRUE)
  qml_linear_obl <- modsem(m_linear, data = data, method = "qml")
  qml_linear_ort <- modsem(m_linear, data = data, method = "qml",
                           orthogonal.x = TRUE, orthogonal.y = TRUE)
})

b.lower <- lms_linear_obl$model$params$bounds$lower
b.upper <- lms_linear_obl$model$params$bounds$upper

testthat::expect_equal(b.lower[["a"]], 0.000)
testthat::expect_equal(b.upper[["a"]], 0.330)
testthat::expect_equal(b.upper[["b"]], 1.000)
testthat::expect_equal(b.lower[["b"]], 0.721)

checkIter <- \(m, n) testthat::expect_equal(m$iterations, n)
checkIter(lms_linear_obl, 2L)
checkIter(lms_linear_ort, 2L)
checkIter(qml_linear_obl, 1L)
checkIter(qml_linear_ort, 1L)

m_nlinear <- '
 X =~ x1 + x2
 Z =~ z1 + z2
 Y1 =~ y1 + y2
 Y2 =~ y3 + y4

 Y1 ~ X + Z + X:Z
 Y2 ~ X + Z + X:Z
'

testthat::expect_no_warning({
  lms_nlinear_obl <- modsem(m_nlinear, data = data, method = "lms")
  lms_nlinear_ort <- modsem(m_nlinear, data = data, method = "lms",
                           orthogonal.x = TRUE, orthogonal.y = TRUE)
  qml_nlinear_obl <- modsem(m_nlinear, data = data, method = "qml")
  qml_nlinear_ort <- modsem(m_nlinear, data = data, method = "qml",
                           orthogonal.x = TRUE, orthogonal.y = TRUE)
})
