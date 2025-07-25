devtools::load_all()

m1 <- modsemify(
                '
y1 ~ x1
x1 ~ x2
x2 ~ x3
x3 ~ x1
                '
)
testthat::expect_warning(trace_path(m1, "x1", "x1"),
                         "non-recursive model \\(infinite loop\\)")
m2 <- modsemify(
                '
y1 ~ x1
x1 ~ x2
x2 ~ x3
x3 ~ x4
x1 ~~ x1
x2 ~~ x2
x3 ~~ x3
x4 ~~ x4
                '
)

testthat::expect_equal(
  unname(trace_path(m2, "x3", "x3")),
  "(x3~~x3 + x3~x4 ^ 2 * x4~~x4)"
)

testthat::expect_equal(
  unname(trace_path(m2, "x4", "x4")),
  "(x4~~x4)"
)

m3 <- modsemify('visual  =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed   =~ x7 + x8 + x9
                visual  ~ speed + textual + speed:textual')

testthat::expect_true(is.na(trace_path(m3, "textual", "textual")))


m4 <- modsemify('
   # Outer Model
   X =~ x1 + x2 +x3
   Y =~ y1 + y2 + y3
   Z =~ z1 + z2 + z3

   # Inner model
   Y ~ X + Z + X:Z
')

testthat::expect_equal(
  unname(trace_path(m4, "X", "X", missing.cov = TRUE)),
  "(X~~X)"
)
