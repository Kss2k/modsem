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

testthat::expect_equal(trace_path(m2, "x3", "x3"), "(x3~~x3 + x3~x4 ^ 2 * x4~~x4)")
testthat::expect_equal(trace_path(m2, "x4", "x4"), "(x4~~x4)")
