devtools::load_all()
library(lavaan)

m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + X:Z
'

est <- sam(m1, oneInt)
parameter_estimates(est)

centered_estimates(est)
standardized_estimates(est)

plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z",
                 vals_z = c(-0.5, 0.5), model = est)
plot_jn(x = "X", z = "Z", y = "Y", model = est)
plot_surface(x = "X", z = "Z", y = "Y", model = est)


m2 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + X:Z

X ~ 1
Z ~ 1
Y ~ 1
x1~0*1
z1~0*1
y1~0*1
'

est_nc <- sam(m2, oneInt)
wrap <- \(expr) (testthat::expect_warning(expr, regex = "Replacing.*"))
wrap(centered_estimates(est_nc))
wrap(standardized_estimates(est_nc))
