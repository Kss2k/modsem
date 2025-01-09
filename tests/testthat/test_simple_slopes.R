devtools::load_all()

m1 <- "
# Outer Model
  X =~ x1
  X =~ x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner model
  Y ~ X + Z + X:Z
"

est1 <- modsem(m1, data = oneInt)
simple_slopes(x = "X", z = "Z", y = "Y", model = est1)
plot_interaction(x = "X", z = "Z", y = "Y",
                 vals_z = c(1, 0), model = est1)
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z",
                 vals_z = c(1, 0), model = est1, ci_type = "prediction")
