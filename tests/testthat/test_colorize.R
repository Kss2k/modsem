devtools::load_all()

set_modsem_colors(positive = "green", 
                  negative = "red",
                  true = "green", 
                  false = "red", 
                  na = "purple")

m1 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
# Inner Model
  Y ~ X + Z + X:Z
"

est <- modsem(m1, data = oneInt)
colorize_output(summary(est))
colorize_output(print(est))

est_lms <- modsem(m1, data = oneInt, method = "lms")
colorize_output(summary(est_lms))

colorize_output(modsem_inspect(est_lms))
