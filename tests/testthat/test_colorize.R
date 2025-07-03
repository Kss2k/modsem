devtools::load_all()

m1 <- "
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
# Inner Model
  Y ~ X + Z + X:Z
"

est <- modsem(m1, data = oneInt)
print(summary(est)) # no colors

# Summary with default colors
set_modsem_colors()
print(summary(est))

# Change colors
set_modsem_colors(numeric.positive = "green", numeric.negative = "red", active = TRUE)
print(summary(est))

# Disable colors
set_modsem_colors(active = FALSE)
print(summary(est, standardized = TRUE))
