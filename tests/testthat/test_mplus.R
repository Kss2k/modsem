devtools::load_all()
m1 <- '
# Outer Model
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2 

# Inner model
  Y ~ X + Z + X:Z
'
run <- tryCatch({
    MplusAutomation::detectMplus()
    TRUE
  },
  error = function(e) FALSE
)
if (run) {
  mplus <- modsem(m1, oneInt, method = "mplus")
  print(summary(mplus))
  plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z", vals_z = c(-0.5, 0.5), model = mplus)
} 
