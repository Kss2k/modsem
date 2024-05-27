devtools::load_all()
m1 <- '
# Outer Model
X =~ x1
X =~x2
Z =~ z1 + z2
Y =~ y1 + y2 


# Inner model
Y ~ X + Z
Y ~ X:Z
'
run <- tryCatch({
    MplusAutomation::detectMplus()
    TRUE
  },
  error = function(e) FALSE
)
if (run) mplus <- modsem(m1, oneInt, method = "mplus")
