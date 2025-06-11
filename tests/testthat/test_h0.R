devtools::load_all()
m1 <- "
 # Outer Model
 X =~ x1 + x2 + x3
 Y =~ y1 + y2 + y3
 Z =~ z1 + z2 + z3

 # Inner model
 Y ~ X + Z + X:Z
"


# modsem_da
est_h1 <- modsem(m1, oneInt, "lms")
est_h0 <- estimate_h0(est_h1, calc.se=FALSE) # std.errors are not needed
compare_fit(est_h1 = est_h1, est_h0 = est_h0)

modsem_inspect(est_h1, what="fit")


# modsem_pi
est_h1 <- modsem(m1, oneInt, method = "dblcent")
est_h0 <- estimate_h0(est_h1, oneInt)

compare_fit(est_h1 = est_h1, est_h0 = est_h0)

est_h1 <- modsem(m1, oneInt, method = "ca")
est_h0 <- estimate_h0(est_h1, oneInt)

compare_fit(est_h1 = est_h1, est_h0 = est_h0)
