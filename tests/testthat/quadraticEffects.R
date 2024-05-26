devtools::load_all()
m1 <- '
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Z =~ z1 + z2 + z3
Y ~ X + Z + X:X
'
methods <- c("rca", "ca", "dblcent", "lms", "qml")
ests <- vector("list", length(methods))
names(ests) <- methods

for (method in methods) {
  ests[[method]] <- modsem(m1, data = oneInt, 
                           method = method, run = TRUE)
}
