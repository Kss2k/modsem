devtools::load_all()
pkgbuild::compile_dll(deb=F, force=TRUE)

m1 <- '
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z + X:Z
'

lms1 <- modsem_gsem(m1, oneInt, nodes = 15)
