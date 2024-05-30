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
  if (method == "ca") {
    ests[[method]] <- modsem(m1, data = oneInt, method = method, match = TRUE)
  } else {
    ests[[method]] <- modsem(m1, data = oneInt, method = method)
  }
}

nlsemModel <- '
ENJ =~ enjoy1 + enjoy2 + enjoy3 + enjoy4 + enjoy5
CAREER =~ career1 + career2 + career3 + career4
SC =~ academic1 + academic2 + academic3 + academic4 + academic5 + academic6
CAREER ~ ENJ + SC + ENJ:ENJ + SC:SC + ENJ:SC
'

ests2 <- vector("list", length(methods))
names(ests2) <- methods
for (method in c("rca", "dblcent", "qml")) {
  ests2[[method]] <- modsem(nlsemModel, data = jordan, method = method)
}
