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
  if (method == "lms") {
    ests[[method]] <- modsem(m1, data = oneInt, method = method, 
                             convergence.abs = 1e-3)
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

est_qml2 <- modsem(nlsemModel, data = jordan, method = "qml", adaptive.quad=TRUE,
                   nodes = 15, mean.observed = FALSE, convergence.rel = 1e-2)
est_rca2 <- modsem(nlsemModel, data = jordan, method = "rca")
est_dblcent2 <- modsem(nlsemModel, data = jordan, method = "dblcent")


std <- standardized_estimates(est_qml2)
covENJENJ_ENJSC <- std[std$lhs == "ENJ:SC" & std$rhs == "ENJ:ENJ" & std$op == "~~", "est"]
testthat::expect_true(covENJENJ_ENJSC > 1) # check that covariances between product terms 
                                              # are handled correctly
testthat::expect_equal(unname(calcVarParTable("CAREER", std)), 1)
