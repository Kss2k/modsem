
underlyingModel <- '

# Outer Model
X =~ 1*x1 + 0.9*x2 + 0.8*x3
x1 ~~ rvar(0.4)*x1
x2 ~~ rvar(0.4)*x2
x3 ~~ rvar(0.4)*x3
x1 ~ 1
x2 ~ 1.2
x3 ~ 0.7

Z =~ 1*z1 + 0.9*z2 + 0.8*z3
z1 ~~ rvar(0.4)*z1
z2 ~~ rvar(0.4)*z2
z3 ~~ rvar(0.4)*z3
z1 ~ 1
z2 ~ 1.2
z3 ~ 0.7

Y =~ 1*y1 + 0.9*y2 + 0.8*y3
y1 ~~ rvar(0.4)*y1
y2 ~~ rvar(0.4)*y2
y3 ~~ rvar(0.4)*y3
z1 ~ 1
z2 ~ 1.2
z3 ~ 0.7

# Inner model
X ~ rvar(1)
Z ~ 0.5*X + rvar(0.5)
Y ~ X + Z + X:Z + rvar(0.5)
'

specifyUnderlyingModel <- function(syntax, N = 100, data = data.frame(id = 1:N)) {
  modEnv$data <- data
  parTable <- modsemify(syntax)

  predictExprs <- parTable[parTable$op == "~", ]
  measureExprs <- parTable[parTable$op == "=~", ]
  covVarExprs <- parTable[parTable$op == "~~", ]

  # Specify latent variables first
  latentVars <- unique(measureExprs$lhs)

  lapply(latentVars,
         FUN = createVarFromSpec,
         parTable = parTable,
         type = "~")

  # Then specify observed variables

}




createVarFromSpec <- function(latentVar, parTable, type = "~") {
  parTableExprs <- parTable[parTable$op == type & parTable$lhs == latentVar, ]

  parTableToExpr(parTableExprs) |>
    rlang::parse_expr() |>
    eval(env = modVarsEnv)
  updateVariablesEnvir()
}




parTableToExpr <- function(parTable) {
  varName <- unique(parTable$lhs)
  if (length(varName) != 1) {
    stop("Expected partable to contain expressions for creating a single variable \n",
         capturePrint(parTable))
  }

  lhs <- paste0("modEnv$data$", varName)
  rhsVec <- c("")
  for (i in 1:nrow(parTable)) {
    rhsVec[[i]] <- parTable$rhs[[i]]

    if (parTable$mod[[i]] != "") {
      rhsVec[[i]] <- paste0(rhsVec[[i]], "*", parTable$mod[[i]])
    }
  }
  rhs <- stringr::str_c(rhsVec, collapse = " + ") |>
    stringr::str_replace_all(":", "*")

  paste(lhs, "<-", rhs)
}



