devtools::load_all()

model <- '
  X =~ x1 + b * x2 + x3
  Z =~ z1 + b * z2 + x3 # let x3 cross-load
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
  b := normal(0.0, 1.2)
'

parTable <- modsemify(model, parentheses.as.string = TRUE)
parTableEx <- expandModsemParTable(parTable)

cat(buildStanSyntaxFromParTable(parTableEx))


m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

parTable <- modsemify(m1)
parTableEx <- expandModsemParTable(parTable)

library(rstan)
stanCode <- buildStanSyntaxFromParTable(parTableEx)
stanModel <- rstan::stan_model(model_code = stanCode)

stanData <- list(
  N = nrow(oneInt),
  IND__x1 = oneInt$x1,
  IND__x2 = oneInt$x2,
  IND__x3 = oneInt$x3,
  IND__z1 = oneInt$z1,
  IND__z2 = oneInt$z2,
  IND__z3 = oneInt$z3,
  IND__y1 = oneInt$y1,
  IND__y2 = oneInt$y2,
  IND__y3 = oneInt$y3,
)

fit <- rstan::sampling(
  object  = stanModel,
  data    = stan_data,
  chains  = chains,
  iter    = iter,
  warmup  = warmup,
  pars    = compiled_model$info$exclude.pars,
  include = FALSE,
  ...
)
