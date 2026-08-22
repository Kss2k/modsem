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
