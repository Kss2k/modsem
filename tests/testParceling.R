test <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y20:mean(y2, y3, y4)
     dem65 =~ y5 + y19:mean(y6, y7, y8)
     dem95 =~ y17 + y18:sum(y11, y12, y13)

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
testIris <- '
mean(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
mean(Sepal.Length, Sepal.Width, Petal.Length)
mean(Sepal.Length, Sepal.Width)

sum(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
sum(Sepal.Length, Sepal.Width, Petal.Length)
sum(Sepal.Length, Sepal.Width)

'

testIrisMean <- '
mean(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
mean(Sepal.Length, Sepal.Width, Petal.Length)
mean(Sepal.Length, Sepal.Width)


'

testIrisSum <- '


sum(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
sum(Sepal.Length, Sepal.Width, Petal.Length)
sum(Sepal.Length, Sepal.Width)

'
computeParcels(testIris, iris)


cat(fixParcelsSyntax(test, getParcelInfo(test)))
cat(fixParcelsSyntax(testIris, getParcelInfo(testIris)))
