devtools::load_all("..")
models <- list(m1 = ' 
               # latent variables 
               ind60 =~ x1 + x2 + x3 
               dem60 =~ y1 + y2 + y3 + y4 
               dem65 =~ y5 + y6 + y7 + y8 
               # regressions
               dem60 ~ ind60 
               dem65 ~ ind60 + dem60 
               dem65 ~ ind60:dem60
               # residual covariances 
               y1 ~~ y5
               y2 ~~ y4 + y6 
               y3 ~~ y7 
               y4 ~~ y8
               y6 ~~ y8
               ',
               m2 = myModel <- ' # regressions
               y1 + y2 ~ f1 + f2 + x1 + x2
               f1 ~ f2 + f3 + f2:f3
               f2 ~ f3 + x1 + x2

               # latent variable definitions 
               f1 =~ y1 + y2 + y3 
               f2 =~ y4 + y5 + y6 
               f3 =~ y7 + y8 + y9 + y10

               # variances and covariances 
               y1 ~~ y1 
               y1 ~~ y2 
               f1 ~~ f2

               # intercepts 
               y1 ~ 1 
               f1 ~ 1
               ',
               m3 = ' visual  =~ x1 + x2 + x3 
               textual =~ x4 + x5 + x6
               speed   =~ x7 + x8 + x9
               visual ~ speed + textual + speed:textual',
               m4 = 'visual  =~ x1 + start(0.8)*x2 + start(1.2)*x3
               textual =~ x4 + start(0.5)*x5 + start(1.0)*x6
               speed   =~ x7 + start(0.7)*x8 + start(1.8)*x9
               visual ~ speed + textual + speed:textual',
               m5 = '# three-factor model
               visual  =~ x1 + x2 + x3
               textual =~ x4 + x5 + x6
               speed   =~ x7 + x8 + x9
               visual ~ speed + textual + 0.1*speed:textual
               # intercepts with fixed values
               x1 + x2 + x3 + x4 ~ 0.5*1',
               m6 = '# three-factor model
               visual  =~ x1 + x2 + x3
               textual =~ x4 + x5 + x6
               speed   =~ x7 + x8 + x9
               visual ~ speed + textual + speed:textual
               # intercepts
               x1 ~ 1
               x2 ~ 1
               x3 ~ 1
               x4 ~ 1
               x5 ~ 1
               x6 ~ 1
               x7 ~ 1
               x8 ~ 1
               x9 ~ 1', 
               m7 = ' # direct effect
               Y ~ c*X
               # mediator
               M ~ a*X
               Y ~ b*M + M:X
               # indirect effect (a*b)
               ab := a*b
               # total effect
               total := c + (a*b)
               '
)


d2 <- vector("list", 12L) 
varNames <- names(d2) <- c("x1", "x2", paste0("y", 1:10))
for (i in varNames) {
  d2[[i]] <- rnorm(500)
}
d2 <- as.data.frame(d2)

set.seed(1234)
X <- rnorm(100)
M <- 0.5*X + rnorm(100)
Y <- 0.7*M + rnorm(100)
d7 <- data.frame(X = X, Y = Y, M = M)

data <- list(d1 = lavaan::PoliticalDemocracy,
             d2 = d2,
             d3 = lavaan::HolzingerSwineford1939, 
             d4 = lavaan::HolzingerSwineford1939, 
             d5 = lavaan::HolzingerSwineford1939, 
             d6 = lavaan::HolzingerSwineford1939, 
             d7 = d7)


estimates <- vector("list", length(models))
for (i in seq_along(estimates)) {
  estimates[[i]]$modsem <- tryCatch({
    est <- modsem(models[[i]], data = data[[i]], estimator = "ML")
    est
  },
  warning = function(e) {
    warning("Warning in modsem model ", i)
    est
  },
  error = function(e) {
    est <- NA
    warning("Error in modsem model ", i) 
  },
  finally = {
    est
  })

}

modsem(models[[3]], data[[3]], estimator = "ML")
