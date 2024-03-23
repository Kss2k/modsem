models <- list(m1 = ' 
               # latent variables 
               ind60 =~ x1 + x2 + x3 
               dem60 =~ y1 + y2 + y3 + y4 
               dem65 =~ y5 + y6 + y7 + y8 
               ind60:dem60 =~ x1:y1 + x2:y2 + x3:mean(y3, y4)
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
               m3 = ' visual  =~ x1 + x2 + x3 
               textual =~ x4 + x5 + x6
               speed   =~ x7 + x8 + x9
               textual:speed =~ x4:x7 + x5:x8 + x6:x9
               visual ~ speed + textual + textual:speed',
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
  d2[[i]] <- stats::rnorm(500)
}
d2 <- as.data.frame(d2)

set.seed(1234)
X <- rnorm(100)
M <- 0.5*X + rnorm(100)
Y <- 0.7*M + rnorm(100)
d7 <- data.frame(X = X, Y = Y, M = M)

data <- list(d1 = lavaan::PoliticalDemocracy,
             d3 = lavaan::HolzingerSwineford1939, 
             d4 = lavaan::HolzingerSwineford1939, 
             d5 = lavaan::HolzingerSwineford1939, 
             d6 = lavaan::HolzingerSwineford1939, 
             d7 = d7)

nativeMethods <- allNativeMethods[!allNativeMethods %in%  c("pind",
                                                           "ca",
                                                           "uca")]
methods <- list(m1 = nativeMethods,
                m2 = nativeMethods,
                m3 = nativeMethods,
                m4 = nativeMethods,
                m5 = nativeMethods,
                m6 = nativeMethods)


estimates <- vector("list", length(models))
for (i in seq_along(estimates)) {
    runMultipleMethods(models[[i]], data = data[[i]], 
                              methods = methods[[i]],
                              estimator = "ML",
    removeFromParTable = "")
}
