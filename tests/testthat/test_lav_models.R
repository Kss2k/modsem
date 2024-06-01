set.seed(123)
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
               x9 ~ 1'
               
)



data <- list(d1 = lavaan::PoliticalDemocracy,
             d3 = lavaan::HolzingerSwineford1939, 
             d4 = lavaan::HolzingerSwineford1939, 
             d5 = lavaan::HolzingerSwineford1939, 
             d6 = lavaan::HolzingerSwineford1939)

nativeMethods <- allNativeMethods[allNativeMethods != "pind"]
methods <- list(m1 = nativeMethods[nativeMethods != "ca"],
                m3 = nativeMethods,
                m4 = nativeMethods,
                m5 = nativeMethods,
                m6 = nativeMethods[nativeMethods != c("uca", "ca")])


estimates <- vector("list", length(models))
for (i in seq_along(estimates)) {
    estimates[[i]] <- runMultipleMethods(models[[i]], data = data[[i]], 
                                         methods = methods[[i]],
                                         estimator = "ML")

}

# testing plot function 
plot_interaction(x = "ind60", z = "dem60", y = "dem65", xz = "ind60:dem60", 
                 vals_z = c(-0.5, 0.5), model = estimates[[1]][["rca"]])
plot_interaction(x = "speed", z = "textual", y = "visual", xz = "speed:textual", 
                 vals_z = c(-0.5, 0.5), model = estimates[[2]][["ca"]])
