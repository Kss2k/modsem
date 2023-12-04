library(mvtnorm)
library(tidyverse)
library(lavaan)
# N <- c(250, 1000)
# sigx1 <- sigx2 <- sigx3 <- 0.4
# sigz1 <- sigz2 <- sigz3 <- 0.4
# sigy1 <- sigy2 <- sigy3 <- 0.4
#
# b0x1 <- 1
# b0x2 <- 1.2
# b0x3 <- 0.9
#
# b0z1 <- 1
# b0z2 <- 1.2
# b0z3 <- 0.9
#
# b0y1 <- 1
# b0y2 <- 1.2
# b0y3 <- 0.9
#
# b1x1 <- 1
# b1x2 <- 0.8
# b1x3 <- 0.9
#
# b1z1 <- 1
# b1z2 <- 0.8
# b1z3 <- 0.9
#
# b1y1 <- 1
# b1y2 <- 0.8
# b1y3 <- 0.9
#
# gammaX <- c(0, 0.2, 0.4, 1)
# gammaZ <- c(0, 0.2, 0.5, 1)
# gammaXZ <- c(0, 0.2, 0.5, 0.9)
#
# zetaY <- c(0.6, 1, 1.5)
#
# corrXZ <- c(0, 0.2, 0.5, 0.7, 0.9)

testModel <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3

  # Inner model
  Y ~ X + Z + X:Z
'


combosParams <- expand.grid(
  N = 1000,
  # Disturbance coefficients
  sigx1 = 0.4,
  sigx2 = 0.4,
  sigx3 = 0.4,
  sigz1 = 0.4,
  sigz2 = 0.4,
  sigz3 = 0.4,
  sigy1 = 0.4,
  sigy2 = 0.4,
  sigy3 = 0.4,
  # Intercepts X
  b0x1 = 1,
  b0x2 = 1.2,
  b0x3 = 0.9,
  # Intercepts Z
  b0z1 = 1,
  b0z2 = 1.2,
  b0z3 = 0.9,
  # Intercepts Y
  b0y1 = 1,
  b0y2 = 1.2,
  b0y3 = 0.9,
  # Loadings X
  b1x1 = 1,
  b1x2 = 0.8,
  b1x3 = 0.9,
  # Loadings Z
  b1z1 = 1,
  b1z2 = 0.8,
  b1z3 = 0.9,
  # Loadings Z
  b1y1 = 1,
  b1y2 = 0.8,
  b1y3 = 0.9,
  # Coefficients Y ~
  gammaX = c(0, 0.2, 0.4, 1),
  gammaZ = c(0, 0.2, 0.5, 1),
  gammaXZ = c(0, 0.2, 0.5, 0.9),
  # Coefficients disturbances Y
  zetaY = c(0.6),
  # Correlation between X and Y
  corrXZ = c(0, 0.2, 0.5, 0.7)
)


createDataParams <- function(paramsRow) {
  sigma <- matrix(c(1, paramsRow$corrXZ,
                    paramsRow$corrXZ, 1), byrow = TRUE, ncol = 2)
  dfXZ <- rmvnorm(paramsRow$N, c(0, 0), sigma) |> as.data.frame()

  # Latent Variable X
  #X <- rnorm(N, 0, 1)
  X <- dfXZ$V1
  x1 <- paramsRow$b1x1 * X + paramsRow$b0x1 + rnorm(paramsRow$N, 0, paramsRow$sigx1)
  x2 <- paramsRow$b1x2 * X + paramsRow$b0x2 + rnorm(paramsRow$N, 0, paramsRow$sigx2)
  x3 <- paramsRow$b1x3 * X + paramsRow$b0x3 + rnorm(paramsRow$N, 0, paramsRow$sigx3)


  # Latent Variable z
  #Z <- rnorm(N, 0, 1)
  Z <- dfXZ$V2
  z1 <- paramsRow$b1z1 * Z + paramsRow$b0z1 + rnorm(paramsRow$N, 0, paramsRow$sigz1)
  z2 <- paramsRow$b1z2 * Z + paramsRow$b0z2 + rnorm(paramsRow$N, 0, paramsRow$sigz2)
  z3 <- paramsRow$b1z3 * Z + paramsRow$b0z3 + rnorm(paramsRow$N, 0, paramsRow$sigz3)



  # Latent Variable Y = X + Z + X:Z
  Y <- paramsRow$gammaX*X + paramsRow$gammaZ*Z + paramsRow$gammaXZ*Z*X + rnorm(paramsRow$N, 0, paramsRow$zetaY)

  #scaledY <- scale(YnoError)*10 + 5

  y1 <- paramsRow$b1y1 * Y + paramsRow$b0y1 + rnorm(paramsRow$N, 0, paramsRow$sigy1)
  y2 <- paramsRow$b1y2 * Y + paramsRow$b0y2 + rnorm(paramsRow$N, 0, paramsRow$sigy2)
  y3 <- paramsRow$b1y3 * Y + paramsRow$b0y3 + rnorm(paramsRow$N, 0, paramsRow$sigy3)

  simData <- data.frame(realY = Y,
                        realX = X,
                        realZ = Z,
                        x1 = x1,
                        x2 = x2,
                        x3 = x3,
                        y1 = y1,
                        y2 = y2,
                        y3 = y3,
                        z1 = z1,
                        z2 = z2,
                        z3 = z3)
  simData
}


getRealParamsRow <- function(paramsRow) {
  coefs <- unlist(paramsRow[c("gammaX", "gammaZ", "gammaXZ")])
  data.frame(rhs = c("X", "Z", "XZ"),
             est = coefs,
             ci.lower = NA,
             ci.upper = NA,
             method = "population")
}


estMethods <- function(model, data, realParams = NULL, corrXZ,
                    methods = c("rca", "uca", "dblcent", "pind"), ...)  {
  est <- lapply(methods,
                FUN = function(method, model, data, ...)
                  modsem(model, data, method = method, ...),
                model = model,
                data = data,
                ... = ...) |>
    lapply(FUN = function(x)
      cbind.data.frame(parameterestimates(x$lavaan), method = attr(x, "method"))) |>
    purrr::list_rbind() |>
    filter(op == "~") |>
    select(c("rhs", "est", "ci.lower", "ci.upper", "method"))

  rbind(est, realParams)  |>
    mutate(corrXZ = corrXZ)
}


simLoop <- function(model, combosParams, ...) {
  estimates <- list()
  for (i in 1:nrow(combosParams)) {
    cat("Iteration", paste0(i, "/", nrow(combosParams)), "\n")
    paramRow <- combosParams[i, ]
    data <- createDataParams(paramRow)
    realParams <- getRealParamsRow(paramRow)
    estimates[[i]] <- estMethods(model, data = data,
                                 realParams = realParams, corrXZ = paramRow[["corrXZ"]], ...)
  }
  estimates
}

estimates <- simLoop(testModel, combosParams, auto.center = "pind",
                     methods = c("rca", "uca", "dblcent", "pind", "ca"))

# labelling simnumber

estDf <-imap(estimates,
                     .f = function(df, idx)
                       mutate(df, simNum = idx)) |>
  list_rbind() |>
  as_tibble()
# calculating bias

#estDf <- mutate(estDf, id = paste(rhs, est, simNum))
estDfWide <- pivot_wider(estDf,
                         names_from = "method",
                         values_from = c("est", "ci.upper", "ci.lower"),
                         values_fill = 0)
# calculating bias
estDfWide[paste0("error_", c(methods, "population"))] <-
  estDfWide[grepl("^est", colnames(estDfWide))] - estDfWide[["est_population"]]

# Converting back into long
estDfLong <- pivot_longer(estDfWide,
                          cols = starts_with(c("est", "ci", "error")),
                          names_to = c(".value", "method"),
                          names_transform = list("method" = "as.factor"),
                          names_pattern = "(ci.upper|ci.lower|est|error)_([a-z]*)")
estDfLong <- mutate(estDfLong, sqrError = error^2) |>
  mutate(sdError = sqrt(sqrError))

## Analysis --------------------------------------------------------------------

  # Bias and sd from population mean for all the gamma coefficients
lapply(estDfWide[grepl("^error", colnames(estDfWide))],
       mean)
lapply(estDfWide[grepl("^sdError", colnames(estDfWide))],
       function(x) sqrt(sum(x^2)/length(x)))

  # Bias and sd from population mean for the interaction terms
lapply(estDfWide[estDfWide$rhs == "XZ", grepl("^error", colnames(estDfWide))],
       mean)
lapply(estDfWide[estDfWide$rhs == "XZ", grepl("^error", colnames(estDfWide))],
       function(x) sqrt(sum(x^2)/length(x)))


  # Trying to predict bias
estDfLong <- mutate(estDfLong, method = relevel(method, "population"))
model <- lm(sdError ~ method, filter(estDfLong, rhs == "XZ"))
coefs <- summary(model)[["coefficients"]] |> as_tibble(rownames = "method")
coefs[c("ci.lower", "ci.upper")] <- confint.lm(model)
ggplot(coefs, aes(x = method, y = Estimate, colour = method,
                      ymin = ci.lower, ymax = ci.upper)) +
  geom_pointrange(position = position_dodge(width = 0.2))

plots <- list()
for (i in sample(1:length(estimates), 20)) {
  plots[[i]] <- ggplot(estimates[[i]], aes(x = rhs, y = est, colour = method,
                             ymin = ci.lower, ymax = ci.upper)) +
    geom_pointrange(position = position_dodge(width = 0.2))
}


#saveRDS(estDfLong, "simulations/256simsTibbleLong031223.rds")



