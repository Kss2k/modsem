library(tidyverse)
library(mvtnorm)
#library(modsem)
library(lavaan)

# TEst data1
N <- 1000
sigx1 <- sigx2 <- sigx3 <- 0.4
sigz1 <- sigz2 <- sigz3 <- 0.4
sigy1 <- sigy2 <- sigy3 <- 0.4

b0x1 <- 1
b0x2 <- 1.2
b0x3 <- 0.9

b0z1 <- 1
b0z2 <- 1.2
b0z3 <- 0.9

b0y1 <- 1
b0y2 <- 1.2
b0y3 <- 0.9

b1x1 <- 1
b1x2 <- 0.8
b1x3 <- 0.9

b1z1 <- 1
b1z2 <- 0.8
b1z3 <- 0.9

b1y1 <- 1
b1y2 <- 0.8
b1y3 <- 0.9

gammaX <- 0.2
gammaZ <- .5
gammaXZ <- .9

zetaY <- 0.6

corrXZ <- 0.7

# Generating X and Y from a normal distribution
sigma <- matrix(c(1, corrXZ, corrXZ, 1), byrow = TRUE, ncol = 2)
dfXZ <- rmvnorm(N, c(0, 0), sigma) |> as.data.frame()

# Latent Variable X
#X <- rnorm(N, 0, 1)
X <- dfXZ$V1
x1 <- b1x1 * X + b0x1 + rnorm(N, 0, sigx1)
x2 <- b1x2 * X + b0x2 + rnorm(N, 0, sigx2)
x3 <- b1x3 * X + b0x3 + rnorm(N, 0, sigx3)


# Latent Variable z
#Z <- rnorm(N, 0, 1)
Z <- dfXZ$V2
z1 <- b1z1 * Z + b0z1 + rnorm(N, 0, sigz1)
z2 <- b1z2 * Z + b0z2 + rnorm(N, 0, sigz2)
z3 <- b1z3 * Z + b0z3 + rnorm(N, 0, sigz3)



# Latent Variable Y = X + Z + X:Z
Y <- gammaX*X + gammaZ*Z + gammaXZ*Z*X + rnorm(N, 0, zetaY)

#scaledY <- scale(YnoError)*10 + 5

y1 <- b1y1 * Y + b0y1 + rnorm(N, 0, sigy1)
y2 <- b1y2 * Y + b0y2 + rnorm(N, 0, sigy2)
y3 <- b1y3 * Y + b0y3 + rnorm(N, 0, sigy3)

simTest <- data.frame(realY = Y,
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
#saveRDS(oneInt, "simTest.rds")
#### Model 1 ####

testModel <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3

  # Inner model
  Y ~ X + Z + X:Z
'

realParams <- data.frame(rhs = c("X", "Z", "XZ"),
                         est = c(gammaX, gammaZ, gammaXZ),
                         ci.lower = NA,
                         ci.upper = NA,
                         method = "population")
methods = c("rca", "uca", "dblcent", "pind", "ca")

simLoop <- function(model, data, realParams = NULL,
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
  rbind(est, realParams)
}

estimates <- simLoop(testModel, simTest, realParams, methods = methods, auto.scale = "pind")
ggplot(estimates, aes(x = rhs, y = est, colour = method,
                      ymin = ci.lower, ymax = ci.upper)) +
  geom_pointrange(position = position_dodge(width = 0.2))


estRca <- modsem(testModel, simTest)
