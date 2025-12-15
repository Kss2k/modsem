devtools::load_all()

n <- 1000
x <- rnorm(n)
z <- rnorm(n)
y <- -1 * x + 1 * z + rnorm(n)


fit <- modsem('y ~ x + z', data = data.frame(x, z, y))
plot_surface(x = "x", z = "z", y = "y", model = fit)

x <- y <- c(1, 2, 3)
z <- t(outer(x, y, FUN = \(x, y) x))
plot_ly(x = ~x, z=~z, y=~y, type = "surface")

