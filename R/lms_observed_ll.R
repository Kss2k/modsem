obsLogLikLms <- function(theta, model, data, P, sign = 1, ...) {
  sum(obsLogLikLms_i(theta, model = model, data = data, P = P, sign = sign))
}


obsLogLikLms_i <- function(theta, model, data, P, sign = 1, ...) {
  modFilled <- fillModel(model = model, theta = theta, method = "lms")

  V <- P$V
  w <- P$w
  N <- nrow(data)
  m <- nrow(V)
  px <- numeric(N)

  for (i in seq_len(m)) {
    z_i     <- V[i, ]
    mu_i    <- muLmsCpp(  model = modFilled, z = z_i)
    sigma_i <- sigmaLmsCpp(model = modFilled, z = z_i)
    dens_i  <- dmvn(data, mean = mu_i, sigma = sigma_i, log = FALSE)
    px <- px + w[i] * dens_i
  }
 
  sign * log(px)
}


gradientObsLogLikLms <- function(theta, model, data, P, sign = -1, epsilon = 1e-4) {
  baseLL <- logLikLms(theta, model = model, data = data, P = P, sign = sign)

  vapply(seq_along(theta), FUN.VALUE = numeric(1L), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (obsLogLikLms(theta, model = model, data = data, P = P, sign = sign) - baseLL) / epsilon
  })
}


# gradient function of logLikLms_i
gradientObsLogLikLms_i <- function(theta, model, data, P, sign = -1, epsilon = 1e-4) {
  baseLL <- logLikLms_i(theta, model, data = data, P = P, sign = sign)

  lapplyMatrix(seq_along(theta), FUN = function(i) {
    theta[[i]] <- theta[[i]] + epsilon
    (logLikLms_i(theta, model, data = data, P = P, sign = sign) - baseLL) / epsilon
  }, FUN.VALUE = numeric(nrow(data)))
}
