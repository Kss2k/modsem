devtools::load_all()
pkgbuild::compile_dll(deb=F, force=TRUE)

m1 <- '
# Outer Model
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z + X:Z
'

lms1 <- modsem_gsem(m1, oneInt, nodes = 15)


rthreshold <- \(k, offset = runif(1, min = -1, max = 1), sigma = 0.35) {
  t <- seq_len(k) - mean(seq_len(k)) + offset
  t <- t + runif(k, min = -sigma, max = sigma)
  c(-Inf, t, Inf)
}


cut_data <- function(data, k = 5, choose = NULL) {
  if (is.null(choose))
    choose <- colnames(data)

  standardize <- \(x) (x - mean(x)) / sd(x)

  thresholds <- list()
  for (var in choose) {
    x <- standardize(data[[var]])
    t <- rthreshold(k)
    y <- cut(x, breaks = t, ordered_result = TRUE)
    y <- as.integer(as.ordered(as.integer(y)))
    min.x <- min(x)
    max.x <- max(x)

    data[[var]]       <- y
    thresholds[[var]] <- t[t >= min.x & t <= max.x]
  }

  list(data = data, thresholds = thresholds)
}


choose <- colnames(oneInt)
set.seed(2837290)
CUTS <- cut_data(oneInt, choose = choose)
oneInt2 <- CUTS$data
lms1 <- modsem_gsem(m1, oneInt2, ordered = choose)
