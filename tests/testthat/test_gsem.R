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

lms1 <- modsem_gsem(m1, oneInt, nodes = 15, n.threads = 5)
#> Analytical Gradient:
#>   lambda11   lambda12   lambda23   lambda24   lambda35   lambda36 
#> -1019.9519  -682.6786 -1024.4444  -642.4775 -1889.3656 -1729.8974 
#>       tau1       tau2       tau3       tau4       tau5       tau6 
#>  -355.7580 -1248.0169  -655.8259  -342.2559 -1234.7537  -654.8358 
#>       tau7       tau8       tau9     theta1    theta11    theta21 
#>  -536.4231 -1392.3630  -888.4010   474.8811   237.6322   444.8463 
#>    theta31    theta41    theta51    theta61    theta71    theta81 
#>   462.7849   234.5910   461.8497   121.0870  -104.0361   -41.6927 
#>       psi6       psi7      psi11      psi16     gamma4     gamma8 
#>  -613.0140  -817.2402  -603.4588 -1730.6229     0.0000 -2258.2403 
#>    gamma12 
#> -2178.9442 
#> Numerical Gradient:
#>    lambda11    lambda12    lambda23    lambda24    lambda35    lambda36 
#> -1019.95174  -682.67844 -1024.44421  -642.47703 -1889.36534 -1729.89705 
#>        tau1        tau2        tau3        tau4        tau5        tau6 
#>  -355.75773 -1248.01671  -655.82568  -342.25573 -1234.75365  -654.83560 
#>        tau7        tau8        tau9      theta1     theta11     theta21 
#>  -536.42292 -1392.36294  -888.40083   474.88120   237.63219   444.84619 
#>     theta31     theta41     theta51     theta61     theta71     theta81 
#>   462.78547   234.59099   461.84978   121.08692  -104.03601   -41.69277 
#>        psi6        psi7       psi11       psi16      gamma4      gamma8 
#>  -613.01384 -1634.47992  -603.45872 -1730.62275 -4516.94053 -2921.11101 
#>     gamma12 
#> -2871.45158 
#> Difference:
#>   lambda11   lambda12   lambda23   lambda24   lambda35   lambda36 
#>     0.0002     0.0002     0.0001     0.0004     0.0003     0.0004 
#>       tau1       tau2       tau3       tau4       tau5       tau6 
#>     0.0003     0.0002     0.0002     0.0002     0.0001     0.0002 
#>       tau7       tau8       tau9     theta1    theta11    theta21 
#>     0.0002     0.0001     0.0002     0.0001     0.0000    -0.0001 
#>    theta31    theta41    theta51    theta61    theta71    theta81 
#>     0.0005     0.0000     0.0001    -0.0001     0.0001    -0.0001 
#>       psi6       psi7      psi11      psi16     gamma4     gamma8 
#>     0.0002  -817.2397     0.0001     0.0001 -4516.9405  -662.8707 
#>    gamma12 
#>  -692.5074 


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
