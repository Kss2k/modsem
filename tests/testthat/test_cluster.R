devtools::load_all()

# no interaction
m1 <- '
level: 1
  fw =~ y1 + y2 + y3
  fw ~ x1 + x2 + x3
level: 2
  fb =~ y1 + y2 + y3
  fb ~ w1 + w2
'

fit <- modsem(model = m1, data = lavaan::Demo.twolevel, cluster = "cluster")

# with interaction
m2 <- '
level: 1
    X =~ x1 + x2 + x3
    Z =~ z1 + z2 + z3
    Y =~ y1 + y2 + y3

    Y ~ X + Z + X:Z

level: 2
    Y =~ y1 + y2 + y3
    Y ~ w1 + w2
'


# Generate data accoring to `m2`
set.seed(123)

# Parameters
n_clusters <- 100
n_per_cluster <- 20
n_total <- n_clusters * n_per_cluster

# Cluster-level predictors
w1 <- rnorm(n_clusters)
w2 <- rnorm(n_clusters)

# Level-2 latent Y
Y_b_lvl2 <- 0.7 * w1 + 0.5 * w2 + rnorm(n_clusters, sd = 0.5)

# Create base data frame
data <- data.frame(cluster = rep(1:n_clusters, each = n_per_cluster))

# Repeat cluster-level variables
data$w1 <- rep(w1, each = n_per_cluster)
data$w2 <- rep(w2, each = n_per_cluster)
data$Y_b_lvl2 <- rep(Y_b_lvl2, each = n_per_cluster)

# Level-1 latent variables
X <- rnorm(n_total)
Z <- rnorm(n_total)
XZ <- X * Z
Y_w <- 0.6 * X + 0.6 * Z + 0.5 * XZ + rnorm(n_total, sd = 0.5)
Y_total <- Y_w + data$Y_b_lvl2

# Generate cluster-level residuals for y1–y3 (random intercepts)
y_between_resid <- matrix(rnorm(n_clusters * 3, sd = 0.3), ncol = 3)  # residuals per cluster
colnames(y_between_resid) <- paste0("y_b_", 1:3)

# Expand residuals to full data
y_b_mat <- y_between_resid[data$cluster, ]

# Generate indicators of Y with both between and within residuals
Y_ind <- matrix(NA, n_total, 3)
colnames(Y_ind) <- paste0("y", 1:3)
for (j in 1:3) {
  within_resid <- rnorm(n_total, sd = 0.3)
  Y_ind[, j] <- Y_total + y_b_mat[, j] + within_resid
}

# Indicators for X and Z (no cluster-level component assumed here)
gen_indicators <- function(latent, sd = 0.3) {
  sapply(1:3, function(i) latent + rnorm(length(latent), sd = sd))
}
X_ind <- gen_indicators(X)
Z_ind <- gen_indicators(Z)
colnames(X_ind) <- paste0("x", 1:3)
colnames(Z_ind) <- paste0("z", 1:3)

# Combine all into final dataset
sim_data <- cbind(data, X_ind, Z_ind, Y_ind)

# Check structure
fit <- modsem(m2, data=sim_data, cluster="cluster")
testthat::expect_warning(summary(fit), regexp = "Comparative fit.*")
standardized_estimates(fit)
testthat::expect_error(standardized_estimates(fit, correction = TRUE),
                       regexp = "Correction of clustered .*not supported!")

# Test robust std.errors
mod <- '
# X1-3 are Level 1 variables
X1 =~ x1 
X2 =~ x2
X3 =~ x3

# W1-2 are Level 2 variables
W1 =~ w1
W2 =~ w2

fw =~ y1 + y2 + y3
fw ~ X1 + X2 + X3 + W1 + W2
'

# Standard errors corrected for clustering
fit.rc <- modsem(mod, lavaan::Demo.twolevel, method = "lms", cluster = "cluster", robust.se = TRUE)
#> Regressions:
#>                  Estimate  Std.Error  z.value  P(>|z|)
#>   fw ~          
#>     X1              0.517      0.027   19.224    0.000
#>     X2              0.384      0.029   13.388    0.000
#>     X3              0.174      0.026    6.677    0.000
#>     W1              0.207      0.089    2.331    0.020
#>     W2              0.158      0.075    2.111    0.035

# Standard errors not corrected for clustering
fit.r <- modsem(mod, lavaan::Demo.twolevel, method = "lms", robust.se = TRUE)
#> Regressions:
#>                  Estimate  Std.Error  z.value  P(>|z|)
#>   fw ~          
#>     X1              0.517      0.029   18.068    0.000
#>     X2              0.384      0.028   13.520    0.000
#>     X3              0.174      0.028    6.286    0.000
#>     W1              0.207      0.031    6.737    0.000
#>     W2              0.158      0.029    5.455    0.000

fit.c <- modsem(mod, lavaan::Demo.twolevel, method = "lms", cluster = "cluster", robust.se = FALSE)
