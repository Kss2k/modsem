devtools::load_all()


m1 <- '
  X =~ x1 + x2
  Z =~ z1 + z2
  Y =~ y1 + y2

  Y ~ X + Z + X:Z
'


rthreshold <- \(k, offset = runif(1, min = -1, max = 1), sigma = 0.35) {
  t <- seq_len(k) - mean(seq_len(k)) + offset
  t <- t + runif(k, min = -sigma, max = sigma)
  c(-Inf, t, Inf)
}


cut_data <- function(data, k = 5, choose = NULL) {
  if (is.null(choose))
    choose <- colnames(data)

  standardize <- \(x) (x - mean(x)) / sd(x)

  for (var in choose) {
    x <- standardize(data[[var]])
    t <- rthreshold(k)
    data[[var]] <- cut(x, breaks = t, ordered_result = TRUE)
  }

  data
}


choose <- colnames(oneInt)
choose <- choose[!grepl("3", choose)]

set.seed(2837290)
data <- cut_data(oneInt, choose = choose)
data <- unique(data[choose])

fit.mod <- modsem(
  m1,
  data,
  method = "lms",
  ordered = choose,
  adaptive = "full",
  integration = "gh",
  nodes = 15
)

fit.mp <- modsem( # just for refereance, we should probably run it manually after the 
                  # syntax has been generated
  m1,
  data,
  method = "lms",
  ordered = choose,
  adaptive = "full",
  integration = "gh",
  nodes = 15
)
