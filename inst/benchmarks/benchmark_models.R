#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    output = "benchmark-results.csv",
    iterations = as.integer(Sys.getenv("MODSEM_BENCH_ITERS", "1"))
  )

  for (a in args) {
    if (grepl("^--output=", a)) {
      opts$output <- sub("^--output=", "", a)
    } else if (grepl("^--iters=", a)) {
      opts$iterations <- suppressWarnings(as.integer(sub("^--iters=", "", a)))
    } else {
      stop("Unknown argument: ", a)
    }
  }

  if (!is.finite(opts$iterations) || opts$iterations < 1L)
    opts$iterations <- 1L

  opts
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("modsem", quietly = TRUE))
    stop("Package 'modsem' is not installed.")
  library(modsem)
})

data("oneInt", package = "modsem")
data("TPB", package = "modsem")

oneInt_model <- "
# Outer Model
  X =~ x1 + lx2 * x2 + lx3 * x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner Model
  Y ~ X + Z

# Constraints
  x1 ~~ vx1*x1
  x2 ~~ vx2*x2
  x3 ~~ vx3*x3

  proj_vx1 := 0.8
  proj_vx2 := (lx2 ^ 2) * 0.8
  proj_vx3 := (lx3 ^ 2) * 0.8

  vx1 == 1 - proj_vx1
  vx2 == 1 - proj_vx2
  vx3 == 1 - proj_vx3
"

tpb_model <- "
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN  =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ PBC:INT
"

benchmarks <- list(
  list(
    name = "oneInt_lms",
    run = function() modsem(
      model.syntax     = oneInt_model,
      data             = oneInt,
      method           = "lms",
      standardize.data = TRUE,
      mean.observed    = FALSE,
      calc.se          = FALSE,
      verbose          = FALSE
    )
  ),
  list(
    name = "oneInt_qml",
    run = function() modsem(
      model.syntax     = oneInt_model,
      data             = oneInt,
      method           = "qml",
      standardize.data = TRUE,
      mean.observed    = FALSE,
      calc.se          = FALSE,
      verbose          = FALSE
    )
  ),
  list(
    name = "TPB_lms",
    run = function() modsem(
      model.syntax = tpb_model,
      data         = TPB,
      method       = "lms",
      calc.se      = FALSE,
      verbose      = FALSE,
      nodes        = 32
    )
  ),
  list(
    name = "TPB_qml",
    run = function() modsem(
      model.syntax = tpb_model,
      data         = TPB,
      method       = "qml",
      calc.se      = FALSE,
      verbose      = FALSE
    )
  )
)

all_timings <- lapply(benchmarks, function(entry) {
  vals <- numeric(args$iterations)
  for (i in seq_len(args$iterations)) {
    cat(sprintf("Running %s (iteration %d/%d)...\n",
                entry$name, i, args$iterations))
    vals[i] <- system.time({
      fit <- entry$run()
      rm(fit)
    })[["elapsed"]]
    invisible(gc())
  }
  data.frame(
    model = entry$name,
    iteration = seq_len(args$iterations),
    elapsed = vals
  )
})

timings <- do.call(rbind, all_timings)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(timings, args$output, row.names = FALSE)

print(timings)

summary_df <- aggregate(elapsed ~ model, data = timings,
  FUN = function(x) c(mean = mean(x), median = stats::median(x), sd = stats::sd(x)))

cat("\nSummary (seconds):\n")
print(summary_df)
