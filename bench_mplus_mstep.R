## How much does ONE M-step iteration cost Mplus, and does it scale with N?
##
## This measures the thing directly rather than inferring it from a regression
## on N. MUITERATIONS caps the iterations Mplus spends inside each M-step, and
## TECH5 shows it running them all -- the derivative sits at ~2e-13 from
## iteration 6 to 10 and it keeps going -- so the budget is honoured exactly
## rather than being cut short by convergence. Then
##
##   T(mu) = startup + EM * (E-step + mu * M-step-iteration)
##   cost of one M-step iteration = (T(mu.long) - T(mu.short))
##                                  / (EM * (mu.long - mu.short))
##
## Everything that does not depend on `mu` -- start-up, the E-step, all fixed
## overhead -- cancels in the difference.
##
## THE DISCRIMINATION. Run it at two sample sizes:
##   * cost per M-step iteration FLAT in N  => the M-step is running off an
##     N-independent sufficient statistic (Bock-Aitkin expected counts on a
##     common node grid), which is the only way ~5 Newton iterations fit inside
##     a 0.2 s M-step when one E-step traversal costs 0.7 s.
##   * cost per M-step iteration PROPORTIONAL to N => each M-step iteration
##     traverses the data, exactly as ours does, and the gap is elsewhere.
##
## Indicators are mixed continuous/ordered for the same reason as in
## bench_mplus_scaling.R: Mplus compresses duplicate rows, and with six
## categorical indicators the distinct-row count saturates as N grows, so
## "N" would stop meaning what the test needs it to mean. The script asserts
## unique rows == N.
MPLUS <- Sys.getenv("MPLUS_BIN", "/opt/mplusdemo/mpdemo")
NS    <- c(500L, 2000L)
NODES <- 15L
EM.ITERS <- 20L
MU.SHORT <- 1L
MU.LONG  <- 11L
PROCESSORS <- 2L

categorise <- function(x, k = 4L) {
  cuts <- stats::quantile(x, probs = seq_len(k - 1L) / k)
  as.integer(cut(x, breaks = c(-Inf, cuts, Inf)))
}
simulate <- function(n, seed) {
  set.seed(seed)
  X <- stats::rnorm(n); Z <- stats::rnorm(n)
  Y <- 0.5 * X + 0.4 * Z + 0.3 * X * Z + stats::rnorm(n, sd = 0.8)
  data.frame(
    x1 = X + stats::rnorm(n, sd = 0.5),
    x2 = categorise(X + stats::rnorm(n, sd = 0.5)),
    z1 = Z + stats::rnorm(n, sd = 0.5),
    z2 = categorise(Z + stats::rnorm(n, sd = 0.5)),
    y1 = categorise(Y + stats::rnorm(n, sd = 0.5)),
    y2 = categorise(Y + stats::rnorm(n, sd = 0.5))
  )
}

dir <- file.path(tempdir(), "mplus_mstep")
dir.create(dir, showWarnings = FALSE, recursive = TRUE)

mplusInput <- function(file, mu) sprintf('TITLE:
M-step cost
DATA:
FILE = "%s";
VARIABLE:
NAMES = x1 x2 z1 z2 y1 y2;
CATEGORICAL = x2 z2 y1 y2;
ANALYSIS:
estimator = ml;
type = random;
algorithm = integration;
processors = %d;
integration = GAUSSHERMITE(%d);
MITERATIONS = %d;
MUITERATIONS = %d;
CONVERGENCE = 0.1D-11;
MODEL:
X BY x1;
X BY x2;
Z BY z1;
Z BY z2;
Y BY y1;
Y BY y2;
Y ON X;
Y ON Z;
Y ON XZ;
XZ | X XWITH Z;
OUTPUT: TECH8;
', file, PROCESSORS, NODES, EM.ITERS, mu)

elapsed <- function(lines) {
  row <- grep("Elapsed Time:", lines, value = TRUE)
  if (!length(row)) return(NA_real_)
  parts <- as.numeric(strsplit(
    trimws(sub(".*Elapsed Time:", "", row[[1]])), ":")[[1]])
  sum(parts * c(3600, 60, 1))
}
iterations <- function(lines) {
  rows <- grep("^\\s+[0-9]+\\s+-?0?\\.[0-9]+D[+-][0-9]+", lines, value = TRUE)
  iter <- suppressWarnings(vapply(strsplit(trimws(rows), "\\s+"),
                                  function(f) as.integer(f[[1]]), integer(1)))
  iter <- iter[!is.na(iter)]
  if (!length(iter)) NA_integer_ else max(iter)
}
run <- function(n, mu, data.file) {
  stem <- sprintf("mu%d_n%d", mu, n)
  writeLines(mplusInput(data.file, mu), file.path(dir, paste0(stem, ".inp")))
  owd <- setwd(dir); on.exit(setwd(owd), add = TRUE)
  live <- suppressWarnings(system2(MPLUS, paste0(stem, ".inp"),
                                   stdout = TRUE, stderr = TRUE))
  setwd(owd)
  out <- file.path(dir, paste0(stem, ".out"))
  lines <- if (file.exists(out)) readLines(out, warn = FALSE) else live
  list(seconds = elapsed(lines), iterations = iterations(lines))
}

results <- data.frame()
for (n in NS) {
  data <- simulate(n, seed = 4000L + n)
  stopifnot(NROW(unique(data)) == n)
  data.file <- sprintf("mstep%d.dat", n)
  utils::write.table(data, file.path(dir, data.file), row.names = FALSE,
                     col.names = FALSE, quote = FALSE, sep = "\t")

  short <- run(n, MU.SHORT, data.file)
  long  <- run(n, MU.LONG,  data.file)

  # Both runs must have done the same number of EM iterations, or the
  # difference is not isolating the M-step.
  same.em <- isTRUE(short$iterations == long$iterations)
  extra <- short$iterations * (MU.LONG - MU.SHORT)
  per.mstep <- if (same.em && extra > 0)
    (long$seconds - short$seconds) / extra else NA_real_

  cat(sprintf("N = %5d   MU=%2d: %3.0fs (%s EM)   MU=%2d: %3.0fs (%s EM)",
              n, MU.SHORT, short$seconds, short$iterations,
              MU.LONG, long$seconds, long$iterations))
  cat(if (same.em) "" else "   [EM counts differ -- unusable]")
  cat(sprintf("\n            one M-step iteration = %s s\n",
              if (is.na(per.mstep)) "--" else sprintf("%.4f", per.mstep)))

  results <- rbind(results, data.frame(
    n = n, em = short$iterations, short = short$seconds, long = long$seconds,
    per.em.iter = short$seconds / max(short$iterations, 1),
    per.mstep.iter = per.mstep))
}

cat("\n=== Mplus M-step iteration cost ===\n")
print(results, row.names = FALSE)

if (sum(!is.na(results$per.mstep.iter)) >= 2L) {
  ratio.n <- results$n[NROW(results)] / results$n[1]
  ratio.t <- results$per.mstep.iter[NROW(results)] / results$per.mstep.iter[1]
  cat(sprintf("\nN grew %.1fx; cost per M-step iteration grew %.2fx\n",
              ratio.n, ratio.t))
  cat(sprintf("  ~1x  => N-independent sufficient statistic\n"))
  cat(sprintf("  ~%.1fx => every M-step iteration traverses the data\n", ratio.n))
}
cat("\nArtefacts in ", dir, "\n", sep = "")
