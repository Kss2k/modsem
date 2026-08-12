## Does Mplus's M-step cost scale with N?
##
## TECH5 shows the measurement M-step converging quadratically -- exact-Hessian
## Newton -- in ~5 iterations, and Mplus runs it to 10. If each of those
## iterations traversed the data, the M-step alone would cost ~10 E-steps. It
## does not, so those Newton iterations must run off a sufficient statistic:
## with a COMMON node grid the measurement objective collapses to
##     Q = sum_l sum_q sum_c  n_lqc * log p(c | m_lq),
##     n_lqc = sum_{i: y_il = c} P_iq
## a Q x C table per indicator, formed once in the E-step. Every Newton
## iteration is then O(Q*C*L) -- INDEPENDENT OF N. (Textbook Bock-Aitkin.)
##
## WHAT THIS SCRIPT CAN AND CANNOT ANSWER. It was written to read the M-step
## off the intercept of time(N) = intercept + slope * N. That does NOT work,
## and the reasoning is recorded here so nobody repeats it: BOTH hypotheses
## predict a near-zero intercept. If each M-step iteration traverses the data
## its cost is O(N*Q) and folds into the SLOPE; if it runs off the sufficient
## statistic its cost is ~1 ms here (about 2100 nodes x 4 categories x 6
## indicators x 10 Newton iterations), invisible against an E-step of 0.35-5.6
## s. Only varying the number of M-step iterations separates them -- that is
## bench_mplus_mstep.R.
##
## What this run does establish, and what it is still worth keeping for:
## whether either engine is super-linear in N, and the constant factor between
## them. Measured: log-log slopes 1.004 (Mplus) and 1.038 (modsem), R^2 0.9999
## and 0.9985, so both are proportional to N with no hidden super-linearity,
## and modsem costs 2.726e-03 s/row against Mplus's 1.405e-03 -- a flat 1.94x.
##
## WHY THE INDICATORS ARE MIXED. Mplus compresses duplicate rows (so do we, as
## of compress.data). With six categorical indicators the row space is finite,
## so raising N stops raising the number of DISTINCT rows and the x-axis of
## this regression quietly saturates. Giving X and Z one continuous indicator
## each makes every row unique with probability one, so unique rows == N and
## compression cannot touch it. The script asserts that.
##
## WHY EACH POINT IS RUN TWICE. Mplus reports only `Elapsed Time`, at ONE
## SECOND resolution, and its TECH8 table here carries no per-iteration column,
## so a single short run cannot be timed. Each configuration is run at two
## iteration budgets and differenced: that removes start-up and makes the
## measured quantity large enough that second-granularity rounding is a small
## relative error. Iteration counts are read back from TECH8 rather than
## assumed, because EM may converge before the budget and silently shorten a
## run -- which would corrupt the difference.
suppressMessages(devtools::load_all(".", quiet = TRUE, compile = FALSE))

MPLUS  <- Sys.getenv("MPLUS_BIN", "/opt/mplusdemo/mpdemo")
NS     <- c(250L, 500L, 1000L, 2000L, 4000L)
NODES  <- 15L          # per dimension; 3 latents
ITERS.SHORT <- 15L
ITERS.LONG  <- 55L
MODSEM.ITERS <- 10L
PROCESSORS <- 2L

## ---- data --------------------------------------------------------------
## x1, z1 continuous; x2, z2, y1, y2 ordered. Y ~ X + Z + X:Z.
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

syntax <- 'X =~ x1 + x2
           Z =~ z1 + z2
           Y =~ y1 + y2
           Y ~ X + Z + X:Z'
ordered.vars <- c("x2", "z2", "y1", "y2")

dir <- file.path(tempdir(), "mplus_scaling")
dir.create(dir, showWarnings = FALSE, recursive = TRUE)

mplusInput <- function(file, iterations) sprintf('TITLE:
N scaling
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
', file, PROCESSORS, NODES, iterations)

mplusElapsed <- function(lines) {
  row <- grep("Elapsed Time:", lines, value = TRUE)
  if (!length(row)) return(NA_real_)
  parts <- as.numeric(strsplit(
    trimws(sub(".*Elapsed Time:", "", row[[1]])), ":")[[1]])
  sum(parts * c(3600, 60, 1))
}

## EM iterations actually performed, counted from the TECH8 table.
mplusIterations <- function(lines) {
  rows <- grep("^\\s+[0-9]+\\s+-?0?\\.[0-9]+D[+-][0-9]+", lines, value = TRUE)
  iter <- suppressWarnings(vapply(strsplit(trimws(rows), "\\s+"),
                                  function(f) as.integer(f[[1]]), integer(1)))
  iter <- iter[!is.na(iter)]
  if (!length(iter)) NA_integer_ else max(iter)
}

runMplus <- function(n, iterations, data.file) {
  stem <- sprintf("scale%d_%d", n, iterations)
  writeLines(mplusInput(data.file, iterations),
             file.path(dir, paste0(stem, ".inp")))
  owd <- setwd(dir)
  on.exit(setwd(owd), add = TRUE)
  live <- suppressWarnings(system2(MPLUS, paste0(stem, ".inp"),
                                   stdout = TRUE, stderr = TRUE))
  setwd(owd)
  out <- file.path(dir, paste0(stem, ".out"))
  lines <- if (file.exists(out)) readLines(out, warn = FALSE) else live
  list(seconds = mplusElapsed(lines), iterations = mplusIterations(lines))
}

haveMplus <- file.exists(MPLUS) || nzchar(Sys.which(MPLUS))
results <- data.frame()

for (n in NS) {
  data <- simulate(n, seed = 4000L + n)

  ## the whole point: no two rows may coincide, or the x-axis saturates
  unique.rows <- NROW(unique(data))
  stopifnot(unique.rows == n)

  data.file <- sprintf("scale%d.dat", n)
  utils::write.table(data, file.path(dir, data.file), row.names = FALSE,
                     col.names = FALSE, quote = FALSE, sep = "\t")

  mplus <- NA_real_
  if (haveMplus) {
    short <- runMplus(n, ITERS.SHORT, data.file)
    long  <- runMplus(n, ITERS.LONG,  data.file)
    gap <- long$iterations - short$iterations
    cat(sprintf("  N = %5d  mplus %3.0fs/%s iters -> %3.0fs/%s iters",
                n, short$seconds, short$iterations, long$seconds,
                long$iterations))
    if (!is.na(gap) && gap > 0) {
      mplus <- (long$seconds - short$seconds) / gap
      cat(sprintf("   = %.4f s/iter\n", mplus))
    } else cat("   [converged early -- unusable]\n")
  }

  t0 <- proc.time()[["elapsed"]]
  invisible(suppressWarnings(suppressMessages(modsem(
    syntax, data, method = "lms", lms.backend = "graph", ordered = ordered.vars,
    adaptive = "full", integration = "gh", nodes = NODES,
    algorithm = "EM", max.iter = 2L, calc.se = FALSE, verbose = FALSE))))
  warm <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  invisible(suppressWarnings(suppressMessages(modsem(
    syntax, data, method = "lms", lms.backend = "graph", ordered = ordered.vars,
    adaptive = "full", integration = "gh", nodes = NODES,
    algorithm = "EM", max.iter = MODSEM.ITERS, calc.se = FALSE,
    verbose = FALSE))))
  full <- proc.time()[["elapsed"]] - t0
  modsem.per.iter <- (full - warm) / (MODSEM.ITERS - 2L)
  cat(sprintf("  N = %5d  modsem %.4f s/iter\n", n, modsem.per.iter))

  results <- rbind(results, data.frame(n = n, unique = unique.rows,
                                       mplus = mplus, modsem = modsem.per.iter))
}

cat("\n=== seconds per EM iteration ===\n")
print(results, row.names = FALSE)

report <- function(label, y) {
  if (all(is.na(y)) || sum(!is.na(y)) < 3L) return(invisible(NULL))
  keep <- !is.na(y)
  fit <- stats::lm(y[keep] ~ results$n[keep])
  intercept <- stats::coef(fit)[[1]]
  slope <- stats::coef(fit)[[2]]
  cat(sprintf("\n%s\n", label))
  cat(sprintf("  time(N) = %.4f + %.3e * N      R^2 = %.4f\n",
              intercept, slope, summary(fit)$r.squared))
  cat(sprintf("  N-independent share at N=2000: %.1f%%\n",
              100 * intercept / (intercept + slope * 2000)))
  cat(sprintf("  log-log slope: %.3f  (1.0 = fully proportional to N)\n",
              stats::coef(stats::lm(log(y[keep]) ~ log(results$n[keep])))[[2]]))
}
report("Mplus", results$mplus)
report("modsem", results$modsem)

cat("\nNOTE: the intercept CANNOT settle the sufficient-statistic question --\n")
cat("see the header. Both hypotheses predict a near-zero intercept. What this\n")
cat("run does establish is whether either engine is super-linear in N, and the\n")
cat("constant factor between them. Use bench_mplus_mstep.R for the M-step.\n")
cat("Artefacts in ", dir, "\n", sep = "")
