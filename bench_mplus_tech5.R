## Default-settings Mplus run on the all-ordered benchmark data, TECH5 + TECH8,
## under both EM and EMA. Nothing is tuned: no MITERATIONS, no MUITERATIONS,
## no CONVERGENCE override -- the point is to see what Mplus does on its own,
## early and at the optimum.
suppressMessages(devtools::load_all(".", quiet = TRUE, compile = FALSE))
MPLUS <- Sys.getenv("MPLUS_BIN", "/opt/mplusdemo/mpdemo")

rthreshold <- function(k, offset = runif(1, -1, 1), sigma = 0.35) {
  t <- seq_len(k) - mean(seq_len(k)) + offset
  c(-Inf, t + runif(k, -sigma, sigma), Inf)
}
cut_data <- function(d, k = 5, choose) {
  z <- function(x) (x - mean(x)) / stats::sd(x)
  for (v in choose) d[[v]] <- cut(z(d[[v]]), breaks = rthreshold(k),
                                  ordered_result = TRUE)
  d
}
vars <- colnames(oneInt)[!grepl("3", colnames(oneInt))]
set.seed(2837290)
dat <- cut_data(oneInt, choose = vars)[vars]

## NOT tempdir(): it is wiped when the R session exits, which destroyed the
## EMA output of an 11-minute run before it could be read.
dir <- file.path(normalizePath("."), "mplus_tech5")
dir.create(dir, showWarnings = FALSE)
codes <- as.data.frame(lapply(dat, as.integer))
utils::write.table(codes, file.path(dir, "ord.dat"), row.names = FALSE,
                   col.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("data: %d rows, %d unique\n", nrow(codes), nrow(unique(codes))))

inp <- function(alg) sprintf('TITLE:
elementary interaction, all ordered, default settings (%s)
DATA:
FILE = "ord.dat";
VARIABLE:
NAMES = %s;
CATEGORICAL = %s;
ANALYSIS:
type = random;
algorithm = integration %s;
MODEL:
X BY x1 x2;
Z BY z1 z2;
Y BY y1 y2;
Y ON X Z;
XZ | X XWITH Z;
Y ON XZ;
OUTPUT: TECH5 TECH8;
', alg, paste(vars, collapse = " "), paste(vars, collapse = " "), alg)

for (alg in c("EM", "EMA")) {
  writeLines(inp(alg), file.path(dir, sprintf("%s.inp", alg)))
  owd <- setwd(dir)
  t0 <- proc.time()[["elapsed"]]
  system2(MPLUS, sprintf("%s.inp", alg), stdout = TRUE, stderr = TRUE)
  cat(sprintf("%s: %.0f s\n", alg, proc.time()[["elapsed"]] - t0))
  setwd(owd)
}
cat("dir:", dir, "\n")
