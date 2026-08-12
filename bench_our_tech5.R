## Our own TECH5. For every EM iteration, report what each M-step block
## actually did: the gradient norm it started from, the objective before and
## after, and -- for the measurement block -- which path was taken and, on the
## Newton path, which backtracking scale was accepted.
##
## The question this is built to answer: is the Newton step being throttled by
## damping (scale 0.5/0.25 accepted), or is a FULL step simply a weaker move
## than nlminb's line search near the optimum? Those have different fixes.
suppressMessages(devtools::load_all(".", quiet = TRUE, compile = FALSE))
args <- commandArgs(trailingOnly = TRUE)
step <- args[[1]]; iters <- as.integer(args[[2]]); alg <- args[[3]]

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

ns <- asNamespace("modsem")
TRACE <- new.env(); TRACE$iter <- 0L; TRACE$rows <- list()

## Report which backtracking scale the Newton path accepted.
newton.orig <- get("lmsGraphNewtonMeasurement", envir = ns)
newton.traced <- function(theta, model, P, measurement, lower, upper,
                          objective, link = "logit") {
  probe <- new.env(); probe$n <- 0L
  wrapped <- function(x) { probe$n <- probe$n + 1L; objective(x) }
  out <- newton.orig(theta, model, P, measurement, lower, upper, wrapped, link)
  TRACE$newton <- if (is.null(out)) "declined"
                  else sprintf("scale=%s", c("1", "0.5", "0.25")[probe$n])
  out
}
unlockBinding("lmsGraphNewtonMeasurement", ns)
assign("lmsGraphNewtonMeasurement", newton.traced, envir = ns)

mstep.orig <- get("mstepLmsGraphEcm", envir = ns)
f <- function(theta, model, P, ...) {
  TRACE$iter <- TRACE$iter + 1L
  TRACE$newton <- "-"
  partition <- lmsGraphEcmPartition(theta, model)
  before <- completeAtEstepLmsGraph(P, sign = -1)
  grad <- gradientCompLogLikLmsGraph(theta, model, P, sign = -1)
  gm <- sqrt(sum(grad[partition$measurement]^2))
  gs <- sqrt(sum(grad[partition$structural]^2))
  out <- mstep.orig(theta, model, P, ...)
  TRACE$rows[[length(TRACE$rows) + 1L]] <- data.frame(
    iter = TRACE$iter, complete.before = before, complete.after = out$objective,
    improvement = before - out$objective, deriv.meas = gm, deriv.struct = gs,
    path = TRACE$newton, mstep.iters = out$iterations)
  out
}
formals(f)$... <- NULL
body(f) <- body(f)   # keep as-is
unlockBinding("mstepLmsGraphEcm", ns)
g <- mstep.orig; formals(g)$measurement.step <- step
assign("mstepLmsGraphEcm", local({ inner <- g; function(theta, model, P, ...) {
  TRACE$iter <- TRACE$iter + 1L
  TRACE$newton <- "-"
  partition <- lmsGraphEcmPartition(theta, model)
  before <- completeAtEstepLmsGraph(P, sign = -1)
  grad <- gradientCompLogLikLmsGraph(theta, model, P, sign = -1)
  out <- inner(theta, model, P, ...)
  TRACE$rows[[length(TRACE$rows) + 1L]] <- data.frame(
    iter = TRACE$iter, complete.before = before,
    complete.after = out$objective, improvement = before - out$objective,
    deriv.meas = sqrt(sum(grad[partition$measurement]^2)),
    deriv.struct = sqrt(sum(grad[partition$structural]^2)),
    path = TRACE$newton, mstep.iters = out$iterations)
  out
} }), envir = ns)

fit <- suppressWarnings(suppressMessages(modsem(
  'X =~ x1 + x2
   Z =~ z1 + z2
   Y =~ y1 + y2
   Y ~ X + Z + X:Z',
  dat, method = "lms", lms.backend = "graph", ordered = vars,
  adaptive = "full", integration = "gh", nodes = 8L,
  algorithm = alg, max.iter = iters, calc.se = FALSE, verbose = FALSE)))

T <- do.call(rbind, TRACE$rows)
cat(sprintf("\n===== %s / %s : logLik %.4f, %d M-steps =====\n",
            step, alg, as.numeric(fit$logLik), nrow(T)))
cat("\n-- first 10 --\n"); print(head(T, 10), row.names = FALSE, digits = 6)
cat("\n-- last 10 --\n");  print(tail(T, 10), row.names = FALSE, digits = 6)
cat("\npath taken:\n"); print(table(T$path))
