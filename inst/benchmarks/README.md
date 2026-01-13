Performance benchmarks
======================

`inst/benchmarks/benchmark_models.R` fits four representative models:
the `oneInt` and `TPB` examples under both LMS and QML estimation. Elapsed
times are written to a CSV so different builds (CRAN, `main`, PR) can be
compared. Run it with any installed version of `modsem`:

```bash
Rscript inst/benchmarks/benchmark_models.R \
  --output=benchmark-results/pr.csv \
  --iters=3
```

Arguments:

- `--output=<file>`: CSV destination (default `benchmark-results.csv`)
- `--iters=<n>`: number of repetitions per model. You can also set
  `MODSEM_BENCH_ITERS` in the environment.

The script prints per-run timing as well as summary statistics, and the CSV
contains one row per model/iteration/method. It respects `R_LIBS_USER`, so you
can point it at different installed versions (CRAN release, `main`, PR build,
…) to compare performance.

GitHub PRs automatically run `.github/workflows/modsem-benchmarks.yml`, which
installs the CRAN release, the current `main` branch, and the PR commit into
separate libraries, runs this script for each, and uploads the resulting CSVs
and summaries as workflow artifacts.
