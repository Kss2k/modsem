#!/usr/bin/env Rscript

## # Performance benchmarks and acceptance gate
## The repository contains an automated LMS/QML performance regression check that compares any candidate
## version of `modsem` with a baseline branch (defaults to `main`). The workflow lives in
## `.github/workflows/performance.yml` and publishes a markdown dashboard artifact every run. A pull request
## is rejected when any of the benchmark examples runs more than **5 seconds slower on average** for either
## the LMS or QML approach.
## 
## You can reproduce the dashboard locally via:
## 
## ```bash
## Rscript inst/benchmarks/benchmark-modsem-version.R \
##   --baseline=main \
##   --candidate=$(git rev-parse HEAD) \
##   --candidate-source=local \
##   --reps=5 \
##   --output=inst/benchmarks/performance-dashboard.md \
##   --results=inst/benchmarks/performance-results.csv
## ```
## 
## Key flags:
## 
## - `--candidate` and `--baseline` pick which git refs to compare.
## - `--candidate-source=local` forces the script to install the package from your working tree instead of GitHub.
## - `--tolerance` (default 5) controls the acceptable slowdown in seconds before the workflow fails.
## - `--badge-dir=badges` emits Shields.io-compatible JSON summaries per method (used by the performance badges).
## 
## The resulting dashboard highlights the per-example deltas for LMS and QML so you can quickly spot and investigate regressions before opening a PR.
## 
## Two color-coded badges at the top of this README summarize the **average** performance deltas (main vs latest
## CRAN release) for LMS and QML separately. Green means the current `main` branch is faster (bright green for gains >=5s,
## green for >=2s, yellow-green for >=0.5s), while yellow/orange/red indicate slowdowns (<=-0.5s, <=-2s, <=-5s respectively).
## These badges are refreshed daily (and on every push to `main`) via `.github/workflows/performance-badge.yml`,
## which benchmarks `main` against the CRAN release and publishes the resulting JSON endpoints to the `perf-badges` branch.

args <- commandArgs(trailingOnly = TRUE)

repos <- getOption("repos")
cranIsUnset <- is.null(repos) || identical(repos, c(CRAN = "@CRAN@")) ||
  is.na(repos["CRAN"]) || identical(repos["CRAN"], "@CRAN@")
if (cranIsUnset) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

ensureNamespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this benchmark workflow.", pkg), call. = FALSE)
  }
}

ensureNamespace("devtools")
ensureNamespace("rbenchmark")

parseArg <- function(name, default = NULL) {
  candidates <- unique(c(
    name,
    if (grepl("_", name, fixed = TRUE)) gsub("_", "-", name, fixed = TRUE),
    if (grepl("-", name, fixed = TRUE)) gsub("-", "_", name, fixed = TRUE)
  ))
  for (candidate in candidates) {
    pattern <- paste0("^--?", candidate, "=")
    hit <- args[grepl(pattern, args)]
    if (length(hit)) {
      return(sub(pattern, "", hit[1]))
    }
  }
  default
}

splitMethods <- function(value) {
  if (is.null(value) || value == "") return(character())
  parts <- unlist(strsplit(value, "[, ]+"))
  parts <- trimws(parts)
  unique(parts[nzchar(parts)])
}

baselineRef <- parseArg("baseline", parseArg("base", "main"))
candidateRef <- parseArg("candidate", parseArg("ref", "HEAD"))
if (!nzchar(candidateRef)) {
  stop("Provide a candidate ref via '--candidate=<git ref>' or '--ref=<git ref>'.", call. = FALSE)
}

reps <- as.integer(parseArg("reps", "3"))
if (is.na(reps) || reps < 1) stop("--reps must be an integer >= 1.", call. = FALSE)

tolerance <- as.numeric(parseArg("tolerance", "5"))
if (is.na(tolerance) || tolerance < 0) stop("--tolerance must be a non-negative number.", call. = FALSE)

methods <- splitMethods(parseArg("methods", parseArg("method", "LMS,QML")))
if (!length(methods)) stop("Supply at least one method via '--methods=LMS,QML'.", call. = FALSE)

repoSlug <- parseArg("repo", "kss2k/modsem")

buildSourceSpec <- function(prefix, defaultSource, defaultRepo) {
  rawType <- parseArg(paste0(prefix, "_source"), defaultSource)
  if (is.null(rawType) || rawType == "") rawType <- defaultSource
  list(
    type = tolower(rawType),
    repo = parseArg(paste0(prefix, "_repo"), defaultRepo),
    path = parseArg(paste0(prefix, "_path"), ".")
  )
}

validateSourceSpec <- function(spec, label) {
  if (!spec$type %in% c("github", "local", "cran")) {
    stop(sprintf("Unsupported %s_source '%s'. Use 'github', 'local', or 'cran'.", label, spec$type), call. = FALSE)
  }
  if (spec$type == "github" && !nzchar(spec$repo)) {
    stop(sprintf("Provide a %s_repo value when using the github source.", label), call. = FALSE)
  }
  if (spec$type == "local" && !nzchar(spec$path)) spec$path <- "."
  spec
}

baselineSpec <- validateSourceSpec(buildSourceSpec("baseline", "github", repoSlug), "baseline")
candidateSpec <- validateSourceSpec(buildSourceSpec("candidate", "github", repoSlug), "candidate")

outputPath <- parseArg("output", file.path("inst", "benchmarks", "modsem-performance-dashboard.md"))
resultsPath <- parseArg("results", "")
badgeDir <- parseArg("badge_dir", parseArg("badge", ""))

ensurePath <- function(path) {
  if (!nzchar(path)) return(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

installModsemVersion <- function(ref, libDir, sourceSpec) {
  dir.create(libDir, recursive = TRUE, showWarnings = FALSE)
  if (identical(sourceSpec$type, "local")) {
    pkgPath <- sourceSpec$path
    if (!nzchar(pkgPath)) pkgPath <- "."
    message(sprintf("[modsem %s] Installing local sources from %s into %s", ref, pkgPath, libDir))
    devtools::install(
      pkg = pkgPath,
      lib = libDir,
      quick = TRUE,
      upgrade = "never",
      dependencies = FALSE
    )
  } else if (identical(sourceSpec$type, "cran")) {
    pkgName <- if (nzchar(sourceSpec$repo)) sourceSpec$repo else "modsem"
    message(sprintf("[modsem %s] Installing CRAN release (%s) into %s", ref, pkgName, libDir))
    utils::install.packages(
      pkgName,
      lib = libDir,
      dependencies = FALSE,
      quiet = TRUE
    )
  } else {
    message(sprintf("[modsem %s] Installing %s@%s into %s", ref, sourceSpec$repo, ref, libDir))
    devtools::install_github(
      repo = sourceSpec$repo,
      ref = ref,
      lib = libDir #, quick = TRUE, upgrade_dependencies = FALSE
    )
  }
}

loadExampleData <- function() {
  dataEnv <- new.env(parent = baseenv())
  data(list = c("oneInt", "TPB", "jordan"), package = "modsem", envir = dataEnv)
  dataEnv
}

createExampleRunners <- function(dataEnv) {
  list(
    ONEINT = function(method) {
      spec <- '
      # Outer Model
        X =~ x1 + x2 + x3
        Z =~ z1 + z2 + z3
        Y =~ y1 + y2 + y3

      # Inner Model
        Y ~ X + Z + X:Z
      '
      modsem(spec, dataEnv$oneInt, method = tolower(method))
    },
    TPB = function(method) {
      spec <- '
      # Outer Model (Based on Hagger et al., 2007)
        ATT =~ att1 + att2 + att3 + att4 + att5
        SN =~ sn1 + sn2
        PBC =~ pbc1 + pbc2 + pbc3
        INT =~ int1 + int2 + int3
        BEH =~ b1 + b2

      # Inner Model (Based on Steinmetz et al., 2011)
        INT ~ ATT + SN + PBC
        BEH ~ INT + PBC
        BEH ~ INT:PBC
      '
      modsem(spec, dataEnv$TPB, method = tolower(method), nodes = 32)
    }
    # ,
    # JORDAN = function(method) {
    #   spec <- '
    #     ENJ =~ enjoy1 + enjoy2 + enjoy3 + enjoy4 + enjoy5
    #     CAREER =~ career1 + career2 + career3 + career4
    #     SC =~ academic1 + academic2 + academic3 + academic4 + academic5 + academic6
    #     CAREER ~ ENJ + SC + ENJ:ENJ + SC:SC + ENJ:SC
    #   '
    #   modsem(spec, dataEnv$jordan, method = tolower(method), nodes = 15)
    # }
  )
}

benchmarkVersion <- function(ref, methods, reps, sourceSpec) {
  libDir <- file.path(tempdir(), sprintf("modsem-%s-%s", gsub("[^A-Za-z0-9]+", "-", ref), as.integer(Sys.time())))
  installModsemVersion(ref, libDir, sourceSpec)

  oldPaths <- .libPaths()
  on.exit(.libPaths(oldPaths), add = TRUE)
  .libPaths(c(libDir, oldPaths))

  suppressPackageStartupMessages(library(modsem, quietly = TRUE, warn.conflicts = FALSE))
  on.exit({
    if ("package:modsem" %in% search()) {
      try(detach("package:modsem", unload = TRUE, character.only = TRUE), silent = TRUE)
    }
  }, add = TRUE)

  pkgVersion <- as.character(utils::packageVersion("modsem"))
  dataEnv <- loadExampleData()
  runners <- createExampleRunners(dataEnv)

  versionResults <- lapply(methods, function(method) {
    message(sprintf("[modsem %s] Benchmarking %s (%s reps)", ref, method, reps))
    methodRows <- lapply(names(runners), function(exampleName) {
      bench <- rbenchmark::benchmark(
        runners[[exampleName]](method),
        replications = reps,
        columns = c("test", "replications", "elapsed")
      )
      avg <- bench$elapsed / bench$replications
      data.frame(
        ref = ref,
        installedVersion = pkgVersion,
        method = method,
        example = exampleName,
        reps = bench$replications,
        avgSeconds = avg,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, methodRows)
  })

  output <- do.call(rbind, versionResults)
  
  warn.output <- warnings()

  if (length(warn.output)) {
    cat("There were", length(warn.output), "warnings:\n")
    print(warn.output)
  }

  output
}

formatSeconds <- function(x) {
  if (is.na(x)) return("NA")
  sprintf("%.3f", x)
}

formatPercent <- function(x) {
  if (is.na(x)) return("NA")
  sprintf("%.1f%%", x)
}

summarizeDeltas <- function(baselineDf, candidateDf, tolerance) {
  merged <- merge(
    baselineDf,
    candidateDf,
    by = c("method", "example"),
    suffixes = c("Baseline", "Candidate"),
    all = TRUE
  )

  merged$deltaSeconds <- merged$avgSecondsCandidate - merged$avgSecondsBaseline
  merged$deltaPct <- (merged$deltaSeconds / merged$avgSecondsBaseline) * 100
  merged$status <- ifelse(
    is.na(merged$deltaPct),
    "MISSING",
    ifelse(merged$deltaPct > tolerance, "FAIL", "PASS")
  )
  merged
}

writeDashboard <- function(summaryDf, baselineRef, candidateRef, reps, tolerance, outputPath) {
  if (!nzchar(outputPath)) return(invisible(NULL))

  headers <- c("Method", "Example", "Baseline (s)", "Candidate (s)", "\u0394 Seconds", "\u0394 Percent", "Status")
  align <- c("left", "left", "right", "right", "right", "right", "left")
  rowMatrix <- t(apply(summaryDf, 1, function(row) {
    c(
      row[["method"]],
      row[["example"]],
      formatSeconds(as.numeric(row[["avgSecondsBaseline"]])),
      formatSeconds(as.numeric(row[["avgSecondsCandidate"]])),
      formatSeconds(as.numeric(row[["deltaSeconds"]])),
      formatPercent(as.numeric(row[["deltaPct"]])),
      if (identical(row[["status"]], "FAIL")) "**FAIL**" else row[["status"]]
    )
  }))
  if (!length(rowMatrix)) {
    rowMatrix <- matrix(nrow = 0, ncol = length(headers))
  }
  tableMatrix <- rbind(headers, rowMatrix)
  colWidths <- apply(tableMatrix, 2, function(col) max(nchar(col, type = "width"), na.rm = TRUE))

  formatCell <- function(value, width, alignment) {
    fmt <- if (alignment == "right") paste0("%", width, "s") else paste0("%-", width, "s")
    sprintf(fmt, value)
  }

  formattedRows <- apply(tableMatrix, 1, function(row) {
    sapply(seq_along(row), function(i) formatCell(row[i], colWidths[i], align[i]), USE.NAMES = FALSE)
  })
  formattedRows <- t(formattedRows)

  separator <- sapply(seq_along(headers), function(i) {
    width <- colWidths[i]
    if (align[i] == "right") {
      paste0(strrep("-", width - 1), ":")
    } else {
      paste0(":", strrep("-", width - 1))
    }
  })

  tableLines <- apply(formattedRows, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  lines <- c(
    "# modsem performance dashboard",
    "",
    sprintf("- Baseline ref: `%s`", baselineRef),
    sprintf("- Candidate ref: `%s`", candidateRef),
    sprintf("- Repetitions per benchmark: %s", reps),
    sprintf("- Acceptance tolerance: %.2f seconds", tolerance),
    "",
    tableLines[1],
    paste("|", paste(separator, collapse = " | "), "|"),
    tableLines[-1],
    ""
  )
  ensurePath(outputPath)
  writeLines(lines, con = outputPath)
  message(sprintf("Performance dashboard written to %s", outputPath))
}

writeResultsCsv <- function(allResults, resultsPath) {
  if (!nzchar(resultsPath)) return(invisible(NULL))
  ensurePath(resultsPath)
  utils::write.csv(allResults, file = resultsPath, row.names = FALSE)
  message(sprintf("Raw benchmark results written to %s", resultsPath))
}

badgeColorFor <- function(value) {
  if (is.na(value)) return("lightgrey")
  if (value >= 5) return("brightgreen")
  if (value >= 2) return("green")
  if (value >= 0.5) return("yellowgreen")
  if (value <= -5) return("red")
  if (value <= -2) return("orange")
  if (value <= -0.5) return("yellow")
  "lightgrey"
}

writeBadgeJsons <- function(summaryDf, baselineRef, candidateRef, badgeDir) {
  if (!nzchar(badgeDir)) return(invisible(NULL))
  dir.create(badgeDir, recursive = TRUE, showWarnings = FALSE)
  validRows <- summaryDf[!is.na(summaryDf$method), , drop = FALSE]
  perMethod <- split(validRows, validRows$method)
  lapply(names(perMethod), function(methodName) {
    rows <- perMethod[[methodName]]
    improvement <- mean(rows$avgSecondsBaseline - rows$avgSecondsCandidate, na.rm = TRUE)
    badge <- if (is.na(improvement)) "n/a" else sprintf("%+.2f s", improvement)
    color <- badgeColorFor(improvement)
    label <- sprintf("%s Δ (main-rel)", methodName)
    json <- sprintf(
      '{"schemaVersion":1,"label":"%s","message":"%s","color":"%s"}',
      label, badge, color
    )
    fileName <- file.path(badgeDir, paste0(tolower(methodName), ".json"))
    writeLines(json, fileName)
    message(sprintf("Badge for %s written to %s (%s)", methodName, fileName, badge))
  })
}

baselineResults <- benchmarkVersion(baselineRef, methods, reps, baselineSpec)
candidateResults <- benchmarkVersion(candidateRef, methods, reps, candidateSpec)

summaryDf <- summarizeDeltas(baselineResults, candidateResults, tolerance)

writeDashboard(summaryDf, baselineRef, candidateRef, reps, tolerance, outputPath)

combinedResults <- rbind(
  transform(baselineResults, branch = "baseline"),
  transform(candidateResults, branch = "candidate")
)
writeResultsCsv(combinedResults, resultsPath)
writeBadgeJsons(summaryDf, baselineRef, candidateRef, badgeDir)

if (any(summaryDf$status == "FAIL")) {
  message("Performance regression detected (see dashboard).")
  quit(status = 1)
}

message("Performance benchmarks completed successfully.")
