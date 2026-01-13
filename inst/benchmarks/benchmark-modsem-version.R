#!/usr/bin/env Rscript

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
  pattern <- paste0("^--?", name, "=")
  hit <- args[grepl(pattern, args)]
  if (!length(hit)) return(default)
  sub(pattern, "", hit[1])
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
      modsem(spec, dataEnv$oneInt, method = method)
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
      modsem(spec, dataEnv$TPB, method = method)
    },
    JORDAN = function(method) {
      spec <- '
        ENJ =~ enjoy1 + enjoy2 + enjoy3 + enjoy4 + enjoy5
        CAREER =~ career1 + career2 + career3 + career4
        SC =~ academic1 + academic2 + academic3 + academic4 + academic5 + academic6
        CAREER ~ ENJ + SC + ENJ:ENJ + SC:SC + ENJ:SC
      '
      modsem(spec, dataEnv$jordan, method = method)
    }
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
        repetitions = reps,
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

  do.call(rbind, versionResults)
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
    is.na(merged$deltaSeconds),
    "MISSING",
    ifelse(merged$deltaSeconds > tolerance, "FAIL", "PASS")
  )
  merged
}

writeDashboard <- function(summaryDf, baselineRef, candidateRef, reps, tolerance, outputPath) {
  if (!nzchar(outputPath)) return(invisible(NULL))

  lines <- c(
    "# modsem performance dashboard",
    "",
    sprintf("- Baseline ref: `%s`", baselineRef),
    sprintf("- Candidate ref: `%s`", candidateRef),
    sprintf("- Repetitions per benchmark: %s", reps),
    sprintf("- Acceptance tolerance: %.2f seconds", tolerance),
    "",
    "| Method | Example | Baseline (s) | Candidate (s) | Δ Seconds | Δ Percent | Status |",
    "| --- | --- | ---: | ---: | ---: | ---: | --- |"
  )

  rowLines <- apply(summaryDf, 1, function(row) {
    baseVal <- formatSeconds(as.numeric(row[["avgSecondsBaseline"]]))
    candVal <- formatSeconds(as.numeric(row[["avgSecondsCandidate"]]))
    deltaVal <- formatSeconds(as.numeric(row[["deltaSeconds"]]))
    pctVal <- formatPercent(as.numeric(row[["deltaPct"]]))
    status <- row[["status"]]
    if (identical(status, "FAIL")) status <- "**FAIL**"
    paste0(
      "| ",
      row[["method"]], " | ",
      row[["example"]], " | ",
      baseVal, " | ",
      candVal, " | ",
      deltaVal, " | ",
      pctVal, " | ",
      status, " |"
    )
  })

  lines <- c(lines, rowLines, "")
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
