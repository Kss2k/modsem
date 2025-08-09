devtools::load_all()

OK_OpenMP <- checkOpenMP_Cpp()
startupMessage <- purrr::quietly(devtools::load_all)()$messages[-1L]
startupMessage <- stringr::str_split(startupMessage, pattern = "\n")[[1L]] |>
  (\(x) (x[x != ""]))() |>
  stringr::str_pad(width = 4) |> paste0(collapse = "\n")


printSessionInfo <- function() {
  info <- sessionInfo()
  info$date <- tryCatch(as.character(Sys.Date()), error = "<NA>")
  info$openmp <- if (OK_OpenMP) "AVAILABLE" else "NOT AVAILABLE"

  show <- (sapply(info, FUN = length) == 1L) & names(info) != "locale"
  lhs <- names(info[show])
  rhs <- unlist(info[show])

  # cut large paths
  isLargePath <- stringr::str_count(rhs, pattern = "/") > 2L
  getstem <- \(path, i) stringr::str_split_i(path, pattern = "/", i = i)

  largePaths <- rhs[isLargePath]
  rhs[isLargePath] <- paste(
    "...",
    # getstem(largePaths, -2L),
    getstem(largePaths, -1L),
    sep = "/"
  )

  flhs <- format(lhs, justify = "left")
  frhs <- format(rhs, justify = "right")

  width <- options("width")[[1L]]
  space <- strrep(" ", max(2L, width - max(nchar(lhs)) - max(nchar(rhs))))
  sep1  <- paste0(strrep("\u2500", width), "\n")
  sep2  <- paste0(strrep("\u2550", width), "\n")

  cat("\n")
  cat(sep1)
  cat("RUNNING TESTS CASES\n")
  cat(sep2)

  for (i in seq_along(lhs))
    cat(flhs[i], space, frhs[i], "\n", sep = "")

  cat(sep1)
  cat("STARTUP MESSAGE:\n")
  cat(startupMessage, "\n")

  cat(sep2, "\n")
}

printSessionInfo()

OpenMP_Pattern <- "OpenMP is not available.*"

if (OK_OpenMP) {
  OK_OpenMP_Message <- !grepl(OpenMP_Pattern, startupMessage)
} else  {
  OK_OpenMP_Message <-  grepl(OpenMP_Pattern, startupMessage)
}
testthat::expect_true(OK_OpenMP_Message)

