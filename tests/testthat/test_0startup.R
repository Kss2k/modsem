devtools::load_all()

OK_OpenMP <- checkOpenMP_Cpp()
startupMessages <- purrr::quietly(devtools::load_all)()$messages[-1L]
startupMessage <- paste0("    ", startupMessages, collapse = "\n")


printSessionInfo <- function() {
  info <- sessionInfo()
  info$date <- tryCatch(as.character(Sys.Date()), error = "<NA>")
  info$openmp <- if (OK_OpenMP) "AVAILABLE" else "NOT AVAILABLE"

  show <- (sapply(info, FUN = length) == 1L) & names(info) != "locale"
  rhs <- names(info[show])
  lhs <- unlist(info[show])

  # cut paths
  largePath <- stringr::str_count(rhs, pattern = "/") > 3L
  browser()
  getstem <- \(path, i) stringr::str_split_i(path, pattern = "/", i = i)
  rhs[ispath] <- paste0(

    stringr::str_split_i(rhs[ispath], pattern = "/", i = -1)
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
  cat(startupMessage)

  cat(sep1, "\n")
}

printSessionInfo()

OpenMP_Pattern <- "OpenMP is not available.*"

if (OK_OpenMP) {
  OK_OpenMP_Message <- !grepl(OpenMP_Pattern, startupMessage)
} else  {
  OK_OpenMP_Message <-  grepl(OpenMP_Pattern, startupMessage)
}
testthat::expect_true(OK_OpenMP_Message)

