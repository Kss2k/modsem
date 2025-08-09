OK_OpenMP <- checkOpenMP_Cpp()
startupMessages <- purrr::quietly(devtools::load_all)()$messages[-1L]
startupMessage <- paste0("    ", startupMessages, collapse = "\n")

printSessionInfo <- function() {
  info <- sessionInfo()
  info$date <- tryCatch(as.character(Sys.Date()), error = "<NA>")
  info$openmp <- if (OK_OpenMP) "AVAILABLE" else "NOT AVAILABLE"

  show <- (sapply(info, FUN = length) == 1L) & names(info) != "locale"
  lhs <- format(names(info[show]), justify = "left")
  rhs <- format(unlist(info[show]), justify = "right")

  width <- options("width")[[1L]]
  space <- strrep(" ", max(2L, width - max(nchar(lhs)) - max(nchar(rhs))))
  sep1  <- paste0(strrep("\u2500", width), "\n")
  sep2  <- paste0(strrep("\u2550", width), "\n")


  cat(sep1)
  cat("RUNNING TESTS CASES\n")
  cat(sep2)

  for (i in seq_along(lhs)) cat(lhs[i], space, rhs[i], "\n", sep = "")

  cat(sep1)
  cat("STARTUP MESSAGE:\n")
  cat(startupMessage)

  cat(sep1)


}

printSessionInfo()

OpenMP_Pattern <- "OpenMP is not available.*"

if (OK_OpenMP) {
  OK_OpenMP_Message <- !grepl(OpenMP_Pattern, startupMessage)
} else  {
  OK_OpenMP_Message <-  grepl(OpenMP_Pattern, startupMessage)
}
testthat::expect_true(OK_OpenMP_Message)

