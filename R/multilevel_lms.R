#' @export
multilevel_lms <- function(model.syntax, data, method = "lms", ...) {
  levelPattern <- "[\s\n]+level([:blank:]*):([:blank:]*)([A-z]|[0-9]+)"
  levelHeaders <- unlist(stringr::str_extract_all(model.syntax, pattern=levelPattern))
  syntaxBlocks <- unlist(stringr::str_split(model.syntax, pattern=levelPattern))

  browser()

}
