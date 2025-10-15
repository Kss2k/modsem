#' @export
multilevel_lms <- function(model.syntax, data, cluster, method = "lms", group = NA, ...) {
  levelPattern <- "[\\s\\n]*level([:blank:]*):([:blank:]*)([A-z]|[0-9]+)"
  levelHeaders <- unlist(stringr::str_extract_all(model.syntax, pattern=levelPattern))
  syntaxBlocks <- unlist(stringr::str_split(model.syntax, pattern=levelPattern))
  syntaxBlocks <- syntaxBlocks[syntaxBlocks!=""]

  parTableL1 <- modsemify(syntaxBlocks[[1L]])
  parTableL2 <- modsemify(syntaxBlocks[[2L]])

  k <- length(unique(data[[cluster]]))

  mapping <- c(`=~` = "m", `~` = "r", `~1` = "i", `~~` = "c")
  # constrain level 1 parameters to unity
  labels <- paste0(parTableL1$lhs, mapping[parTableL1$op], parTableL1$rhs)
  parTableL1$mod <- ifelse(parTableL1$mod == "", yes = labels, no = parTableL1$mod)

  multigroupSyntax <- parTableToSyntax(rbind(parTableL1, parTableL2))
  
  modsem(model.syntax = multigroupSyntax, data = data, method = method,
         group = cluster, multilevel.free.pars = TRUE, ...)
}
