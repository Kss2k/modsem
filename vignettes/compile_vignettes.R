devtools::install()

setwd("vignettes")

files <- dir()
origfiles <- files[grepl("\\.Rmd\\.orig$", files)]

for (origfile in origfiles) {
  rmdfile <- stringr::str_remove(origfile, pattern = "\\.orig$")
  knitr::knit(input = origfile, output = rmdfile)
}