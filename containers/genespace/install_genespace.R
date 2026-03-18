cran_repo <- "https://cloud.r-project.org"

cran_pkgs <- c(
  "abind","askpass","base64enc","BH","bit","bit64","brew","brio","bslib",
  "cachem","callr","cli","clipr","commonmark","cpp11","crayon","credentials",
  "curl","data.table","dbscan","desc","devtools","diffobj","digest","downlit",
  "dplyr","ellipsis","evaluate","fansi","farver","fastmap","fontawesome",
  "formatR","fs","futile.logger","futile.options","generics","getopt","gert",
  "ggplot2","gh","gitcreds","glue","gtable","highr","hms","htmltools",
  "htmlwidgets","httpuv","httr","httr2","igraph","ini","isoband","jquerylib",
  "jsonlite","knitr","labeling","lambda.r","later","lifecycle","magrittr",
  "memoise","mime","miniUI","openssl","optparse","pillar","pkgbuild",
  "pkgconfig","pkgdown","pkgload","praise","prettyunits","processx","profvis",
  "progress","promises","ps","purrr","R.methodsS3","R.oo","R.utils","R6",
  "ragg","rappdirs","rcmdcheck","RColorBrewer","Rcpp","RCurl","readr",
  "remotes","restfulr","rjson","rlang","rmarkdown","roxygen2","rprojroot",
  "rstudioapi","rversions","S7","sass","scales","sessioninfo","shiny","snow",
  "sourcetools","stringi","stringr","sys","systemfonts","testthat",
  "textshaping","tibble","tidyselect","tinytex","tzdb","urlchecker","usethis",
  "utf8","vctrs","viridisLite","vroom","waldo","whisker","withr","xfun","XML",
  "xml2","xopen","xtable","yaml","zip"
)

bioc_pkgs <- c(
  "Biobase","BiocGenerics","BiocIO","BiocParallel","Biostrings",
  "DelayedArray","GenomeInfoDb","GenomeInfoDbData","GenomicAlignments",
  "GenomicRanges","IRanges","MatrixGenerics","Rhtslib","Rsamtools",
  "rtracklayer","S4Arrays","S4Vectors","SparseArray","SummarizedExperiment",
  "UCSC.utils","XVector"
)

install_if_missing <- function(pkgs, install_fun) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install_fun(missing)
}

install_if_missing(cran_pkgs, function(pkgs) {
  install.packages(pkgs, repos = cran_repo, Ncpus = max(1, parallel::detectCores() - 1))
})

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = cran_repo)
}

install_if_missing(bioc_pkgs, function(pkgs) {
  BiocManager::install(pkgs, ask = FALSE, update = FALSE, Ncpus = max(1, parallel::detectCores() - 1))
})

if (!requireNamespace("GENESPACE", quietly = TRUE)) {
  remotes::install_github("jtlovell/GENESPACE@1.3.1", upgrade = "never")
}

cat("GENESPACE version:", as.character(utils::packageVersion("GENESPACE")), "\n")
