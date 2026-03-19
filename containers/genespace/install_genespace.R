cran_repo <- "https://cloud.r-project.org"
options(repos = c(CRAN = cran_repo))

# Install required Bioconductor packages if not successfully installed by conda
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Biostrings is not installed in the conda environment.")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  stop("rtracklayer is not installed in the conda environment.")
}

# Install remotes to install GENESPACE from github
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "jtlovell/GENESPACE@7561036be821e333c1dfb52461ccec7222e95582",
  upgrade = "never",
  dependencies = TRUE
)
