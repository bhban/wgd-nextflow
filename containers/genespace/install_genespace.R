cran_repo <- "https://cloud.r-project.org"

# Install BiocManager before installing GENESPACE
install.packages("BiocManager", repos = cran_repo)
BiocManager::install(version = "3.21", ask = FALSE)

# Install required BiocManager packages
BiocManager::install(
  c("Biostrings", "rtracklayer"),
  ask = FALSE,
  update = FALSE
)

# Need devtools to install GENESPACE. Must install remotes first
install.packages("remotes", repos = cran_repo)
install.packages("devtools", repos = cran_repo)

devtools::install_github("jtlovell/GENESPACE@7561036be821e333c1dfb52461ccec7222e95582", upgrade = 'never')
