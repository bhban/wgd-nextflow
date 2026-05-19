# Dockerfile for rediploidisation species-tree plotting with ggtree
#
# Required R packages:
#   optparse
#   ape
#   ggplot2
#   ggtree
#   readr
#   dplyr
#   stringr
#   tidyr
#
# ggtree is installed from Bioconductor via BiocManager.

FROM rocker/r-ver:4.4.2

LABEL maintainer="Brianna Banting"
LABEL description="R/ggtree container for rediploidisation species-tree plotting"

ENV DEBIAN_FRONTEND=noninteractive

# System libraries needed by common tidyverse/Bioconductor packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    zlib1g-dev \
    pandoc \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages first.
RUN Rscript -e "install.packages(c( \
    'BiocManager', \
    'optparse', \
    'ape', \
    'ggplot2', \
    'readr', \
    'dplyr', \
    'stringr', \
    'tidyr' \
), repos = 'https://cloud.r-project.org')"

# Install ggtree from Bioconductor.
RUN Rscript -e "BiocManager::install(c('ggtree'), ask = FALSE, update = FALSE)"

# Quick import test so the image build fails early if anything is missing.
RUN Rscript -e "suppressPackageStartupMessages({ \
    library(optparse); \
    library(ape); \
    library(ggplot2); \
    library(ggtree); \
    library(readr); \
    library(dplyr); \
    library(stringr); \
    library(tidyr) \
}); cat('R species-tree plotting environment OK\n')"

WORKDIR /work

CMD ["Rscript", "--version"]
