FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL maintainer="Brianna Banting"
LABEL description="R/ggtree container for rediploidisation species-tree plotting"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
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
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('remotes', 'BiocManager'), repos = 'https://cloud.r-project.org')"

RUN Rscript -e "remotes::install_version('ggplot2', version = '3.5.2', repos = 'https://cloud.r-project.org')"

RUN Rscript -e "install.packages(c( \
    'optparse', \
    'ape', \
    'readr', \
    'dplyr', \
    'stringr', \
    'tidyr' \
), repos = 'https://cloud.r-project.org')"

RUN Rscript -e "BiocManager::install('ggtree', ask = FALSE, update = FALSE); \
    if (!requireNamespace('ggtree', quietly = TRUE)) stop('ggtree did not install')"

RUN Rscript -e "suppressPackageStartupMessages({ \
    library(optparse); \
    library(ape); \
    library(ggplot2); \
    library(ggtree); \
    library(readr); \
    library(dplyr); \
    library(stringr); \
    library(tidyr) \
}); \
cat('ggplot2=', as.character(packageVersion('ggplot2')), '\n'); \
cat('ggtree=', as.character(packageVersion('ggtree')), '\n'); \
cat('R species-tree plotting environment OK\n')"

WORKDIR /work

CMD ["Rscript", "--version"]
