FROM continuumio/miniconda3:24.7.1-0

USER root

COPY install_genespace.R /tmp/install_genespace.R

RUN conda install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    r-base=4.5.1 \
    git \
    make \
    orthofinder=2.5.5 \
    diamond=2.1.13 \
    blast=2.17.0 \
    mafft=7.526 \
    iqtree=3.0.1 \
    raxml-ng=1.2.2 \
    raxml=8.2.13 \
    fasttree=2.2.0 \
    fastme=2.1.6.3 \
    muscle=3.8.1551 \
    mmseqs2 \
    mcscanx=1.0.0 \
    mcl \
    ete3=3.1.3 \
    entrez-direct \
    aria2 \
    openmpi \
    gsl \
    biopython=1.85 \
    numpy=2.3.3 \
    scipy=1.16.2 \
    scikit-learn=1.7.2 \
    joblib \
    lxml \
    packaging \
    pip \
    r-remotes \
    r-devtools \
    r-data.table \
    r-ggplot2 \
    r-igraph \
    r-rmarkdown \
    r-knitr \
    r-testthat \
    r-yaml \
    r-readr \
    r-dplyr \
    r-optparse \
    r-getopt \
    r-progress \
    r-vroom \
    r-bit \
    r-bit64 \
    r-hms \
    r-tidyselect \
    r-tzdb \
    r-rlang \
    r-vctrs \
    r-lifecycle \
    bioconductor-biobase \
    bioconductor-biocgenerics \
    bioconductor-biocio \
    bioconductor-biocparallel \
    bioconductor-biostrings \
    bioconductor-delayedarray \
    bioconductor-genomeinfodb \
    bioconductor-genomeinfodbdata \
    bioconductor-genomicalignments \
    bioconductor-genomicranges \
    bioconductor-iranges \
    bioconductor-matrixgenerics \
    bioconductor-rhtslib \
    bioconductor-rsamtools \
    bioconductor-rtracklayer \
    bioconductor-s4arrays \
    bioconductor-s4vectors \
    bioconductor-sparsearray \
    bioconductor-summarizedexperiment \
    bioconductor-ucsc.utils \
    bioconductor-xvector \
    && conda clean --all --yes

RUN Rscript /tmp/install_genespace.R && rm -f /tmp/install_genespace.R

ENV PATH=/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /work
