FROM continuumio/miniconda3:24.7.1-0

LABEL software="GENESPACE" \
      git_commit="7561036be821e333c1dfb52461ccec7222e95582" \
      org.opencontainers.image.source="https://github.com/bhban/wgd-nextflow"

USER root

COPY install_genespace.R /tmp/install_genespace.R
COPY ofInBlk_engine_patched.R /tmp/ofInBlk_engine_patched.R

RUN conda create -y -n genespace --override-channels \
    -c conda-forge -c bioconda \
    python=3.11.13 \
    numpy=2.3.3 \
    scipy=1.16.2 \
    r-base=4.5 \
    orthofinder=2.5.5 \
    mcscanx=1.0.0 \
    pandas \
    bioconductor-biostrings \
    bioconductor-rtracklayer \
    r-remotes \
    r-igraph \
    r-yaml \
    r-optparse \
    r-dplyr \
    r-readr \
    r-data.table \
    xz \
    && conda clean --all --yes

RUN conda run -n genespace Rscript /tmp/install_genespace.R && \
    rm -f /tmp/install_genespace.R /tmp/ofInBlk_engine_patched.R

ENV PATH=/opt/conda/envs/genespace/bin:/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /work
