FROM continuumio/miniconda3:24.7.1-0

USER root

COPY install_genespace.R /tmp/install_genespace.R

RUN conda create -y -n genespace --override-channels \
    -c conda-forge -c bioconda \
    python=3.11 \
    r-base=4.5 \
    orthofinder=2.5.5 \
    mcscanx=1.0.0 \
    pandas \
    bioconductor-biostrings \
    bioconductor-rtracklayer \
    r-remotes \
    r-igraph \
    xz \
    && conda clean --all --yes

RUN conda run -n genespace Rscript /tmp/install_genespace.R && \
    rm -f /tmp/install_genespace.R

ENV PATH=/opt/conda/envs/genespace/bin:/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /work
