FROM continuumio/miniconda3:24.7.1-0

LABEL software="redip" \
      org.opencontainers.image.source="https://github.com/bhban/wgd-nextflow"

USER root

RUN apt-get update && apt-get install -y \
    bash \
    coreutils \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y -c conda-forge -c bioconda \
    python=3.11 \
    ete3=3.1.3 \
    pandas \
    numpy \
    biopython \
    circos=0.69.9 \
    && conda clean -afy

ENV PATH="/opt/conda/bin:${PATH}"

WORKDIR /work

CMD ["/bin/bash"]
