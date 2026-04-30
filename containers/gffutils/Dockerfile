FROM mambaorg/micromamba:2.0.5

LABEL software="gffutils" \
      container_version="1.0" \
      org.opencontainers.image.source="https://github.com/bhban/wgd-nextflow"

USER root

WORKDIR /work

RUN micromamba install -y -n base \
      -c conda-forge \
      -c bioconda \
      python=3.14 \
      pandas \
      numpy \
      biopython \
      pyfaidx \
      gffutils \
      gffread \
      agat \
      bedtools \
      samtools \
      seqkit \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

WORKDIR /work
