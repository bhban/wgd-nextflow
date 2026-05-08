FROM mambaorg/micromamba:1.5.8

LABEL software="MAFFT" \
      version="7.525" \
      org.opencontainers.image.source="https://github.com/bhban/wgd-nextflow"

USER root

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      mafft=7.525 \
      procps-ng \
      coreutils \
      bash \
    && micromamba clean -a -y

ENV PATH="/opt/conda/bin:${PATH}"

SHELL ["/bin/bash", "-c"]

CMD ["mafft", "--help"]
