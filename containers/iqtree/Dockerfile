FROM continuumio/miniconda3:24.7.1-0

LABEL software="IQ-TREE" \
      org.opencontainers.image.source="https://github.com/bhban/wgd-nextflow"
      
USER root

RUN conda create -y -n iqtree --override-channels \
    -c conda-forge -c bioconda \
    python=3.14 \
    iqtree=3.0.1 \
    pandas=3.0.1 \
    && conda clean --all --yes

ENV PATH=/opt/conda/envs/iqtree/bin:/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /work
