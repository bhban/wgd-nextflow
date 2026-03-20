FROM continuumio/miniconda3:24.7.1-0

LABEL software="ANNEVO" \
      git_commit="b07cb8e344e4c49165c8f818d9bf7708a495240b"
      org.opencontainers.image.source="https://github.com/bhban/wgd-nextflow"
      
USER root

RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# clone repo
RUN git clone https://github.com/xjtu-omics/ANNEVO.git \
    && cd ANNEVO \
    && git checkout b07cb8e344e4c49165c8f818d9bf7708a495240b

WORKDIR /opt/ANNEVO

# create environment from YAML
RUN conda env create -f ANNEVO.yml -n annevo \
    && conda clean --all --yes

ENV PATH=/opt/conda/envs/annevo/bin:/opt/conda/bin:$PATH

WORKDIR /work
