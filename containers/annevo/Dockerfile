FROM continuumio/miniconda3:24.7.1-0

USER root

WORKDIR /opt

RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/xjtu-omics/ANNEVO.git

WORKDIR /opt/ANNEVO

RUN conda create -y -n annevo \
    python=3.10 \
    && conda clean --all --yes

RUN conda install -y -n annevo \
    -c pytorch -c nvidia \
    pytorch==2.1.0 \
    torchvision==0.16.0 \
    torchaudio==2.1.0 \
    pytorch-cuda=12.1 \
    && conda clean --all --yes

RUN /opt/conda/bin/conda run -n annevo pip install --no-cache-dir -r requirements.txt

ENV PATH=/opt/conda/envs/annevo/bin:/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /work
