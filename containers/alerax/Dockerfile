FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    ca-certificates \
    git \
    build-essential \
    cmake \
    make \
    pkg-config \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN git clone --recursive https://github.com/BenoitMorel/AleRax.git

WORKDIR /opt/AleRax

RUN ./install.sh

RUN cp build/bin/alerax /usr/local/bin/alerax

WORKDIR /work
