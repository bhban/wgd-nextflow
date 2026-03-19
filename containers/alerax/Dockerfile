FROM ubuntu:24.04

LABEL software="AleRax" \
      version="1.4.1" \
      git_commit="775d91d14c06d14ba5afcc855f58bd1ab520c48c"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    ca-certificates \
    git \
    build-essential \
    cmake \
    make \
    pkg-config \
    libgomp1 \
    libstdc++6 \
    openmpi-bin \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN git clone --recursive https://github.com/BenoitMorel/AleRax.git \
    && cd AleRax \
    && git checkout 775d91d14c06d14ba5afcc855f58bd1ab520c48c \
    && git submodule update --init --recursive \
    && ./install.sh

ENV PATH="/opt/AleRax/build/bin:${PATH}"

WORKDIR /work
