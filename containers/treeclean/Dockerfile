FROM condaforge/mambaforge:24.3.0-0

LABEL org.opencontainers.image.title="treeclean"
LABEL org.opencontainers.image.description="Tree cleaning utilities with TreeShrink and ETE3"
LABEL org.opencontainers.image.authors="Brianna Banting"
LABEL org.opencontainers.image.source="https://github.com/bhban/wgd-nextflow"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN mamba create -y -n treeclean \
      -c conda-forge \
      -c smirarab \
      python=3.9 \
      treeshrink=1.3.9 \
      numpy \
      scipy \
      pandas \
      sed \
      coreutils \
      bash \
    && mamba clean -a -y

RUN conda run -n treeclean pip install --no-cache-dir ete3==3.1.3

ENV PATH=/opt/conda/envs/treeclean/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN run_treeshrink.py -h >/dev/null \
    && python - <<'PY'
from ete3 import Tree
t = Tree("(A:0.1,B:0.2);", format=1)
assert len(t) == 2
print("treeclean container OK")
PY
