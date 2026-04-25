#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/../apptainer"

mkdir -p "$OUT_DIR"

apptainer pull "$OUT_DIR/genespace_1.1.sif" docker://ghcr.io/bhban/genespace:1.1
apptainer pull "$OUT_DIR/macse_1.0.sif" docker://ghcr.io/bhban/macse:1.0
apptainer pull "$OUT_DIR/iqtree_1.0.sif" docker://ghcr.io/bhban/iqtree:1.0
apptainer pull "$OUT_DIR/alerax_1.0.sif" docker://ghcr.io/bhban/alerax:1.0
apptainer pull "$OUT_DIR/annevo_1.0.sif" docker://ghcr.io/bhban/annevo:1.0
apptainer pull "$OUT_DIR/gffutils_1.0.sif" docker://ghcr.io/bhban/gffutils:1.0
