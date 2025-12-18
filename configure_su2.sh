#!/usr/bin/env bash
set -euo pipefail

# Always run from the repo root (the directory where this script lives)
cd "$(dirname "$0")"

# Avoid cross-branch contamination
rm -rf build

meson setup build \
  --prefix=/Users/reumschuessel/SU2 \
  -Dwith-mpi=enabled \
  -Denable-autodiff=true \
  -Denable-directdiff=true \
  -Denable-pywrapper=true \
  -Denable-tests=true \
  -Dwith-omp=false

./ninja -C build install 