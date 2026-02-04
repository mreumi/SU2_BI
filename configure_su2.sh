#!/usr/bin/env bash
set -euo pipefail

# Always run from the repo root (the directory where this script lives)
cd "$(dirname "$0")"

# Avoid cross-branch contamination
rm -rf build

# Make sure to adjust the python path to your su2 environment
cat > native-python.ini <<'EOF'
[binaries]
python = '/Users/reumschuessel/anaconda3/envs/su2_IP_env/bin/python'
EOF

meson setup build \
  --prefix=/Users/reumschuessel/SU2 \
  --native-file native-python.ini \
  -Dwith-mpi=enabled \
  -Denable-autodiff=true \
  -Denable-directdiff=true \
  -Denable-pywrapper=true \
  -Denable-tests=true \
  -Dwith-omp=false \
  -Denable-tecio=false

./ninja -C build install 

# after reboot:
# export SU2_RUN=/Users/reumschuessel/SU2
# export PATH=$SU2_RUN:$PATH
# export PYTHONPATH=$SU2_RUN:$PYTHONPATH
# export PYTHONPATH="/Users/reumschuessel/SU2/bin:$PYTHONPATH"
