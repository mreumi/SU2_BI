#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

if [ ! -d subprojects/cantera ]; then
  git clone -b feature_differentiation \
    https://github.com/bigfooted/cantera.git \
    subprojects/cantera
else
  echo "Cantera already present."
fi

