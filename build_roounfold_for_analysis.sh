#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOUNFOLD_SRC_REF=${ROOUNFOLD_SRC_REF:-/home/llr/cms/kalipoliti/gitRepos/RooUnfoldOld}
ROOUNFOLD_DIR=${SCRIPT_DIR}/RooUnfold_build_test

if [ ! -d "${ROOUNFOLD_SRC_REF}" ]; then
  echo "ERROR: reference RooUnfold source not found: ${ROOUNFOLD_SRC_REF}"
  echo "Set ROOUNFOLD_SRC_REF=/path/to/RooUnfold source and rerun."
  exit 1
fi

unset ROOT_INCLUDE_PATH
unset CPLUS_INCLUDE_PATH
unset CPATH
# Avoid mixing an already-loaded CMSSW runtime with the LCG view.
unset LD_LIBRARY_PATH
unset PYTHONPATH

# LCG setup scripts are not nounset-safe.
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a/x86_64-el9-gcc14-opt/setup.sh
set -u

mkdir -p "${ROOUNFOLD_DIR}"
if [ ! -d "${ROOUNFOLD_DIR}/src" ]; then
  echo "Copying RooUnfold source from ${ROOUNFOLD_SRC_REF}"
  rsync -a --exclude build "${ROOUNFOLD_SRC_REF}/" "${ROOUNFOLD_DIR}/src/"
fi

mkdir -p "${ROOUNFOLD_DIR}/build"
cd "${ROOUNFOLD_DIR}/build"
cmake ../src
cmake --build . --target RooUnfold -j4

missing=$(ldd libRooUnfold.so | grep 'not found' || true)
if [ -n "${missing}" ]; then
  echo "ERROR: unresolved RooUnfold dependencies:"
  echo "${missing}"
  exit 1
fi

echo "Built RooUnfold successfully: ${ROOUNFOLD_DIR}/build/libRooUnfold.so"
echo "Use ROOUNFOLD_DIR=${ROOUNFOLD_DIR} in scripts, or run scripts from this analysis directory."
