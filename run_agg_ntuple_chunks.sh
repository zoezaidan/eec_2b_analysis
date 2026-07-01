#!/bin/bash

set -u

WORK=/home/llr/cms/mnguyen/eec_2b_analysis
INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_recalJP_chunks_20260617
OUT_BASE=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks_20260701
LOG_DIR=${OUT_BASE}/logs

mkdir -p "${LOG_DIR}"

cd "${WORK}" || exit 1

for i in $(seq 0 9); do
  block=$(printf "000%d" "${i}")
  input="${INPUT_DIR}/merged_block_${block}_Pythia8_recalJP_20260617.root"
  outdir="${OUT_BASE}/block_${block}"
  mkdir -p "${outdir}"

  root -l -b -q "create_files_for_template_fit.cpp+(3,2,80,200,2,1,true,true,0.868,false,false,true,0,-1,\"${input}\",\"${outdir}\")" \
    > "${LOG_DIR}/block_${block}.log" 2>&1 &

  echo "submitted block ${block}, pid $!"
done

wait
echo "all chunk jobs finished"
