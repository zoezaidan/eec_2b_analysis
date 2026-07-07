#!/bin/bash

set -u

WORK=./
INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_recalJP_chunks_20260617
OUT_BASE=$mydata/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks_20260701
LOG_DIR=${OUT_BASE}/logs

mkdir -p "${LOG_DIR}"

cd "${WORK}" || exit 1

COMPILE_LOG="${LOG_DIR}/compile.log"
echo "compiling create_files_for_template_fit.cpp once before launching parallel jobs"

# Avoid stale or half-written ACLiC products from a previous parallel launch.
rm -f create_files_for_template_fit_cpp.d \
      create_files_for_template_fit_cpp.so \
      create_files_for_template_fit_cpp_ACLiC_dict.* \
      create_files_for_template_fit_cpp_ACLiC_map.*

root -l -b -q -e '.L create_files_for_template_fit.cpp++' > "${COMPILE_LOG}" 2>&1
compile_status=$?
if [ ${compile_status} -ne 0 ]; then
  echo "compile failed; see ${COMPILE_LOG}"
  tail -n 80 "${COMPILE_LOG}"
  exit ${compile_status}
fi

echo "compile finished; launching chunk jobs"

for i in $(seq 0 9); do
  block=$(printf "000%d" "${i}")
  input="${INPUT_DIR}/merged_block_${block}_Pythia8_recalJP_20260617.root"
  outdir="${OUT_BASE}/block_${block}"
  mkdir -p "${outdir}"

  nice -n 10 root -l -b -q -e "gSystem->Load(\"${WORK}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,2,80,200,2,1,true,true,0.868,true,true,true,0,-1,\"${input}\",\"${outdir}\")" \
       > "${LOG_DIR}/block_${block}.log" 2>&1 &

  echo "submitted block ${block}, pid $!"
done

wait
echo "all chunk jobs finished"
