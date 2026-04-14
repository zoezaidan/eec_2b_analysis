#!/bin/bash
source config.sh
# Merge per-job ROOT files into a single output using hadd.
#
# Usage:
#   ./merge.sh <output_hist_base> <btag|nobtag> <n_jobs>
#
# Example (bjet MC, btag, 20 jobs):
#   ./merge.sh template_for_fit_histos_3D_bjet_test btag 20


## arguments added to config.sh
## can be used for fasts test (when njobs differs from the total output)
#OUTPUT_HIST=$1   # e.g. template_for_fit_histos_3D_bjet_test
#LABEL=$2         # btag or nobtag
#N_JOBS=$3


OUTPUT_FOLDER="/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/condor/condor_output"
MERGED="${OUTPUT_FOLDER}/${OUTPUT_HIST}_${LABEL}.root"

FILES=""
for i in $(seq 0 $((N_JOBS - 1))); do
  f="${OUTPUT_FOLDER}/${OUTPUT_HIST}_${TRIGGER_LABEL}_${LABEL}_job${i}.root"
  if [[ ! -f "$f" ]]; then
    echo "WARNING: missing $f"
  else
    FILES="$FILES $f"
  fi
done

if [[ -z "$FILES" ]]; then
  echo "No job output files found, aborting."
  exit 1
fi

echo "Merging into ${MERGED} ..."
hadd -f "${MERGED}" ${FILES}
echo "Done: ${MERGED}"
