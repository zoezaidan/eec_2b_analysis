#!/bin/bash
# Merge per-job ROOT files into a single output using hadd.
#
# Usage:
#   ./merge.sh <output_hist_base> <btag|nobtag> <n_jobs>
#
# Example (bjet MC, btag, 20 jobs):
#   ./merge.sh template_for_fit_histos_3D_bjet_test btag 20

OUTPUT_HIST=$1   # e.g. template_for_fit_histos_3D_bjet_test
LABEL=$2         # btag or nobtag
N_JOBS=$3

if [[ -z "$OUTPUT_HIST" || -z "$LABEL" || -z "$N_JOBS" ]]; then
  echo "Usage: $0 <output_hist_base> <btag|nobtag> <n_jobs>"
  exit 1
fi

OUTPUT_FOLDER="/data_CMS/cms/zaidan/analysis_lise"
MERGED="${OUTPUT_FOLDER}/${OUTPUT_HIST}_${LABEL}.root"

FILES=""
for i in $(seq 0 $((N_JOBS - 1))); do
  f="${OUTPUT_FOLDER}/${OUTPUT_HIST}_${LABEL}_job${i}.root"
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
