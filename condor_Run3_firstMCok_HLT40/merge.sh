#!/bin/bash
# Merge per-job ROOT files into a single output using hadd.
#
# Usage:
#   ./merge.sh <output_hist_base> <btag|nobtag> <n_jobs>
#
# Example (bjet MC, btag, 20 jobs):
#   ./merge.sh template_for_fit_histos_3D_bjet_test btag 20

## output: Run3_WP90_template_for_fit_histos_3D_qcd_file_0_job0_btag.root
# Run3_WP90_template_for_fit_histos_3D_qcd_file_$1_job$3_$2
## Merge per dataset
OUTPUT_HIST=$1   # e.g. Run3_WP90_template_for_fit_histos_3D_qcd
LABEL=$2         # btag or nobtag
N_JOBS=$3	 # 100    

OUTPUT_FOLDER="/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3/JobResult"
MERGED_OUTPUT_FOLDER="/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3/MergedJobResult"
MERGED="${MERGED_OUTPUT_FOLDER}/MergedAllFiles_${OUTPUT_HIST}_${LABEL}.root"

FILES=""
for FILEINDEX in $(seq 0 9); do
	for i in $(seq 0 $((N_JOBS - 1))); do
  		f="${OUTPUT_FOLDER}/${OUTPUT_HIST}_file_${FILEINDEX}_job${i}_${LABEL}.root"
  		if [[ ! -f "$f" ]]; then
    			echo "WARNING: missing $f"
  		else
    			FILES="$FILES $f"
  		fi
	done
done


if [[ -z "$FILES" ]]; then
  echo "No job output files found, aborting."
  exit 1
fi

echo "Merging into ${MERGED} ..."
hadd -f "${MERGED}" ${FILES}
echo "Done: ${MERGED}"
