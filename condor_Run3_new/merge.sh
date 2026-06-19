#!/bin/bash
# Merge per-job ROOT files into a single output using hadd.
#
# Usage:
#   ./merge.sh <output_hist_base> <btag|nobtag> <n_jobs>
#
# Example (qcd MC, btag, 100 jobs):
#   ./merge.sh Run3_secondbinsplitting_June_WP090_template_for_fit_histos_3D_qcd btag 100
# Example for Run3 data:
#   ./merge.sh Run3_secondbinsplitting_June_WP090_template_for_fit_histos_3D_alltrgData btag 100 

## output: Run3_secondbinsplitting_June_WP090_template_for_fit_histos_3D_qcd_file_9_btag_job99.root
# Run3_secondbinsplitting_June_WP090_template_for_fit_histos_3D_qcd_file_fileindex_$2_job$3
## Merge per dataset
OUTPUT_HIST=$1   # e.g. Run3_WP90_template_for_fit_histos_3D_qcd
LABEL=$2         # btag or nobtag
N_JOBS=$3	 # 100    

OUTPUT_FOLDER="/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3_new/JobResult"
MERGED_OUTPUT_FOLDER="/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3_new/MergedJobResult"
MERGED="${MERGED_OUTPUT_FOLDER}/MergedAllFiles_${OUTPUT_HIST}_${LABEL}.root"

FILES=""

#for FILEINDEX in $(seq 0 9); do
for FILEINDEX in $(seq 0 3); do
	for i in $(seq 0 $((N_JOBS - 1))); do
  		f="${OUTPUT_FOLDER}/${OUTPUT_HIST}_file_${FILEINDEX}_${LABEL}_job${i}.root"
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
