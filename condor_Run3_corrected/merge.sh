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


## Output.root 
#DATA:  Run3_secondbinsplitting_WP0872_template_for_fit_histos_3D_data _(fileindex) _$2 _job$3.root
#MC:    Run3_secondbinsplitting_June_WP0872_template_for_fit_histos_3D_qcd _(fileindex) _$2 _job$3.root

## Merge per dataset
OUTPUT_HIST=$1   # e.g. Run3_secondbinsplitting_WP0872_template_for_fit_histos_3D_data or Run3_secondbinsplitting_WP0872_template_for_fit_histos_3D_qcd
LABEL=$2         # btag or nobtag
N_JOBS=$3	 			 # 100    

OUTPUT_FOLDER="/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3_corrected/JobResult"
MERGED_OUTPUT_FOLDER="/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3_corrected/MergedJobResult"
MERGED="${MERGED_OUTPUT_FOLDER}/MergedAllFiles_${OUTPUT_HIST}_${LABEL}.root"

FILES=""

# For qcdMC
#for FILEINDEX in $(seq 0 9); do
# For data : 0 - 4
for FILEINDEX in $(seq 0 9); do
	for i in $(seq 0 $((N_JOBS - 1))); do
  		f="${OUTPUT_FOLDER}/${OUTPUT_HIST}_${FILEINDEX}_${LABEL}_job${i}.root"
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
