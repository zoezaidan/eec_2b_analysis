#!/bin/bash
# Generates one shell script per job in exec_jobs/, then submits.
#
# Usage:
#   ./generate_jobs.sh <n_jobs> <dataType> <pT_low> <pT_high> <n> <btag> <isMC>
#
# Example (HighEG data, btag, 100 jobs):
#   ./generate_jobs.sh 100 0 80 140 1 1 0 [file index from 0-9]

# to choose MC run 3: ./generate_jobs.sh 100 2 80 200 1 1 0 0
### To run this script: condor_submit submit.sub
### MAke sure this script is chmod +x generate_jobs.sh

JOBID=$1
N_JOBS=$2
DATATYPE=$3
PT_LOW=$4
PT_HIGH=$5
N=$6
BTAG=$7
ISMC=$8
# Only for Run3 MC 
FILEINDEX=$9 


if [[ -z "$N_JOBS" || -z "$DATATYPE" ]]; then
  echo "Usage: $0 <n_jobs> <dataType> <pT_low> <pT_high> <n> <btag> <isMC> <FILEINDEX>"
  echo "  dataType: -1=LowEG data  0=HighEG data  1=bjet MC  2=dijet MC"
  exit 1
fi

# shatat path 
CONDORDIR=/grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3

mkdir -p ${CONDORDIR}/logfiles

cd ${CONDORDIR}
unset CMSSW_BASE CMSSW_VERSION CMSSW_RELEASE_BASE ROOT_INCLUDE_PATH
## Try compile first
#root -b -q -e 'gSystem->CompileMacro("run_condor_job.C","kf"); exit(0);'
## If use precompile --> use c+ not c++
root -b -q -l "run_condor_job_Run3.C++(${JOBID}, ${N_JOBS}, ${DATATYPE}, ${PT_LOW}, ${PT_HIGH}, ${N}, ${BTAG}, ${ISMC}, ${FILEINDEX})"

