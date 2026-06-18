#!/bin/bash
# Generates one shell script per job in exec_jobs/, then submits.
#
# Usage:
#   ./generate_jobs.sh <n_jobs> <dataType> <pT_low> <pT_high> <n> <btag> <isMC>
#
# Example (HighEG data, btag, 100 jobs):
#   ./generate_jobs.sh 100 0 80 140 1 1 0

N_JOBS=$1
DATATYPE=$2
PT_LOW=$3
PT_HIGH=$4
N=$5
BTAG=$6
ISMC=$7

if [[ -z "$N_JOBS" || -z "$DATATYPE" ]]; then
  echo "Usage: $0 <n_jobs> <dataType> <pT_low> <pT_high> <n> <btag> <isMC>"
  echo "  dataType: -1=LowEG data  0=HighEG data  1=bjet MC  2=dijet MC"
  exit 1
fi

CONDORDIR=/grid_mnt/vol_home/llr/cms/zaidan/analysis_lise/updates_2002/condor
MACRODIR=/grid_mnt/vol_home/llr/cms/zaidan/analysis_lise/updates_2002

mkdir -p ${CONDORDIR}/exec_jobs
mkdir -p ${CONDORDIR}/logfiles

# Pre-compile the shared library ONCE here before submitting.
# This way all jobs reuse the cached .so instead of 100 jobs
# simultaneously trying to compile and filling the disk.
echo "Pre-compiling run_condor_job.C (this may take a minute)..."
unset CMSSW_BASE CMSSW_VERSION CMSSW_RELEASE_BASE ROOT_INCLUDE_PATH
cd ${CONDORDIR}
root -b -q -e 'gSystem->CompileMacro("run_condor_job.C","kf"); exit(0);'
if [[ $? -ne 0 ]]; then
  echo "ERROR: pre-compilation failed, aborting."
  exit 1
fi
echo "Pre-compilation done."

cd ${CONDORDIR}
echo "Generating ${N_JOBS} job scripts..."

for i in $(seq 0 $((N_JOBS - 1))); do
  SCRIPT=${CONDORDIR}/exec_jobs/job_${i}.sh
  cat > ${SCRIPT} << EOF
#!/bin/bash
# Unset CMSSW so ACLiC doesn't repick up CMSSW headers
unset CMSSW_BASE CMSSW_VERSION CMSSW_RELEASE_BASE ROOT_INCLUDE_PATH
cd ${CONDORDIR}
# Use .C+ (not .C++) so ROOT reuses the already-compiled .so
root -b -q -l "run_condor_job.C+(${i}, ${N_JOBS}, ${DATATYPE}, ${PT_LOW}, ${PT_HIGH}, ${N}, ${BTAG}, ${ISMC})"
EOF
  chmod +x ${SCRIPT}
done

echo "Done. Updating queue count in submit.sub to ${N_JOBS}..."
sed -i "s/^queue .*/queue ${N_JOBS}/" ${CONDORDIR}/submit.sub

echo "Submitting..."
condor_submit submit.sub
