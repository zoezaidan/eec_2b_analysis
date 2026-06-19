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

CONDORDIR=/grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/test_condor_oldroot
MACRODIR=/grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis

mkdir -p ${CONDORDIR}/exec_jobs
mkdir -p ${CONDORDIR}/logfiles


echo "Generating ${N_JOBS} job scripts..."

for i in $(seq 0 $((N_JOBS - 1))); do
  SCRIPT=${CONDORDIR}/exec_jobs/job_${i}.sh
  cat > ${SCRIPT} << EOF
#!/bin/bash

# ===============================
# CMSSW + ROOT ENV (inside Singularity)
# ===============================

source /cvmfs/cms.cern.ch/cmsset_default.sh	
	# Go to your CMSSW release area
	#cd /home/llr/cms/shatat/CMSAnalysis/CMSSW_10_6_48/src
export SCRAM_ARCH=slc7_amd64_gcc700
# Create CMSSW release locally on worker node
scram project CMSSW CMSSW_10_6_48
cd CMSSW_10_6_48/src


# Set runtime environment
# Activate CMSSW release
eval \`scram runtime -sh\`

echo "=== CMSSW environment loaded for job ${i} ==="
echo "===== Environment checks ====="
echo "SCRAM_ARCH = \$SCRAM_ARCH"
echo "ROOT executable:"
which root

echo "ROOT version:"
root-config --version
echo "=============================="


cd ${CONDORDIR}


# Run analysis macro (compiled on the fly per job)
root -b -l -q "run_condor_job.C+(${i}, ${N_JOBS}, ${DATATYPE}, ${PT_LOW}, ${PT_HIGH}, ${N}, ${BTAG}, ${ISMC})"

EOF

  chmod +x ${SCRIPT}
done

echo "Done. Updating queue count in submit.sub..."

sed -i "s/^queue .*/queue ${N_JOBS}/" ${CONDORDIR}/submit.sub

echo "Submitting jobs..."
condor_submit ${CONDORDIR}/submit.sub

