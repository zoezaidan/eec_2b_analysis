### instead of  #!/bin/bash use the next one when run condor  #!/usr/bin/env bash 

#!/bin/bash
source config.sh # instead of add argumnets at runtime 


# Generates one shell script per job in exec_jobs/, then submits.
#
# Usage:
#   ./generate_jobs.sh <n_jobs> <dataType> <pT_low> <pT_high> <n> <btag> <isMC>
#
# Example (HighEG data, btag, 100 jobs):
#   ./generate_jobs.sh 100 0 80 140 1 1 0

# path to the run code for each job 
CONDORDIR=/grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/condor/Run_bjet_HLT80
# path for run_condor_job.C
MACRODIR=/grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/condor/common_condor_header 
 


# make directories for the work 
mkdir -p ${CONDORDIR}/exec_jobs
mkdir -p ${CONDORDIR}/logfiles

# Pre-compile the shared library ONCE here before submitting.
# This way all jobs reuse the cached .so instead of 100 jobs
# simultaneously trying to compile and filling the disk.
echo "Pre-compiling run_condor_job.C (this may take a minute)..."
unset CMSSW_BASE CMSSW_VERSION CMSSW_RELEASE_BASE ROOT_INCLUDE_PATH

# info about root 
which root
root --version

# cd ${CONDORDIR} # important 
cd ${MACRODIR} 
root -b -q -e 'gSystem->CompileMacro("run_condor_job.C","kf"); exit(0);' # to compile only once and get .so compilation file (used for all jobs if exist), kf:k keep compilation, f: force frech compilation. exit(0): to exit root cleanly.
if [[ $? -ne 0 ]]; then # check exist status from previous command. exit 0 is success, 1 is fail.
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
root -b -q -l "${MACRODIR}/run_condor_job.C+(${i}, ${N_JOBS}, ${DATATYPE}, ${PT_LOW}, ${PT_HIGH}, ${N}, ${BTAG}, ${ISMC}, ${HLT80})"
EOF
  chmod +x ${SCRIPT}
done


echo "Done. Updating queue count in submit.sub to ${N_JOBS}..."
# update submit.sub file directly: find where line start by "queue" and replace by  queue ${N_JOBS}  -> this is needed for condor_submit how many#jobs
sed -i "s/^queue .*/queue ${N_JOBS}/" ${CONDORDIR}/submit.sub 

echo "Submitting..."
# launch the jobs at condor ! 
condor_submit submit.sub
