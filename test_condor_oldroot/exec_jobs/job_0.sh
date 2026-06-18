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
eval `scram runtime -sh`

echo "=== CMSSW environment loaded for job 0 ==="
echo "===== Environment checks ====="
echo "SCRAM_ARCH = $SCRAM_ARCH"
echo "ROOT executable:"
which root

echo "ROOT version:"
root-config --version
echo "=============================="


cd /grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/test_condor_oldroot


# Run analysis macro (compiled on the fly per job)
root -b -l -q "run_condor_job.C+(0, 1, 0, 80, 200, 1, 1, 0)"

