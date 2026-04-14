#!/bin/bash
# Unset CMSSW so ACLiC doesn't repick up CMSSW headers
unset CMSSW_BASE CMSSW_VERSION CMSSW_RELEASE_BASE ROOT_INCLUDE_PATH
cd /grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/condor/Run_data_HLT80
# Use .C+ (not .C++) so ROOT reuses the already-compiled .so
root -b -q -l "/grid_mnt/vol_home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/condor/common_condor_header/run_condor_job.C+(79, 100, 0, 80, 140, 1, 1, 0, 1)"
