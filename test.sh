cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis

unset ROOT_INCLUDE_PATH
unset CPLUS_INCLUDE_PATH
unset CPATH
unset LD_LIBRARY_PATH
unset PYTHONPATH

set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a/x86_64-el9-gcc14-opt/setup.sh
set -u

which gcc
which root

root -l -b -q 'compile_create_files_roounfold.C("RooUnfold_build_test/src/src","RooUnfold_build_test/build")'
