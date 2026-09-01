#!/bin/bash
set -u


SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=${SCRIPT_DIR}
source "${WORK}/setup_roounfold_env.sh"


# ============================ SET THE SAMPLE HERE ============================
# Two flags; input dir, output dir, chunk filenames and the macro's sample tag all
# follow from them.
SAMPLE=qcd          # qcd | bjet
GENERATOR=herwig    # pythia | herwig

# b-tag production of the chunks. Data runs over negTagFix, so use negTagFix here too;
# negTag is the older production. Only Pythia8 QCD has negTag/negTagFix chunks.
INPUT_TAG=UParTV2   # UParTV2 | negTagFix | negTag

case "${SAMPLE}_${GENERATOR}" in
  qcd_pythia)   IN_SUBDIR="QCD/Pythia8_${INPUT_TAG}_chunks";  OUT_SUBDIR=QCD        ;;
  bjet_pythia)  IN_SUBDIR="bJet/Pythia8_${INPUT_TAG}_chunks"; OUT_SUBDIR=bJet       ;;
  qcd_herwig)   IN_SUBDIR="QCDHerwig/${INPUT_TAG}_chunks";    OUT_SUBDIR=QCDHerwig  ;;
  bjet_herwig)  IN_SUBDIR="bJetHerwig/${INPUT_TAG}_chunks";   OUT_SUBDIR=bJetHerwig ;;
  *) echo "unknown SAMPLE/GENERATOR '${SAMPLE}/${GENERATOR}': use qcd|bjet and pythia|herwig"; exit 1 ;;
esac

# Pythia8 chunks carry the generator in the filename, the Herwig ones do not:
#   merged_block_0000_Pythia8_UParTV2.root  vs  merged_block_0000_UParTV2.root
FILE_TAG=$([ "${GENERATOR}" = pythia ] && echo "Pythia8_${INPUT_TAG}" || echo "${INPUT_TAG}")

INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/${IN_SUBDIR}
OUT_BASE=$mydata/bJetAggRun3/PPRef2024/${OUT_SUBDIR}/agg_ntuple_chunks


SAMPLE_TAG=${SAMPLE}

OUT_TAG=upartv2
OUT_SUFFIX="${OUT_TAG:+_${OUT_TAG}}"

case "${INPUT_TAG}" in
  UParTV2)   BTAG_WP=0.712 ;;
  negTagFix) BTAG_WP=0.868 ;;
  negTag)    BTAG_WP=0.868 ;;
  *) echo "unknown INPUT_TAG '${INPUT_TAG}': set BTAG_WP for it explicitly"; exit 1 ;;
esac

echo "sample ${SAMPLE} ${GENERATOR} (${INPUT_TAG}), b-tag WP ${BTAG_WP}, sample tag ${SAMPLE_TAG}"
echo "input  ${INPUT_DIR}"
echo "output ${OUT_BASE}"
# =============================================================================


LOG_DIR=${OUT_BASE}/logs

ACLIC_BUILD_DIR=$mydata/bJetAggRun3/PPRef2024/build

mkdir -p "${LOG_DIR}"

cd "${WORK}" || exit 1

COMPILE_LOG="${LOG_DIR}/compile.log"
echo "compiling create_files_for_template_fit.cpp once before launching parallel jobs"

# Drop half-written ACLiC products from a previous parallel launch.
rm -f "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp.d \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp.so \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp_ACLiC_dict.* \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp_ACLiC_map.*

## compilation of create_files_for_template_fit.cpp
root -l -b <<EOF
.L create_files_for_template_fit.cpp++
.q
EOF
# Check compilation status 
compile_status=$?
if [ ${compile_status} -ne 0 ]; then
  echo "compile failed; see ${COMPILE_LOG}"
  tail -n 80 "${COMPILE_LOG}"
  exit ${compile_status}
fi

echo "compile finished; launching chunk jobs"



shopt -s nullglob
inputs=("${INPUT_DIR}"/merged_block_*_${FILE_TAG}.root)
shopt -u nullglob

if [ ${#inputs[@]} -eq 0 ]; then
  echo "no chunks matching ${INPUT_DIR}/merged_block_*_${FILE_TAG}.root"
  exit 1
fi
echo "found ${#inputs[@]} chunks"


for input in "${inputs[@]}"; do

  base=$(basename "${input}")
  block=${base#merged_block_}
  block=${block%%_*}

  outdir="${OUT_BASE}/block_${block}"
  stagedir="${outdir}/.stage${OUT_SUFFIX}"
  rm -rf "${stagedir}"
  mkdir -p "${stagedir}"

  (
  
  # Remeber: create_files_for_template_fit(Int_t RunN = 3, Float_t pT_low = 80, Float_t etaCut = 2, Int_t n = 1,bool btag = true, bool isMC = true, Double_t btagWP = 0.712, bool makeTemplates = true, bool createRmatrix = true, bool makeAggNtuple = true, Long64_t ev_first = 0, Long64_t ev_last = -1, const char* inputFileOverride = "", const char* outputFolderOverride = "", const char* sampleTag = "")
  # MC vs data is the isMC argument; bjet vs qcd is SAMPLE_TAG (last argument).

  nice -n 10 root -l -b -q -e "gSystem->Load(\"${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,80,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"${input}\",\"${stagedir}\",\"${SAMPLE_TAG}\")"

    for f in "${stagedir}"/*.root; do
      [ -e "${f}" ] || continue
      base=$(basename "${f}")
      mv "${f}" "${outdir}/${base%.root}${OUT_SUFFIX}.root"
    done
    rmdir "${stagedir}" 2>/dev/null
  ) > "${LOG_DIR}/block_${block}${OUT_SUFFIX}.log" 2>&1 &


  echo "submitted block ${block}, pid $!"
done

wait
echo "all chunk jobs finished"
