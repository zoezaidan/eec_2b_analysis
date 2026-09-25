#!/bin/bash
# Publish the momentum-balance (pT balance) plots to the CERNBox website, under
# https://zzaidanc.web.cern.ch/pt_balance/
#
# Needs a Kerberos ticket FIRST -- this cannot be done from an automated session:
#     kinit zzaidanc@CERN.CH
#     ./sync_pt_balance_to_eos.sh
#
# Unlike sync_plots_to_eos.sh this copies ONLY .png/.pdf. The template-fit folders also
# contain multi-MB .root files, and a web area is not where they belong.
set -u

CERN_USER=zzaidanc
DEST_ROOT=/eos/user/${CERN_USER:0:1}/${CERN_USER}/www/pt_balance
R=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results

# <source dir>:<subfolder on the website>
#
# ⚠️ The observable was renamed z -> B on 2026-09-22, so these source folders are the NEW
# "_B" ones and do not exist until the chain has been re-run end to end (step 1 into
# OUT_TAG=upartv2_B, then the fit and the unfolding). Until then every entry below is
# "skip (missing)" and the website keeps whatever the "_z" run last published. The old
# "_z" folders are deliberately not listed -- publishing both would put two spellings of
# the same measurement side by side on the site.
MAP=(
  "$R/B_first_look_both_pythia_noeecw_upartv2_B:dNdB"
  "$R/B_first_look_both_pythia_upartv2_B:eec_weighted"
  "$R/TemplateFit_Run3/TemplateFits_both_pythia_noeecw_B_upartv2:templatefit_dNdB"
  "$R/TemplateFit_Run3/TemplateFits_both_pythia_B_upartv2:templatefit_eec_weighted"
  "$R/momentum_balance_study_both_pythia_upartv2:binning_study"
  # The unfolded dN/dB with its band -- data_vs_gen, systematics_curves. Matrix inversion is
  # the unfolding of record, hence the matrix_inversion/ prefix; the folder carries
  # "_sfupartoff" because no UParT SF exists in B.
  "$R/matrix_inversion/unfolding_both_pythia_B_sfupartoff_noeecw_upartv2:unfolded_dNdB"
  "$R/B_first_look_qcd_pythia_noeecw_upartv2_B:dNdB_qcd_only"
  "$R/B_first_look_bjet_pythia_noeecw_upartv2_B:dNdB_bjet_only"
)

if ! klist -s 2>/dev/null; then
  echo "ERROR: no Kerberos ticket. Run:  kinit ${CERN_USER}@CERN.CH"
  exit 1
fi

for entry in "${MAP[@]}"; do
  SRC="${entry%%:*}"; SUB="${entry##*:}"
  if [ ! -d "$SRC" ]; then echo "skip (missing): $SRC"; continue; fi
  n=$(ls "$SRC"/*.png "$SRC"/*.pdf 2>/dev/null | wc -l)
  if [ "$n" -eq 0 ]; then echo "skip (no images): $SRC"; continue; fi
  echo "--> $SUB  ($n images)"
  ssh "${CERN_USER}@lxplus.cern.ch" "mkdir -p ${DEST_ROOT}/${SUB}" || exit 1
  # --no-perms/--no-group: EOS fuse rejects rsync's ownership calls.
  rsync -rltvz --no-perms --no-group \
        --include='*/' --include='*.png' --include='*.pdf' --exclude='*' \
        "$SRC/" "${CERN_USER}@lxplus.cern.ch:${DEST_ROOT}/${SUB}/" || exit 1
done

echo
echo "done: https://${CERN_USER}.web.cern.ch/pt_balance/"
echo "(directory listing is off by default -- copy CERN's index.php into each folder to browse as a gallery)"
