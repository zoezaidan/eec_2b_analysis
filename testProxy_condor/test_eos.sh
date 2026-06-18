#!/bin/bash

## test line 
export X509_USER_PROXY=~/.t3/proxy.cert

echo "=== Identity ==="
whoami

echo "=== EOS test ==="
eos root://eos.grif.fr ls /eos/grif/cms/llr/store/user/mnguyen/bJetAggRun3/QCD_pThat-15to1200_TuneCP5_5p36TeV_pythia8/bJetAgg_2024PPRef_QCD/0/
