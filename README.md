# eec_2b_analysis

EEC (energy-energy correlator) between the two b-jets in double-b events, PPRef2024.

**`workflow/` is the live code.** Everything at the top level and in `Run3/`,
`TemplateFit_Run3/`, `Unfolding/` is an older copy — read those, do not run them.

Run everything on LLR (`ssh llruicms01.in2p3.fr`), from
`/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow`. See `CLAUDE.md` for the
sshfs mount and the LCG/Condor environment, and `workflow/readme_workflow.md` for the
detailed reference behind every step below.

## Setup

```bash
source setup_roounfold_env.sh          # LCG 106a + the private RooUnfold build
root -l -b -q -e ".L create_files_for_template_fit.cpp++"   # compile the .so once
```

## The chain

### 1. MC ntuples — `run_agg_ntuple_chunks.sh`

- Runs `create_files_for_template_fit.cpp` over the MC chunks, one background ROOT job per block, producing the response matrices, gen histograms and aggregated B-hadron ntuples.
- Edit `SAMPLE` (`qcd`|`bjet`) and `GENERATOR` (`pythia`|`herwig`) at the top, or override them on the command line.
- Takes ~1 h and prints nothing while it runs; never run two copies at once.

```bash
./run_agg_ntuple_chunks.sh
SAMPLE=qcd GENERATOR=herwig ./run_agg_ntuple_chunks.sh
TRACK_EFF_UNC=true ./run_agg_ntuple_chunks.sh     # tracking-efficiency systematic
```

### 2. Data ntuples — `make_hardprobes_condor_scripts.sh`

- Writes the 50 per-job Condor scripts that run the same macro over HardProbes data; usually nothing to edit.
- `condor.submit` is the submit file they are launched with.

```bash
./make_hardprobes_condor_scripts.sh
condor_submit condor.submit
condor_q <cluster> -af:h JobStatus RemoteHost    # 1 idle, 2 running, 5 held
```

### 3. Template fit — `template_fit.cpp`

- Fits the 2B signal and the effective background to data and writes the signal fraction per ΔR bin, plus the nominal and varied (`var0B_2`, `var0B_0`) fit files.

```bash
root -l -b -q 'template_fit.cpp("both","pythia")'
```

### 4. Unfolding — `apply_unfolding_2d.C`

- Sums the per-block files of steps 1–2, applies the signal fraction and the UParT SF to data, unfolds and applies the purity and efficiency corrections.
- **Matrix inversion is the default** (`unfoldBayes = false`), for both observables: it has no regularisation, so there is no per-observable iteration count to justify (Bayesian needed 7 for ΔR and 21–26 for B). Pass `true` for the Bayesian result.
- Nothing to edit — every knob is an argument: `(SAMPLE, UNFOLDING_GENERATOR, test_mode, unfoldBayes, scan_niter, TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION, SFUPART_VARIATION, EEC_WEIGHT_OFF, OBSERVABLE, NITER)`; `test_mode` is `0` full-MC closure, `1` split test, `2` data; `OBSERVABLE` is `dr` or `B` (the momentum balance). The three b/non-b fraction observables (`fnb`, `fb`, `lnfb`) were commented out on 2026-09-22; see `readme_workflow.md`.
- Any observable other than `dr` forces `SFUPART_VARIATION = "off"` (the UParT SF is binned in reco ΔR and has no equivalent) and `EEC_WEIGHT_OFF = true` (it is a yield: dN/dB). Both rules live in `result_paths.h`, one function each.
- One results folder per argument combination, under `/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/`. **Matrix-inversion results go in their own subdirectory**, `results/matrix_inversion/unfolding_*`, so a whole MI production — every variation, both observables, plots and logs — is one directory to sync or compare. Bayesian keeps the flat layout it always had, so nothing that exists today moved.

```bash
root -l -b -q 'apply_unfolding_2d.C("both","pythia")'                  # nominal
root -l -b -q 'apply_unfolding_2d.C("both","herwig")'                  # detector response
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",true)'               # tracking eff
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_2")'    # mistag up
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_0")'    # mistag down
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","jpcalib_hf")'  # UParT SF
root -l -b -q 'apply_unfolding_2d.C("both","pythia",0)'                # closure test
# the momentum balance (dN/dB): nominal, then its two mistag variations
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","off",true,"B")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_2","off",true,"B")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_0","off",true,"B")'
```

### 5. Systematics and final plots — `apply_weights_and_systematics.C`

- Reads one step-4 result per booked variation, builds the systematic band, and writes the final ROOT file and the two plots you show; a missing variation is reported and skipped.
- Outputs go into the nominal folder of that observable: `final_*_with_systematics.root`, `data_vs_gen_*`, `systematics_curves_*`, `sfupart_before_after_*`.
- Run it **once per observable**: `OBSERVABLE` is the last argument and picks both the results it reads and the variation list. ΔR gets the full set (template fit, detector response, tracking, mistag pair, UParT SF stat + calibration); B gets the mistag pair, because the Herwig and track-drop MC have no B axis in them yet — the other entries are in the file, commented, with what each needs.
- To add a systematic: one entry in the `variations` vector, and its tag in `result_paths.h`.

```bash
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'   # EEC(ΔR)
root -l -b -q 'apply_weights_and_systematics.C("both","pythia",2,false,"h_data_fully_corrected_2D","pythia",2,false,"B")'   # dN/dB
```

## Other macros

| File | What it does |
|---|---|
| `plot_purity_efficiency_response.cpp` | Plots the response matrix, purity and efficiency from the step-1 files. Also holds `momentum_balance_mc_study()`, the purity/stability study that chose the B binning. |
| `template_fit.cpp` → `plot_z_first_look()` | The pre-fit step for any observable (last argument): flavour templates, m_2B separation per bin, purity/efficiency, data vs MC, and the raw yields. |
| `apply_weights_and_systematics.C` → `plot_systematics_summary()`, `plot_final_generator_band()` | Plots every Herwig variation as a relative shift against the Pythia nominal (needs both nominals from step 5). |
| `plotNice_UParT_roc.C` / `plotNice_bdt_roc.C` | ROC curves for the UParT b-tagger and for the BDT. |
| `merge_RMatrixTH2D_qcd.C` / `_bjet.C` / `_bjet_qcd_merged.C` | Merge per-block response matrices into one file (step 4 sums in memory, so this is only for inspection). |
| `Help_Functions.h` → `compare_root_contents()` | Bin-by-bin comparison of two ROOT files, to prove a rerun reproduces exactly. |
| `sync_plots_to_eos.sh` | Copies plots from `/data_CMS` to the CERNBox website (needs `kinit zzaidanc@CERN.CH`). |

## Headers

| File | What it holds |
|---|---|
| `result_paths.h` | The one definition of every result path and variation tag; add a variation here and both macros follow. |
| `binning_histos_small.h` | The ΔR, jet-pT and mass binning. |
| `central_selections.h` | The jet and event selections. |
| `tracking_efficiency_syst.h` | The hashed 3% track drop used by the tracking systematic. |
| `tTree.h` | The input ntuple branches. |
| `Help_Functions.h`, `Draw_EEC.h` | Plotting helpers and the EEC drawing used by the template fit. |
| `rootlogon.C` | Sends ACLiC build products to `/data_CMS`; loaded automatically by `root` in this folder. |
