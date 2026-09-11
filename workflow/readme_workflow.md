# workflow — running the EEC 2b chain
mount; `/Users/zoezaidan/llruicms01/...` == `/home/llr/cms/zaidan/...`.
Remote work dir: `/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow`

## Environment

```bash
source setup_roounfold_env.sh   # LCG 106a + the private RooUnfold build
```

## The chain

| Step | Script | Produces |
|---|---|---|
| 1. MC chunks | `run_agg_ntuple_chunks.sh` | per-block `RMatrix_*`, `*MCGEN*`, `AggBHadronNtuple_*` |
| 2. Data chunks | `make_hardprobes_condor_scripts.sh` + `condor_submit` | per-block data templates |
| 3. Unfolding | `apply_unfolding_2d.C` | unfolded EEC + closure/scan plots |
| 4. Systematics | `apply_weights_and_systematics.C` | final result + band, `data_vs_gen_*`, `systematics_curves_*` |

Steps 1 and 2 both run the same macro, `create_files_for_template_fit.cpp`. Step 3 reads
what they wrote, and step 4 reads one step-3 output per systematic variation.

---

## 1. MC — `run_agg_ntuple_chunks.sh`

**Change (usually only these two):**

| Line | Knob | Values |
|---|---|---|
| 13 | `SAMPLE` | `qcd` \| `bjet` |
| 14 | `GENERATOR` | `pythia` \| `herwig` |

Everything else follows: input dir, output dir, chunk filename pattern, the macro's sample
tag, and the list of blocks (taken from what is on disk — 10 qcd/pythia, 9 bjet/pythia,
8 qcd/herwig, 9 bjet/herwig).

Occasionally:

| Line | Knob | Note |
|---|---|---|
| 18 | `INPUT_TAG` | `UParTV2` \| `negTagFix` \| `negTag`. Also picks `BTAG_WP`. Only Pythia8 QCD has the `negTag*` chunks. |
| 38 | `OUT_TAG` | suffix on every output file; change it to avoid overwriting an earlier run |
| 56 | `TRACK_EFF_UNC` | `false` nominal, `true` the tracking-efficiency systematic — see below |

**Run:**

```bash
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && ./run_agg_ntuple_chunks.sh'
```

It prints the resolved config first — check it before walking away:

```
sample bjet herwig (UParTV2), b-tag WP 0.712, sample tag bjet
input  /data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/bJetHerwig/UParTV2_chunks
output /data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/bJetHerwig/agg_ntuple_chunks
found 9 chunks
```

It compiles the macro once, then runs one `nice`d ROOT job per block in parallel and waits.
Outputs land in `<OUT_BASE>/block_NNNN/`, logs in `<OUT_BASE>/logs/block_NNNN_<OUT_TAG>.log`.

**It takes ~50–60 min and looks hung the whole time.** The macro prints
`Processing events [0, N)` once and the staging files only grow at autoflush boundaries. To
check it is alive, watch the counters rise:

```bash
for p in $(pgrep -u zaidan -f "root.exe.*create_files_for_template_fit"); do
  awk '{print $14+$15}' /proc/$p/stat; grep ^read_bytes /proc/$p/io
done
```

Do **not** kill and restart (no resume — each block starts from event 0), and **never run
two copies at once**: it does `rm -rf` on the staging dirs and wipes the shared ACLiC build.

### Tracking-efficiency systematic (`TRK_DROP`)

Track reconstruction in the jet core is mismodelled in simulation. To propagate that,
set `TRACK_EFF_UNC=true`: 3% of the **reconstructed** tracks in MC are thrown away during the
B reconstruction (`makeSvtxs_withBDT`, before any vertex is built), so the SVs, m_2B, ΔR and the EEC weight are all built from a jet
that is genuinely missing tracks. The gen-level b hadrons (`refTrk*`) are untouched — the
variation has to move reco while truth stays put, which is what makes the response matrix
and the corrections change. Data is never varied (the macro guards on `isMC`).

**It is reproducible by construction.** The number that decides a track's fate is not
drawn from a stream, it is *hashed* from `(TRK_SEED, entry number, track index)`. So the
same track always gets the same verdict, no matter how many blocks the sample is split
into, in what order they run, or how many times you rerun. A seeded `std::mt19937`
consumed as a stream would **not** give this: it only reproduces if the tracks are
consumed in exactly the same order by exactly one process, which the parallel block jobs
violate. See `tracking_efficiency_syst.h` for the argument in full.

Verified on `block_0000` of Pythia8 QCD: two runs of the same varied job give bit-identical
histogram contents (46 + 13 histograms, max |Δ| = 0), and `[0,50k) + [50k,100k)` summed
equals a single `[0,100k)` run up to 1e-16 rounding.

Output files are tagged `_trkdrop030`, so a variation run cannot overwrite the nominal and
both can sit in the same `block_NNNN/` directory. Nominal filenames are unchanged.

**The full procedure:**

```bash
# 1. produce the varied MC — TRACK_EFF_UNC=true in run_agg_ntuple_chunks.sh, once per
#    SAMPLE/GENERATOR you use in the nominal (qcd and bjet if you unfold with "both")
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && ./run_agg_ntuple_chunks.sh'

# 2. unfold the SAME data with it (last argument), which lands in its own results folder
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",true)'

# 3. the band: the "tracking_eff" entry in apply_weights_and_systematics.C picks the
#    varied result up automatically and adds |nominal - varied| to the quadrature sum
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'
```

⚠️ The signal fraction from the template fit is **not** re-derived on the varied
templates. `template_fit.cpp` reads pre-merged MCGEN files (the Pythia ones from Afnan's
area, under a different naming scheme), so pointing it at a varied production is a
separate job. The template fit is already booked as its own systematic axis
(`TF_GENERATOR`); decide whether you also want it inside the tracking one.

### Light-jet mistagging systematic (`TF_VARIATION`)

The **0B template** is the jets that passed the b tag with no gen b hadron in them
(`jtNbHad == 0`) — mistagged light and charm jets. It is not a free component of the
template fit: the fit has two PDFs, 2B (signal) and an *effective* background, and 0B is
folded into that background at the ratio MC predicts,

```
c' = 0B / (0B + 1B)          effective bkg = (1-c') x 1B_shape + c' x 0B_shape
```

so nothing in the data constrains the mistag rate. It reaches the measurement only through
the **shape** of the effective background, and what moves is the fitted signal fraction.

The systematic rebuilds that shape with the 0B admixture scaled by ×2 and ×0 — a
deliberately conservative ±100% on the mistag rate, which brackets the measured
light-flavour mistag scale factors (typically 20–50% at a tight WP) by a wide margin and
needs no external SF input. `template_fit.cpp` already runs all of them: every fit
directory holds `nominal_`, `var0B_2_` and `var0B_0_` copies of
`Run3_TemplateFits_histos_3d_80_inf.root`, each with its own `h_sig_fraction_fit`.

**So there is no MC to produce and no refit — only two more unfoldings:**

```bash
# 1. the two varied signal fractions, each into its own results folder
#    (_tf0Bx2 / _tf0Bx0)
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"var0B_2")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"var0B_0")'

# 2. the band: the "mistag_0B_up" / "mistag_0B_down" entries pick them up automatically
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'
```

The two entries share `Variation::group = "mistag_0B"`, so the total takes their **per-bin
envelope once** instead of adding both in quadrature — they are two sides of one
uncertainty, not two sources. Any other up/down pair can be combined the same way by giving
it a shared group.

`TF_VARIATION` only makes sense with `test_mode = 2`: the signal fraction is applied to data
alone, so the macro refuses the closure modes rather than write a folder that reads as a
result and is not one.

**Size, measured on `both`/`pythia` (Pythia templates, UParT WP 0.712):** the shift on the
signal fraction is **≤ 2.4 %, and under 1 % in most bins** — smaller than the fit's own
statistical error everywhere except the first two ΔR bins. It is the smallest of the
systematics booked so far. The largest shifts sit at the two ends: the lowest-ΔR bins, where
the sample is background-dominated (signal fraction 0.12–0.14) so the background shape
matters most, and the highest-ΔR bins, where `c'` itself is largest (up to 0.18).

There is a guard in `template_fit.cpp`: `eff_bkg1B = 1 - eff_bkg0B`, so `c' > 0.5` would
make the ×2 variation build a background PDF with **negative** 1B weight. The normalisation
check further down would still pass — the two weights sum to 1 by construction — and
`RooHistPdf` would silently clip the negative bins. It cannot happen at the current WP
(`2c' ≤ 0.37`), but a looser b tag or coarser ΔR bins can reach it, so the value is now
clipped with a warning that says the variation is truncated in that bin.

⚠️ `c'` is computed from the **dijet (qcd) integrals only** — `h_bb`/`h_b`/`h_nob`, never
the `_bjet` ones — while the shapes being mixed are qcd+bjet when `also_bjet`. That is
deliberate (qcd is the sample whose flavour composition matches the data; the bjet sample is
there for template statistics and would bias `c'` downwards if it entered the ratio), but it
is worth knowing when reading the numbers.

## 2. Data — `make_hardprobes_condor_scripts.sh`

**Usually change nothing.** It only ever runs HardProbes data, so `SAMPLE_TAG=data`
(line 19) is fixed. Occasionally line 15 `INPUT_TAG` or line 22 `OUT_TAG`.

⚠️ The `BTAG_WP` case here is duplicated from `run_agg_ntuple_chunks.sh` — if you change a
working point, change it in **both** files or data and MC come out at different WPs.

**Run — two steps, and the `.so` must exist first:**

```bash
# 0. compile once (run_agg_ntuple_chunks.sh also does this); jobs exit 3 if it is missing
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && \
  source setup_roounfold_env.sh && root -l -b -q -e ".L create_files_for_template_fit.cpp++"'

# 1. generate the 50 job scripts
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && ./make_hardprobes_condor_scripts.sh'

# 2. submit
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && condor_submit condor.submit'

# 3. monitor  (JobStatus: 1 idle, 2 running, 5 held; gone from the queue = finished)
ssh llruicms01.in2p3.fr 'condor_q <cluster> -af:h JobStatus RemoteHost'
```

Read the job logs locally through the mount rather than over SSH.

## 3. Unfolding — `apply_unfolding_2d.C`

**Nothing to edit — everything is a macro argument:**

| # | Argument | Values | Default |
|---|---|---|---|
| 1 | `SAMPLE` | `qcd` \| `bjet` \| `both` | `qcd` |
| 2 | `GENERATOR` | `pythia` \| `herwig` | `pythia` |
| 3 | `test_mode` | `0` full-MC closure, `1` split test, `2` data | `2` |
| 4 | `unfoldBayes` | `true` Bayesian, `false` matrix inversion | `true` |
| 5 | `scan_niter` | `true` scan niter 1..100, pick the best by refolding GoF | `true` |

A wrong `SAMPLE`/`GENERATOR` string is rejected with a message instead of silently falling
back to the default.

**Run** (on LLR, after `source setup_roounfold_env.sh`):

```bash
root -l -b -q 'apply_unfolding_2d.C("qcd","pythia")'      # data, Bayesian, scan
root -l -b -q 'apply_unfolding_2d.C("both","herwig")'     # qcd+bjet Herwig
root -l -b -q 'apply_unfolding_2d.C("bjet","herwig",0)'   # full-MC closure
```

In one line from the Mac (note the nested quoting):

```bash
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && \
  source setup_roounfold_env.sh && root -l -b -q "apply_unfolding_2d.C(\"both\",\"herwig\")"'
```

Results go to `/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/unfolding_<SAMPLE>_<GENERATOR>_upartv2/`
— one directory per flag combination — with a timestamped `.log` per run (stdout is
redirected there, so the terminal stays quiet).

The macro reads the **per-block** files of step 1 directly and sums them in memory; there is
no merge step to run first. `SAMPLE=both` appends the bjet blocks to the qcd ones.

### How "both" is summed

- **count histograms** and the **`RooUnfoldResponse`** are added (`TH1::Add`,
  `RooUnfoldResponse::Add`) — the same operation as `merge_RMatrixTH2D_bjet_qcd_merged.C`;
- **purity/efficiency ratios are never added.** They are recomputed from the summed counts.
  Adding the N stored ratios gives N× the true value (17 Herwig blocks peak at ~16.9).

Sanity check: with `test_mode=0` the closure must be exact — refolding `chi2/ndf = 0.000,
p = 1.000`. If it is not, the response and the corrections are not from the same events.

Physics caveat: `both` adds a b-enriched sample onto inclusive QCD, which double-counts
b-jets unless the sample weights already handle the overlap. Comparing `qcd` vs `both`
shows how much it moves.

The **template fit** file (signal fraction, needed only for `test_mode=2`) is not produced
here — it is read from Afnan's results area and is Pythia-only, so `GENERATOR` does not
affect it. `SAMPLE=both` picks the `qcdbjet` fit, otherwise the `qcd` fit.

---

## 4. Systematics and the final plots — `apply_weights_and_systematics.C`

Reads one step-3 result per variation, builds the band, and writes the two plots you
actually show. **Nothing to edit for a normal run** — the variations are the
`std::vector<Variation> variations` list near the top of the macro, and each entry names the
unfolding run it needs.

```bash
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'
```

A variation whose result file is missing is reported with the exact `apply_unfolding_2d.C`
call that would produce it, and skipped — the nominal alone still works, you just get a
smaller band.

**Outputs**, all into the *nominal* result folder
(`unfolding_both_pythia_upartv2/`):

| File | What it is |
|---|---|
| `data_vs_gen_<label>.{pdf,png}` | the money plot: unfolded data with stat bars and the systematic box, against particle-level MC |
| `sfupart_before_after_<label>.{pdf,png}` | the UParT SF on its own: the result before and after it, and their ratio. Only produced while the display-only `sfupart_off` entry is booked |
| `systematics_curves_<label>.{pdf,png}` | **per-source breakdown** — one curve per variation, and below it each source's fractional pull against the total band |
| `final_<label>_with_systematics.root` | `h_nominal`, `h_syst_total_up`/`_down`, `h_syst_total`, `h_result_stat_syst`, `g_syst_band` (asymmetric), `h_syst_*`/`h_delta_*`/`h_var_*` per source, and `h_gen_*` |

### Where the results are, and which of them are results

The results area gets one folder per flag combination, and **most of them are inputs to the
band, not results**. `apply_weights_and_systematics.C` prints a `FINAL OUTPUTS` block at the
end of every run naming the deliverables, so this never has to be worked out by hand.

All under `/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/`:

| Folder | What it is |
|---|---|
| `unfolding_both_pythia_upartv2/` | **the result.** Nominal run, and where all four deliverables are written |
| `unfolding_both_herwig_upartv2/` | input — "Detector response" |
| `unfolding_both_pythia_tfherwig_upartv2/` | input — "MC template modeling" |
| `unfolding_both_pythia_trkdrop030_upartv2/` | input — "Tracking efficiency" |
| `unfolding_both_pythia_tf0Bx2/x0_upartv2/` | input — "Light jet mistagging" pair |
| `unfolding_both_pythia_sfupartjphf/qqup/qqdn_upartv2/` | input — UParT SF systematics (currently commented out of the band) |
| `unfolding_both_pythia_sfupartoff_upartv2/` | comparison — the before-SF run, display only |

The four deliverables, all in the nominal folder:

| File | What it is |
|---|---|
| `final_<label>_with_systematics.root` | the result: central values, stat, σ⁺/σ⁻, per-source histograms |
| `data_vs_gen_<label>.{pdf,png}` | the measurement vs particle-level MC |
| `systematics_curves_<label>.{pdf,png}` | per-source breakdown |
| `sfupart_before_after_<label>.{pdf,png}` | the UParT SF on its own |

⚠️ Three folders predate the current naming and are **not** produced by this code any more:
`unfolding_upartv2/`, `unfolding_both_pythia_respherwig_upartv2/`,
`unfolding_both_herwig_tfherwig_upartv2/`. Nothing reads them. Left in place rather than
deleted, but do not mistake them for current results.

⚠️ `/data_CMS` is **not** on the sshfs mount, so none of this opens from the Mac. Use
`sync_plots_to_eos.sh` or `scp`.

### Paths and variation tags live in `result_paths.h`

Every result path and every variation tag has **one** definition, in `result_paths.h`, which
both `apply_unfolding_2d.C` (which writes them) and `apply_weights_and_systematics.C` (which
reads them) include.

This is not tidiness — it is a bug fix. The two macros used to carry private copies of the
same mapping, and on 2026-09-11 a new variation (`"off"`) was added to one copy only. The
reader's copy did not recognise it, returned `""`, and `""` **is** the nominal tag — so it
loaded the nominal result and compared it against itself. Nothing errored; the column just
read `-0`, which looks like "the correction does nothing" rather than "wrong file".

Three guards now, in increasing order of how early they catch it:

1. **Validation up front.** `apply_weights_and_systematics.C` checks every booked
   variation's sample/generator/tfVariation/sfupartVariation against `isKnown*()` before
   opening anything, and refuses to run naming the bad entry. It also rejects two variations
   sharing a name, which would silently overwrite each other's histograms.
2. **Unknown ⇒ poison, never nominal.** An unrecognised value returns `kUnknownTag()`
   (`_UNKNOWN_VARIATION`), not `""`, so the path cannot alias onto the nominal — it simply
   fails to open, with an error naming the bad value.
3. **The writer checks the reader's path.** After each unfolding, `apply_unfolding_2d.C`
   rebuilds the path `apply_weights_and_systematics.C` will look for and warns if it does not
   exist — so writer/reader drift surfaces at write time, not days later.

**To add a variation, add it in one place: `result_paths.h`.** Add the value to the relevant
`isKnown*()` and its tag to the matching `*VarTag()`, and both macros follow.

### The UParT b-tag efficiency scale factor

The b-tag efficiency SF for the UParT selection on jets with **at least two reconstructed
secondary vertices**, from `lifetime_jp_sfb_agg_dr_wp712_eec_systematics.root` in this
folder, binned in reco R_BB over the same 9 ΔR bins as the analysis.

**It is applied inside `apply_unfolding_2d.C`, to the reco-level data, before unfolding** —
right after the signal fraction multiplies `h_data_after_fit`. That is the point of it: the
SF is measured in *reco* R_BB, so it belongs on a distribution that is still in reco R_BB,
and the response matrix then carries it through the migration. Applying it to the unfolded
result instead would put a reco-binned correction on particle-level bins and silently ignore
that migration.

It is **divided** in: the chain corrects data with the efficiency measured in MC, and the
true efficiency is `ε_data = SF · ε_MC`, so `data/ε_MC` has to be divided by SF again. Data
only — it never touches MC, and the closure modes are refused. It runs 0.933 at low ΔR to
1.149 at bin 8, so it is a real reshaping, not a wash.

The central SF is a **correction, not an option**: every data run applies it, and an untagged
results path means "with the nominal SF", not "without one".

⚠️ **The KEYS and the histogram names differ in that file.** Key `h_SFb_dr_central` holds a
histogram whose internal name is `h_SFb_dr_wp712_eec_central`. `TFile::Get()` takes the
**key**, so the internal name returns null. `upartSFHist()` uses key names.

`SFUPART_VARIATION` (9th argument) picks which SF to divide in, and tags the results folder.
The name is namespaced to UParT on purpose — more scale factors are coming, and each should
get its own `sfupart*`-style plumbing rather than sharing one generic `SF_` knob:

| `SFUPART_VARIATION` | Key | Tag | Treatment in the band |
|---|---|---|---|
| `nominal` | `h_SFb_dr_central` | — | the correction itself |
| `jpcalib_hf` | `h_SFb_dr_jpcalib_hf` | `_sfupartjphf` | one alternative calibration → one-sided, symmetrised |
| `qqrate_up` | `h_SFb_dr_qqrate_up` | `_sfupartqqup` | grouped pair with ↓, envelope once, stays asymmetric |
| `qqrate_down` | `h_SFb_dr_qqrate_down` | `_sfupartqqdn` | ditto |
| `statup` | `h_SFb_dr_central` **+ its bin error** | `_sfupartstatup` | grouped pair with ↓ — the calibration's statistical precision |
| `statdn` | `h_SFb_dr_central` **− its bin error** | `_sfupartstatdn` | ditto |
| `off` | — | `_sfupartoff` | no SF at all; display-only comparison, never in the band |

The `statup`/`statdn` shift is **coherent** across ΔR — every bin moves together. That is
deliberate: the analysis's last two ΔR bins share one calibration bin, so their errors really
are 100% correlated, and a coherent shift gets that right while staying conservative
elsewhere.

```bash
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"nominal","jpcalib_hf")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"nominal","qqrate_up")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"nominal","qqrate_down")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"nominal","statup")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"nominal","statdn")'
```

`SFUPART_VARIATION = "off"` applies **no** SF at all. It is not a systematic — it exists so
the same macro produces the before-SF result (tag `_sfupartoff`) for a like-for-like
comparison. It is booked in `apply_weights_and_systematics.C` as a **display-only** entry
(`Variation::inBand = false`): drawn and printed, but contributing nothing to σ⁺/σ⁻ and
nothing to the tag, because it is a correction, not an uncertainty. It drives two things —
the "Before UParT SF" curve on `systematics_curves_*`, and the dedicated
`sfupart_before_after_*` plot (the two curves and their ratio, nothing else). Comment that
one entry out and both disappear.

⚠️ `sfupartVarTag()` is duplicated in **both** macros and they must agree. An unrecognised
value falls through to `""`, which resolves to the **nominal** result rather than failing —
so a typo shows up as a suspiciously exact zero difference, not an error.

The SF's **statistical** error is **not** folded into the data's error. Applying the SF is a
pure rescale — content and the data's statistical error multiplied by the same 1/SF, so the
relative error is unchanged and the error bar stays a *data* error bar — and the SF's own
precision is booked as the `statup`/`statdn` systematic instead. Doing both would count it
twice.

`h_SFb_dr_syst_up`/`_syst_down` are the quadrature of the two calibration systematics and are
**not** used: booking `jpcalib_hf` and `qqrate_up`/`down` separately gives the same total
while keeping each visible as its own column and curve.

⚠️ Because the SF now lives inside the unfolding, **every** unfolding run changes when it
changes — the nominal and all six systematic variations have to be re-run together, not just
the SF ones.

### MC-derived corrections: bin-by-bin rescaling

Purity (before unfolding) and the three post-unfolding corrections — reconstruction
efficiency, 2SV+b-tag efficiency, EEC weight — are applied by `applyCorrection()` as a plain
**bin-by-bin rescale**: content and error multiplied by the *same* factor, so the data's
**relative** error is unchanged.

Deliberately **not** `TH2::Multiply`/`Divide`. Those propagate the correction's own error
into the result, which folds the MC's limited statistics into the *data's* statistical error.
The convention here is to treat every MC-derived correction as exact. The helper also runs
over underflow/overflow, checks the binning explicitly instead of failing silently, and
reports any bin where a zero correction emptied non-empty data.

⚠️ **The MC statistical uncertainty is therefore not counted anywhere.** It is out of the
stat error and no systematic is booked for it — a deliberate choice (2026-09-11), not an
oversight. It is *not* negligible: removing it cut the statistical error by 13–22%, which
backs out to an MC contribution roughly as large as the data's own. Revisit if the error
budget is ever questioned.

⚠️ Corrections applied at **fill time** — the EEC weight and sample weights inside
`create_files_for_template_fit.cpp` — cannot get this treatment. They are baked into the
histogram contents with their fluctuations already folded in, so there is no separable factor
left to rescale, and their MC statistics stay mixed into the data's error.

### The band is asymmetric

σ⁺ and σ⁻ are computed and stored separately, because a genuine up/down variation need not
move the result by the same amount each way. Which side a source feeds depends on whether it
is *really* two-sided:

- **Grouped entries** (a shared `Variation::group` — an up/down pair) are treated **signed**.
  The member that pushes the result up sets σ⁺, the one that pushes it down sets σ⁻, and the
  source contributes its largest shift on each side once.
- **Ungrouped entries** have one alternative (a Pythia↔Herwig swap, the 3% track drop). One
  alternative tells you the *size* of the difference, not its direction — nothing says the
  truth lies on the far side of the nominal — so these are **symmetrised**: |Δ| into both
  sides.

With only the mistag pair genuinely two-sided today, the two bands come out equal to ~3
significant figures. **That is expected** — the machinery is there for the asymmetric sources
still to come, and a new one is asymmetric the moment it is booked as a group.

If both members of a pair happen to push the same way in a bin, that source would contribute
nothing to one side, which claims zero uncertainty there. Too aggressive, so such a group
falls back to its envelope on both sides and prints a `NOTE:` naming itself.

`h_syst_total` is still written as the **larger** of the two sides, for anything that wants
one number per bin (`h_result_stat_syst`'s bin errors, and any script predating the split).

⚠️ **Two sign conventions, deliberately.** The printed table shows `variation − nominal`, the
shift of the *result*, so its sign matches which of σ⁺/σ⁻ the source feeds. The stored
`h_delta_*` keep the original `nominal − variation`, because `plot_systematics_summary.C`
reads them and flips the sign itself. Same numbers, opposite sign.

### The per-source breakdown plot

`systematics_curves_*` is where you read **what each systematic does on its own**. The
bottom pad is the content: `variation / nominal − 1` for every source, one line each, over
the shaded total band — so a source's line reaching the edge of the band means it dominates
that bin, and a line sitting on zero means it does nothing there. The top pad is the same
variations as absolute EEC curves, which mostly shows that they all sit on top of the
nominal; it is a reference, not the message.

The pads are split 55/45 rather than the usual 70/30 for exactly that reason — the pulls are
a few percent and need the room. Up/down pairs appear as **two separate lines** here (not
the envelope that enters the total), which is what you want for debugging: it shows whether
a pair is genuinely two-sided or whether one side does all the work.

Per-bin numbers for the same thing are in the printed table — one signed column per source —
and in `final_*_with_systematics.root` as `h_delta_<name>` (signed shift), `h_syst_<name>`
(the symmetrised `|shift|` that enters the band) and `h_var_<name>` (that variation's own
EEC curve).

### The data-vs-gen plot

Top pad: red = the measurement (unfolded data, stat error bars, systematic as a shaded box
per point), blue = particle-level MC of the nominal generator, purple dashed = the other
generator when its unfolding exists. Bottom pad: **MC / Data**, with the data's own
uncertainty drawn around 1 — the systematic as the wide shaded box, the statistical error as
the narrow hatched one, so the two are separable by eye.

Upper left is a **"Systematic uncertainties" tag** listing the sources in the band. It is
built from the `sysLabel` field of the entries that *actually entered* the band, deduplicated
so an up/down pair is listed once — a source whose result file was missing and got skipped
cannot be claimed on the plot. Leave `sysLabel` empty to keep an entry out of the tag while
it still counts in the total. The second particle-level curve is off by default; flip
`show_alt_gen` to draw the other generator next to the nominal one.

The gen curve is `h_mc_true_2D`, which in data mode is `hgenjet_2b_all`: the particle-level
reference the corrected data is meant to land on. It is only re-normalised to unit area,
never weighted — the scale factors in `Weights` correct *data* for detector effects and have
no meaning at particle level.

⚠️ **`h_mc_true_2D`'s `Write()` was disabled until 2026-09-10**, so any result unfolded
before then does not carry it and the plot is skipped with a message. Re-run
`apply_unfolding_2d.C` for that configuration to get it.

### Not booked yet

The band is deliberately these five for now. Known gaps, recorded so they are not
rediscovered from scratch:

| Missing | Status |
|---|---|
| Template **fit range** | `fitRange0to8` is already produced by `template_fit.cpp` on every run (it is inside the `ivar < 4` loop) and consumed nowhere. Now that `TF_VARIATION` exists it is two more unfoldings to book — the cheapest of these by far. `FITRANGE_0_7` is in the enum but the loop stops at 4, so it is not even produced. |
| Unfolding **regularisation** | Nothing varies the Bayesian iteration count; `scan_niter=false` pins it at 4. |
| **JES / JER** | Absent from the whole chain. Moves jets across the pT bin edges, so it does not fully cancel in a unit-area shape. |
| b-tag **efficiency** SF | Only the light-jet *mistag* is covered. `Weights::sf_file` is empty, so no SF is applied at all. A dR-flat efficiency cancels in the normalisation; a dR-dependent one does not. |

⚠️ **Double counting between two that ARE booked.** "Detector response" (`unfolding_model`)
swaps the *entire* unfolding MC to the other generator, which includes the EEC weight
correction — so it already contains "MC modeling of EEC weight". Booking both counts that
difference twice. Read it as the deliberately conservative option, or drop
`unfolding_model` and book the corrections individually.

### Adding a systematic

One entry in the `variations` vector; everything else — the quadrature total, the printed
table, both plots, the written histograms — iterates that vector. The fields are the
`apply_unfolding_2d.C` arguments of the run that produced it:

```cpp
{ "name", SAMPLE, generator, tfGenerator, "legend text", trackEffUnc, tfVariation, group, "Tag name" }
```

`group` is the only non-obvious one. Leave it `""` for an independent source and its
symmetrised `|nominal - variation|` enters the quadrature sum on its own. Give **two**
entries the *same* non-empty group when they are the two sides of one uncertainty (an
up/down pair): the total then takes their per-bin **envelope once** instead of adding both,
which would inflate that source by √2. `mistag_0B_up` / `mistag_0B_down` are the worked
example.

---

## Macro arguments (`create_files_for_template_fit.cpp`)

```cpp
create_files_for_template_fit(RunN, pT_low, etaCut, n, btag, isMC, btagWP,
                              makeTemplates, createRmatrix, makeAggNtuple,
                              ev_first, ev_last, inputFileOverride,
                              outputFolderOverride, sampleTag,
                              track_eff_unc, trkSeed)
```

`track_eff_unc` defaults to `false` (nominal, a no-op) and `trkSeed` to 20260908; the run
script sets them from `TRACK_EFF_UNC` / `TRK_SEED`. The 3% itself is
`TrkEffSyst::kDropFraction`.

There is no `dataType` any more: **`isMC`** is data vs MC, and **`sampleTag`**
(`data`/`bjet`/`qcd`, set by the run scripts) is what names the output histograms
`..._3D_<sampleTag>_f`. `pT_high` is also gone — the upper pT edge lives in
`jtpt_binsVector` in `binning_histos_small.h`.

## Conventions

- Compile the `.so` **once**; jobs `gSystem->Load()` it. Parallel ACLiC compilation in the
  same directory produces half-written artifacts.
- ACLiC products, ROOT files, plots and logs are gitignored — do not commit them.
- `*~`, `#*#` (Emacs) and `._*` (macOS/sshfs) files litter the tree; ignore them.
- Unwanted analysis code is commented out with `/* ---- disabled (kept for reference) ---- */`
  rather than deleted.
