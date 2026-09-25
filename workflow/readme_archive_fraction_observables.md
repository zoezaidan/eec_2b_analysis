# Archived: the b/non-b momentum fractions `fnb`, `fb`, `lnfb`

Moved out of the workflow readme (now the repo-root `README.md`) on 2026-09-24, verbatim. None of the three is measured:
they were disabled on 2026-09-22 and the measurement is EEC(ΔR) and dN/dB only. Their code
is still in `observables.h`, `binning_histos_small.h`, `result_paths.h` and
`create_files_for_template_fit.cpp`, inside `DISABLED 2026-09-22` blocks. Paths and
commands below are as they were when written. The balance was still called `z` then, and
re-enabling any of these needs a fresh production, not just an uncomment.

### The `fnb` observable — DISABLED 2026-09-22 (archived)

> ⚠️ **Not measured any more.** `fnb`, `fb` and `lnfb` were commented out on 2026-09-22 —
> the measurement is ΔR (EEC) and the momentum balance `B` only. The code for all three is
> still in `observables.h`, `binning_histos_small.h`, `result_paths.h` and
> `create_files_for_template_fit.cpp`, inside `DISABLED 2026-09-22` blocks, so re-enabling
> one is four uncomments. Everything below is kept as the record of what was measured and
> what it cost; the numbers in it are still the numbers those productions gave.
>
> ⚠️ Re-enabling also needs a **fresh production**: their productions on disk predate the
> `z` → `B` rename and name the balance axis `_z`, which this code no longer reads.


`OBSERVABLE = "fnb"` measures how much of the jet is **not** in the two B hadrons:

```
fnb = 1/(1 + pT_b/pT_nonb) = pT_nonb / (pT_b + pT_nonb)          in [0, 1]

  pT_b     pT(B1) + pT(B2), the two aggregated B hadrons (scalar sum, as B uses)
  pT_nonb  RECO  jet tracks, pT > 1, whose BDT score fails the b-track cut
           GEN   jet refTracks, pT > 1, with refTrkSta < 100 (not from a B decay)
```

At reco that is exactly the complement of the B reconstruction — the aggregation absorbs
every no-SV track that passes the cut into the nearer vertex, so nothing else is left over.
It is computed inside `makeSvtxs_withBDT`, the only place that knows which tracks were kept.
No BDT cut is applied on the non-b side: picking the non-b component with a b-likeness score
would bias it, and gen has no BDT to match.

Like B it is a **yield** (`dN/dfnb`, EEC weight off) and carries **no UParT SF** — both from
the same `result_paths.h` functions.

⚠️ **Reco is strongly biased up against gen**: mean 0.264 vs 0.122 on the same jets (400k
events, Pythia8 QCD block_0000). The BDT rejects genuine B-decay tracks, which moves pT from
the b side to the non-b side — the same effect that puts reco SV pT at 0.81 of the gen B pT.
So fnb migrates much more than dr or B, and its response matrix is strongly off-diagonal.
Run the equivalent of `momentum_balance_mc_study.C` on fnb before trusting an unfolded result.

**Produced 2026-09-17** — MC (qcd + bjet), data (50 Condor jobs), fit, unfolding and band.
The commands, for the record or to redo it:

```bash
# 1. MC, both samples. OUT_TAG must match obsProdTag("fnb") in result_paths.h
SAMPLE=qcd  EEC_WEIGHT_OFF=true OUT_TAG=upartv2_3obs ./run_agg_ntuple_chunks.sh
SAMPLE=bjet EEC_WEIGHT_OFF=true OUT_TAG=upartv2_3obs ./run_agg_ntuple_chunks.sh
hadd -f <QCD|bJet>/agg_ntuple_chunks/..._fMCGEN_noeecw_upartv2_3obs.root <same>/block_*/<same>
# 2. data (the .so must NOT be recompiled while these run -- do the MC first)
EEC_WEIGHT_OFF=true OUT_TAG=upartv2_3obs ./make_hardprobes_condor_scripts.sh && condor_submit condor.submit
hadd -f HardProbes/agg_template_chunks/..._data_fMCGEN_noeecw_upartv2_3obs.root <same>/HardProbes*/block_*/<same>
# 3. fit, unfolding (nominal + the mistag pair), band
root -l -b -q 'template_fit.cpp("both","pythia",false,true,"fnb")'
for v in nominal var0B_2 var0B_0; do
  root -l -b -q "apply_unfolding_2d.C(\"both\",\"pythia\",2,false,false,\"pythia\",false,\"$v\",\"off\",true,\"fnb\")"
done
root -l -b -q 'apply_weights_and_systematics.C("both","pythia",2,false,"h_data_fully_corrected_2D","pythia",2,false,"fnb")'
```

**Seeing the whole chain for one observable.** Every step draws itself, so an observable can
be read end to end. For fnb:

| Step | Macro | Where the plots land |
|---|---|---|
| 1. raw yields, pre-fit | `plot_z_first_look.C(...,"fnb")` | `results/fnb_first_look_both_pythia_noeecw_upartv2_3obs/` |
| 2. template fit | `template_fit.cpp(...,"fnb")` | `results/TemplateFit_Run3/TemplateFits_both_pythia_noeecw_fnb_upartv2/` |
| 3. unfolding | `apply_unfolding_2d.C(...,"fnb")` | `results/matrix_inversion/unfolding_both_pythia_fnb_sfupartoff_noeecw_upartv2/` |
| 4. systematics | `apply_weights_and_systematics.C(...,"fnb")` | the same folder as step 3 |

`plot_z_first_look.C` keeps its name but takes the **observable as its last argument** and
works for any of them — the five plots it writes (flavour templates, m_2B separation per bin,
purity/efficiency, data vs MC, raw yields) are named and foldered after that observable, so B
lands exactly where it always did. It refuses a production tag that does not belong to the
observable asked for, checked against `obsProdTag()`.

**What the pre-fit look says about fnb** (both/pythia, 80 < pT < 100):

- the flavour templates **do** separate: 0b sits at high fnb (33% of it above 0.3, against
  13% for 2b), so the signal fraction genuinely varies along this axis
- m_2B still separates signal from background **inside every one of the six bins**, which is
  what the fit needs
- matching purity 0.74–0.82 and efficiency 0.72–0.76, comparable to ΔR's 0.71–0.81
- at reco, data/MC runs 0.60 → 1.37 across the range: data sits at systematically higher fnb
  than Pythia before any correction
- 927k raw jets in the lower pT bin, 0.2–0.4% Poisson per bin

**First result** (uniform 0.2 binning, set 2026-09-18). Full-MC closure is exact (refolding
chi2/ndf = 0.000, p = 1), so the chain itself is sound on this axis.

| bin | dN/df_nb | stat | syst | raw jets |
|---|---|---|---|---|
| [0.0, 0.2] | 1.999 | 0.6% | 0.9% | 645,084 |
| [0.2, 0.4] | 1.927 | 0.9% | 1.0% | 593,821 |
| [0.4, 0.6] | 0.726 | 2.8% | 0.3% | 343,623 |
| [0.6, 0.8] | 0.324 | 6.9% | 3.9% | 161,523 |
| [0.8, 1.0] | 0.024 | **122%** | 113% | 29,036 |

⚠️ **The top bin is consistent with zero and unusable**, and the one below it is marginal.
It is NOT a statistics problem -- 29k raw jets sit in it. It is that particle-level fnb has
~0.2% of its weight above 0.8, so the inversion is asked to recover an almost empty truth bin
from a reco bin filled mostly by migration, and unregularised inversion turns that into a
±122% error. The shape-driven binning this replaced ({0, .05, .10, .20, .30, .50, 1.0}) gave
1.4-2.7% stat in every one of its six bins.

Quote four bins, merge [0.6, 1.0] into one, or go back to a shape-driven binning -- but do
not show the last bin as it stands.

**Data and Pythia disagree strongly in shape**, and they already do so at reco level, before
any fit or correction (data/MC 0.60 -> 1.37 across the range in the pre-fit look). So the
disagreement is not made by the unfolding; the corrections change its size, not its
existence. Read it against the migration warning above -- the reco definition (BDT-rejected
tracks) and the gen one (tracks not from a B decay) agree only as far as the BDT separates
them, and reco sits at mean 0.264 against gen 0.122.

That production fills **all three** observables in one pass, so it also gives B (and dr) a
production that a track-drop and Herwig counterpart could be built from — which is exactly
what the B band is missing today.

### The `fb` observable — pT_b/(pT_b + pT_nonb), linear — DISABLED 2026-09-22 (archived)

`OBSERVABLE = "fb"` is the b momentum fraction on a plain linear axis (requested 2026-09-21).
Uniform 0.2 bins, production tag `_upartv2_fb`, yield (`dN/dfb`), no UParT SF.

It is the SAME quantity as `fnb` and `lnfb`: fb = 1 − fnb, lnfb = ln(1 + fb). All three are
filled in one pass. On the matching uniform edges `fb` and `fnb` use, **fb bin k is fnb bin
6−k exactly** — verified bit-for-bit on the gen, purity and efficiency histograms of one
production. That mirror is a free closure check: if an fb result is ever not the reverse of
the corresponding fnb one, one of them is wrong.

⚠️ **`fb` FAILS its split test, under both methods**, and in both cases bin 1 alone is the
cause: ratio 0.524 ± 28% (MI, 3.2σ) and 0.733 ± 9.7% (Bayesian, 2.7σ). Bins 2–5 are fine and
pass comfortably on their own. This is the mirror of the fnb top bin already recorded above as
unusable at ±122%, so it is a property of that corner of the jet — almost no truth, a reco bin
filled mostly by migration — not of the new code. Unlike `lnfb`, where the weak bins were only
imprecise, here the bias is significant: **do not show fb bin 1.**

Raw response κ = 1701; the matrix actually inverted has κ = 4.89.

### The `lnfb` observable — ln(1 + pT_b/(pT_b + pT_nonb)) — DISABLED 2026-09-22 (archived)

`OBSERVABLE = "lnfb"` is the same split as `fnb`, monotonically transformed:

```
lnfb = ln(1 + pT_b/(pT_b + pT_nonb))      in (0, ln2 = 0.69315]
```

That formula is evaluated as written, straight from the two momenta — **not** through a named
fraction and **not** through a `1/(1 + ratio)` rewrite (settled 2026-09-21). The quantity is
the complement of `fnb`, but nothing in the code routes through `fnb` to get it. The labels
on every plot say the formula rather than `f_b`, so a plot is readable without this file.

⚠️ Renaming that cost nothing: the code already evaluated the direct form, so **no value
changed and no production was redone** — only the axis labels, which the plot macros read from
`ObsDef` at draw time, not from the stored histogram titles. (The stored titles inside the
production `.root` files still carry the old text; nothing reads them.) For the record, the
`1/(1 + pT_nonb/pT_b)` form is not even numerically identical — it differs by one ulp in about
31% of cells — which is a small argument for the direct form on its own.

Same pT_b / pT_nonb definitions as `fnb`, so it carries no new information — what changes is
where the binning puts the resolution. Like B and fnb it is a yield (`dN/dlnfb`, EEC weight
off) with **no UParT SF**, both forced from `result_paths.h`.

**Binning is uniform as of 2026-09-21** — five bins of ln2/5 = 0.138629 over the whole range,
replacing the shape-driven `{0, .45, .55, .62, .66, .69315}`. The binning is **baked into the
MC production**, so it got its own tag: `obsProdTag("lnfb") = "_upartv2_lnfblin"`, and the
old-edge results stay reproducible from `_upartv2_lnfb`.

⚠️ **The distribution is not uniform on this axis**, and the closure runs show what that
costs. Gen weight per bin spans a factor ~500 from bin 1 to bin 5:

| bin | edges | gen dN/dlnfb | full-closure stat | split-test stat | split unfolded/truth |
|---|---|---|---|---|---|
| 1 | [0.000, 0.139] | 0.0112 | 17.2% | **30.8%** | **0.82** |
| 2 | [0.139, 0.277] | 0.1105 | 7.0% | 11.5% | **0.87** |
| 3 | [0.277, 0.416] | 0.3236 | 4.5% | 6.1% | 1.04 |
| 4 | [0.416, 0.555] | 0.9489 | 2.2% | 3.2% | 0.96 |
| 5 | [0.555, 0.693] | 5.8193 | 0.4% | 0.5% | 1.01 |

Bins 1–2 are the mirror of fnb's unusable top bin: an almost-empty truth bin recovered from a
reco bin filled mostly by migration. They are not *biased* at a significant level (0.82 sits
0.6σ from 1, 0.87 sits 1.2σ), but at 31% and 12% statistical precision they say nothing.
**Quote bins 3–5; merge 1–2 or drop them.**

The uniform edges are the better-conditioned choice, which is the one thing that improved:
condition number **4.41**, against 9.39 on the shape-driven edges.

⚠️ **Stale siblings in the results folder.** The results-folder tag carries the observable but
**not** `obsProdTag()`, so `unfolding_both_pythia_lnfb_sfupartoff_noeecw_upartv2/` holds the
new closure/split files next to a `_MI_data_` result and a `final_*_with_systematics.root`
from 2026-09-18 that were built on the **old** edges. Nothing was overwritten — the filenames
differ only by test mode — but the data result in that folder does not belong to the binning
the closure files in it were made with. Re-run data + band on `_lnfblin` before reading them
together.
