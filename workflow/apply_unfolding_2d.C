#include "binning_histos_small.h"
#include "result_paths.h"   // the ONE definition of every result path and variation tag
#include <algorithm>
#include <vector>

namespace ROCColor {
  Color_t blue()   { return TColor::GetColor("#4C72B0"); }
  Color_t red()    { return TColor::GetColor("#C44E52"); }
  Color_t green()  { return TColor::GetColor("#4F8F52"); }
  Color_t purple() { return TColor::GetColor("#8C6BB1"); }
  Color_t orange() { return TColor::GetColor("#C97430"); }
  // Neutral highlight for the ratio-pad gap bands: a hue used by none of the five curves above,
  // so a filled region never reads as one of the plotted quantities. Drawn at two alphas.
  Color_t teal()   { return TColor::GetColor("#1F8A8A"); }
}

template <typename T>
T *getOrWarn(TFile *f, const TString &name)
{
    T *obj = dynamic_cast<T *>(f->Get(name.Data()));
    if (!obj)
        std::cerr << "ERROR: '" << name << "' (" << T::Class_Name() << ") not found in "
                  << f->GetName() << std::endl;
    return obj;
}

TFile *openOrWarn(const TString &name)
{
    TFile *f = TFile::Open(name);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << name << std::endl;
        return nullptr;
    }
    return f;
}

std::vector<TFile *> openAllOrWarn(const std::vector<TString> &names)
{
    std::vector<TFile *> files;
    for (const TString &n : names) {
        TFile *f = openOrWarn(n);
        if (!f) {
            for (TFile *g : files) g->Close();
            files.clear();
            return files;
        }
        files.push_back(f);
    }
    return files;
}

// ---- Summing inputs across files ----------------------------------------------------
// One list covers both levels of summing: the per-block outputs of run_agg_ntuple_chunks.sh,
// and -- with sample "both" -- the bjet and qcd samples stacked on top of that. Same
// operation as merge_RMatrixTH2D_bjet_qcd_merged.C, done in memory so no pre-merged file
// is needed.
//
// ONLY count histograms and the response may be summed like this. The pre-divided ratios
// (h_full_purity_tf, h_full_efficiency_tf, ...) must NOT be: adding N of them gives N x the
// true value. They are recomputed from the summed counts by ratioFromCounts() below.
template <typename T>
T *sumOverFiles(const std::vector<TFile *> &files, const TString &name)
{
    T *out = nullptr;
    for (TFile *f : files) {
        T *h = getOrWarn<T>(f, name);
        if (!h) { delete out; return nullptr; }
        if (!out) {
            out = (T *) h->Clone(name + "_sum");
            out->SetDirectory(nullptr); // outlives the input files
        } else {
            out->Add(h);
        }
    }
    if (!out) std::cerr << "ERROR: no input files to sum '" << name << "' from" << std::endl;
    return out;
}

// RooUnfoldResponse::Add() adds the migration matrix together with the truth and measured
// projections, which is what keeps the efficiency built into the response consistent.
RooUnfoldResponse *sumResponseOverFiles(const std::vector<TFile *> &files, const TString &name)
{
    RooUnfoldResponse *out = nullptr;
    for (TFile *f : files) {
        RooUnfoldResponse *r = getOrWarn<RooUnfoldResponse>(f, name);
        if (!r) { delete out; return nullptr; }
        if (!out) {
            out = new RooUnfoldResponse(*r);
            out->SetName(r->GetName());
            out->SetTitle(r->GetTitle());
        } else {
            out->Add(*r);
        }
    }
    if (!out) std::cerr << "ERROR: no input files to sum '" << name << "' from" << std::endl;
    return out;
}

// ---- Sample selection ---------------------------------------------------------------
// sample "qcd" | "bjet" | "both", generator "pythia" | "herwig" -> the agg_ntuple_chunks
// directories written by run_agg_ntuple_chunks.sh for those flags.
std::vector<TString> samplesIn(const TString &sample)
{
    if (sample == "both") return {"qcd", "bjet"};
    return {sample};
}

TString sampleSubdir(const TString &sample, const TString &generator)
{
    const bool herwig = (generator == "herwig");
    if (sample == "qcd")  return herwig ? "QCDHerwig"  : "QCD";
    if (sample == "bjet") return herwig ? "bJetHerwig" : "bJet";
    std::cerr << "ERROR: unknown sample '" << sample << "' (use qcd | bjet | both)" << std::endl;
    return "";
}

// Which of template_fit.cpp's variations supplies the signal fraction, and every path tag,
// now live in result_paths.h -- shared with apply_weights_and_systematics.C so the writer
// and the reader cannot disagree. varNames[] in Help_Functions.h are the same strings.

// The template fit that supplies the signal fraction, the one generator-dependent input
// that is NOT part of the unfolding MC. Returns "" when that fit does not exist, so both
// the driver and apply_unfolding() can refuse before doing any work.
TString templateFitFile(const TString &sample, const TString &tfGenerator,
                        bool track_eff_unc = false,
                        const TString &tfVariation = "nominal")
{
    // Written by template_fit.cpp(SAMPLE, GENERATOR), one directory per combination.
    // Its SAMPLE is qcd | both only -- there is no bjet-only fit, since h3D_0b exists in
    // the qcd sample alone -- so a bjet unfolding takes the qcd fit.
    const TString tf_sample = (sample == "both") ? "both" : "qcd";
    // With track_eff_unc the fit itself is redone on the 3%-track-drop templates
    // (template_fit.cpp(SAMPLE, GENERATOR, true)), so the signal fraction moves with the
    // tracking efficiency instead of being held at its nominal value. The data being
    // fitted is the same either way.
    const TString trk_tag = track_eff_unc ? "_trkdrop030" : "";
    return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/TemplateFits_"
         + tf_sample + "_" + tfGenerator + trk_tag + "_upartv2/"
         + tfVariation + "_Run3_TemplateFits_histos_3d_80_inf.root";

    /* ---- disabled (kept for reference): Afnan's fits, Pythia8 only ----
    // These were the nominal until template_fit.cpp gained its own flags. Mixing them with
    // a Herwig fit produced here would confound generator with producer and b-tag WP
    // naming (btagWP712 vs btagWP0712), so both generators now come from the same macro.
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/TemplateFits_btagWP712_qcdbjet_upartv2/nominal_Run3_TemplateFits_histos_3d_80_inf.root" // qcd + bjet
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/TemplateFits_btagWP712_qcd_upartv2/nominal_Run3_TemplateFits_histos_3d_80_inf.root"     // qcd
    "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFits_Run3_minHLT60_LinearBin_upartv2/nominal_Run3_TemplateFits_histos_3d_80_inf.root"            // Zoe, older
    ---- */
}

// MC systematic-variation tag, as create_files_for_template_fit.cpp writes it into the
// per-block INPUT filenames. A different namespace from the results-path tags -- this one
// names MC productions, those name result folders -- but it has to be the same literal, so
// it delegates to trkTag() in result_paths.h rather than repeating "_trkdrop030".
// (That literal tracks TrkEffSyst::kDropFraction in tracking_efficiency_syst.h.)
// It sits between the MCGEN part of the name and the run script's OUT_TAG, so a variation
// run looks for <varTag>_upartv2.
TString mcVarTag(bool track_eff_unc)
{
    return trkTag(track_eff_unc);
}

// One per-block file, or the top-level merged one when block < 0. prefix is "RMatrix_" or
// "" (the MCGEN templates), mid is "" or "MCGEN".
TString aggChunkFile(const TString &base, const TString &subdir, int block,
                     const TString &prefix, const TString &tag, const TString &mid,
                     const TString &btagTag, const TString &outSuffix)
{
    TString path = base + "/" + subdir + "/agg_ntuple_chunks/";
    if (block >= 0) path += Form("block_%04d/", block);
    return path + prefix + "Run3_" + btagTag + "_template_for_fit_histos_3D_" + tag + "_f"
                + mid + outSuffix + ".root";
}

// Every per-block file that exists, for one or both samples. Block counts differ per sample
// (10 qcd/pythia, 9 bjet/pythia, 8 qcd/herwig, 9 bjet/herwig), so take what is on disk
// rather than a fixed range -- same rule as the glob in run_agg_ntuple_chunks.sh.
std::vector<TString> aggChunkFiles(const TString &base, const TString &sample,
                                   const TString &generator, const TString &prefix,
                                   const TString &mid, const TString &btagTag,
                                   const TString &outSuffix, int max_blocks = 50)
{
    std::vector<TString> files;
    for (const TString &s : samplesIn(sample)) {
        const TString subdir = sampleSubdir(s, generator);
        if (subdir.Length() == 0) return {};
        int found = 0;
        for (int b = 0; b < max_blocks; ++b) {
            TString f = aggChunkFile(base, subdir, b, prefix, s, mid, btagTag, outSuffix);
            if (gSystem->AccessPathName(f) == kFALSE) { files.push_back(f); ++found; }
        }
        std::cout << "   " << subdir << ": " << found << " block files" << std::endl;
        if (found == 0)
            std::cerr << "ERROR: no block files under " << base << "/" << subdir
                      << "/agg_ntuple_chunks/ -- run run_agg_ntuple_chunks.sh for this sample"
                      << std::endl;
    }
    return files;
}


struct RefoldGof {
    double chi2   = 0.;
    int    ndf    = 0;
    double pvalue = -1.;
    double chi2ndf() const { return ndf > 0 ? chi2 / ndf : -1.; }
};

RefoldGof refoldChi2(const TH2D *h_refolded, const TH2D *h_input,
                     int ix_min, int ix_max, int iy_min, int iy_max)
{
    RefoldGof g;
    if (!h_refolded || !h_input) return g;
    for (int iy = iy_min; iy <= iy_max; ++iy) {
        for (int ix = ix_min; ix <= ix_max; ++ix) {
            const double d  = h_input->GetBinContent(ix, iy);
            const double ed = h_input->GetBinError(ix, iy);
            const double r  = h_refolded->GetBinContent(ix, iy);
            if (d == 0. || ed <= 0.) continue;
            const double pull = (r - d) / ed;
            g.chi2 += pull * pull;
            ++g.ndf;
        }
    }
    // ndf = number of bins compared. Refolding is not a fit to the input distribution,
    // so there are no fitted parameters to subtract.
    if (g.ndf > 0) g.pvalue = TMath::Prob(g.chi2, g.ndf);
    return g;
}

bool normalizeToUnitArea(TH1D *h)
{
    if (!h) return false;
    const Double_t integral = h->Integral();
    if (integral <= 0.) {
        std::cerr << "WARNING: " << h->GetName() << " has integral " << integral
                  << ", skipping normalisation" << std::endl;
        return false;
    }
    h->Scale(1. / integral, "width");
    return true;
}


void apply_unfolding(TString &label, TString &folder, bool btag, Int_t n, TString pT_selection,
                     int test_mode, bool unfoldBayes, bool scan_niter,
                     TString sample = "qcd", TString generator = "pythia",
                     TString tfGenerator = "pythia", bool track_eff_unc = false,
                     TString tfVariation = "nominal", TString sfupartVariation = "nominal")
{
    // Unfolding options. test_mode and unfoldBayes are passed in as arguments:
    //   0 = FULL-MC closure : h3D_bb, full-sample corrections, truth h_full_efficiency_denominator_tf.
    //       Same events in and out, so a technical (not statistically independent) closure.
    //   1 = SPLIT test      : h3D_pseudodata_bb (even half), odd-half corrections,
    //       truth h_pseudodata_truth_tf (even half) -> real, independent closure.
    //   2 = DATA            : h3D_data + template-fit signal fraction, full-sample corrections.
    //   unfoldBayes: true = Bayesian, false = matrix inversion.
    //   scan_niter:  true = scan niter 1..100 and pick the optimal one from the refolding
    //       goodness-of-fit test; false = a single unfolding at niter = 4. Bayesian only.
    const bool split_test       = (test_mode == 1);  // -> pseudo (odd-half) corrections + even-half truth
    const bool is_data          = (test_mode == 2);
    const bool multiply_sigfrac = is_data;           // only real data needs the bb template fit
    // Matching purity = (reco_pass && gen_pass) / (reco_pass), both already restricted to
    // true-2b jets. It is NOT a bb-purity, so h3D_pseudodata_bb still needs it: that input
    // contains jets whose gen jet fails gen_pass, and the response has no truth row for them.
    bool apply_purity = true;


    const Color_t blue = ROCColor::blue();
    const Color_t red = ROCColor::red();
    const Color_t green = ROCColor::green();
    const Color_t purple = ROCColor::purple();
    const Color_t orange = ROCColor::orange();
    const Color_t teal = ROCColor::teal();

    // ---- MC inputs follow the sample/generator flags ---------------------------------
    // The response, the corrections and the MC truth all come from the same file list, as
    // they must: one set of events, one pT binning, or the full-MC closure cannot be 1.
    // sample "both" appends the bjet blocks to the qcd ones and everything is summed.
    // Switching generator switches all of it -- that is the unfolding-model systematic.
    const TString mc_base   = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
    const TString btag_tag  = "btagWP0712"; // follows BTAG_WP in the run scripts (%04d of 1000*WP)
    // OUT_TAG in the run scripts, preceded by the systematic-variation tag the macro
    // writes. track_eff_unc = false is the nominal and gives exactly "_upartv2" as before;
    // true reads the 3%-tracking-efficiency production instead. The response, the
    // corrections, the truth and (in modes 0/1) the input all come from that same list,
    // so the whole unfolding moves with the variation -- which is the point.
    const TString out_tag   = mcVarTag(track_eff_unc) + "_upartv2";

    std::cout << "Sample: " << sample << ", unfolding generator: " << generator
              << ", template-fit generator: " << tfGenerator << std::endl;

    std::vector<TString> filenames_response =
        aggChunkFiles(mc_base, sample, generator, "RMatrix_", "", btag_tag, out_tag);
    if (filenames_response.empty()) return;

    /* ---- disabled (kept for reference): response-matrix-only variation ----
    // Took ONLY the migration matrix from the other generator and left the corrections,
    // the truth and the input alone, isolating the response (the closure is then not
    // expected to be 1 -- that difference IS the uncertainty). Driven by a
    // responseGenerator argument, with a _resp<gen> tag on the folder and filenames.
    std::vector<TString> filenames_corrections =
        aggChunkFiles(mc_base, sample, generator, "RMatrix_", "", btag_tag, out_tag);
    std::vector<TString> filenames_response = filenames_corrections;
    if (responseGenerator != generator)
        filenames_response =
            aggChunkFiles(mc_base, sample, responseGenerator, "RMatrix_", "", btag_tag, out_tag);
    // ... and the response was then read from its own open-file list, everything else
    // from the corrections list.
    ---- */

    /* ---- disabled (kept for reference): pre-merged single response files ----
    // Afnan's merged files, Pythia8 only. The bjet+qcd one was produced by
    // merge_RMatrixTH2D_bjet_qcd_merged.C; sample="both" now does that sum in memory.
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_qcd_f_80_9999_2_merged.root"
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_f_80_9999_2_merged.root"
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_qcd_merged.root"
    ---- */

    //--  Signal fraction from the template fit. This is the OTHER generator-dependent
    // input and it is deliberately its own flag: tfGenerator, independent of the
    // unfolding generator above. The sample still picks qcd vs qcd+bjet templates.
    TString filename_template_fit = templateFitFile(sample, tfGenerator, track_eff_unc, tfVariation);

    if (multiply_sigfrac) {
        if (filename_template_fit.Length() == 0) {
            std::cerr << "ERROR: no template fit available for tfGenerator '" << tfGenerator
                      << "' -- produce it and set the path in apply_unfolding_2d.C" << std::endl;
            return;
        }
        std::cout << "Using template file (" << tfGenerator << "): "
                  << filename_template_fit << std::endl;
    }

    // Dataset ntuples to unfold: data or MC (test mode decision).
    // Modes 0/1 unfold MC (h3D_bb / h3D_pseudodata_bb) from the MCGEN files of the same
    // sample as the response, so they are summed the same way.
    // Mode 2 unfolds real data, which is one file and independent of the flags.
    std::vector<TString> filenames_data;
    if (is_data) {
        filenames_data.push_back("/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks/Run3_btagWP712_template_for_fit_histos_3D_data_f_80_9999_2MCGEN.root");
        // filenames_data.push_back("/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks/Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN_upartv2.root"); // Zoe
    } else {
        filenames_data = aggChunkFiles(mc_base, sample, generator, "", "MCGEN", btag_tag, out_tag);
        if (filenames_data.empty()) return;
    }

    std::cout << "Getting data from " << filenames_data.size() << " file(s), first: "
              << filenames_data.front() << std::endl;
    //Select central pT bin
    int ibin_pt = 2;
    //Print options
    std::cout << "Options:"
              << "\n\tunfoldBayes:" << unfoldBayes
              << "\n\tibin_pt:" << ibin_pt
              << "\n\ttest_mode:" << test_mode
              << "\n\tmultiply_sigfrac:" << multiply_sigfrac
              << "\n\tsplit_test:" << split_test
              << "\n\tapply_purity:" << apply_purity
              << "\n\tscan_niter:" << scan_niter
              << std::endl;



    label += unfoldBayes ? "_bayesian" : "_MI";
    if(test_mode==0)  label += "_full_closure";
    if(test_mode==1)  label += "_split_test";
    if(test_mode==2)  label += "_data";
    if(!apply_purity) label += "_nopurity";
    TString fout_name;
    fout_name = folder + "histos_" + label + "_after_unfolding_2D.root";

    // ---------------- Plotting setup ------------
    gSystem->Load("libRooUnfold.so");
    gStyle->SetErrorX(0.5);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    TGaxis::SetMaxDigits(3);   // "x10^-3" on the axis instead of 0.000123 on every label

    const Float_t  font_scale = 1200. / 800.;
    const Style_t  font_code  = 43;
    const Float_t  label_size = 15. * font_scale;
    const Float_t  title_size = 15. * font_scale;
    const Float_t  legend_size = 15. * font_scale;
    // One offset for every axis title, so they all sit the same distance from their numbers.
    const Float_t  title_offset = 3.0; // 1.0 not enough 

    // ---- Grab response matrix + corrections
    std::cout << "Getting response + corrections from " << filenames_response.size()
              << " file(s), first: " << filenames_response.front() << std::endl;
    std::vector<TFile *> fin_unfolding = openAllOrWarn(filenames_response);
    if (fin_unfolding.empty()) return;

    // hadd (and any summing merger) sums the per-block RATIO histograms across blocks,
    // giving N_blocks x the true value (e.g. hEECweightEff ~ 10 x 1.67 ~ 17). So DON'T read
    // the pre-divided ratios; recompute every correction here from the COUNT histograms,
    // which sum correctly. That applies to summing bjet onto qcd as well.
    // opt "b" = binomial (num is a subset of den), "" = normal errors.
    auto ratioFromCounts = [&](const char* num, const char* den, const char* opt) -> TH2D* {
        TH2D* hn = sumOverFiles<TH2D>(fin_unfolding, num);
        TH2D* hd = sumOverFiles<TH2D>(fin_unfolding, den);
        if (!hn || !hd) return nullptr;
        TH2D* r = (TH2D*) hn->Clone(Form("ratio_%s_over_%s", num, den));
        r->SetDirectory(0);
        r->Divide(hn, hd, 1., 1., opt);
        return r;
    };

    // ----------- Grab data -----------
    std::vector<TFile *> fin_data = openAllOrWarn(filenames_data);
    if (fin_data.empty()) return;
    TString histname = (test_mode == 0) ? "h3D_bb"
                     : (test_mode == 1) ? "h3D_pseudodata_bb"
                     :                     "h3D_data";
    std::cout << "Using input histogram: " << histname << std::endl;
    TH3D *h_data_reco_3D_in = sumOverFiles<TH3D>(fin_data, histname);
    if (!h_data_reco_3D_in) {
        std::cerr << "       re-run create_files_for_template_fit.cpp to produce it." << std::endl;
        return;
    }
    TH3D *h_data_reco_3D = (TH3D*) h_data_reco_3D_in->Clone("h_data_reco_3D"); //att
    TH2D *h_data_reco = (TH2D*)h_data_reco_3D->Project3D("zy");
    TH2D *h_data_after_fit = (TH2D*) h_data_reco->Clone("h_data_after_fit");



    // Multiply histograms by signal fraction

    if (multiply_sigfrac) {
        std::cout << "\t---->Multiplying by signal fraction" << std::endl;
        // Grab signal fraction from template fit
        TString fname_fit = filename_template_fit;
        std::cout << "Getting signal fraction from " << fname_fit << std::endl;
        TFile *fin_fit = openOrWarn(fname_fit);
        if (!fin_fit) return;
        TH2D *h_sig_fraction = getOrWarn<TH2D>(fin_fit, "h_sig_fraction_fit");
        if (!h_sig_fraction) return;
        h_data_after_fit->Multiply(h_sig_fraction);

    }
    else {
        std::cout << "\t---->Not multiplying by signal fraction" << std::endl;
    }

    // ---- UParT b-tag efficiency scale factor, on the RECO-level data ------------------
    // Applied HERE, before unfolding, and not to the unfolded result: the SF is measured in
    // reco R_BB bins, so it belongs on the distribution that is still in reco R_BB. The
    // response matrix then carries it through the migration like any other feature of the
    // data. Applying it afterwards would put a reco-binned correction on particle-level
    // bins and silently ignore that migration.
    //
    // It is DIVIDED in: the chain corrects the data with the efficiency measured in MC, and
    // the true efficiency is eps_data = SF * eps_MC, so data/eps_MC has to be further
    // divided by SF. Data only -- there is no data in the closure modes, and the SF must
    // never touch MC.
    if (is_data && sfupartVariation == "off") {
        std::cout << "\t---->UParT SF NOT applied (SFUPART_VARIATION=off): this is the BEFORE-SF "
                  << "result, for comparison only" << std::endl;
    }
    else if (is_data) {
        const TString sf_hist = sfupartHist(sfupartVariation);
        std::cout << "\t---->Dividing reco data by the UParT SF (" << sfupartVariation
                  << " -> " << sf_hist << ")" << std::endl;
        TFile *fin_sf = openOrWarn(sfupartFile());
        if (!fin_sf) return;
        TH1D *h_sf = getOrWarn<TH1D>(fin_sf, sf_hist);
        if (!h_sf) return;

        const int nx = h_data_after_fit->GetNbinsX();   // x is dr, y is pt
        if (h_sf->GetNbinsX() != nx) {
            std::cerr << "ERROR: UParT SF has " << h_sf->GetNbinsX() << " dr bins but the "
                      << "reco data has " << nx << std::endl;
            return;
        }

        // Bin by bin: read the content, multiply by the INVERSE of the SF in that dr bin,
        // write it back. The table prints the reco content before and after for the pT bin
        // the analysis uses, so the correction can be checked by eye against the SF column.
        printf("\n  UParT SF applied bin by bin (reco level, pT bin %d)\n", ibin_pt);
        printf("  %3s %13s %9s %9s %14s %14s %8s\n",
               "bin", "dr range", "SF", "1/SF", "before", "after", "ratio");
        printf("  ---------------------------------------------------------------------------\n");

        for (int ix = 1; ix <= nx; ++ix) {
            const double sf  = h_sf->GetBinContent(ix);
            const double esf = h_sf->GetBinError(ix);
            if (sf == 0.) {
                std::cerr << "ERROR: UParT SF is zero in dr bin " << ix << std::endl;
                return;
            }
            const double inv_sf = 1. / sf;       // the number every bin is multiplied by
            const double rel_sf = esf / sf;

            const double before = h_data_after_fit->GetBinContent(ix, ibin_pt);

            for (int iy = 1; iy <= h_data_after_fit->GetNbinsY(); ++iy) {
                const double c = h_data_after_fit->GetBinContent(ix, iy);
                const double e = h_data_after_fit->GetBinError(ix, iy);
                const double c_new = c * inv_sf;
                // The SF's own statistical error rides along into the data's error, so it
                // ends up in the unfolded result's statistical error. It is uncorrelated
                // between dr bins (each is a separate calibration measurement), which is
                // what a bin-by-bin combination assumes.
                const double e_new = (c != 0.)
                    ? std::fabs(c_new) * std::sqrt((e / c) * (e / c) + rel_sf * rel_sf)
                    : e * inv_sf;
                h_data_after_fit->SetBinContent(ix, iy, c_new);
                h_data_after_fit->SetBinError(ix, iy, e_new);
            }

            const double after = h_data_after_fit->GetBinContent(ix, ibin_pt);
            printf("  %3d [%5.2f,%5.2f] %9.5f %9.5f %14.6g %14.6g %8.4f\n", ix,
                   h_sf->GetXaxis()->GetBinLowEdge(ix), h_sf->GetXaxis()->GetBinUpEdge(ix),
                   sf, inv_sf, before, after, (before != 0. ? after / before : 0.));
        }
        printf("\n");
        fin_sf->Close();
    }



     //Dimension of the response matrix
    int dim = bins_pt*bins_dr;
    int ibin_dr_min = 1;
    int ibin_dr_max = bins_dr;

    std::cout << "dim = " << dim << std::endl;
    std::cout << "bins_pt = " << bins_pt << std::endl;
    std::cout << "bins_dr = " << bins_dr << std::endl;

    // Note: Result = unfold(raw * purity) * 1 / (efficiency)
    //       fakes are negligible

    


    TH2D *h_full_purity = nullptr;
    TH2D *h_full_efficiency = nullptr;
    TH2D *h_mc_reco = nullptr;
    RooUnfoldResponse *response = nullptr;
    TH2D *h_mc_true_no_eff = nullptr;

    // The reco-level MC comparison must match what the data curve has had done to it:
    // purity-corrected data -> matched reco (numerator); uncorrected data -> all reco (denominator).
    if (split_test){
        std::cout << "\t----> Doing split test" << std::endl;
        // --- 
        // h_full_purity = getOrWarn<TH2D>(fin_unfolding, "h_full_pseudo_purity_tf"); // wrong ratio due to merged files 
        // h_full_efficiency = getOrWarn<TH2D>(fin_unfolding, "h_full_pseudo_efficiency_tf"); // wrong ratio due to merged files 
        // --- 

        h_full_purity = ratioFromCounts("h_full_pseudo_purity_numerator_tf", "h_full_pseudo_purity_denominator_tf", "h_pseudo_purity");
        h_full_efficiency = ratioFromCounts("h_full_pseudo_efficiency_numerator_tf", "h_full_pseudo_efficiency_denominator_tf", "h_pseudo_efficiency");
        
        h_mc_reco = sumOverFiles<TH2D>(fin_unfolding, apply_purity ? "h_full_pseudo_purity_numerator_tf"
                                                                  : "h_full_pseudo_purity_denominator_tf");
        response = sumResponseOverFiles(fin_unfolding, "response_tf_pseudo_full");
        h_mc_true_no_eff = sumOverFiles<TH2D>(fin_unfolding, "h_full_pseudo_efficiency_numerator_tf");
    }
    else {
        // ---  
        // h_full_purity = getOrWarn<TH2D>(fin_unfolding, "h_full_purity_tf"); // Wrong ratio due to merged files 
        // h_full_efficiency = getOrWarn<TH2D>(fin_unfolding, "h_full_efficiency_tf"); // wrong ratio due to merged files 
        // --- 

        h_full_purity     = ratioFromCounts("h_full_purity_numerator_tf", "h_full_purity_denominator_tf", "h_purity");
        h_full_efficiency = ratioFromCounts("h_full_efficiency_numerator_tf", "h_full_efficiency_denominator_tf", "h_efficiency");

        h_mc_reco = sumOverFiles<TH2D>(fin_unfolding, apply_purity ? "h_full_purity_numerator_tf"
                                                                  : "h_full_purity_denominator_tf");
        response = sumResponseOverFiles(fin_unfolding, "response_tf_full");
        h_mc_true_no_eff = sumOverFiles<TH2D>(fin_unfolding, "h_full_efficiency_numerator_tf");
    }
    if (!h_full_purity || !h_full_efficiency || !h_mc_reco || !response) return;



    // ---- Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    std::cout << "\t---->Condition number nominal = " << cond_number
              << std::endl;

    // ---- Grab the truth level MC ----
    // Must be the SAME production as the response files: the truth, the response and the
    // purity/efficiency corrections have to come from one set of events with one pT
    // binning, otherwise the full-MC closure (test_mode 0) cannot be exactly 1.
    std::vector<TFile *> &fin_response_truth = fin_unfolding;
    std::cout << "Getting truth from the same " << fin_response_truth.size()
              << " response file(s)" << std::endl;
    TH2D *h_mc_true = nullptr;
    if (test_mode == 1) {
        // Split test: the even half's gen distribution -- same gen_pass gate and same w_reco
        // weight the correction chain outputs. Independent from the odd-half corrections.
        h_mc_true = sumOverFiles<TH2D>(fin_response_truth, "h_pseudodata_truth_tf");
    } else if (test_mode == 0) {
        // Full-MC closure: the full-sample all-gen distribution the chain recovers by construction.
        h_mc_true = sumOverFiles<TH2D>(fin_response_truth, "h_full_efficiency_denominator_tf");
    } else {
        // Data: gen reference. After the combined SV-reco + b-tag correction below, the
        // result is at the "all true 2b" level, so compare against hgenjet_2b_all (the
        // combined-efficiency denominator), not hgenjet_2b_passbtag (2SV+btag level).
        h_mc_true = sumOverFiles<TH2D>(fin_response_truth, "hgenjet_2b_all");
    }
    if (!h_mc_true) return;


    //------- Apply purity correction
    TH2D *h_data_purity_corrected = (TH2D *) h_data_after_fit->Clone("h_data_purity_corrected");
    if (apply_purity) {
        std::cout << "\t---->Multiplying data by purity" << std::endl;
        h_data_purity_corrected->Multiply(h_full_purity);
    } else {
        std::cout << "\t---->NOT multiplying data by purity" << std::endl;
    }

    //double underflowY = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1, 0, 0);
    //double overflowY  = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1,
                          //h_data_purity_corrected->GetNbinsY()+1, h_data_purity_corrected->GetNbinsY()+1);
    //std::cout << "underflowY=" << underflowY << " overflowY=" << overflowY << std::endl;

    // ---- Scan the Bayesian iteration count, one PDF page per value (only when scan_niter).
    // ---- disabled (kept for reference): iteration scan output in the LLR home dir ----
    // TString home_llr = "/home/llr/cms/zaidan/";
    // label already encodes the option flags (_bayesian / _split_test / _nopurity), so toggling
    // those bools writes to distinct files instead of overwriting the previous run's outputs.
    TString iter_pdf = folder + "bayesian_unfolding_iterations_" + label + ".pdf";  // multi-page PDF of all iterations
    TString png_dir  = folder + "bayesian_unfolding_iterations_" + label + "_png/"; // one PNG per iteration lands here
    if (scan_niter) gSystem->mkdir(png_dir, kTRUE);                     // create the PNG folder if it does not exist
    const int niter_min = 0;
    const int niter_max = scan_niter ? 49 : 0;   // scan niter = 1..niter_max+1; one pass when not scanning
    const Int_t prev_ignore = gErrorIgnoreLevel;
    if (scan_niter) gErrorIgnoreLevel = kError;  // silence "Replacing existing histogram" spam from re-cloning each pass

    // Refolding goodness-of-fit vs iteration (see refoldChi2 above). Filled on every pass;
    // only used for the regularisation plot when scan_niter is on.
    // The optimal iteration is the FIRST one whose refolding p-value clears this threshold:
    // earlier iterations are still non-closing, later ones only overfit the input further.
    const double pvalue_threshold = 0.05;
    std::vector<double> v_niter, v_pvalue_all, v_pvalue_pt, v_chi2ndf_all, v_chi2ndf_pt;

    for (int iter = niter_min; iter <= niter_max; ++iter) {

    // Single non-scan pass keeps the nominal 7 iterations; the scan walks 1..100.
    const Int_t niter_now = scan_niter ? iter + 1 : 7;

    // ---- Unfold
    std::cout << "\t---->Unfolding (niter = " << niter_now << ")" << std::endl;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    TH2D *h_data_unfolded = nullptr;
    TMatrixD covariance_matrix_before_unfolding(dim,dim);
    TMatrixD covariance_matrix_after_unfolding(dim,dim);
    if (unfoldBayes) {
        RooUnfoldBayes unfold(response, h_data_purity_corrected, niter_now);
        // Clone before `unfold` leaves scope: some RooUnfold versions hand back a cached
        // histogram that the unfolder owns and deletes with itself.
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment)->Clone("h_data_unfolded");
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    } else {
        RooUnfoldInvert unfold(response, h_data_purity_corrected);
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment)->Clone("h_data_unfolded");
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    }

    // ---- Fold back
    std::cout << "\t---->Refolding" << std::endl;
    TH2D *h_data_refolded = (TH2D *) response->ApplyToTruth(h_data_unfolded, "h_data_refolded");

    // ---- Refolding goodness of fit, against the distribution that went INTO the unfolding
    // (h_data_purity_corrected -- the same object handed to RooUnfoldBayes above), not the
    // raw h_data_reco: the refolded distribution lives at the purity-corrected level too.
    // "all" = the whole 2D reco space that was unfolded; "pt" = the pT bin actually plotted.
    const RefoldGof gof_all = refoldChi2(h_data_refolded, h_data_purity_corrected,
                                         ibin_dr_min, ibin_dr_max, 1, bins_pt);
    const RefoldGof gof_pt  = refoldChi2(h_data_refolded, h_data_purity_corrected,
                                         ibin_dr_min, ibin_dr_max, ibin_pt, ibin_pt);
    std::cout << Form("\t---->Refolding GoF  niter = %3d | 2D: chi2/ndf = %8.2f/%3d = %6.3f, "
                      "p = %8.3e | pT bin %d: chi2/ndf = %8.2f/%3d = %6.3f, p = %8.3e",
                      niter_now,
                      gof_all.chi2, gof_all.ndf, gof_all.chi2ndf(), gof_all.pvalue,
                      ibin_pt, gof_pt.chi2, gof_pt.ndf, gof_pt.chi2ndf(), gof_pt.pvalue)
              << std::endl;
    v_niter.push_back(niter_now);
    v_pvalue_all.push_back(gof_all.pvalue);
    v_pvalue_pt.push_back(gof_pt.pvalue);
    v_chi2ndf_all.push_back(gof_all.chi2ndf());
    v_chi2ndf_pt.push_back(gof_pt.chi2ndf());

    // ---- Apply efficiency correction
    std::cout << "\t---->Dividing by recostruction efficiency" << std::endl;
    TH2D *h_data_efficiency_corrected = (TH2D *) h_data_unfolded->Clone("h_data_efficiency_corrected");
    h_data_efficiency_corrected->Divide(h_full_efficiency);

    // ---- Final corrections
    TH2D *h_data_fully_corrected = (TH2D *) h_data_efficiency_corrected->Clone("h_data_fully_corrected");

    // ---- Combined SV-reconstruction + b-tag efficiency correction (data only) ----
    // Peels BOTH the 2-SV reconstruction and the b-tag selection off the truth level,
    // taking the result from "2b jets reconstructed with 2 SV and tagged" to "all true 2b".
    // Orthogonal to the reconstruction-kinematic efficiency (h_full_efficiency) applied above.
    // Skipped in closure modes (0/1), whose truth targets live at the 2SV+btag level.
    if (test_mode == 2) {
        TH2D *h_svbtag_eff = ratioFromCounts("hgenjet_2b_reco_btag", "hgenjet_2b_all", "b");
        if (h_svbtag_eff) {
            std::cout << "\t---->Dividing by combined SV-reco + b-tag efficiency" << std::endl;
            // systematic hook: scale h_svbtag_eff by the CMS b-tag SF map here, then vary SF +/-.
            h_data_fully_corrected->Divide(h_svbtag_eff);
        }
        // EEC-weight correction: convert the reco-EEC-weighted result to gen-EEC-weighted.
        // r_eec = sum(eec_gen)/sum(eec_reco) over reconstructed 2b jets -> MULTIPLY.
        TH2D *h_eec_weight_eff = ratioFromCounts("hgenjet_2b_reco_btag", "hgenjet_2b_passbtag", "");
        if (h_eec_weight_eff) {
            std::cout << "\t---->Multiplying by EEC-weight (reco->gen) correction" << std::endl;
            h_data_fully_corrected->Multiply(h_eec_weight_eff);
        }
    }


    // ---- Correction-stage plot: raw EEC (top) + incremental ratios (bottom), data only ----
    // Top pad: EEC(dr) after unfolding with each correction stacked on (absolute).
    // Bottom pad: each correction isolated as the fractional effect (ratio - 1); 0 = no effect.
    // Same main/ratio layout, palette and fonts as the bottomline plot.
    if (test_mode == 2) {
        // normalise_corr_stages = true  -> every curve to unit area; the ratio pad shows the
        //   SHAPE effect of each correction. Writes ..._norm.{pdf,png}.
        // false -> absolute EEC; the ratio pad shows the true SIZE of each correction, and the
        //   gen MC is scaled to the data integral as a shape reference. Writes unsuffixed files.
        // Currently OFF, because normalising divides out exactly what this plot is meant to show.
        // To restore it set this true AND uncomment the normalizeToUnitArea loop below.
        const bool normalise_corr_stages = false;
        std::cout << "\t---->Making correction-stage EEC + ratio plot ("
                  << (normalise_corr_stages ? "normalised" : "absolute") << ")" << std::endl;
        TH2D *h_btag_c = ratioFromCounts("hgenjet_2b_passbtag",  "hgenjet_2b",          "b"); // b-tag only, after 2SV
        TH2D *h_svb_c  = ratioFromCounts("hgenjet_2b_reco_btag", "hgenjet_2b_all",      "b"); // 2SV + b-tag combined
        TH2D *h_eec_c  = ratioFromCounts("hgenjet_2b_reco_btag", "hgenjet_2b_passbtag", "");  // EEC weight

        // 2SV alone. No count pair gives it directly (hgenjet_2b / hgenjet_2b_all mixes in the
        // reco/gen EEC-weight ratio, different weights), so factorise the combined efficiency:
        //   eff(2SV) = eff(2SV + b-tag) / eff(b-tag) = h_svb_c / h_btag_c.
        // The errors treat the two ratios as independent when they are not; fine for a diagnostic.
        TH2D *h_2sv_c = nullptr;
        if (h_svb_c && h_btag_c) {
            h_2sv_c = (TH2D*) h_svb_c->Clone("ratio_2sv_only");
            h_2sv_c->SetDirectory(0);
            h_2sv_c->Divide(h_btag_c);
        }

        // Stage order: after unfolding -> + 2SV -> + 2SV + b-tag -> + EEC weight. The combined
        // stage is the same histogram whichever cut is peeled first (h_svb_c covers both), so
        // reordering only changes which single correction is applied on its own in between.
        TH2D *h2_base   = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_base");
        TH2D *h2_2sv    = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_2sv");
        TH2D *h2_svbtag = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_svbtag");
        TH2D *h2_full   = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_full");
        if (h_2sv_c) h2_2sv->Divide(h_2sv_c);
        if (h_svb_c)  h2_svbtag->Divide(h_svb_c);
        if (h_svb_c)  h2_full->Divide(h_svb_c);
        if (h_eec_c)  h2_full->Multiply(h_eec_c);

        TH1D *p_base   = h2_base  ->ProjectionX(Form("p_cmp_base_pt%d",   ibin_pt), ibin_pt, ibin_pt);
        TH1D *p_2sv    = h2_2sv   ->ProjectionX(Form("p_cmp_2sv_pt%d",    ibin_pt), ibin_pt, ibin_pt);
        TH1D *p_svbtag = h2_svbtag->ProjectionX(Form("p_cmp_svbtag_pt%d", ibin_pt), ibin_pt, ibin_pt);
        TH1D *p_full   = h2_full  ->ProjectionX(Form("p_cmp_full_pt%d",   ibin_pt), ibin_pt, ibin_pt);
        // Gen MC at the level the fully corrected result targets -- for test_mode 2 that is
        // hgenjet_2b_all, "all true 2b". This is the curve p_full is meant to land on, so the
        // plot shows not just what each correction does but whether they arrive somewhere sane.
        TH1D *p_true   = h_mc_true->ProjectionX(Form("p_cmp_true_pt%d",   ibin_pt), ibin_pt, ibin_pt);

        int b1 = 1, b2 = p_base->GetNbinsX();
        while (b1 <= p_base->GetNbinsX() && p_base->GetBinContent(b1) <= 0.) ++b1;
        while (b2 >= 1                   && p_base->GetBinContent(b2) <= 0.) --b2;

        // Whatever the mode, do this BEFORE the ratios below so they inherit it.
        for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true}) h->GetXaxis()->SetRange(b1, b2);
        if (normalise_corr_stages) {
            // Unit-area normalisation of every stage -- commented out, see normalise_corr_stages above.
            //for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true}) normalizeToUnitArea(h);
        } else {
            // The gen MC cannot be overlaid absolutely: it carries an MC event weight with no
            // luminosity scaling, so its integral has no relation to the data's. Put it on the
            // data's scale instead -- shape reference only, and the legend says so.
            const double I_data_cmp = p_full->Integral(), I_true_cmp = p_true->Integral();
            if (I_data_cmp > 0. && I_true_cmp > 0.) p_true->Scale(I_data_cmp / I_true_cmp);
        }

        // ratio of each stage to the PREVIOUS one, isolating that single correction:
        //   2SV = 2sv/after-unfolding, +b-tag = svbtag/2sv, +EEC = full/svbtag.  1 = no change.
        TH1D *r_2sv  = (TH1D*) p_2sv   ->Clone(Form("r_2sv_pt%d",  ibin_pt));  r_2sv ->Divide(p_base);
        TH1D *r_btag = (TH1D*) p_svbtag->Clone(Form("r_btag_pt%d", ibin_pt));  r_btag->Divide(p_2sv);
        TH1D *r_eec  = (TH1D*) p_full  ->Clone(Form("r_eec_pt%d",  ibin_pt));  r_eec ->Divide(p_svbtag);

        // ROCColor palette from the bottomline plot; each step keeps its colour top & bottom.
        // Red is the measurement (p_full, the fully corrected result) and blue the gen MC it is
        // compared to, as everywhere else; the intermediate stages take the remaining colours,
        // in the order they are applied (first added correction purple, second orange).
        p_base  ->SetLineColor(green);  p_base  ->SetMarkerColor(green);  p_base  ->SetMarkerStyle(kOpenCircle);
        p_2sv   ->SetLineColor(purple); p_2sv   ->SetMarkerColor(purple); p_2sv   ->SetMarkerStyle(kFullCircle);
        p_svbtag->SetLineColor(orange); p_svbtag->SetMarkerColor(orange); p_svbtag->SetMarkerStyle(kFullSquare);
        p_full  ->SetLineColor(red);    p_full  ->SetMarkerColor(red);    p_full  ->SetMarkerStyle(kFullTriangleUp);
        p_true  ->SetLineColor(blue);   p_true  ->SetMarkerColor(blue);   p_true  ->SetMarkerStyle(kOpenTriangleUp);
        r_2sv ->SetLineColor(purple); r_2sv ->SetMarkerColor(purple); r_2sv ->SetMarkerStyle(kFullCircle);
        r_btag->SetLineColor(orange); r_btag->SetMarkerColor(orange); r_btag->SetMarkerStyle(kFullSquare);
        r_eec ->SetLineColor(red);    r_eec ->SetMarkerColor(red);    r_eec ->SetMarkerStyle(kFullTriangleUp);
        for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true, r_2sv, r_btag, r_eec}) {
            h->SetStats(0); h->SetLineWidth(2); h->SetMarkerSize(1.2); h->GetXaxis()->SetRange(b1, b2);
        }

        double ymax = 0.;
        for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true})
            for (int i = b1; i <= b2; ++i) ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
        double rmin = 1., rmax = 1.;
        for (TH1D* h : {r_2sv, r_btag, r_eec})
            for (int i = b1; i <= b2; ++i) { double v = h->GetBinContent(i); if (v <= 0.) continue; rmin = std::min(rmin, v); rmax = std::max(rmax, v); }
        double rpad = 0.15 * (rmax - rmin) + 1e-6;

        Float_t ptlo = h_data_reco->GetYaxis()->GetBinLowEdge(ibin_pt);
        Float_t pthi = h_data_reco->GetYaxis()->GetBinUpEdge(ibin_pt);

        // canvas + two pads (same main/ratio split as the bottomline plot)
        TCanvas *c_cmp = new TCanvas(Form("c_corr_stages_pt%d", ibin_pt), "", 800, 800);
        TPad *pad_main  = new TPad(Form("pad_cmp_main_%d",  ibin_pt), "", 0., 0.3, 1., 1.);
        TPad *pad_ratio = new TPad(Form("pad_cmp_ratio_%d", ibin_pt), "", 0., 0.,  1., 0.3);
        for (TPad *p : {pad_main, pad_ratio}) { p->SetTicks(1, 0); p->SetFillColor(0); }
        // pad_main ->SetMargin(0.13, 0.05, 0.00, 0.08);
        // pad_ratio->SetMargin(0.13, 0.05, 0.32, 0.00);
        pad_main->SetMargin(0.1, 0.1, 0.0, 0.1);
        pad_ratio->SetMargin(0.1, 0.1, 0.23, 0.0);


        c_cmp->cd(); pad_main->Draw(); pad_ratio->Draw();

        // ----- top: raw EEC(dr) with each correction -----
        pad_main->cd();
        p_base->SetTitle("");
        p_base->GetYaxis()->SetRangeUser(0., ymax * 1.6);
        p_base->GetYaxis()->SetTitle(normalise_corr_stages ? "normalised EEC(#Delta r)"
                                                          : "EEC(#Delta r)");
        p_base->GetYaxis()->CenterTitle(true);
        p_base->GetYaxis()->SetTitleFont(font_code); p_base->GetYaxis()->SetTitleSize(title_size); p_base->GetYaxis()->SetTitleOffset(1.5);
        p_base->GetYaxis()->SetLabelFont(font_code); p_base->GetYaxis()->SetLabelSize(label_size);
        // The y range starts at 0 and this pad's bottom edge is shared with the ratio pad, so the
        // "0" label sits exactly on the join and gets clipped in half. Hide it (label size 0).
        p_base->GetYaxis()->ChangeLabel(1, -1, 0.);
        p_base->GetXaxis()->SetLabelSize(0); p_base->GetXaxis()->SetTitleSize(0);   // x belongs to the ratio pad
        // Every stage, "after unfolding" included, as points (X0 = no horizontal / bin-width
        // error bars), so the starting point reads as one of the stages rather than as a frame.
        p_base->SetLineWidth(2);
        // Gen MC as an outline with its error band, like the bottomline plot, so it reads as the
        // reference rather than as another correction stage.
        p_base->Draw("PE X0"); p_true->Draw("HIST E SAME");
        p_2sv->Draw("PE X0 SAME"); p_svbtag->Draw("PE X0 SAME"); p_full->Draw("PE X0 SAME");

        TLegend *lg = new TLegend(0.17, 0.60, 0.60, 0.85);
        lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetMargin(0.15);
        lg->SetTextFont(font_code); lg->SetTextSize(legend_size * 0.8);
        lg->SetHeader(pthi == jtpt_max ? Form("p_{T}^{jet} > %.0f GeV", ptlo)
                                       : Form("%.0f < p_{T}^{jet} < %.0f GeV", ptlo, pthi));
        lg->AddEntry(p_base,   "after unfolding",            "pe1");
        lg->AddEntry(p_2sv,    "+ 2SV eff.",                 "pe1");
        lg->AddEntry(p_svbtag, "+ 2SV + b-tag eff.",         "pe1");
        lg->AddEntry(p_full,   "+ 2SV + b-tag + EEC weight", "pe1");
        lg->AddEntry(p_true,   normalise_corr_stages ? "Gen MC" : "Gen MC (scaled to data)", "l");
        lg->Draw();

        // CMS Internal above the top axis (outside the frame)
        TLatex cms; cms.SetNDC();
        cms.SetTextFont(62); cms.SetTextSize(0.042); cms.DrawLatex(0.13, 0.945, "CMS");
        cms.SetTextFont(52); cms.SetTextSize(0.034); cms.DrawLatex(0.235, 0.945, "Internal");
        cms.SetTextFont(42); cms.SetTextSize(0.034); cms.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");
        pad_main->RedrawAxis();

        // ----- bottom: incremental ratios (ratio - 1), colours matched to the added stage -----
        pad_ratio->cd();
        // r_2sv is drawn first, so it owns this pad's axes: style them on it and nowhere else.
        r_2sv->SetTitle("");
        r_2sv->GetYaxis()->SetRangeUser(rmin - rpad, rmax + rpad);
        // Short enough to fit the ratio pad's height: "ratio to prev. stage" runs past the pad
        // edge and gets clipped mid-word.
        r_2sv->GetYaxis()->SetTitle("ratio to prev.");
        r_2sv->GetYaxis()->CenterTitle(true);
        r_2sv->GetYaxis()->SetTitleFont(font_code); r_2sv->GetYaxis()->SetTitleSize(title_size); r_2sv->GetYaxis()->SetTitleOffset(1.5);
        r_2sv->GetYaxis()->SetLabelFont(font_code); r_2sv->GetYaxis()->SetLabelSize(label_size);
        r_2sv->GetYaxis()->SetNdivisions(505);
        r_2sv->GetXaxis()->SetTitle("#Delta r");
        r_2sv->GetXaxis()->CenterTitle(true);

        // Precision-43 offsets scale off the (short) ratio-pad height: 3.2 put the title clean
        // off the bottom of the canvas, so no x title was drawn at all. 1.3 lands it under the labels.
        r_2sv->GetXaxis()->SetTitleFont(font_code); r_2sv->GetXaxis()->SetTitleSize(title_size); r_2sv->GetXaxis()->SetTitleOffset(title_offset);
        r_2sv->GetXaxis()->SetLabelFont(font_code); r_2sv->GetXaxis()->SetLabelSize(label_size);
        r_2sv->Draw("PE X0"); r_btag->Draw("PE X0 SAME"); r_eec->Draw("PE X0 SAME");

        // Unity reference. A cloned histogram drawn with "HIST L" joins BIN CENTRES, so it stopped
        // half a bin short of each end of the axis; a TLine spans the plotted range edge to edge.
        TLine *ref1 = new TLine(r_2sv->GetXaxis()->GetBinLowEdge(b1), 1.,
                                r_2sv->GetXaxis()->GetBinUpEdge(b2),  1.);
        ref1->SetLineColor(kGray+2); ref1->SetLineStyle(2); ref1->SetLineWidth(1);
        ref1->Draw("same");

        // labels for the three isolated corrections (colours match the added stage above)
        TLegend *lg_r = new TLegend(0.15, 0.84, 0.95, 0.99);
        lg_r->SetNColumns(3);
        lg_r->SetFillStyle(0); lg_r->SetBorderSize(0); lg_r->SetMargin(0.12);
        lg_r->SetTextFont(font_code); lg_r->SetTextSize(legend_size * 0.8);
        lg_r->AddEntry(r_2sv,  "2SV",  "pe1");
        lg_r->AddEntry(r_btag, "+b-tag", "pe1");
        lg_r->AddEntry(r_eec,  "+EEC weight", "pe1");
        lg_r->Draw();
        pad_ratio->RedrawAxis();

        // Suffix so the normalised and absolute versions do not overwrite each other.
        TString cmp_stem = folder + "correction_stages_" + label + Form("_pt%d", ibin_pt)
                         + (normalise_corr_stages ? "_norm" : "");
        c_cmp->Print(cmp_stem + ".pdf");
        c_cmp->Print(cmp_stem + ".png");
    }


    // ---- DEBUG: reco SV pt vs gen B pt, to explain the size of the EEC-weight correction ----
    if (test_mode == 2) {
        TH2D *h_svpt = sumOverFiles<TH2D>(fin_unfolding, "h_svpt_vs_bpt");
        TH1D *h_ptr  = sumOverFiles<TH1D>(fin_unfolding, "h_ptratio");
        if (h_svpt) {
            std::cout << "\t---->Making DEBUG reco-vs-gen pt plot" << std::endl;
            gStyle->SetPalette(kViridis);
            TCanvas *c_dbg = new TCanvas("c_debug_svpt", "", 850, 700);
            c_dbg->SetLeftMargin(0.13); c_dbg->SetRightMargin(0.18); c_dbg->SetTicks(1, 1);
            h_svpt->SetStats(0);
            h_svpt->SetTitle(";gen partial-B p_{T} [GeV];reco SV p_{T} [GeV];weighted entries");
            h_svpt->GetXaxis()->SetTitleFont(font_code); h_svpt->GetXaxis()->SetTitleSize(title_size); h_svpt->GetXaxis()->SetTitleOffset(1.2);
            h_svpt->GetXaxis()->SetLabelFont(font_code); h_svpt->GetXaxis()->SetLabelSize(label_size);
            h_svpt->GetYaxis()->SetTitleFont(font_code); h_svpt->GetYaxis()->SetTitleSize(title_size); h_svpt->GetYaxis()->SetTitleOffset(1.5);
            h_svpt->GetYaxis()->SetLabelFont(font_code); h_svpt->GetYaxis()->SetLabelSize(label_size);
            h_svpt->GetZaxis()->SetTitleFont(font_code); h_svpt->GetZaxis()->SetTitleSize(title_size); h_svpt->GetZaxis()->SetTitleOffset(1.4);
            h_svpt->GetZaxis()->SetLabelFont(font_code); h_svpt->GetZaxis()->SetLabelSize(label_size);
            h_svpt->Draw("COLZ");
            // mean reco/gen pt on an opaque white box (readable over the colour map)
            if (h_ptr) {
                TLegend *lg_d = new TLegend(0.16, 0.80, 0.50, 0.87);
                lg_d->SetFillColor(kWhite); lg_d->SetFillStyle(1001); lg_d->SetBorderSize(1);
                lg_d->SetTextFont(font_code); lg_d->SetTextSize(legend_size * 0.75);
                lg_d->AddEntry((TObject*)0, Form("#LTp_{T}^{reco}/p_{T}^{gen}#GT = %.2f", h_ptr->GetMean()), "");
                lg_d->Draw();
            }
            // CMS Internal above the top axis
            TLatex cmsd; cmsd.SetNDC();
            cmsd.SetTextFont(62); cmsd.SetTextSize(0.042); cmsd.DrawLatex(0.13, 0.945, "CMS");
            cmsd.SetTextFont(52); cmsd.SetTextSize(0.034); cmsd.DrawLatex(0.235, 0.945, "Internal");
            c_dbg->Print(folder + "debug_svpt_vs_bpt_" + label + ".pdf");
            c_dbg->Print(folder + "debug_svpt_vs_bpt_" + label + ".png");
        }
    }


    // ---- Graphical bottomline test

    std::cout << "Performing graphical bottomline test" << std::endl;
    TH1D *h_mc_reco_2D = h_mc_reco->ProjectionX("h_mc_reco_2D", ibin_pt, ibin_pt);
    TH1D *h_mc_true_2D = h_mc_true->ProjectionX("h_mc_true_2D", ibin_pt, ibin_pt);
    TH1D *h_data_purity_corrected_2D = h_data_purity_corrected->ProjectionX("h_data_purity_corrected_2D", ibin_pt, ibin_pt);
    TH1D *h_data_fully_corrected_2D = h_data_fully_corrected->ProjectionX("h_data_fully_corrected_2D", ibin_pt, ibin_pt);
    TH1D *h_data_refolded_2D = h_data_refolded->ProjectionX("h_data_refolded_2D", ibin_pt, ibin_pt);
    TH1D *h_data_after_fit_2D = h_data_after_fit->ProjectionX("h_data_after_fit_2D", ibin_pt, ibin_pt);
    TH1D *h_data_unfolded_2D = h_data_unfolded->ProjectionX("h_data_unfolded_2D", ibin_pt, ibin_pt);

    // A Divide() across incompatible binnings returns silently and leaves the ratio empty.
    if (h_mc_true_2D->GetNbinsX() != h_data_fully_corrected_2D->GetNbinsX()) {
        std::cerr << "ERROR: truth has " << h_mc_true_2D->GetNbinsX() << " dr bins but the unfolded result has "
                  << h_data_fully_corrected_2D->GetNbinsX() << std::endl;
        return;
    }

    if (true) {

        double ymax = 0.;
        for (auto h : {
                    h_mc_reco_2D,
                    h_mc_true_2D,
                    h_data_purity_corrected_2D,
                    h_data_fully_corrected_2D,
                    h_data_refolded_2D,
                    h_data_after_fit_2D,
                    h_data_unfolded_2D
                    }) {
                        if (!h) continue;
                        h->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
                        normalizeToUnitArea(h);
                    }
        // The y range has to clear the error bars of the curves that are actually drawn, not
        // just their central values, and not the ones left out of the plot.
        for (auto h : {h_mc_reco_2D, h_mc_true_2D, h_data_purity_corrected_2D, h_data_fully_corrected_2D}) {
            for (int i = ibin_dr_min; i <= ibin_dr_max; ++i)
                ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
        }

        Float_t pt_min_plot = h_data_reco->GetYaxis()->GetBinLowEdge(ibin_pt);
        Float_t pt_max_plot = h_data_reco->GetYaxis()->GetBinUpEdge(ibin_pt);

        TLegend *leg = new TLegend(0.15, 0.55, 0.55, 0.80);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetMargin(0.15);
        leg->SetTextFont(font_code);
        leg->SetTextSize(legend_size);
        // The last pT bin is open-ended (jets above jtpt_max are folded into it by
        // jtpt_fill()), so quote it as a threshold rather than a closed range.
        leg->SetHeader(pt_max_plot == jtpt_max
                       ? Form("p_{T}^{jet} > %.0f GeV", pt_min_plot)
                       : Form("%.0f < p_{T}^{jet} < %.0f GeV", pt_min_plot, pt_max_plot));

        // Name every curve by what it is and, in the split test, which half it came from.
        TString lbl_reco     = split_test ? "Reco pseudodata"  : "Reco data";
        TString lbl_mc_reco  = split_test ? (apply_purity ? "Reco MC, matched"
                                                          : "Reco MC, all reco")
                                          : "Reco MC";
        TString lbl_unfolded = split_test ? "Unfolded pseudodata"  : "Unfolded data";
        // Same gen distribution, two disjoint samples. Only the half differs.
        TString lbl_mc_true  = split_test ? "Gen MC"  : "Gen MC";

        // Only the first histogram drawn owns the pad's axes; style them there and nowhere else.
        h_data_purity_corrected_2D->SetTitle("");//Data, " + label + " " + " response matrix");
        h_data_purity_corrected_2D->SetStats(0);
        h_data_purity_corrected_2D->SetMarkerColor(purple);
        h_data_purity_corrected_2D->SetLineColor(purple);
        h_data_purity_corrected_2D->SetMarkerStyle(kFullCircle);
        h_data_purity_corrected_2D->SetMarkerSize(1);
        // 1.75x headroom keeps the legend clear of the EEC peak.
        h_data_purity_corrected_2D->GetYaxis()->SetRangeUser(0., ymax*1.75);
        // TLatex (#Delta), not TMathText (\Delta\mbox{r}): TMathText axis titles do not get
        // placed reliably — the x-axis title was silently dropped altogether.
        h_data_purity_corrected_2D->GetYaxis()->SetTitle("EEC(#Delta r)");
        h_data_purity_corrected_2D->GetYaxis()->CenterTitle(true);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleFont(font_code);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleSize(title_size);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleOffset(1.5);
        h_data_purity_corrected_2D->GetYaxis()->SetLabelFont(font_code);
        h_data_purity_corrected_2D->GetYaxis()->SetLabelSize(label_size);
        h_data_purity_corrected_2D->GetXaxis()->SetLabelSize(0); // x labels belong to the ratio pad
        h_data_purity_corrected_2D->GetXaxis()->SetTitleSize(0);
        leg->AddEntry(h_data_purity_corrected_2D, lbl_reco, "pe1");


        //h_data_after_fit_2D->SetMarkerColor(orange);
        //h_data_after_fit_2D->SetLineColor(orange);
        //h_data_after_fit_2D->SetMarkerStyle(kFullTriangleUp);
        //h_data_after_fit_2D->SetMarkerSize(1);
        //leg->AddEntry(h_data_after_fit_2D, "Before purity correction", "pe1");

        //h_data_unfolded_2D->SetMarkerColor(green);
        //h_data_unfolded_2D->SetLineColor(green);
        //h_data_unfolded_2D->SetMarkerStyle(kFullTriangleUp);
        //h_data_unfolded_2D->SetMarkerSize(1);
        //leg->AddEntry(h_data_unfolded_2D, "After unfolding (no efficiency correction)", "pe1");

        h_mc_reco_2D->SetStats(0);
        h_mc_reco_2D->SetMarkerColor(green);
        h_mc_reco_2D->SetLineColor(green);
        h_mc_reco_2D->SetMarkerStyle(kFullTriangleUp);
        h_mc_reco_2D->SetMarkerSize(1);
        leg->AddEntry(h_mc_reco_2D, lbl_mc_reco, "pe1");

        // Red: this is the measurement. Reserved for it across every plot in the analysis.
        h_data_fully_corrected_2D->SetMarkerColor(red);
        h_data_fully_corrected_2D->SetLineColor(red);
        h_data_fully_corrected_2D->SetMarkerStyle(kOpenCross);
        h_data_fully_corrected_2D->SetMarkerSize(1);
        leg->AddEntry(h_data_fully_corrected_2D, lbl_unfolded, "pe1");

        // Blue: the particle-level MC the measurement is compared to (the "theory" curve).
        h_mc_true_2D->SetMarkerColor(blue);
        h_mc_true_2D->SetLineColor(blue);
        h_mc_true_2D->SetMarkerStyle(kOpenTriangleUp);
        h_mc_true_2D->SetMarkerSize(1);
        leg->AddEntry(h_mc_true_2D, lbl_mc_true, "lp");   // line + open-triangle marker

        // Reco-level data, like h_data_purity_corrected_2D above, so it shares that curve's
        // purple and is told apart by the marker; blue now belongs to the gen MC.
        h_data_refolded_2D->SetMarkerColor(purple);
        h_data_refolded_2D->SetLineColor(purple);
        h_data_refolded_2D->SetMarkerStyle(kFullCross);
        h_data_refolded_2D->SetMarkerSize(1);


      
        
        TCanvas *c_unfold = new TCanvas("c_unfold", "", 800, 800);

        TPad *pad_main  = new TPad("pad_main",  "", 0., 0.3, 1., 1.);   // placeholder coords;
        TPad *pad_ratio = new TPad("pad_ratio", "", 0., 0.,  1., 0.3);  // applySquareLayout owns them
        
    
        
        
        for (TPad *p : {pad_main, pad_ratio}) {
            p->SetTicks(1, 0);   // mirrored ticks on all four sides
            p->SetFillColor(0);
        
        }
        //pad_ratio->SetLogx();
        //pad_main->SetLogx();

        pad_main->SetMargin(0.1, 0.1, 0.0, 0.1);
        pad_ratio->SetMargin(0.1, 0.1, 0.23, 0.0);

    
        // Attach the pads to the canvas before drawing into them, so gPad and the pad's
        // primitive list agree at every Draw().
        c_unfold->cd();
        pad_main->Draw();
        pad_ratio->Draw();

        pad_main->cd();
        h_data_purity_corrected_2D->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
        h_data_purity_corrected_2D->Draw("PE X0");

        //h_data_unfolded_2D->Draw("pe1 same");
        h_mc_reco_2D->Draw("PE X0 same");
        h_data_fully_corrected_2D->Draw("PE X0 same");
        //h_data_after_fit_2D->Draw("pe1 same");
        h_mc_true_2D->Draw("hist E same");

        leg->Draw();
        TLatex *test_info_text = new TLatex;
        
        test_info_text->SetNDC();
        test_info_text->SetTextFont(font_code);
        test_info_text->SetTextSize(label_size);
        test_info_text->DrawLatex(0.57, 0.75,
            test_mode == 0 ? "MC full-sample closure" :
            test_mode == 1 ? "MC split closure test"  :
                             "Data Run 3");
        test_info_text->DrawLatex(0.57, 0.70, "2D unfolding");
        if(unfoldBayes) {
            test_info_text->DrawLatex(0.57, 0.65, "Bayesian unfolding");
            test_info_text->DrawLatex(0.57, 0.60, Form("N iter = %d", niter_now));
            // Refolding GoF of THIS iteration, so every scan page carries its own number.
            // Data only: in the closure modes the input is (mode 0) or comes from (mode 1) the
            // same MC as the response, so this number is a technical check, not a measured GoF,
            // and quoting it on the plot invites reading it as one. It is still printed to stdout
            // for both modes, and still drives the iteration scan.
            if (is_data) {
                test_info_text->DrawLatex(0.57, 0.55,
                    Form("refold #chi^{2}/ndf = %.2f", gof_all.chi2ndf()));
                test_info_text->DrawLatex(0.57, 0.50, Form("p = %.3f", gof_all.pvalue));
                test_info_text->DrawLatex(0.57, 0.45, Form("Response matrix CN = %.1f", cond_number));
            }
            else {
                test_info_text->DrawLatex(0.57, 0.55, Form("Response matrix CN = %.1f", cond_number));
            }
        }
        else {
            test_info_text->DrawLatex(0.57, 0.65, "Matrix inversion unfolding");
            test_info_text->DrawLatex(0.57, 0.60, Form("Response matrix CN = %.1f", cond_number));
        }
        //drawHeader();

        // CMS Internal above the top axis (outside the frame)
        TLatex cms_bl; cms_bl.SetNDC();
        cms_bl.SetTextFont(62); cms_bl.SetTextSize(0.045); cms_bl.DrawLatex(0.10, 0.945, "CMS"); // cms_bl.SetTextSize(0.042); 
        cms_bl.SetTextFont(52); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.205, 0.945, "Internal"); // cms_bl.SetTextSize(0.034); 
        cms_bl.SetTextFont(42); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV"); // cms_bl.SetTextSize(0.034);

        pad_main->RedrawAxis();

        
        TH1D *h_data_mc_reco_ratio = (TH1D *) h_data_purity_corrected_2D->Clone("h_data_mc_reco_ratio");
        
        h_data_mc_reco_ratio->SetTitle("");
        h_data_mc_reco_ratio->SetStats(0);
        h_data_mc_reco_ratio->Divide(h_mc_reco_2D);
        h_data_mc_reco_ratio->SetMarkerStyle(kFullCircle);
        h_data_mc_reco_ratio->SetMarkerColor(purple);
        h_data_mc_reco_ratio->SetLineColor(purple);
        h_data_mc_reco_ratio->SetMarkerSize(1.1);

        TH1D *h_data_mc_true_ratio = (TH1D *) h_data_fully_corrected_2D->Clone("h_data_mc_true_ratio");
        h_data_mc_true_ratio->Divide(h_mc_true_2D);
        h_data_mc_true_ratio->SetMarkerStyle(kOpenCross);
        h_data_mc_true_ratio->SetMarkerColor(red);
        h_data_mc_true_ratio->SetLineColor(red);
        h_data_mc_true_ratio->SetMarkerSize(1.3);

        TH1D *h_mc_gen_reco_ratio = (TH1D *) h_mc_true_2D->Clone("h_mc_gen_reco_ratio");


        h_mc_gen_reco_ratio->Divide(h_mc_reco_2D);
        // Each ratio carries the colour of its NUMERATOR in the main pad: this one is the gen MC,
        // so blue, matching the gen MC curve above.
        h_mc_gen_reco_ratio->SetLineColor(blue);
        h_mc_gen_reco_ratio->SetLineWidth(1);
        h_mc_gen_reco_ratio->SetMarkerSize(0);
                h_mc_gen_reco_ratio->GetYaxis()->SetTitle("ratio");
                h_mc_gen_reco_ratio->GetXaxis()->SetTitle("#Delta r");
                        h_mc_gen_reco_ratio->GetYaxis()->CenterTitle(true);
                        h_mc_gen_reco_ratio->GetXaxis()->CenterTitle(true);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleSize(0.2); 
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleOffset(0.28);

        h_mc_gen_reco_ratio->GetXaxis()->SetLabelSize(0.06);
        h_mc_gen_reco_ratio->GetYaxis()->SetLabelSize(0.1);
        h_mc_gen_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleOffset(0.24);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleSize(0.2);
        h_mc_gen_reco_ratio->GetYaxis()->SetLimits(0.5, 2);
        h_mc_gen_reco_ratio->SetLineStyle(2);
        h_mc_gen_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);

        // Match the ratio-pad axes to the main pad: pixel (precision-43) fonts with the same
        // title size, label size and title offset, so every axis name is the same size and the
        // same distance from its numbers. Overrides the fractional sizes set just above.
        for (TAxis *ax : { h_mc_gen_reco_ratio->GetXaxis(), h_mc_gen_reco_ratio->GetYaxis() }) {
            ax->SetTitleFont(font_code);
            ax->SetTitleSize(title_size);
            ax->SetTitleOffset(title_offset); //title_offset = 1.0 default   not enough for the unfolding plot 
            ax->SetLabelFont(font_code);
            ax->SetLabelSize(label_size);
        }


        // Orange, not red: its numerator is the unfolded data, but red is already drawn in this pad
        // as the unfolded/gen-MC points, and a dashed red line next to red crosses reads as the same
        // quantity. Green and purple are taken by the reco MC and the reco data, so orange it is.
        TH1D *h_unfolded_reco_ratio = (TH1D *) h_data_fully_corrected_2D->Clone("h_unfolded_reco_ratio");
        h_unfolded_reco_ratio->Divide(h_data_purity_corrected_2D);
        h_unfolded_reco_ratio->SetLineColor(orange);
        h_unfolded_reco_ratio->SetLineWidth(1);
        h_unfolded_reco_ratio->SetMarkerSize(0);
        h_unfolded_reco_ratio->SetLineStyle(2);
        

        h_data_mc_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);
        h_data_mc_reco_ratio->GetYaxis()->SetTitle("ratio");
        h_data_mc_reco_ratio->GetYaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleFont(font_code);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleSize(title_size);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleOffset(0.05);
        h_data_mc_reco_ratio->GetYaxis()->SetLabelFont(font_code);
        h_data_mc_reco_ratio->GetYaxis()->SetLabelSize(label_size);
        h_data_mc_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_data_mc_reco_ratio->GetXaxis()->SetTitle("#Delta r");
        h_data_mc_reco_ratio->GetXaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleFont(font_code);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleSize(title_size);
        
        // With precision-3 fonts this offset scales off the (short) pad height: 3.2 pushed the
        // title clean off the canvas. 1.35 lands it just under the tick labels.
        h_data_mc_reco_ratio->GetXaxis()->SetTitleOffset(0.05);
        h_data_mc_reco_ratio->GetXaxis()->SetLabelFont(12);
        h_data_mc_reco_ratio->GetXaxis()->SetLabelSize(16);
        h_data_mc_reco_ratio->GetXaxis();
       

        // -- Add style also for h_mc_gen_reco_ratio (drawn first in ratio)
        h_mc_gen_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitle("ratio");
        h_mc_gen_reco_ratio->GetYaxis()->CenterTitle(true);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleFont(font_code);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleSize(title_size);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleOffset(1.5);
        h_mc_gen_reco_ratio->GetYaxis()->SetLabelFont(font_code);
        h_mc_gen_reco_ratio->GetYaxis()->SetLabelSize(label_size);
        h_mc_gen_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitle("#Delta r");
        h_mc_gen_reco_ratio->GetXaxis()->CenterTitle(true);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleFont(font_code);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleSize(title_size);

       

        

        // The unity line and the gap bands span the axis actually plotted, not the full dr
        // range: the two disagree as soon as ibin_dr_min/ibin_dr_max restrict the view.
        const Double_t x_lo = h_mc_reco_2D->GetXaxis()->GetBinLowEdge(ibin_dr_min);
        const Double_t x_hi = h_mc_reco_2D->GetXaxis()->GetBinUpEdge(ibin_dr_max);

        pad_ratio->cd();
        //h_data_mc_reco_ratio->Draw("axis");   // axes first, so the bands land behind everything

      
        TLine *line = new TLine(x_lo, 1., x_hi, 1.);
        line->SetLineWidth(1.);
        line->SetLineStyle(kSolid);   // solid, so the unity line reads apart from the dashed blue/orange curves
        line->SetLineColor(kGray+2);

        // "Before" band: for each dr bin, fill the full bin width between the gen MC/reco MC line
        // (blue) and the unfolded/reco line (orange) -- the size of the correction. Stepped closed
        // polygon: forward along blue at the bin edges, back along orange, so each bin is filled
        // edge to edge. Drawn behind, translucent light grey.
        const int nb_band = ibin_dr_max - ibin_dr_min + 1;
        TGraph *band_before_after = new TGraph(4 * nb_band);
        int ip = 0;
        for (int i = 0; i < nb_band; ++i) {              // top edge: along blue, left -> right
            const int bin = ibin_dr_min + i;
            const double xlo = h_mc_gen_reco_ratio->GetXaxis()->GetBinLowEdge(bin);
            const double xhi = h_mc_gen_reco_ratio->GetXaxis()->GetBinUpEdge(bin);
            const double yr  = h_mc_gen_reco_ratio->GetBinContent(bin);
            band_before_after->SetPoint(ip++, xlo, yr);
            band_before_after->SetPoint(ip++, xhi, yr);
        }
        for (int i = nb_band - 1; i >= 0; --i) {         // bottom edge: along orange, right -> left
            const int bin = ibin_dr_min + i;
            const double xlo = h_unfolded_reco_ratio->GetXaxis()->GetBinLowEdge(bin);
            const double xhi = h_unfolded_reco_ratio->GetXaxis()->GetBinUpEdge(bin);
            const double yg  = h_unfolded_reco_ratio->GetBinContent(bin);
            band_before_after->SetPoint(ip++, xhi, yg);
            band_before_after->SetPoint(ip++, xlo, yg);
        }
        // Solid light grey, NOT SetFillColorAlpha: the PDF backend does not render alpha
        // transparency (it shows on screen in XQuartz but disappears in the PDF), whereas a
        // solid fill renders identically everywhere. Curves are redrawn on top, so opaque is fine.
        band_before_after->SetFillColor(TColor::GetColor(224, 224, 224));
        band_before_after->SetFillStyle(1001);
        band_before_after->SetLineWidth(0);

        h_mc_gen_reco_ratio->Draw("hist");     // axes + blue line
        band_before_after->Draw("f same");     // band behind the curves
        h_mc_gen_reco_ratio->Draw("hist same");// redraw blue crisp over the band
        line->Draw("same");

        h_unfolded_reco_ratio->Draw("hist same");
        h_data_mc_reco_ratio->Draw("PE X0 same");
        h_data_mc_true_ratio->Draw("PE X0 same");

        // Curves grouped by family (points = closure, lines = size of correction), then the two
        // band swatches, so the legend reads top-down as the story: correction, closure, and the
        // two gaps that connect them.
        TLegend *leg_ratio = new TLegend(0.16, 0.66, 0.96, 0.985);
        leg_ratio->SetFillStyle(0);
        leg_ratio->SetBorderSize(0);
        leg_ratio->SetTextFont(font_code);
        leg_ratio->SetTextSize(legend_size - 2.);
        leg_ratio->SetNColumns(2);
        leg_ratio->AddEntry(h_data_mc_reco_ratio, "reco data/reco MC", "pe1");
        leg_ratio->AddEntry(h_mc_gen_reco_ratio, "gen MC/reco MC", "l");
        leg_ratio->AddEntry(h_data_mc_true_ratio, "unfolded data/gen MC", "pe1");
        leg_ratio->AddEntry(h_unfolded_reco_ratio, "unfolded/reco", "l");
        leg_ratio->Draw();
        pad_ratio->RedrawAxis();
        c_unfold->cd();
        c_unfold->Update();
        
        h_data_mc_reco_ratio->GetXaxis();
        c_unfold->Modified();
        c_unfold->Update();
       
        if (scan_niter) {
            // One page per iteration: "(" opens the PDF on the first page, ")" closes it on the last.
            TString iter_page = iter_pdf;
            if (iter == niter_min)      iter_page += "(";
            else if (iter == niter_max) iter_page += ")";
            c_unfold->Print(iter_page, "pdf");

            // Also drop a standalone PNG per iteration into png_dir (zero-padded so they sort 001..100).
            c_unfold->Print(Form("%siteration_%03d.png", png_dir.Data(), iter + 1), "png");
        } else {
            // Single run: the usual one-off bottomline plot next to the other outputs.
            TString plot_stem = folder + "unfolding_plot_" + label + "_bottomline_test_eec_2D";
            c_unfold->Print(plot_stem + ".pdf");
            c_unfold->Print(plot_stem + ".png");
        }

        // Only delete the canvas when another pass will recreate it (avoids the name clash during
        // the scan). On the final/only pass keep it, so the window stays open to view.
        if (iter < niter_max) delete c_unfold;


    }

    // Write the corrected histograms once, from the final iteration's result.
    if (iter == niter_max) {
    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    h_data_fully_corrected->SetName("h_data_unfolded_corrected");
    h_data_fully_corrected->Write();
    h_data_after_fit->SetName("h_data");
    h_data_after_fit->Write();
    h_mc_reco->SetName("h_mc_reco");
    h_mc_reco->Write();

    //Write 2D histograms for plotting
    h_mc_reco_2D->Write();
    // The particle-level gen reference -- hgenjet_2b_all in data mode. Re-enabled so
    // apply_weights_and_systematics.C can draw it against the final result: it is the only
    // curve in this file that lives at gen level, and without it the theory comparison
    // would have to re-open and re-sum every per-block MC file itself. Already restricted
    // to [ibin_dr_min, ibin_dr_max] and unit-normalised by the bottomline block above,
    // exactly like h_data_fully_corrected_2D next to it.
    h_mc_true_2D->Write();
    h_data_purity_corrected_2D->Write();
    h_data_fully_corrected_2D->Write();
    h_data_refolded_2D->Write();

    fout->Close();
    delete fout;
    }

    }   // end loop over iteration counts
    gErrorIgnoreLevel = prev_ignore;

    // ---- Regularisation: pick the number of iterations from the refolding test ------------
    // p-value of refolded vs input against the iteration number. It rises as the nonclosure
    // shrinks, then saturates. The optimal iteration is the smallest niter with
    // p >= pvalue_threshold; going further only feeds fluctuations back through the prior.
    if (unfoldBayes && v_niter.size() > 1) {

        auto firstAbove = [&](const std::vector<double> &p) -> int {
            for (size_t i = 0; i < p.size(); ++i)
                if (p[i] >= pvalue_threshold) return (int) v_niter[i];
            return -1;   // never reaches the threshold within the scanned range
        };
        const int best_all = firstAbove(v_pvalue_all);
        const int best_pt  = firstAbove(v_pvalue_pt);

        // Secondary diagnostic: with a well-conditioned response the p-value can clear the
        // threshold at the first iteration, and chi2/ndf still shows where the refolding starts
        // over-describing the input (chi2/ndf well below 1).
        auto closestToOne = [&](const std::vector<double> &c) -> int {
            int best = -1; double dmin = 1e30;
            for (size_t i = 0; i < c.size(); ++i) {
                if (c[i] <= 0.) continue;
                const double d = std::fabs(c[i] - 1.);
                if (d < dmin) { dmin = d; best = (int) v_niter[i]; }
            }
            return best;
        };
        const int c2one_all = closestToOne(v_chi2ndf_all);

        std::cout << "\n================ Refolding goodness-of-fit scan ================\n"
                  << "  input   : h_data_purity_corrected (reconstructed level)\n"
                  << "  refolded: response->ApplyToTruth(unfolded)\n"
                  << Form("  criterion: smallest niter with p >= %.3f\n", pvalue_threshold);
        std::cout << Form("  %6s  %12s %10s   %12s %10s\n",
                          "niter", "chi2/ndf 2D", "p 2D", "chi2/ndf pT", "p pT");
        for (size_t i = 0; i < v_niter.size(); ++i)
            std::cout << Form("  %6.0f  %12.3f %10.3e   %12.3f %10.3e\n",
                              v_niter[i], v_chi2ndf_all[i], v_pvalue_all[i],
                              v_chi2ndf_pt[i], v_pvalue_pt[i]);
        std::cout << Form("  ----> optimal niter (full 2D reco space) = %d\n", best_all)
                  << Form("  ----> optimal niter (pT bin %d only)      = %d\n", ibin_pt, best_pt)
                  << Form("  ----> cross-check, chi2/ndf closest to 1 = %d\n", c2one_all)
                  << "================================================================\n"
                  << std::endl;

        const int np = (int) v_niter.size();
        TGraph *g_p_all    = new TGraph(np, v_niter.data(), v_pvalue_all.data());
        TGraph *g_p_pt     = new TGraph(np, v_niter.data(), v_pvalue_pt.data());
        TGraph *g_c2_all   = new TGraph(np, v_niter.data(), v_chi2ndf_all.data());
        TGraph *g_c2_pt    = new TGraph(np, v_niter.data(), v_chi2ndf_pt.data());
        g_p_all ->SetLineColor(blue);   g_p_all ->SetMarkerColor(blue);   g_p_all ->SetMarkerStyle(kFullCircle);
        g_p_pt  ->SetLineColor(orange); g_p_pt  ->SetMarkerColor(orange); g_p_pt  ->SetMarkerStyle(kOpenSquare);
        g_c2_all->SetLineColor(blue);   g_c2_all->SetMarkerColor(blue);   g_c2_all->SetMarkerStyle(kFullCircle);
        g_c2_pt ->SetLineColor(orange); g_c2_pt ->SetMarkerColor(orange); g_c2_pt ->SetMarkerStyle(kOpenSquare);
        for (TGraph *g : {g_p_all, g_p_pt, g_c2_all, g_c2_pt}) { g->SetLineWidth(2); g->SetMarkerSize(0.9); }

        const double x_lo_it = v_niter.front() - 0.5;
        const double x_hi_it = v_niter.back()  + 0.5;

        TCanvas *c_reg = new TCanvas("c_regularisation", "", 800, 800);
        TPad *pad_p  = new TPad("pad_reg_p",  "", 0., 0.42, 1., 1.);
        TPad *pad_c2 = new TPad("pad_reg_c2", "", 0., 0.,   1., 0.42);
        for (TPad *p : {pad_p, pad_c2}) { p->SetTicks(1, 1); p->SetFillColor(0); }
        pad_p ->SetMargin(0.13, 0.05, 0.00, 0.08);
        pad_c2->SetMargin(0.13, 0.05, 0.22, 0.00);
        c_reg->cd(); pad_p->Draw(); pad_c2->Draw();

        // ---- top: p-value vs iteration, full 2D reco space only (the space that was unfolded,
        // and the one the selected iteration is taken from). Log y: the first iterations sit
        // many orders of magnitude below the threshold, and the approach to it has to be readable.
        pad_p->cd();
        //pad_p->SetLogy();
        double p_min = 1.;
        for (double p : v_pvalue_all) if (p > 0.) p_min = std::min(p_min, p);
        // One whole decade below the smallest p, so (a) the first iterations stay inside the
        // frame -- with 20 iterations the first p can be ~1e-17 -- and (b) the lowest labelled
        // decade sits one step above the frame edge instead of being clipped by it.
        const double y_lo_p = std::max(std::pow(10., std::floor(std::log10(p_min)) - 1.), 1e-30);
        const double y_hi_p = 1.05;  // p is bounded by 1; just enough headroom to clear the curve
        TH1F *fr_p = pad_p->DrawFrame(x_lo_it, y_lo_p, x_hi_it, y_hi_p);
        fr_p->GetYaxis()->SetTitle("p-value");
        fr_p->GetYaxis()->CenterTitle(true);
        fr_p->GetYaxis()->SetTitleFont(font_code); fr_p->GetYaxis()->SetTitleSize(title_size); fr_p->GetYaxis()->SetTitleOffset(1.5);
        fr_p->GetYaxis()->SetLabelFont(font_code); fr_p->GetYaxis()->SetLabelSize(label_size);
        // Label only every few decades: over ~10 decades one label per power of ten is unreadable.
        // On a log axis the primary-division count is how many decades ROOT puts between labels.
        fr_p->GetYaxis()->SetNdivisions(506);
        // The two pads share an edge, and this pad's LOWEST y label sits exactly on it, so it
        // gets clipped in half (and would collide with the bottom pad's highest label anyway).
        // Hide it: size 0. labNum = 1 is the first (lowest) label.
        fr_p->GetYaxis()->ChangeLabel(1, -1, 0.);
        fr_p->GetXaxis()->SetLabelSize(0); fr_p->GetXaxis()->SetTitleSize(0);  // x belongs to the bottom pad
        g_p_all->Draw("PL SAME");

        TLine *l_thr = new TLine(x_lo_it, pvalue_threshold, x_hi_it, pvalue_threshold);
        l_thr->SetLineColor(kGray + 2); l_thr->SetLineStyle(2); l_thr->Draw("same");
        // Vertical marker at the chosen iteration. TLine takes user coordinates, which on a
        // log pad are the values themselves (not their logs), so the frame range is what to use.
        if (best_all > 0) {
            TLine *l_best = new TLine(best_all, y_lo_p, best_all, y_hi_p);
            l_best->SetLineColor(red); l_best->SetLineStyle(7); l_best->SetLineWidth(2);
            l_best->Draw("same");
        }
        // Label the threshold line where it sits rather than adding a legend. User coordinates
        // so the label tracks the line; it goes above the line on linear y, below on log y.
        TLatex tx_thr;
        tx_thr.SetTextFont(font_code); tx_thr.SetTextSize(label_size); tx_thr.SetTextColor(kGray + 2);
        tx_thr.DrawLatex(x_hi_it - 0.25 * (x_hi_it - x_lo_it),
                         pad_p->GetLogy() ? pvalue_threshold / 8.
                                          : pvalue_threshold + 0.04 * (y_hi_p - y_lo_p),
                         Form("p = %.2f", pvalue_threshold));

        TLatex tx_reg; tx_reg.SetNDC();
        tx_reg.SetTextFont(font_code); tx_reg.SetTextSize(label_size);
        // Right of the rise, in the empty band under the saturated curve: on the left the text
        // sits on top of the part of the curve the plot is about.
        tx_reg.DrawLatex(0.45, 0.30,
            test_mode == 0 ? "MC full-sample closure" :
            test_mode == 1 ? "MC split closure test"  : "Data Run 3");
        tx_reg.DrawLatex(0.45, 0.24, "D'Agostini, early stopping");
        tx_reg.DrawLatex(0.45, 0.18, best_all > 0 ? Form("optimal N iter = %d", best_all)
                                                  : "no iteration reaches the threshold");
        TLatex cms_reg; cms_reg.SetNDC();
        cms_reg.SetTextFont(62); cms_reg.SetTextSize(0.042); cms_reg.DrawLatex(0.13, 0.945, "CMS");
        cms_reg.SetTextFont(52); cms_reg.SetTextSize(0.034); cms_reg.DrawLatex(0.205, 0.945, "Internal");
        cms_reg.SetTextFont(42); cms_reg.SetTextSize(0.034); cms_reg.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");
        pad_p->RedrawAxis();

        // ---- bottom: the same information as chi2/ndf, where 1 is the target
        pad_c2->cd();
        pad_c2->SetLogy();
        double c2_max = 1., c2_min = 1.;
        for (double c : v_chi2ndf_all) if (c > 0.) { c2_max = std::max(c2_max, c); c2_min = std::min(c2_min, c); }
        TH1F *fr_c2 = pad_c2->DrawFrame(x_lo_it, c2_min * 0.5, x_hi_it, c2_max * 2.);
        fr_c2->GetYaxis()->SetTitle("#chi^{2}/ndf");
        fr_c2->GetYaxis()->CenterTitle(true);
        fr_c2->GetYaxis()->SetTitleFont(font_code); fr_c2->GetYaxis()->SetTitleSize(title_size); fr_c2->GetYaxis()->SetTitleOffset(1.5);
        fr_c2->GetYaxis()->SetLabelFont(font_code); fr_c2->GetYaxis()->SetLabelSize(label_size);
        // Same shared edge from the other side: hide this pad's HIGHEST y label (labNum = -1).
        fr_c2->GetYaxis()->ChangeLabel(-1, -1, 0.);
        fr_c2->GetXaxis()->SetTitle("number of iterations");
        fr_c2->GetXaxis()->CenterTitle(true);
        // Offset in precision-43 fonts scales off the (short) pad height: anything near 2 puts
        // the title off the canvas entirely. 1.3 lands it just below the tick labels.
        fr_c2->GetXaxis()->SetTitleFont(font_code); fr_c2->GetXaxis()->SetTitleSize(title_size); fr_c2->GetXaxis()->SetTitleOffset(1.3);
        fr_c2->GetXaxis()->SetLabelFont(font_code); fr_c2->GetXaxis()->SetLabelSize(label_size);
        g_c2_all->Draw("PL SAME");
        TLine *l_one = new TLine(x_lo_it, 1., x_hi_it, 1.);
        l_one->SetLineColor(kGray + 2); l_one->SetLineStyle(2); l_one->Draw("same");
        pad_c2->RedrawAxis();

        // Next to the per-iteration scan outputs, which this plot summarises.
        TString reg_stem = folder + "refolding_pvalue_vs_iteration_" + label;
        c_reg->Print(reg_stem + ".pdf");
        c_reg->Print(reg_stem + ".png");

        // Keep the scan itself, so the choice can be re-plotted without re-unfolding.
        TFile *f_reg = new TFile(reg_stem + ".root", "recreate");
        g_p_all ->Write("g_pvalue_vs_niter_2D");
        g_p_pt  ->Write(Form("g_pvalue_vs_niter_pt%d", ibin_pt));
        g_c2_all->Write("g_chi2ndf_vs_niter_2D");
        g_c2_pt ->Write(Form("g_chi2ndf_vs_niter_pt%d", ibin_pt));
        f_reg->Close();
        delete f_reg;
        std::cout << "Wrote regularisation scan to " << reg_stem << ".{pdf,png,root}" << std::endl;
    }

    // gApplication -> Terminate(0);
}

// SAMPLE:      qcd | bjet | both   ("both" sums the bjet response, corrections and truth
//              onto the qcd ones -- counts and RooUnfoldResponse are added, the
//              purity/efficiency ratios are recomputed from the summed counts)
// There are two independent generator choices, one per half of the chain:
//
// UNFOLDING_GENERATOR: pythia | herwig
//              Everything the unfolding uses: the migration matrix, the purity and
//              efficiency corrections, the truth reference, h_mc_reco, the hgenjet_2b_*
//              corrections, and in modes 0/1 the MC being unfolded. Read from the
//              per-block output of run_agg_ntuple_chunks.sh, so no pre-merged file.
// TF_GENERATOR: pythia | herwig
//              The template fit that supplies the signal fraction (data mode only).
//              Independent of UNFOLDING_GENERATOR: only a Pythia fit exists today, so
//              a Herwig unfolding still takes its signal fraction from the Pythia fit.
//              Setting it to herwig stops with a message until that fit is produced.
// test_mode:   0 = full-MC closure, 1 = split test, 2 = data
// unfoldBayes: true = Bayesian, false = matrix inversion
// scan_niter:  true = scan niter = 1..100 and pick the optimal one from the refolding
//              goodness-of-fit test; false = one unfolding at the nominal niter = 4.
//              Bayesian only -- matrix inversion has no iterations to scan.
// TRACK_EFF_UNC: false (nominal) | true
//              Which MC production to unfold with. true reads the files written with 3%
//              of the reconstructed tracks thrown away during the B reconstruction
//              (TRACK_EFF_UNC=true in run_agg_ntuple_chunks.sh) -- the tracking-efficiency
//              systematic. It varies BOTH halves of the MC: the response/corrections AND
//              the template fit that supplies the signal fraction, so the propagation is
//              complete -- which means template_fit.cpp(SAMPLE, TF_GENERATOR, true) has to
//              have been run first. Only the MC changes; the data being unfolded and the
//              data being fitted are the same files either way.
//              Results go to their own folder, so the nominal is never overwritten.
// TF_VARIATION: nominal | var0B_2 | var0B_0
//              Which of template_fit.cpp's variations supplies the signal fraction (data
//              mode only). var0B_2 / var0B_0 are the LIGHT-JET MISTAG variations: the 0B
//              template -- the jets with no gen b hadron that passed the b tag, i.e.
//              mistagged light and charm -- has its contribution to the background PDF
//              doubled or removed, and the fit is redone. Both fits already exist in every
//              template-fit directory (template_fit.cpp writes one file per variation), so
//              this costs one more unfolding and no refit. Results go to their own folder.
//              Unlike TRACK_EFF_UNC this varies ONLY the signal fraction: the response and
//              the corrections are the nominal ones, which is right -- the mistag rate
//              affects what fraction of the data is signal, not how signal migrates.
//
// e.g.  root -l -b -q 'apply_unfolding_2d.C("both","pythia")'   // nominal
//       root -l -b -q 'apply_unfolding_2d.C("both","herwig")'   // Herwig unfolding MC
//       root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",true)'
//                                                             // tracking-efficiency variation
//       root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,true,false,"pythia",false,"var0B_2")'
//                                                             // light-jet mistag up
// SFUPART_VARIATION: nominal | jpcalib_hf | qqrate_up | qqrate_down
//              Which UParT b-tag efficiency scale factor to divide into the RECO-level data
//              (data mode only). The SF is measured in reco R_BB bins, so it is applied
//              before unfolding and the response matrix carries it through the migration.
//              nominal is the central SF -- it is a CORRECTION, applied in every run, not an
//              opt-in. The other three are its calibration systematics: jpcalib_hf swaps the
//              JP heavy-flavour calibration, qqrate_up/down move the qq rate by +/-25%.
//              Each lands in its own results folder.
void apply_unfolding_2d(TString SAMPLE = "both", TString UNFOLDING_GENERATOR = "pythia",
                        int test_mode = 2, bool unfoldBayes = true, bool scan_niter = false,
                        TString TF_GENERATOR = "herwig", bool TRACK_EFF_UNC = false,
                        TString TF_VARIATION = "nominal",
                        TString SFUPART_VARIATION = "nominal"){

    // Kept so the old positional name still reads naturally below.
    TString GENERATOR = UNFOLDING_GENERATOR;

    // A typo must not fall through to the default sample: check before reading anything.
    if (SAMPLE != "qcd" && SAMPLE != "bjet" && SAMPLE != "both") {
        std::cerr << "ERROR: unknown SAMPLE '" << SAMPLE << "' (use qcd | bjet | both)" << std::endl;
        return;
    }
    if (GENERATOR != "pythia" && GENERATOR != "herwig") {
        std::cerr << "ERROR: unknown UNFOLDING_GENERATOR '" << GENERATOR
                  << "' (use pythia | herwig)" << std::endl;
        return;
    }
    if (TF_GENERATOR != "pythia" && TF_GENERATOR != "herwig") {
        std::cerr << "ERROR: unknown TF_GENERATOR '" << TF_GENERATOR
                  << "' (use pythia | herwig)" << std::endl;
        return;
    }
    if (TF_VARIATION != "nominal" && TF_VARIATION != "var0B_2" && TF_VARIATION != "var0B_0") {
        std::cerr << "ERROR: unknown TF_VARIATION '" << TF_VARIATION
                  << "' (use nominal | var0B_2 | var0B_0)" << std::endl;
        return;
    }
    // A varied signal fraction only ever multiplies real data, so in the closure modes it
    // would change nothing while still writing to a separate folder -- which reads as a
    // result and is not one.
    if (TF_VARIATION != "nominal" && test_mode != 2) {
        std::cerr << "ERROR: TF_VARIATION '" << TF_VARIATION << "' needs test_mode 2 (data): "
                  << "the signal fraction is only applied to data" << std::endl;
        return;
    }
    if (sfupartHist(SFUPART_VARIATION).Length() == 0) {
        std::cerr << "ERROR: unknown SFUPART_VARIATION '" << SFUPART_VARIATION
                  << "' (use nominal | jpcalib_hf | qqrate_up | qqrate_down)" << std::endl;
        return;
    }
    // Same reasoning as TF_VARIATION: the SF is only ever applied to data.
    if (SFUPART_VARIATION != "nominal" && SFUPART_VARIATION != "off" && test_mode != 2) {
        std::cerr << "ERROR: SFUPART_VARIATION '" << SFUPART_VARIATION << "' needs test_mode 2 (data): "
                  << "the UParT SF is only applied to data" << std::endl;
        return;
    }
    // Refuse before creating anything, as above: a missing SF file would otherwise surface
    // only in the log, after the results folder exists.
    if (test_mode == 2 && SFUPART_VARIATION != "off") {
        if (gSystem->AccessPathName(sfupartFile())) {
            std::cerr << "ERROR: missing UParT SF file " << sfupartFile() << std::endl;
            return;
        }
    }

    // Refuse before creating anything: past this point stdout goes to the log file, so a
    // failure inside apply_unfolding() would leave an empty results folder and a silent
    // terminal.
    if (test_mode == 2) {
        const TString tf_file = templateFitFile(SAMPLE, TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION);
        if (tf_file.Length() == 0 || gSystem->AccessPathName(tf_file)) {
            std::cerr << "ERROR: no template fit for TF_GENERATOR '" << TF_GENERATOR
                      << "'" << (TRACK_EFF_UNC ? " with TRACK_EFF_UNC" : "")
                      << ", variation '" << TF_VARIATION << "'"
                      << ": " << tf_file << "\n       run template_fit.cpp(\""
                      << ((SAMPLE == "both") ? "both" : "qcd") << "\",\"" << TF_GENERATOR
                      << "\"" << (TRACK_EFF_UNC ? ",true" : "") << ") first" << std::endl;
            return;
        }
    }

    // Tag only what actually differs from the nominal, so nominal names are unchanged.
    // Both built by result_paths.h, the same functions apply_weights_and_systematics.C uses
    // to FIND them. "" tags mean nominal, so nominal paths are exactly what they always
    // were; a variation lands in its own directory next to the nominal, never over it.
    TString dataset = SAMPLE + "_" + GENERATOR
                    + variationTag(TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION, SFUPART_VARIATION);

    TString folder = resultFolder(SAMPLE, GENERATOR, TF_GENERATOR, TRACK_EFF_UNC,
                                  TF_VARIATION, SFUPART_VARIATION);

    // An unrecognised variation name poisons the tag rather than silently becoming the
    // nominal (see result_paths.h). Catch it here, before a directory is created.
    if (folder.Contains(kUnknownTag())) {
        std::cerr << "ERROR: refusing to run -- one of the variation names is not recognised."
                  << std::endl;
        return;
    }

    gSystem->mkdir(folder, kTRUE);
    TString pT_selection = "80_inf";
    bool btag = true;
    Int_t n = 1;

    // ---- Log file ----
    TDatime now;
    TString timestamp = Form(
        "%04d-%02d-%02d_%02d-%02d-%02d",
        now.GetYear(),
        now.GetMonth(),
        now.GetDay(),
        now.GetHour(),
        now.GetMinute(),
        now.GetSecond()
    );
    TString label;
    label += unfoldBayes ? "_bayesian" : "_MI";
    if(test_mode==0)  label += "_full_closure";
    if(test_mode==1)  label += "_split_test";
    if(test_mode==2)  label += "_data";
    if(scan_niter)    label += "_iterative";
    TString logfile = folder +  "unfolding_" + dataset + label + "_" +  timestamp + ".log";
    gSystem->RedirectOutput(logfile);  // "a" = append
	apply_unfolding(dataset, folder, btag, n, pT_selection, test_mode, unfoldBayes, scan_niter,
	                SAMPLE, GENERATOR, TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION, SFUPART_VARIATION);

    // Restore terminal output
    gSystem->RedirectOutput(nullptr);

    // Close the loop between writer and reader: resultFile() is the path
    // apply_weights_and_systematics.C will look for. Check here, at write time, that it
    // actually resolves -- otherwise a mismatch only shows up much later as a variation
    // that is "missing" or, worse, as one that silently resolved to something else.
    const TString expected = resultFile(SAMPLE, GENERATOR, test_mode, unfoldBayes,
                                        TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION,
                                        SFUPART_VARIATION);
    if (gSystem->AccessPathName(expected)) {
        std::cerr << "WARNING: the run finished but its result is not where the systematics "
                  << "macro will look for it:\n         " << expected
                  << "\n         Check " << logfile << " -- the unfolding may have returned "
                  << "early, or the naming in result_paths.h has drifted from what "
                  << "apply_unfolding() writes." << std::endl;
    } else {
        std::cout << "Result: " << expected << std::endl;
    }

}
