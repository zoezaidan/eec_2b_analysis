// Post-processing of the unfolded result.
//
// Takes the output of apply_unfolding_2d.C, applies the final weights, builds a
// systematic band from the variations (Pythia vs Herwig to start with) and writes the
// result carrying stat and syst errors separately plus their quadrature sum.
//
// Prerequisite: run apply_unfolding_2d.C once per (SAMPLE, GENERATOR) you want in here.
// A variation whose file is missing is reported and skipped, so the nominal alone still
// works -- you just get stat errors until the Herwig run exists.

#include "result_paths.h"   // the ONE definition of every result path and variation tag
#include <vector>
#include <cmath>
#include <map>

// ---- Paths ---------------------------------------------------------------------------
// resultFolder/resultLabel/resultFile and every variation tag now come from
// result_paths.h, which apply_unfolding_2d.C includes too. They used to be duplicated here
// and drifted -- see the note at the top of that header for what that cost.

// A colour blended towards white, i.e. what "alpha over a white canvas" looks like.
// SetFillColorAlpha() is honoured in the PDF but ignored by ROOT's PNG backend, which
// paints the fill solid; blending gives the same light tone in both.
Color_t blendWithWhite(const char *hex, double alpha)
{
    TColor *c = gROOT->GetColor(TColor::GetColor(hex));
    if (!c) return kWhite;
    Float_t r, g, b;
    c->GetRGB(r, g, b);
    return (Color_t) TColor::GetColor(Float_t(1. - alpha + alpha * r),
                                      Float_t(1. - alpha + alpha * g),
                                      Float_t(1. - alpha + alpha * b));
}

// Generator name as it should appear on a plot.
TString prettyGen(const TString &g)
{
    if (g == "herwig") return "Herwig";
    if (g == "pythia") return "Pythia";
    return g;
}

// A correction-level variation: no second unfolding, the alternative generator's correction
// is applied to the final result and the result re-normalised. inverse = the correction is
// DIVIDED into the data (an efficiency); false = it MULTIPLIES the data (the EEC weight).
struct CorrVariation {
    TString name;     // printout column and histogram name (no spaces)
    TString label;    // legend text
    TString effType;  // plotNice_UParT_efficiency effType, for the "run this first" message
    TString file;     // that macro's output .root
    bool    inverse;
    TString sysLabel; // name in the "systematics applied" tag on the data-vs-gen plot
};

// One systematic variation: the same unfolding redone with something changed.
struct Variation {
    TString name;        // printout column and output histogram name (no spaces)
    TString sample;
    TString generator;   // UNFOLDING_GENERATOR of that run
    TString tfGenerator; // TF_GENERATOR of that run ("pythia" unless the fit is varied)
    TString label;       // legend text
    bool    trackEffUnc; // TRACK_EFF_UNC of that run (false unless the MC itself is varied)
    TString tfVariation = "nominal"; // TF_VARIATION of that run (the template-fit variation)
    TString sfupartVariation = "nominal"; // SFUPART_VARIATION of that run (the UParT SF variation)
    // Entries sharing a non-empty group are the two SIDES OF ONE uncertainty (an up/down
    // pair). Their envelope -- the larger of the two |shifts| in each bin -- enters the
    // quadrature total once. An empty group means an independent source that enters on its
    // own, which is how every entry behaved before groups existed.
    TString group = "";
    // Name in the "systematics applied" tag on the data-vs-gen plot. The two sides of a
    // pair carry the SAME sysLabel and are listed once. Leave it empty to keep an entry
    // out of the tag while it still counts in the band.
    TString sysLabel = "";
    // false = DISPLAY ONLY. The curve is loaded, drawn and printed, but contributes nothing
    // to sigma+/sigma- and nothing to the tag. For showing the effect of a correction --
    // "here is the result without the UParT SF" -- which is a comparison, not a systematic.
    bool    inBand = true;
};

// ---- Weights = scale factors ---------------------------------------------------------
// Applied to the nominal AND to every variation, before any difference is taken: a shape
// systematic has to be compared after the same weighting, otherwise a pure normalisation
// difference between the samples leaks into the shape band.
//
// This is the GENERIC per-bin scale-factor hook, deliberately not named for any one SF --
// the UParT SF has its own sfupart* plumbing in apply_unfolding_2d.C (it is applied at reco
// level, before unfolding), and this stays free for a future SF that genuinely belongs on
// the final result.
//
// sf_file/sf_hist is the per-bin scale factor histogram, multiplied bin by bin (it must
// have the same binning as the result). global_scale is a single number for a flat SF.
struct Weights {
    TString sf_file             = "";    // per-bin scale factors, "" = off
    TString sf_hist             = "";
    bool    sf_errors           = true;  // false = treat the SF as exact (ignore its errors)
    // An EFFICIENCY scale factor is DIVIDED into the result, it does not multiply it. The
    // chain corrected the data with the efficiency measured in MC, and the true efficiency
    // is eps_data = SF * eps_MC, so
    //     corrected = raw / eps_MC   ->   true = raw / (SF * eps_MC) = corrected / SF.
    // Same convention as CorrVariation::inverse elsewhere in this file. Set false for a
    // factor that genuinely multiplies (a weight, a luminosity ratio).
    bool    sf_divide           = true;
    double  global_scale        = 1.0;   // flat scale factor / luminosity / prescale; 1 = off
    bool    normalise_unit_area = true;  // Scale(1/integral, "width") -- per bin width, unit area
};

bool applyWeights(TH1D *h, const Weights &w)
{
    if (!h) return false;

    if (w.sf_file.Length() > 0) {
        TFile *fw = TFile::Open(w.sf_file);
        if (!fw || fw->IsZombie()) {
            std::cerr << "ERROR: cannot open scale-factor file " << w.sf_file << std::endl;
            return false;
        }
        TH1D *hsf = dynamic_cast<TH1D *>(fw->Get(w.sf_hist));
        if (!hsf) {
            std::cerr << "ERROR: '" << w.sf_hist << "' not in " << w.sf_file << std::endl;
            fw->Close();
            return false;
        }
        if (hsf->GetNbinsX() != h->GetNbinsX()) {
            std::cerr << "ERROR: scale factors have " << hsf->GetNbinsX() << " bins, result has "
                      << h->GetNbinsX() << std::endl;
            fw->Close();
            return false;
        }
        // TH1::Multiply/Divide propagate the SF errors into the result's bin errors, i.e.
        // into its STATISTICAL error. Zero them to treat the SF as exact. (The UParT SF no
        // longer comes through here -- it is applied at reco level in apply_unfolding_2d.C,
        // which propagates its statistical error the same way.)
        TH1D *hsf_use = (TH1D *) hsf->Clone("sf_use");
        hsf_use->SetDirectory(nullptr);
        if (!w.sf_errors)
            for (int i = 1; i <= hsf_use->GetNbinsX(); ++i) hsf_use->SetBinError(i, 0.);
        for (int i = 1; i <= hsf_use->GetNbinsX(); ++i) {
            if (hsf_use->GetBinContent(i) != 0.) continue;
            std::cerr << "ERROR: scale factor is zero in bin " << i << " of " << w.sf_hist
                      << " -- refusing to " << (w.sf_divide ? "divide" : "multiply") << std::endl;
            delete hsf_use;
            fw->Close();
            return false;
        }
        if (w.sf_divide) h->Divide(hsf_use);
        else             h->Multiply(hsf_use);
        delete hsf_use;
        fw->Close();
    }

    if (w.global_scale != 1.0) h->Scale(w.global_scale);

    if (w.normalise_unit_area) {
        const double integral = h->Integral();
        if (integral <= 0.) {
            std::cerr << "WARNING: " << h->GetName() << " has integral " << integral
                      << ", not normalising" << std::endl;
            return false;
        }
        h->Scale(1. / integral, "width");
    }
    return true;
}

// Read one result histogram, detached from its file, with the weights already applied.
TH1D *loadWeighted(const TString &filename, const TString &histname,
                   const TString &newname, const Weights &w)
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "MISSING: " << filename << std::endl;
        return nullptr;
    }
    TH1D *h_in = dynamic_cast<TH1D *>(f->Get(histname));
    if (!h_in) {
        std::cerr << "ERROR: '" << histname << "' not found in " << filename << std::endl;
        f->Close();
        return nullptr;
    }
    TH1D *h = (TH1D *) h_in->Clone(newname);
    h->SetDirectory(nullptr);
    f->Close();
    if (!applyWeights(h, w)) return nullptr;
    return h;
}

// =====================================================================================
// GENERATOR / TF_GENERATOR are the UNFOLDING_GENERATOR and TF_GENERATOR of the NOMINAL
// unfolding run; the variations below say which of them to flip.
void apply_weights_and_systematics(TString SAMPLE = "both", TString GENERATOR = "pythia",
                                   int test_mode = 2, bool unfoldBayes = true,
                                   TString histname = "h_data_fully_corrected_2D",
                                   TString TF_GENERATOR = "pythia",
                                   int ibin_pt = 2)
{
    if (SAMPLE != "qcd" && SAMPLE != "bjet" && SAMPLE != "both") {
        std::cerr << "ERROR: unknown SAMPLE '" << SAMPLE << "' (use qcd | bjet | both)" << std::endl;
        return;
    }
    if (GENERATOR != "pythia" && GENERATOR != "herwig") {
        std::cerr << "ERROR: unknown GENERATOR '" << GENERATOR << "' (use pythia | herwig)" << std::endl;
        return;
    }
    // A generator systematic only means something when both unfoldings ran on the SAME
    // input, which is only true for data. In the closure modes each run unfolds its own MC,
    // so the difference is not a systematic on a measurement.
    if (test_mode != 2)
        std::cerr << "WARNING: test_mode " << test_mode << " unfolds MC, not data -- the "
                  << "variation differences below are NOT a systematic on a measurement."
                  << std::endl;

    // =========================== SET SCALE FACTORS HERE ===============================
    // The UParT b-tag efficiency SF is NOT applied here any more. It is measured in reco
    // R_BB bins, so it is divided into the RECO-level data inside apply_unfolding_2d.C,
    // before unfolding -- see SFUPART_VARIATION there. Applying it here would put a reco-binned
    // correction on particle-level bins and ignore the migration between them.
    // sf_file/sf_hist stay available for a genuinely final-state, per-bin factor.
    Weights w;
    w.sf_file             = "";    // per-bin scale factors: file and histogram name
    w.sf_hist             = "";
    w.sf_errors           = true;  // propagate the SF errors into the stat error
    w.sf_divide           = true;  // an efficiency SF divides into the data (see Weights)
    w.global_scale        = 1.0;   // flat scale factor; 1 = off
    w.normalise_unit_area = true;  // the EEC is a shape: unit area, divided by bin width
    // ==================================================================================

    // ============================= SET VARIATIONS HERE ================================
    // Each entry is re-read, re-weighted and compared bin by bin against the nominal.
    //
    // Nominal is (UNFOLDING_GENERATOR, TF_GENERATOR) = (pythia, pythia). Each variation
    // flips ONE of the two and is subtracted from it:
    //   template_fit    (pythia, herwig) -- the signal fraction from the Herwig fit
    //   unfolding_model (herwig, pythia) -- the whole unfolding MC: migration matrix AND
    //                                       purity, efficiency, truth, hgenjet_2b_*
    // They are different sources, so the symmetrised |nominal - variation| shifts are added
    // in quadrature below. (If you ever want them treated as two sides of ONE uncertainty,
    // give them a shared Variation::group -- the total loop then takes their envelope once,
    // which is how the light-jet mistag pair below is combined.)
    const TString other_gen = (GENERATOR == "pythia") ? "herwig" : "pythia";
    std::vector<Variation> variations = {
        { "template_fit",    SAMPLE, GENERATOR, other_gen,    "TF " + prettyGen(other_gen),        false,
          "nominal", "nominal", "", "MC template modeling" },
        { "unfolding_model", SAMPLE, other_gen, TF_GENERATOR, "Unfolding " + prettyGen(other_gen), false,
          "nominal", "nominal", "", "Detector response" },
        // Tracking efficiency: the same unfolding redone on MC in which 3% of the
        // reconstructed tracks were thrown away during the substructure extraction. That
        // value covers the residual data/MC difference for tracking in the jet core seen
        // in the jet energy scale determination. It is an independent source from the two
        // above -- a detector effect, not a generator choice -- so its symmetrised
        // |nominal - variation| enters the quadrature sum on its own.
        // Produce it with TRACK_EFF_UNC=true in run_agg_ntuple_chunks.sh, then
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, true)
        { "tracking_eff",    SAMPLE, GENERATOR, TF_GENERATOR, "Tracking eff. (-3% tracks)", true,
          "nominal", "nominal", "", "Tracking efficiency" },
        // Light-jet mistagging: the two sides of ONE uncertainty, so they share a group and
        // their envelope counts once in the total (see Variation::group).
        //
        // The 0B template is the jets that passed the b tag with NO gen b hadron in them --
        // mistagged light and charm. It is not a free component of the fit: it is folded
        // into the background PDF at the MC-predicted ratio c' = 0B/(0B+1B), so the mistag
        // rate reaches the measurement only through the SHAPE of the effective background,
        // and the fitted signal fraction is what moves. The two variations rebuild that
        // shape with the 0B admixture doubled and removed -- a conservative +/-100% on the
        // mistag rate, which brackets the measured light-flavour mistag scale factors
        // (20-50% at a tight WP) by a wide margin and needs no external SF input.
        //
        // These need no new MC and no refit: template_fit.cpp already writes one fit per
        // variation into every template-fit directory. Produce each with
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, false, "var0B_2")
        // and the same with "var0B_0".
        { "mistag_0B_up",   SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 2)", false,
          "var0B_2", "nominal", "mistag_0B", "Light jet mistagging" },
        { "mistag_0B_down", SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 0)", false,
          "var0B_0", "nominal", "mistag_0B", "Light jet mistagging" },
        // DISPLAY ONLY (inBand = false): the same unfolding with no UParT SF applied, so
        // the plot shows the result before and after the correction. It is a correction,
        // not a systematic, so it must not enter the band -- hence the trailing false.
        // Produce it with
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, false,
        //                        "nominal", "off")
        { "sfupart_off", SAMPLE, GENERATOR, TF_GENERATOR, "Before UParT SF", false,
          "nominal", "off", "", "", false },

/* ---- disabled (kept for reference): UParT SF calibration systematics ----
        // Step-by-step: the SF CORRECTION is applied (at reco level, inside the unfolding),
        // but its two calibration systematics are held out of the band for now. Uncomment to
        // put them back -- the three unfolding runs they need already exist
        // (_sfjphf, _sfqqup, _sfqqdn) and nothing else has to change.
        //
        //   JP HF swap : one alternative calibration -> one-sided, ungrouped, symmetrised.
        //   qq rate    : a genuine +/-25% pair straddling the central -> grouped, so the
        //                total takes its envelope once and it stays asymmetric.
        // Produce each with
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, false,
        //                        "nominal", "<jpcalib_hf|qqrate_up|qqrate_down>")
        { "sfupart_jphf",   SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (JP HF swap)", false,
          "nominal", "jpcalib_hf",  "", "UParT SF: JP HF swap" },
        { "sfupart_qq_up",  SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (qq rate +25%)", false,
          "nominal", "qqrate_up",   "sfupart_qqrate", "UParT SF: qq rate #pm25%" },
        { "sfupart_qq_down",SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (qq rate #minus25%)", false,
          "nominal", "qqrate_down", "sfupart_qqrate", "UParT SF: qq rate #pm25%" },
---- end disabled: UParT SF systematics ---- */
        // { "sample", (SAMPLE == "both") ? "qcd" : "both", GENERATOR, TF_GENERATOR, "qcd only", false }, // MC composition
    };

    // ---- Correction-level variations (no second unfolding; see the loop below) --------
    // NOTE ON DOUBLE COUNTING: "unfolding_model" above already swaps the whole unfolding MC
    // to the other generator, which includes these very corrections. Keeping both in the
    // quadrature total counts the correction difference twice. Use these when you want the
    // corrections booked on their own -- and then drop "unfolding_model", or read it as the
    // deliberately conservative option.
    const TString corr_dir = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/";
    const TString corr_tag = Form("_%s_pt%d.root", SAMPLE.Data(), ibin_pt);
    std::vector<CorrVariation> corrVariations = {
        /* ---- disabled (kept for reference): 2SV+UParT efficiency variation ----
        // Dropped on request 2026-09-09: not used in the result. Commenting the entry
        // out removes it from the quadrature total, the printout, both plots and the
        // written histograms in one go -- everything downstream iterates this vector.
        // Re-enable by uncommenting; nothing else has to change.
        { "eff_2sv_btag", "2SV+UParT eff. " + prettyGen(other_gen), "2sv_btag",
          corr_dir + "UParT_efficiency_2sv_btag" + corr_tag, true  },
        ---- end disabled: eff_2sv_btag ---- */
        { "eec_weight",   "EEC weight "     + prettyGen(other_gen), "eec_weight",
          corr_dir + "UParT_efficiency_eec_weight" + corr_tag,  false, "MC modeling of EEC weight" },
    };
    // ==================================================================================

    // ---- Validate every booked variation name BEFORE touching a file ------------------
    // An unrecognised name poisons its tag (result_paths.h) instead of aliasing onto the
    // nominal, so it would show up as a missing file. Checking up front turns that into one
    // clear message naming the offending entry, rather than a variation quietly dropping
    // out of the band.
    {
        bool ok = true;
        for (const Variation &v : variations) {
            if (!isKnownSample(v.sample)) {
                std::cerr << "ERROR: variation '" << v.name << "' has unknown sample '"
                          << v.sample << "'" << std::endl; ok = false;
            }
            if (!isKnownGenerator(v.generator) || !isKnownGenerator(v.tfGenerator)) {
                std::cerr << "ERROR: variation '" << v.name << "' has an unknown generator"
                          << std::endl; ok = false;
            }
            if (!isKnownTfVariation(v.tfVariation)) {
                std::cerr << "ERROR: variation '" << v.name << "' has unknown tfVariation '"
                          << v.tfVariation << "'" << std::endl; ok = false;
            }
            if (!isKnownSfupartVariation(v.sfupartVariation)) {
                std::cerr << "ERROR: variation '" << v.name << "' has unknown sfupartVariation '"
                          << v.sfupartVariation << "'" << std::endl; ok = false;
            }
        }
        // Two entries writing to the same name would silently overwrite each other's
        // histograms in the output file.
        for (size_t a = 0; a < variations.size(); ++a)
            for (size_t b = a + 1; b < variations.size(); ++b)
                if (variations[a].name == variations[b].name) {
                    std::cerr << "ERROR: two variations are both named '"
                              << variations[a].name << "'" << std::endl; ok = false;
                }
        if (!ok) {
            std::cerr << "Refusing to run: fix the variation list above." << std::endl;
            return;
        }
    }

    // The nominal always has trackEffUnc = false: every variation above is measured
    // against the unvaried MC.
    const TString label       = resultLabel(SAMPLE, GENERATOR, test_mode, unfoldBayes, TF_GENERATOR, false);
    const TString nominalFile = resultFile(SAMPLE, GENERATOR, test_mode, unfoldBayes, TF_GENERATOR, false);

    std::cout << "Nominal : " << nominalFile << std::endl;
    std::cout << "Histogram: " << histname << std::endl;

    TH1D *h_nom = loadWeighted(nominalFile, histname, "h_nominal", w);
    if (!h_nom) {
        std::cerr << "ERROR: no nominal result -- run apply_unfolding_2d.C(\"" << SAMPLE
                  << "\",\"" << GENERATOR << "\") first" << std::endl;
        return;
    }

    const int nbins = h_nom->GetNbinsX();

    // Per-variation shifts, kept so each one can be plotted and debugged: h_shifts holds
    // the symmetrised |delta| that enters the band, h_deltas the signed one, h_vars the
    // variation curve itself (for the comparison plot).
    std::vector<TH1D *> h_shifts;
    std::vector<TH1D *> h_deltas;
    std::vector<TH1D *> h_vars;
    std::vector<TString> shift_names;
    std::vector<TString> shift_labels;
    // Parallel to h_shifts: the group id each entry belongs to ("" = independent source).
    std::vector<TString> shift_groups;
    // Parallel to h_shifts: false = drawn and printed but kept out of the band.
    std::vector<bool> shift_in_band;
    // Particle-level MC curves, filled by the data-vs-gen plot below and written out with
    // the rest so the comparison can be redrawn from the final file alone.
    std::vector<TH1D *> h_gen_curves;
    std::vector<TString> gen_curve_names;
    // Display names for the "systematics applied" tag on the data-vs-gen plot, in booking
    // order and deduplicated -- an up/down pair shares one sysLabel and is listed once.
    // Built from the entries that actually made it into the band, so the tag cannot claim a
    // source whose result file was missing and got skipped.
    std::vector<TString> syst_tag_labels;
    auto addTagLabel = [&](const TString &lbl) {
        if (!lbl.Length()) return;
        for (const TString &t : syst_tag_labels) if (t == lbl) return;
        syst_tag_labels.push_back(lbl);
    };

    for (const Variation &v : variations) {
        const TString varFile =
            resultFile(v.sample, v.generator, test_mode, unfoldBayes, v.tfGenerator,
                       v.trackEffUnc, v.tfVariation, v.sfupartVariation);
        std::cout << "Variation '" << v.name << "' : " << varFile << std::endl;

        TH1D *h_var = loadWeighted(varFile, histname, "h_var_" + v.name, w);
        if (!h_var) {
            std::cerr << "   -> skipped (run apply_unfolding_2d.C(\"" << v.sample << "\",\""
                      << v.generator << "\"," << test_mode << "," << (unfoldBayes ? "true" : "false")
                      << ",false,\"" << v.tfGenerator << "\","
                      << (v.trackEffUnc ? "true" : "false")
                      << ",\"" << v.tfVariation << "\",\"" << v.sfupartVariation << "\""
                      << ") to include it)" << std::endl;
            continue;
        }
        if (h_var->GetNbinsX() != nbins) {
            std::cerr << "   -> skipped: " << h_var->GetNbinsX() << " bins vs " << nbins
                      << " in the nominal" << std::endl;
            continue;
        }

        // Symmetrised one-sided variation: the band is +/- |nominal - variation|. The signed
        // shift is kept too, so the printout shows which way each variation pulls.
        TH1D *h_shift  = (TH1D *) h_nom->Clone("h_syst_" + v.name);
        TH1D *h_signed = (TH1D *) h_nom->Clone("h_delta_" + v.name);
        h_shift->SetDirectory(nullptr);  h_shift->Reset();
        h_signed->SetDirectory(nullptr); h_signed->Reset();
        for (int i = 1; i <= nbins; ++i) {
            const double d = h_nom->GetBinContent(i) - h_var->GetBinContent(i);
            h_signed->SetBinContent(i, d);
            h_shift->SetBinContent(i, std::fabs(d));
        }

        h_shifts.push_back(h_shift);
        h_deltas.push_back(h_signed);
        h_vars.push_back(h_var);
        shift_names.push_back(v.name);
        shift_labels.push_back(v.label.Length() ? v.label : v.name);
        shift_groups.push_back(v.group);
        shift_in_band.push_back(v.inBand);
        if (v.inBand) addTagLabel(v.sysLabel);
        else std::cout << "   -> display only, not in the band" << std::endl;
    }

    // ---- Correction-level variations, evaluated on the final result -------------------
    // These do not need another unfolding run: the correction enters the result as a
    // per-dr factor, so swapping the generator's correction for the other one is a bin-by-bin
    // rescale of the nominal, followed by RE-NORMALISING to unit area.
    //
    // The re-normalisation is the point. The EEC is a shape, so a correction ratio that is
    // flat in dr cancels completely and only its SHAPE contributes. Comparing the raw
    // correction curves (as plotNice_UParT_efficiency draws them) overstates the effect --
    // e.g. the 2SV+UParT efficiency differs by 3-4%, but almost all of that is a constant
    // offset that self-normalisation removes.
    //
    //   efficiency (divided into the data) : varied = nominal * eff_nominalGen / eff_otherGen
    //   EEC weight (multiplied)            : varied = nominal * w_otherGen  / w_nominalGen
    //
    // Inputs are the histograms plotNice_UParT_efficiency() writes; run it first for each
    // correction, at the SAME pT bin as the unfolding.
    for (const CorrVariation &cv : corrVariations) {
        TFile *fc = TFile::Open(cv.file);
        if (!fc || fc->IsZombie()) {
            std::cerr << "Correction variation '" << cv.name << "': MISSING " << cv.file
                      << "\n   -> skipped (run plotNice_UParT_efficiency(\"" << SAMPLE << "\","
                      << ibin_pt << ",\"" << cv.effType << "\") to include it)" << std::endl;
            continue;
        }
        TH1D *h_gen_nom = dynamic_cast<TH1D *>(fc->Get(GENERATOR == "pythia" ? "h_eff_pythia"
                                                                             : "h_eff_herwig"));
        TH1D *h_gen_alt = dynamic_cast<TH1D *>(fc->Get(GENERATOR == "pythia" ? "h_eff_herwig"
                                                                             : "h_eff_pythia"));
        if (!h_gen_nom || !h_gen_alt) {
            std::cerr << "ERROR: efficiency histograms not found in " << cv.file << std::endl;
            fc->Close();
            continue;
        }
        if (h_gen_nom->GetNbinsX() != nbins) {
            std::cerr << "Correction variation '" << cv.name << "': " << h_gen_nom->GetNbinsX()
                      << " bins vs " << nbins << " in the result -- skipped" << std::endl;
            fc->Close();
            continue;
        }

        TH1D *h_var = (TH1D *) h_nom->Clone("h_var_" + cv.name);
        h_var->SetDirectory(nullptr);
        for (int i = 1; i <= nbins; ++i) {
            const double c_nom = h_gen_nom->GetBinContent(i);
            const double c_alt = h_gen_alt->GetBinContent(i);
            if (c_nom == 0. || c_alt == 0.) continue;
            const double scale = cv.inverse ? (c_nom / c_alt)   // correction divided into data
                                            : (c_alt / c_nom);  // correction multiplied in
            h_var->SetBinContent(i, h_var->GetBinContent(i) * scale);
            h_var->SetBinError(i, h_var->GetBinError(i) * scale);
        }
        fc->Close();

        // Self-normalise, exactly as the nominal is normalised.
        if (w.normalise_unit_area) {
            const double integral = h_var->Integral();
            if (integral > 0.) h_var->Scale(1. / integral, "width");
        }

        TH1D *h_shift  = (TH1D *) h_nom->Clone("h_syst_" + cv.name);
        TH1D *h_signed = (TH1D *) h_nom->Clone("h_delta_" + cv.name);
        h_shift->SetDirectory(nullptr);  h_shift->Reset();
        h_signed->SetDirectory(nullptr); h_signed->Reset();
        for (int i = 1; i <= nbins; ++i) {
            const double d = h_nom->GetBinContent(i) - h_var->GetBinContent(i);
            h_signed->SetBinContent(i, d);
            h_shift->SetBinContent(i, std::fabs(d));
        }

        std::cout << "Correction variation '" << cv.name << "' : " << cv.file << std::endl;
        h_shifts.push_back(h_shift);
        h_deltas.push_back(h_signed);
        h_vars.push_back(h_var);
        shift_names.push_back(cv.name);
        shift_labels.push_back(cv.label);
        shift_groups.push_back("");   // correction variations are independent sources
        shift_in_band.push_back(true);
        addTagLabel(cv.sysLabel);
    }

    if (h_shifts.empty())
        std::cerr << "WARNING: no variation available -- the output carries stat errors only"
                  << std::endl;

    // ---- Total systematic: an ASYMMETRIC quadrature sum -----------------------------
    // Two bands, sigma+ and sigma-, because a genuine up/down variation does not have to
    // move the result by the same amount in each direction. Which side a source feeds
    // depends on whether it is really two-sided:
    //
    //   GROUPED entries (a shared, non-empty Variation::group) are the two sides of ONE
    //     uncertainty. Each member is signed: the ones that push the result UP set sigma+,
    //     the ones that push it DOWN set sigma-, and the source contributes its largest
    //     shift on each side -- once, never both members to the same side.
    //
    //   UNGROUPED entries have a single alternative (a Pythia<->Herwig swap, the 3% track
    //     drop). One alternative tells you the SIZE of the difference, not its direction --
    //     nothing says the truth lies on the far side of the nominal -- so these are
    //     symmetrised: |delta| into BOTH sides, exactly as before.
    //
    // With only the mistag pair genuinely two-sided the two bands come out all but equal;
    // the point of the machinery is the asymmetric sources still to be added.
    //
    // Degenerate case: if both members of a pair happen to push the SAME way, one side
    // would get nothing from that source, which claims an uncertainty of zero there. That
    // is too aggressive, so such a group falls back to its envelope on both sides and says
    // so.
    TH1D *h_syst_up = (TH1D *) h_nom->Clone("h_syst_total_up");
    TH1D *h_syst_dn = (TH1D *) h_nom->Clone("h_syst_total_down");
    for (TH1D *h : {h_syst_up, h_syst_dn}) { h->SetDirectory(nullptr); h->Reset(); }
    std::vector<TString> one_sided_groups;
    for (int i = 1; i <= nbins; ++i) {
        double s2p = 0., s2m = 0.;
        std::map<TString, double> grp_up, grp_dn;
        for (size_t k = 0; k < h_shifts.size(); ++k) {
            if (!shift_in_band[k]) continue;   // display-only comparison, not a systematic
            // h_deltas is (nominal - variation); the shift OF THE RESULT is the negative.
            const double d = -h_deltas[k]->GetBinContent(i);
            if (shift_groups[k].Length()) {
                double &u = grp_up[shift_groups[k]];   // value-initialised to 0 on insert
                double &m = grp_dn[shift_groups[k]];
                if (d > 0. && d  > u) u =  d;
                if (d < 0. && -d > m) m = -d;
            } else {
                s2p += d * d;
                s2m += d * d;
            }
        }
        for (auto &g : grp_up) {
            double up = g.second, dn = grp_dn[g.first];
            if (up == 0. || dn == 0.) {
                // Both sides of the pair went the same way in this bin: symmetrise it.
                const double env = std::max(up, dn);
                up = dn = env;
                bool seen = false;
                for (const TString &t : one_sided_groups) if (t == g.first) seen = true;
                if (!seen) one_sided_groups.push_back(g.first);
            }
            s2p += up * up;
            s2m += dn * dn;
        }
        h_syst_up->SetBinContent(i, std::sqrt(s2p));
        h_syst_dn->SetBinContent(i, std::sqrt(s2m));
    }
    for (const TString &g : one_sided_groups)
        std::cerr << "NOTE: both sides of pair '" << g << "' shift the result the same way "
                  << "in at least one bin -- symmetrised there rather than claiming zero "
                  << "uncertainty on the other side." << std::endl;

    // The symmetrised total, for anything that wants one number per bin (the TH1 bin errors
    // below, and any downstream script that predates the split). The larger side, so it is
    // never the optimistic one.
    TH1D *h_syst = (TH1D *) h_nom->Clone("h_syst_total");
    h_syst->SetDirectory(nullptr);
    h_syst->Reset();
    for (int i = 1; i <= nbins; ++i)
        h_syst->SetBinContent(i, std::max(h_syst_up->GetBinContent(i),
                                          h_syst_dn->GetBinContent(i)));

    // ---- Result with stat + syst ----
    // h_total keeps the central values with the quadrature error in the bin errors; the
    // graph is what you draw as the systematic box around the points.
    TH1D *h_total = (TH1D *) h_nom->Clone("h_result_stat_syst");
    h_total->SetDirectory(nullptr);
    // Asymmetric: this is the band you draw, and the two sides are independent numbers.
    TGraphAsymmErrors *g_syst = new TGraphAsymmErrors(nbins);
    g_syst->SetName("g_syst_band");

    // One signed-shift column per variation -- the shift OF THE RESULT (variation - nominal),
    // so the sign in the table is the direction that source actually pushes the point, and
    // it matches which of syst+/syst- the source feeds.
    printf("\n bin |        value |     stat");
    for (size_t k = 0; k < shift_names.size(); ++k)
        printf(" | %15.15s", shift_in_band[k] ? shift_names[k].Data()
                                              : ("*" + shift_names[k]).Data());
    printf(" |    syst+ |    syst- |    total |  syst/value\n");
    printf("-----+--------------+---------");
    for (size_t k = 0; k < shift_names.size(); ++k) printf("-+----------------");
    printf("-+----------+----------+----------+------------\n");

    for (int i = 1; i <= nbins; ++i) {
        const double v     = h_nom->GetBinContent(i);
        const double stat  = h_nom->GetBinError(i);
        const double sup   = h_syst_up->GetBinContent(i);
        const double sdn   = h_syst_dn->GetBinContent(i);
        const double syst  = h_syst->GetBinContent(i);           // the larger side
        const double tot   = std::sqrt(stat * stat + syst * syst);

        h_total->SetBinError(i, tot);
        g_syst->SetPoint(i - 1, h_nom->GetBinCenter(i), v);
        g_syst->SetPointError(i - 1, 0.5 * h_nom->GetBinWidth(i), 0.5 * h_nom->GetBinWidth(i),
                              sdn, sup);   // low, high

        printf(" %3d | %12.5g | %8.4g", i, v, stat);
        for (TH1D *h : h_deltas) printf(" | %+15.4g", -h->GetBinContent(i));
        printf(" | %8.4g | %8.4g | %8.4g | %10.3f\n",
               sup, sdn, tot, (v != 0. ? syst / std::fabs(v) : 0.));
    }
    {
        bool any_display_only = false;
        for (bool b : shift_in_band) if (!b) any_display_only = true;
        if (any_display_only)
            printf("\n * = shown for comparison only; NOT included in syst+/syst-.\n");
    }

    // ---- Plot: every variation's EEC curve against the nominal ----------------------
    // Top pad: the EEC itself, nominal with its stat errors plus one curve per variation.
    // Bottom pad: variation/nominal - 1, i.e. the fractional pull of each systematic, with
    // the total band drawn around zero so each source can be read against the total.
    if (!h_vars.empty()) {
        // Palette convention: red = the measurement, the rest for the variations.
        const Color_t col_nom  = (Color_t) TColor::GetColor("#C44E52");
        // Light enough to sit behind the curves, dark enough to actually read. 0.10 is
        // effectively white on paper; raise or lower this one number to taste.
        const Color_t col_band = blendWithWhite("#C44E52", 0.2); // systematic band
        // Red stays the measurement; the variations take blue, green, purple, then the
        // spares. k % size() only reads as distinct while size() >= the number of variations
        // booked, so this list has to grow with the band -- at six, the three UParT SF
        // curves wrapped onto the same colours as the generator swaps.
        const std::vector<Color_t> col_var = { (Color_t) TColor::GetColor("#4C72B0"),   // blue
                                               (Color_t) TColor::GetColor("#4F8F52"),   // green
                                               (Color_t) TColor::GetColor("#8C6BB1"),   // purple
                                               (Color_t) TColor::GetColor("#C97430"),   // orange
                                               (Color_t) TColor::GetColor("#3D9E9E"),   // teal
                                               (Color_t) TColor::GetColor("#937860"),   // brown
                                               (Color_t) TColor::GetColor("#CC79A7"),   // pink
                                               (Color_t) TColor::GetColor("#7F7F2E"),   // olive
                                               (Color_t) TColor::GetColor("#56638A") }; // slate
        if (h_vars.size() > col_var.size())
            std::cerr << "NOTE: " << h_vars.size() << " variations but only " << col_var.size()
                      << " colours -- some curves repeat a colour and are told apart only by "
                      << "line style. Add entries to col_var." << std::endl;
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetLegendBorderSize(0);
        gStyle->SetLegendFillColor(0);

        // Same layout and fonts as the bottomline plot in apply_unfolding_2d.C.
        const Float_t font_scale   = 1200. / 800.;
        const Style_t font_code    = 43;
        const Float_t label_size   = 15. * font_scale;
        const Float_t title_size   = 15. * font_scale;
        const Float_t legend_size  = 15. * font_scale;
        const Float_t title_offset = 3.0;

        TCanvas *c_sys = new TCanvas("c_syst_curves", "", 800, 800);
        // 45/55 rather than the usual 30/70: on THIS plot the ratio pad is the content --
        // the per-source pulls are a few percent and all six curves sit on top of each
        // other in the top pad, so squeezing them into 0.3 of the canvas hides the answer.
        TPad *p_top = new TPad("p_top", "", 0., 0.45, 1., 1.);
        TPad *p_bot = new TPad("p_bot", "", 0., 0.,   1., 0.45);
        for (TPad *p : {p_top, p_bot}) {
            p->SetTicks(1, 0);
            p->SetFillColor(0);
        }
        p_top->SetMargin(0.1, 0.1, 0.0,  0.1);   // left, right, bottom, top
        p_bot->SetMargin(0.1, 0.1, 0.16, 0.0);
        c_sys->cd();
        p_top->Draw(); p_bot->Draw();

        // ---- top ----
        p_top->cd();
        TH1D *h_nom_draw = (TH1D *) h_nom->Clone("h_nom_draw");
        h_nom_draw->SetDirectory(nullptr);
        double ymax = 0.;
        for (int i = 1; i <= nbins; ++i)
            ymax = std::max(ymax, h_nom_draw->GetBinContent(i)
                                + h_nom_draw->GetBinError(i) + h_syst_up->GetBinContent(i));
        for (TH1D *h : h_vars)
            for (int i = 1; i <= nbins; ++i) ymax = std::max(ymax, h->GetBinContent(i));
        h_nom_draw->SetTitle("");
        h_nom_draw->GetYaxis()->SetRangeUser(0., ymax * 1.75);
        // TLatex (#Delta), not TMathText: TMathText axis titles do not render here.
        h_nom_draw->GetYaxis()->SetTitle("EEC(#Delta r)");
        h_nom_draw->GetYaxis()->CenterTitle(true);
        h_nom_draw->GetYaxis()->SetTitleFont(font_code); h_nom_draw->GetYaxis()->SetTitleSize(title_size);
        h_nom_draw->GetYaxis()->SetTitleOffset(1.5); // as the bottomline plot's main pad
        h_nom_draw->GetYaxis()->SetLabelFont(font_code); h_nom_draw->GetYaxis()->SetLabelSize(label_size);
        h_nom_draw->GetXaxis()->SetTitleSize(0);
        h_nom_draw->GetXaxis()->SetLabelSize(0);
        h_nom_draw->SetLineColor(col_nom);   h_nom_draw->SetMarkerColor(col_nom);
        h_nom_draw->SetMarkerStyle(kFullCircle); h_nom_draw->SetMarkerSize(1);
        h_nom_draw->SetLineWidth(2);
        h_nom_draw->Draw("AXIS");

        // Total systematic as a light red band on the nominal, drawn first so the points
        // and the variation curves stay readable on top of it. Same colour as the band in
        // the ratio pad below, so the two read as one quantity.
        g_syst->SetFillColor(col_band);
        g_syst->SetLineWidth(0);
        g_syst->Draw("2 same");

        // Two columns across the top: eight rows in one column reach down into the peak
        // whichever side they are parked on, and halving the row count clears it.
        const int n_rows = (2 + (int) h_vars.size() + 1) / 2;  // nominal + band + variations
        TLegend *leg = new TLegend(0.30, 0.87 - 0.0425 * n_rows, 0.95, 0.87);
        leg->SetNColumns(2);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetMargin(0.15);
        leg->SetTextFont(font_code);
        leg->SetTextSize(legend_size * 0.8);   // two columns: the long labels have half the width
        leg->AddEntry(h_nom_draw, "nominal (stat)", "pe1");
        leg->AddEntry(g_syst, "quadrature systematic", "f");

        for (size_t k = 0; k < h_vars.size(); ++k) {
            TH1D *h = h_vars[k];
            h->SetLineColor(col_var[k % col_var.size()]);
            h->SetMarkerColor(col_var[k % col_var.size()]);
            h->SetLineWidth(2);
            h->SetLineStyle(2 + (int) k);
            h->Draw("HIST same");
            leg->AddEntry(h, shift_labels[k], "l");
        }
        h_nom_draw->Draw("PE X0 same");
        leg->Draw();

        // Same banner as the bottomline plot.
        TLatex cms_bl;
        cms_bl.SetNDC();
        cms_bl.SetTextFont(62); cms_bl.SetTextSize(0.045); cms_bl.DrawLatex(0.10, 0.945, "CMS");
        cms_bl.SetTextFont(52); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.205, 0.945, "Internal");
        cms_bl.SetTextFont(42); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

        // ---- bottom: fractional pull of each variation ----
        p_bot->cd();
        TH1D *h_frame = (TH1D *) h_nom->Clone("h_ratio_frame");
        h_frame->SetDirectory(nullptr);
        h_frame->Reset();
        double rmax = 0.;
        for (int i = 1; i <= nbins; ++i) {
            const double v = h_nom->GetBinContent(i);
            if (v == 0.) continue;
            for (TH1D *h : h_vars)
                rmax = std::max(rmax, std::fabs(h->GetBinContent(i) / v - 1.));
            rmax = std::max(rmax, std::max(h_syst_up->GetBinContent(i),
                                           h_syst_dn->GetBinContent(i)) / std::fabs(v));
        }
        if (rmax <= 0.) rmax = 0.05;
        h_frame->SetTitle("");
        h_frame->GetYaxis()->SetRangeUser(-1.6 * rmax, 1.6 * rmax);
        h_frame->GetXaxis()->SetTitle("#Delta r");
        // Short enough to fit the pad height -- "variation / nominal - 1" ran off the end
        // and lost its "1", which made the panel look like it plotted a ratio, not a shift.
        h_frame->GetYaxis()->SetTitle("Variation / nominal #minus 1");
        h_frame->GetXaxis()->CenterTitle(true);
        h_frame->GetYaxis()->CenterTitle(true);
        h_frame->GetYaxis()->SetNdivisions(505);
        for (TAxis *ax : { h_frame->GetXaxis(), h_frame->GetYaxis() }) {
            ax->SetTitleFont(font_code);
            ax->SetTitleSize(title_size);
            ax->SetLabelFont(font_code);
            ax->SetLabelSize(label_size);
        }
        // title_offset (3.0) is what the bottomline ratio pad uses, but with these margins
        // it puts both titles outside the pad, so they are set to what actually renders.
        h_frame->GetXaxis()->SetTitleOffset(1.1);
        h_frame->GetYaxis()->SetTitleOffset(1.5);
        h_frame->Draw("AXIS");

        // Total systematic as a light red band around zero, for scale.
        TGraphAsymmErrors *g_band = new TGraphAsymmErrors(nbins);
        g_band->SetName("g_syst_band_rel");
        for (int i = 1; i <= nbins; ++i) {
            const double v  = h_nom->GetBinContent(i);
            const double hw = 0.5 * h_nom->GetBinWidth(i);
            g_band->SetPoint(i - 1, h_nom->GetBinCenter(i), 0.);
            g_band->SetPointError(i - 1, hw, hw,
                                  (v != 0.) ? h_syst_dn->GetBinContent(i) / std::fabs(v) : 0.,
                                  (v != 0.) ? h_syst_up->GetBinContent(i) / std::fabs(v) : 0.);
        }
        g_band->SetFillColor(col_band);
        g_band->SetLineWidth(0);
        g_band->Draw("2 same");

        for (size_t k = 0; k < h_vars.size(); ++k) {
            TH1D *h_r = (TH1D *) h_vars[k]->Clone("h_relshift_" + shift_names[k]);
            h_r->SetDirectory(nullptr);
            for (int i = 1; i <= nbins; ++i) {
                const double v = h_nom->GetBinContent(i);
                h_r->SetBinContent(i, (v != 0.) ? h_vars[k]->GetBinContent(i) / v - 1. : 0.);
                h_r->SetBinError(i, 0.);
            }
            h_r->SetLineColor(col_var[k % col_var.size()]);
            h_r->SetLineWidth(2);
            h_r->SetLineStyle(2 + (int) k);
            h_r->Draw("HIST same");
        }

        TLine *l0 = new TLine(h_frame->GetXaxis()->GetXmin(), 0.,
                              h_frame->GetXaxis()->GetXmax(), 0.);
        l0->SetLineStyle(2);
        l0->Draw();
        p_bot->RedrawAxis();

        const TString plot_stem = resultFolder(SAMPLE, GENERATOR, TF_GENERATOR)
                                + "systematics_curves_" + label;
        c_sys->Print(plot_stem + ".pdf");
        c_sys->Print(plot_stem + ".png");
        std::cout << "\nWrote " << plot_stem << ".{pdf,png}" << std::endl;
    }

    // ---- Plot: the final measurement against particle-level MC ----------------------
    // The money plot. Red = the measurement (unfolded data, stat bars + systematic box),
    // blue = particle-level Pythia, purple = particle-level Herwig when that unfolding
    // exists. h_mc_true_2D is hgenjet_2b_all in data mode: the gen reference the corrected
    // data is supposed to land on, written by apply_unfolding_2d.C already restricted to
    // the drawn dr range and unit-normalised, exactly like the result next to it.
    //
    // The gen curves are only RE-normalised here, never weighted: the scale factors in w
    // correct data for detector effects and have no meaning at particle level.
    {
        auto loadGen = [&](const TString &file, const TString &newname) -> TH1D * {
            TFile *f = TFile::Open(file);
            if (!f || f->IsZombie()) { if (f) delete f; return nullptr; }
            TH1D *h = dynamic_cast<TH1D *>(f->Get("h_mc_true_2D"));
            if (!h) { f->Close(); return nullptr; }
            TH1D *out = (TH1D *) h->Clone(newname);
            out->SetDirectory(nullptr);
            f->Close();
            if (out->GetNbinsX() != nbins) {
                std::cerr << "Gen MC '" << newname << "': " << out->GetNbinsX() << " bins vs "
                          << nbins << " in the result -- skipped" << std::endl;
                delete out;
                return nullptr;
            }
            if (w.normalise_unit_area) {
                const double integral = out->Integral();
                if (integral > 0.) out->Scale(1. / integral, "width");
            }
            return out;
        };

        // Draw the OTHER generator's particle level alongside the nominal one? Off: the
        // plot is the measurement against one theory curve. Flip to true to put Herwig
        // back next to Pythia -- nothing else has to change, every step below is guarded
        // on h_gen_alt being non-null.
        const bool show_alt_gen = false;

        TH1D *h_gen_nom = loadGen(nominalFile, "h_gen_" + GENERATOR);
        TH1D *h_gen_alt = show_alt_gen
                        ? loadGen(resultFile(SAMPLE, other_gen, test_mode, unfoldBayes,
                                             TF_GENERATOR, false), "h_gen_" + other_gen)
                        : nullptr;

        if (!h_gen_nom) {
            std::cerr << "\nNo particle-level MC in " << nominalFile
                      << " (h_mc_true_2D) -- skipping the data-vs-gen plot.\n"
                      << "   That histogram's Write() was disabled until 2026-09-10, so a "
                      << "result unfolded before then does not carry it.\n"
                      << "   Re-run apply_unfolding_2d.C(\"" << SAMPLE << "\",\"" << GENERATOR
                      << "\") to get it." << std::endl;
        } else {
            const Color_t col_data = (Color_t) TColor::GetColor("#C44E52");  // measurement
            const Color_t col_pl   = (Color_t) TColor::GetColor("#4C72B0");  // particle-level MC
            const Color_t col_pl2  = (Color_t) TColor::GetColor("#8C6BB1");  // second PL MC
            const Color_t col_box  = blendWithWhite("#C44E52", 0.2);
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            gStyle->SetLegendBorderSize(0);
            gStyle->SetLegendFillColor(0);

            const Float_t font_scale  = 1200. / 800.;
            const Style_t font_code   = 43;
            const Float_t label_size  = 15. * font_scale;
            const Float_t title_size  = 15. * font_scale;
            const Float_t legend_size = 15. * font_scale;

            TCanvas *c_gen = new TCanvas("c_data_vs_gen", "", 800, 800);
            TPad *q_top = new TPad("q_top", "", 0., 0.3, 1., 1.);
            TPad *q_bot = new TPad("q_bot", "", 0., 0.,  1., 0.3);
            for (TPad *q : {q_top, q_bot}) { q->SetTicks(1, 0); q->SetFillColor(0); }
            q_top->SetMargin(0.1, 0.1, 0.0,  0.1);
            q_bot->SetMargin(0.1, 0.1, 0.23, 0.0);
            c_gen->cd();
            q_top->Draw(); q_bot->Draw();

            // ---- top: EEC ----
            q_top->cd();
            TH1D *h_meas = (TH1D *) h_nom->Clone("h_meas_draw");
            h_meas->SetDirectory(nullptr);
            double ymax = 0.;
            for (int i = 1; i <= nbins; ++i)
                ymax = std::max(ymax, h_nom->GetBinContent(i) + h_nom->GetBinError(i)
                                    + h_syst_up->GetBinContent(i));
            for (TH1D *h : {h_gen_nom, h_gen_alt})
                if (h) for (int i = 1; i <= nbins; ++i)
                    ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
            h_meas->SetTitle("");
            // Headroom: the legend sits upper right, the systematics tag upper left, and
            // both have to clear the peak. 1.75 is what fits the tag's rows above it.
            h_meas->GetYaxis()->SetRangeUser(0., ymax * 1.75);
            h_meas->GetYaxis()->SetTitle("EEC(#Delta r)");
            h_meas->GetYaxis()->CenterTitle(true);
            h_meas->GetYaxis()->SetTitleFont(font_code); h_meas->GetYaxis()->SetTitleSize(title_size);
            h_meas->GetYaxis()->SetTitleOffset(1.5);
            h_meas->GetYaxis()->SetLabelFont(font_code); h_meas->GetYaxis()->SetLabelSize(label_size);
            h_meas->GetXaxis()->SetTitleSize(0);
            h_meas->GetXaxis()->SetLabelSize(0);
            h_meas->SetLineColor(col_data);   h_meas->SetMarkerColor(col_data);
            h_meas->SetMarkerStyle(kFullCircle); h_meas->SetMarkerSize(1);
            h_meas->SetLineWidth(2);
            h_meas->Draw("AXIS");

            // Systematic as a box on each point, drawn under the MC so the curves stay read
            g_syst->SetFillColor(col_box);
            g_syst->SetLineWidth(0);
            g_syst->Draw("2 same");

            // Upper RIGHT, unlike the systematics plot: the EEC peaks at small dr and the
            // left half of this pad is where the curves and the biggest systematic box are.
            // Sized by entry count -- a TLegend spreads its rows over the whole box, so a
            // fixed height leaves three entries floating with too much air between them.
            const int n_leg = 3 + (h_gen_alt ? 1 : 0);
            TLegend *leg_g = new TLegend(0.55, 0.87 - 0.0625 * n_leg, 0.93, 0.87);
            leg_g->SetFillStyle(0);
            leg_g->SetBorderSize(0);
            leg_g->SetMargin(0.15);
            leg_g->SetTextFont(font_code);
            leg_g->SetTextSize(legend_size);
            leg_g->AddEntry(h_meas, "Data, unfolded (stat)", "pe1");
            leg_g->AddEntry(g_syst, "Systematic uncertainty", "f");

            h_gen_nom->SetLineColor(col_pl); h_gen_nom->SetMarkerColor(col_pl);
            h_gen_nom->SetLineWidth(2); h_gen_nom->SetLineStyle(1);
            h_gen_nom->Draw("HIST same");
            leg_g->AddEntry(h_gen_nom, prettyGen(GENERATOR) + " (particle level)", "l");
            if (h_gen_alt) {
                h_gen_alt->SetLineColor(col_pl2); h_gen_alt->SetMarkerColor(col_pl2);
                h_gen_alt->SetLineWidth(2); h_gen_alt->SetLineStyle(2);
                h_gen_alt->Draw("HIST same");
                leg_g->AddEntry(h_gen_alt, prettyGen(other_gen) + " (particle level)", "l");
            }
            h_meas->Draw("PE X0 same");
            leg_g->Draw();

            TLatex cms_g;
            cms_g.SetNDC();
            cms_g.SetTextFont(62); cms_g.SetTextSize(0.045); cms_g.DrawLatex(0.10, 0.945, "CMS");
            cms_g.SetTextFont(52); cms_g.SetTextSize(0.037); cms_g.DrawLatex(0.205, 0.945, "Internal");
            cms_g.SetTextFont(42); cms_g.SetTextSize(0.037); cms_g.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

            // ---- "systematics applied" tag, upper left ----
            // Listed from syst_tag_labels, which is built from the entries that actually
            // entered the band -- so a source whose result file was missing and got skipped
            // cannot be claimed here.
            if (!syst_tag_labels.empty()) {
                TLatex tag;
                tag.SetNDC();
                tag.SetTextFont(42);
                tag.SetTextSize(0.030);
                tag.SetTextAlign(13);           // left, top
                double y_tag = 0.855;
                tag.SetTextFont(62);
                tag.DrawLatex(0.145, y_tag, "Systematic uncertainties");
                tag.SetTextFont(42);
                for (const TString &t : syst_tag_labels) {
                    y_tag -= 0.042;
                    tag.DrawLatex(0.160, y_tag, "#bullet " + t);
                }
            }

            // ---- bottom: MC / data ----
            // Ratio to the MEASUREMENT, so the data uncertainty shows as a band around 1 and
            // each MC reads directly against it. Each ratio keeps its numerator's colour.
            q_bot->cd();
            TH1D *q_frame = (TH1D *) h_nom->Clone("h_gen_ratio_frame");
            q_frame->SetDirectory(nullptr);
            q_frame->Reset();
            double qmax = 0.;
            for (int i = 1; i <= nbins; ++i) {
                const double v = h_nom->GetBinContent(i);
                if (v == 0.) continue;
                for (TH1D *h : {h_gen_nom, h_gen_alt})
                    if (h) qmax = std::max(qmax, std::fabs(h->GetBinContent(i) / v - 1.));
                qmax = std::max(qmax,
                                (std::max(h_syst_up->GetBinContent(i), h_syst_dn->GetBinContent(i))
                                 + h_nom->GetBinError(i)) / std::fabs(v));
            }
            if (qmax <= 0.) qmax = 0.05;
            q_frame->SetTitle("");
            // 1.7, not 1.5: with one MC curve qmax is set by the bin-1 systematic alone and
            // 1.5 puts the frame edge just under a tick, so ROOT draws the top label half
            // outside the pad. The extra headroom keeps every label inside.
            q_frame->GetYaxis()->SetRangeUser(1. - 1.7 * qmax, 1. + 1.7 * qmax);
            q_frame->GetXaxis()->SetTitle("#Delta r");
            q_frame->GetYaxis()->SetTitle("MC / Data");
            q_frame->GetXaxis()->CenterTitle(true);
            q_frame->GetYaxis()->CenterTitle(true);
            q_frame->GetYaxis()->SetNdivisions(505);
            for (TAxis *ax : { q_frame->GetXaxis(), q_frame->GetYaxis() }) {
                ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
                ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
            }
            q_frame->GetXaxis()->SetTitleOffset(1.1);
            q_frame->GetYaxis()->SetTitleOffset(1.5);
            q_frame->Draw("AXIS");

            // Data uncertainty around 1: the systematic as the wide light box, the
            // statistical error as a narrower, more saturated box drawn on top of it. Both
            // solid -- a hatched fill at this width renders as a few stray ticks and reads
            // as nothing at all.
            TGraphAsymmErrors *g_dsyst = new TGraphAsymmErrors(nbins);
            TGraphAsymmErrors *g_dstat = new TGraphAsymmErrors(nbins);
            g_dsyst->SetName("g_data_syst_rel");
            g_dstat->SetName("g_data_stat_rel");
            for (int i = 1; i <= nbins; ++i) {
                const double v = h_nom->GetBinContent(i);
                const double x = h_nom->GetBinCenter(i);
                const double hw = 0.5 * h_nom->GetBinWidth(i);
                const double st = (v != 0.) ? h_nom->GetBinError(i) / std::fabs(v) : 0.;
                g_dsyst->SetPoint(i - 1, x, 1.);
                g_dstat->SetPoint(i - 1, x, 1.);
                // The systematic box is asymmetric about 1; the statistical one is not.
                g_dsyst->SetPointError(i - 1, hw, hw,
                                       (v != 0.) ? h_syst_dn->GetBinContent(i) / std::fabs(v) : 0.,
                                       (v != 0.) ? h_syst_up->GetBinContent(i) / std::fabs(v) : 0.);
                g_dstat->SetPointError(i - 1, hw * 0.30, hw * 0.30, st, st);
            }
            g_dsyst->SetFillColor(col_box);                        g_dsyst->SetLineWidth(0);
            g_dsyst->Draw("2 same");
            g_dstat->SetFillColor(blendWithWhite("#C44E52", 0.5)); g_dstat->SetLineWidth(0);
            g_dstat->Draw("2 same");

            for (TH1D *h : {h_gen_nom, h_gen_alt}) {
                if (!h) continue;
                TH1D *h_q = (TH1D *) h->Clone(TString(h->GetName()) + "_ratio");
                h_q->SetDirectory(nullptr);
                for (int i = 1; i <= nbins; ++i) {
                    const double v = h_nom->GetBinContent(i);
                    h_q->SetBinContent(i, (v != 0.) ? h->GetBinContent(i) / v : 0.);
                    h_q->SetBinError(i, 0.);
                }
                h_q->SetLineColor(h->GetLineColor());
                h_q->SetLineWidth(2);
                h_q->SetLineStyle(h->GetLineStyle());
                h_q->Draw("HIST same");
            }

            TLine *l1 = new TLine(q_frame->GetXaxis()->GetXmin(), 1.,
                                  q_frame->GetXaxis()->GetXmax(), 1.);
            l1->SetLineStyle(2);
            l1->Draw();
            q_bot->RedrawAxis();

            const TString gen_stem = resultFolder(SAMPLE, GENERATOR, TF_GENERATOR)
                                   + "data_vs_gen_" + label;
            c_gen->Print(gen_stem + ".pdf");
            c_gen->Print(gen_stem + ".png");
            std::cout << "\nWrote " << gen_stem << ".{pdf,png}" << std::endl;

            // Kept for the write block below.
            h_gen_curves.push_back(h_gen_nom);
            gen_curve_names.push_back("h_gen_" + GENERATOR);
            if (h_gen_alt) {
                h_gen_curves.push_back(h_gen_alt);
                gen_curve_names.push_back("h_gen_" + other_gen);
            }
        }
    }

    // ---- Plot: the UParT SF, before and after, and nothing else ----------------------
    // The systematics plot carries the before-SF curve too, but buried among every other
    // variation. This one shows the correction on its own: the two curves and their ratio,
    // which is what you look at to judge the SF itself.
    //
    // Driven by the display-only "sfupart_off" entry, so it appears only when that entry is
    // booked and its unfolding exists. Comment the entry out and this plot stops being
    // produced, along with the curve on the systematics plot.
    {
        int k_off = -1;
        for (size_t k = 0; k < shift_names.size(); ++k)
            if (shift_names[k] == "sfupart_off") k_off = (int) k;

        if (k_off >= 0) {
            TH1D *h_before = h_vars[k_off];   // no SF applied
            TH1D *h_after  = h_nom;           // the nominal: SF applied

            // Palette: red is the measurement, and the result WITH the SF is the
            // measurement. Before-SF is an intermediate correction stage -> green.
            const Color_t col_after  = (Color_t) TColor::GetColor("#C44E52");
            const Color_t col_before = (Color_t) TColor::GetColor("#4F8F52");
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            gStyle->SetLegendBorderSize(0);
            gStyle->SetLegendFillColor(0);

            const Float_t font_scale  = 1200. / 800.;
            const Style_t font_code   = 43;
            const Float_t label_size  = 15. * font_scale;
            const Float_t title_size  = 15. * font_scale;
            const Float_t legend_size = 15. * font_scale;

            TCanvas *c_sf = new TCanvas("c_sfupart_before_after", "", 800, 800);
            TPad *s_top = new TPad("s_top", "", 0., 0.35, 1., 1.);
            TPad *s_bot = new TPad("s_bot", "", 0., 0.,   1., 0.35);
            for (TPad *q : {s_top, s_bot}) { q->SetTicks(1, 0); q->SetFillColor(0); }
            s_top->SetMargin(0.1, 0.1, 0.0,  0.1);
            s_bot->SetMargin(0.1, 0.1, 0.20, 0.0);
            c_sf->cd();
            s_top->Draw(); s_bot->Draw();

            // ---- top: the two EEC curves ----
            s_top->cd();
            TH1D *h_a = (TH1D *) h_after->Clone("h_sfupart_after_draw");
            TH1D *h_b = (TH1D *) h_before->Clone("h_sfupart_before_draw");
            for (TH1D *h : {h_a, h_b}) h->SetDirectory(nullptr);

            double ymax = 0.;
            for (int i = 1; i <= nbins; ++i)
                ymax = std::max(ymax, std::max(h_a->GetBinContent(i) + h_a->GetBinError(i),
                                               h_b->GetBinContent(i) + h_b->GetBinError(i)));
            h_a->SetTitle("");
            h_a->GetYaxis()->SetRangeUser(0., ymax * 1.45);
            h_a->GetYaxis()->SetTitle("EEC(#Delta r)");
            h_a->GetYaxis()->CenterTitle(true);
            h_a->GetYaxis()->SetTitleFont(font_code); h_a->GetYaxis()->SetTitleSize(title_size);
            h_a->GetYaxis()->SetTitleOffset(1.5);
            h_a->GetYaxis()->SetLabelFont(font_code); h_a->GetYaxis()->SetLabelSize(label_size);
            h_a->GetXaxis()->SetTitleSize(0);
            h_a->GetXaxis()->SetLabelSize(0);

            h_b->SetLineColor(col_before); h_b->SetMarkerColor(col_before);
            h_b->SetLineWidth(2); h_b->SetLineStyle(2);
            h_a->SetLineColor(col_after);  h_a->SetMarkerColor(col_after);
            h_a->SetMarkerStyle(kFullCircle); h_a->SetMarkerSize(1); h_a->SetLineWidth(2);

            h_a->Draw("AXIS");
            h_b->Draw("HIST same");
            h_a->Draw("PE X0 same");

            TLegend *leg_sf = new TLegend(0.52, 0.68, 0.93, 0.85);
            leg_sf->SetFillStyle(0);
            leg_sf->SetBorderSize(0);
            leg_sf->SetMargin(0.15);
            leg_sf->SetTextFont(font_code);
            leg_sf->SetTextSize(legend_size);
            leg_sf->AddEntry(h_b, "Before UParT SF", "l");
            leg_sf->AddEntry(h_a, "After UParT SF",  "pe1");
            leg_sf->Draw();

            TLatex cms_sf;
            cms_sf.SetNDC();
            cms_sf.SetTextFont(62); cms_sf.SetTextSize(0.045); cms_sf.DrawLatex(0.10, 0.945, "CMS");
            cms_sf.SetTextFont(52); cms_sf.SetTextSize(0.037); cms_sf.DrawLatex(0.205, 0.945, "Internal");
            cms_sf.SetTextFont(42); cms_sf.SetTextSize(0.037); cms_sf.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

            // ---- bottom: after / before ----
            // The effect of the SF, on its own. Not 1/SF: the unit-area renormalisation and
            // the migration through the unfolding both change it, which is the whole reason
            // the SF is applied at reco level rather than here.
            s_bot->cd();
            TH1D *h_r = (TH1D *) h_a->Clone("h_sfupart_ratio");
            h_r->SetDirectory(nullptr);
            double rlo = 1., rhi = 1.;
            for (int i = 1; i <= nbins; ++i) {
                const double b = h_b->GetBinContent(i);
                const double r = (b != 0.) ? h_a->GetBinContent(i) / b : 1.;
                h_r->SetBinContent(i, r);
                h_r->SetBinError(i, 0.);
                rlo = std::min(rlo, r);
                rhi = std::max(rhi, r);
            }
            const double pad = std::max(0.02, 0.25 * (rhi - rlo));
            h_r->SetTitle("");
            h_r->GetYaxis()->SetRangeUser(rlo - pad, rhi + pad);
            h_r->GetXaxis()->SetTitle("#Delta r");
            h_r->GetYaxis()->SetTitle("After / Before");
            h_r->GetXaxis()->CenterTitle(true);
            h_r->GetYaxis()->CenterTitle(true);
            h_r->GetYaxis()->SetNdivisions(505);
            for (TAxis *ax : { h_r->GetXaxis(), h_r->GetYaxis() }) {
                ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
                ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
            }
            h_r->GetXaxis()->SetTitleOffset(1.0);
            h_r->GetYaxis()->SetTitleOffset(1.5);
            h_r->SetLineColor(col_after);
            h_r->SetLineWidth(2);
            h_r->SetMarkerSize(0);
            h_r->Draw("HIST");

            TLine *l_one = new TLine(h_r->GetXaxis()->GetXmin(), 1.,
                                     h_r->GetXaxis()->GetXmax(), 1.);
            l_one->SetLineStyle(2);
            l_one->Draw();
            s_bot->RedrawAxis();

            const TString sfupart_stem = resultFolder(SAMPLE, GENERATOR, TF_GENERATOR)
                                  + "sfupart_before_after_" + label;
            c_sf->Print(sfupart_stem + ".pdf");
            c_sf->Print(sfupart_stem + ".png");
            std::cout << "\nWrote " << sfupart_stem << ".{pdf,png}" << std::endl;
        }
    }

    // ---- Write ----
    const TString fout_name = resultFolder(SAMPLE, GENERATOR, TF_GENERATOR)
                            + "final_" + label + "_with_systematics.root";
    TFile fout(fout_name, "RECREATE");
    h_nom->Write("h_nominal");          // central values, stat errors
    h_syst_up->Write("h_syst_total_up");   // absolute systematic per bin, upward side
    h_syst_dn->Write("h_syst_total_down"); // ... and downward side
    h_syst->Write("h_syst_total");         // the larger of the two, for one-number consumers
    h_total->Write("h_result_stat_syst");  // central values, stat (+) the larger syst side
    g_syst->Write("g_syst_band");          // TGraphAsymmErrors: the band to draw
    for (size_t k = 0; k < h_shifts.size(); ++k) {
        h_shifts[k]->Write("h_syst_"  + shift_names[k]);   // symmetrised |delta|
        // NOTE the sign: the STORED convention is (nominal - variation), unchanged, because
        // plot_systematics_summary.C reads these and flips the sign itself. The printed
        // table above shows the negative of this -- the shift of the RESULT -- so that its
        // sign matches which of syst+/syst- the source feeds. Same numbers, opposite sign.
        h_deltas[k]->Write("h_delta_" + shift_names[k]);   // signed nominal - variation
        h_vars[k]  ->Write("h_var_"   + shift_names[k]);   // the variation's own EEC curve
    }
    for (size_t k = 0; k < h_gen_curves.size(); ++k)
        h_gen_curves[k]->Write(gen_curve_names[k]);        // particle-level MC, unit area
    fout.Close();

    std::cout << "\nWrote " << fout_name << std::endl;

    // ---- What of all this is actually a result --------------------------------------
    // The results area accumulates one folder per flag combination, and most of them are
    // INPUTS to the band, not results. Print the short list, so the deliverables can be
    // found without reading the code or guessing from folder names.
    {
        const TString out = resultFolder(SAMPLE, GENERATOR, TF_GENERATOR);
        printf("\n");
        printf("================================ FINAL OUTPUTS ================================\n");
        printf("All in %s\n\n", out.Data());
        printf("  %-46s %s\n", TString("final_" + label + "_with_systematics.root").Data(),
                                "the result: central values, stat, syst+/-");
        printf("  %-46s %s\n", TString("data_vs_gen_" + label + ".pdf").Data(),
                                "the measurement vs particle-level MC");
        printf("  %-46s %s\n", TString("systematics_curves_" + label + ".pdf").Data(),
                                "per-source breakdown");
        bool has_sf_cmp = false;
        for (const TString &n : shift_names) if (n == "sfupart_off") has_sf_cmp = true;
        if (has_sf_cmp)
            printf("  %-46s %s\n", TString("sfupart_before_after_" + label + ".pdf").Data(),
                                    "the UParT SF on its own (comparison)");
        printf("\nEvery OTHER unfolding_* folder is an INPUT to the band -- one systematic\n");
        printf("variation each -- not a result. The variation list above names them.\n");
        printf("===============================================================================\n");
    }
}
