#ifndef RESULT_PATHS_H
#define RESULT_PATHS_H

// Where unfolding results live, and how a systematic variation is spelled in a path.
//
// WHY THIS FILE EXISTS
// --------------------
// apply_unfolding_2d.C WRITES these paths and apply_weights_and_systematics.C READS them.
// The two used to build them from their own private copies of the same mapping functions,
// which is exactly the kind of duplication that fails silently: on 2026-09-11 a new
// variation ("off", the no-SF run) was added to one copy and not the other. The reader's
// copy did not recognise it, fell through to "", and resolved to the NOMINAL folder -- so
// the macro compared the nominal against itself and printed a column of -0. Nothing
// errored. It looked like "the scale factor does nothing" rather than "I opened the wrong
// file".
//
// Both macros now include this header, so there is ONE definition of every mapping and the
// two cannot drift. Two further guards make a mistake loud instead of silent:
//
//   1. An unrecognised variation returns kUnknownTag, not "". The path then fails to open
//      with a message naming the bad value, instead of aliasing onto the nominal.
//   2. isKnown*() lets a caller validate its whole list of variations up front, before any
//      file is touched.
//
// If you add a variation, add it in ONE place -- here -- and both macros follow.

#include "TString.h"
#include <iostream>

// A tag no path will ever contain, returned for a variation name that is not recognised.
// The point is that it must NOT be "": an empty tag is the nominal, so an unknown value
// silently becoming nominal is the bug this header exists to prevent.
inline const TString &kUnknownTag()
{
    static const TString t = "_UNKNOWN_VARIATION";
    return t;
}

inline void warnUnknown(const char *what, const TString &value, const char *allowed)
{
    std::cerr << "ERROR: unknown " << what << " '" << value << "' (use " << allowed << "). "
              << "Paths built from it carry " << kUnknownTag() << " and will not open -- "
              << "this is deliberate, so the value cannot alias onto the nominal result."
              << std::endl;
}

// ---- Known values -------------------------------------------------------------------
// The canonical lists. Validate against these before doing any work.
inline bool isKnownGenerator(const TString &g)
{
    return g == "pythia" || g == "herwig";
}

inline bool isKnownSample(const TString &s)
{
    return s == "qcd" || s == "bjet" || s == "both";
}

// Which of template_fit.cpp's fits supplies the signal fraction. These are the varNames[]
// strings of Help_Functions.h, which are also the prefixes on its per-variation output
// files, so a fit file can be found from the variation name alone.
inline bool isKnownTfVariation(const TString &v)
{
    return v == "nominal" || v == "var0B_2" || v == "var0B_0";
}

// Which UParT b-tag efficiency scale factor is divided into the reco-level data.
// "off" applies no SF at all -- a comparison, not a systematic.
inline bool isKnownSfupartVariation(const TString &v)
{
    return v == "nominal" || v == "off"
        || v == "jpcalib_hf" || v == "qqrate_up" || v == "qqrate_down"
        || v == "statup"     || v == "statdn";
}

// ---- Path tags ----------------------------------------------------------------------
// "" means nominal in every one of these, so that nominal paths carry no tag at all and
// are unchanged from before any variation existed.

// TF_GENERATOR: which generator's template fit supplied the signal fraction.
inline TString tfTag(const TString &tfGenerator)
{
    if (!isKnownGenerator(tfGenerator)) {
        warnUnknown("TF_GENERATOR", tfGenerator, "pythia | herwig");
        return kUnknownTag();
    }
    return (tfGenerator != "pythia") ? ("_tf" + tfGenerator) : TString("");
}

// TRACK_EFF_UNC: the 3% track-drop MC production. The literal must stay in step with
// TrkEffSyst::kDropFraction in tracking_efficiency_syst.h.
inline TString trkTag(bool trackEffUnc)
{
    return trackEffUnc ? TString("_trkdrop030") : TString("");
}

// TF_VARIATION: the template-fit variation supplying the signal fraction.
inline TString tfVarTag(const TString &tfVariation)
{
    if (!isKnownTfVariation(tfVariation)) {
        warnUnknown("TF_VARIATION", tfVariation, "nominal | var0B_2 | var0B_0");
        return kUnknownTag();
    }
    if (tfVariation == "var0B_2") return "_tf0Bx2";
    if (tfVariation == "var0B_0") return "_tf0Bx0";
    return "";
}

// SFUPART_VARIATION: the UParT SF variation. Namespaced to UParT because more scale
// factors are coming and each should get its own mapping rather than share a generic one.
inline TString sfupartVarTag(const TString &sfupartVariation)
{
    if (!isKnownSfupartVariation(sfupartVariation)) {
        warnUnknown("SFUPART_VARIATION", sfupartVariation,
                    "nominal | off | jpcalib_hf | qqrate_up | qqrate_down | statup | statdn");
        return kUnknownTag();
    }
    if (sfupartVariation == "off")         return "_sfupartoff";
    if (sfupartVariation == "jpcalib_hf")  return "_sfupartjphf";
    if (sfupartVariation == "qqrate_up")   return "_sfupartqqup";
    if (sfupartVariation == "qqrate_down") return "_sfupartqqdn";
    if (sfupartVariation == "statup")      return "_sfupartstatup";
    if (sfupartVariation == "statdn")      return "_sfupartstatdn";
    return "";
}

// ---- The UParT scale factor itself ---------------------------------------------------
inline TString sfupartFile()
{
    // "finalbins": source bins 8-9 merged into one output bin, so this file has 8 dr bins
    // against the analysis's 9. Its last bin spans 0.35-0.45 and therefore applies to BOTH
    // of the analysis's last two bins -- the application maps by bin centre so that happens
    // on its own. The merge also tamed the last bins: 1.149 +/- 0.020 and 1.086 +/- 0.029
    // became a single 1.010 +/- 0.017.
    // The 9-bin predecessor (..._eec_systematics.root, no "finalbins") is superseded.
    return "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow/"
           "lifetime_jp_sfb_agg_dr_wp712_eec_finalbins_systematics.root";
}

// NOTE: these are the KEY names in that file. The histograms' own internal names differ
// (key "h_SFb_dr_central" holds a histogram named "h_SFb_dr_wp712_eec_central"), and
// TFile::Get() takes the key -- the internal name returns null.
// "off" returns the marker "(none)": no SF is applied, the caller checks for it.
inline TString sfupartHist(const TString &sfupartVariation)
{
    if (!isKnownSfupartVariation(sfupartVariation)) {
        warnUnknown("SFUPART_VARIATION", sfupartVariation,
                    "nominal | off | jpcalib_hf | qqrate_up | qqrate_down | statup | statdn");
        return "";
    }
    if (sfupartVariation == "off")         return "(none)";
    // statup/statdn shift the CENTRAL SF by +/- its own bin error, so they read the same
    // histogram; the shift is applied where the SF is used, not here. Coherent across all
    // dr bins, because the analysis's last two bins share ONE calibration bin and their
    // errors are therefore 100% correlated -- a coherent shift gets that right, and is the
    // conservative choice elsewhere.
    if (sfupartVariation == "nominal"
     || sfupartVariation == "statup"
     || sfupartVariation == "statdn")      return "h_SFb_dr_central";
    if (sfupartVariation == "jpcalib_hf")  return "h_SFb_dr_jpcalib_hf";
    if (sfupartVariation == "qqrate_up")   return "h_SFb_dr_qqrate_up";
    if (sfupartVariation == "qqrate_down") return "h_SFb_dr_qqrate_down";
    return "";
}

// ---- Result paths --------------------------------------------------------------------
// The combined tag every path carries. One function, so a folder and the label inside it
// can never disagree about which variation they are.
inline TString variationTag(const TString &tfGenerator, bool trackEffUnc,
                            const TString &tfVariation, const TString &sfupartVariation)
{
    return tfTag(tfGenerator) + trkTag(trackEffUnc)
         + tfVarTag(tfVariation) + sfupartVarTag(sfupartVariation);
}

inline TString resultFolder(const TString &sample, const TString &generator,
                            const TString &tfGenerator = "pythia",
                            bool trackEffUnc = false,
                            const TString &tfVariation = "nominal",
                            const TString &sfupartVariation = "nominal")
{
    return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/unfolding_"
         + sample + "_" + generator
         + variationTag(tfGenerator, trackEffUnc, tfVariation, sfupartVariation)
         + "_upartv2/";
}

// <sample>_<generator><variation tags>_<bayesian|MI>_<mode>.
// This is the "dataset" string apply_unfolding_2d.C names its output histograms file with,
// so the reader reconstructs exactly what the writer wrote.
inline TString resultLabel(const TString &sample, const TString &generator,
                           int test_mode, bool unfoldBayes,
                           const TString &tfGenerator = "pythia",
                           bool trackEffUnc = false,
                           const TString &tfVariation = "nominal",
                           const TString &sfupartVariation = "nominal")
{
    TString label = sample + "_" + generator
                  + variationTag(tfGenerator, trackEffUnc, tfVariation, sfupartVariation);
    label += unfoldBayes ? "_bayesian" : "_MI";
    if (test_mode == 0) label += "_full_closure";
    if (test_mode == 1) label += "_split_test";
    if (test_mode == 2) label += "_data";
    return label;
}

inline TString resultFile(const TString &sample, const TString &generator,
                          int test_mode, bool unfoldBayes,
                          const TString &tfGenerator = "pythia",
                          bool trackEffUnc = false,
                          const TString &tfVariation = "nominal",
                          const TString &sfupartVariation = "nominal")
{
    return resultFolder(sample, generator, tfGenerator, trackEffUnc, tfVariation,
                        sfupartVariation)
         + "histos_" + resultLabel(sample, generator, test_mode, unfoldBayes, tfGenerator,
                                   trackEffUnc, tfVariation, sfupartVariation)
         + "_after_unfolding_2D.root";
}

#endif // RESULT_PATHS_H
