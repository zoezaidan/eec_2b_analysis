// plot_z_first_look.C
//
// First look at a measured observable, straight out of step 1 of the chain -- BEFORE any
// template fit, signal fraction, unfolding or correction. This is the pre-fit step of the
// analysis, the one that answers "is this observable measurable at all".
//
// ⚠️ THE FILE NAME IS HISTORICAL. The momentum balance it was written for is called B, not
// z, since 2026-09-22; the macro keeps its name because ROOT requires the entry function to
// match the file, and renaming both would break every command line in readme_workflow.md
// and every bookmark to the plots it writes. Read "z" in this filename as "the balance".
//
// It takes the observable as its LAST argument:
//
//   root -l -b -q 'plot_z_first_look.C("both","pythia","noeecw_upartv2_B","B")'
//
// The tag must be a production that actually filled that observable's axis -- checked
// against obsProdTag(), not assumed. ⚠️ Pre-rename productions ("_zfirst", "_3obs",
// "_lnfblin", "_fb") name the balance axis "_z" and will not be found under "B".
//
// WHAT THIS IS, AND WHAT IT IS NOT
// --------------------------------
// Nothing here is a measurement. The data curve is the raw EEC-weighted 2b-tagged yield per
// B bin: it still contains the mistagged and single-b background that the template fit
// exists to subtract, it has not been unfolded, and no efficiency, purity or EEC-weight
// correction has been applied. The MC curves are reco-level too, so data and MC are
// compared like for like -- but neither is a particle-level number.
//
// ⚠️ And unlike the dR chain, NO UParT scale factor is applied anywhere here: that
// calibration is binned in reco dR and has no B equivalent (2026-09-15 decision). The dR
// result has a correction that this does not.
//
// The point is to answer, before the template fit is written:
//
//   1. do the flavour templates SEPARATE in B? (plot 1)
//      If the 2b, 1b and 0b B shapes were identical, a fit per B bin would be pointless --
//      the B dependence of the signal fraction would be flat and could be taken from dR.
//   2. does m_2B still separate signal from background INSIDE each B bin? (plot 2)
//      This is the one that decides whether the fit will work: the fit discriminates on
//      m_2B, sliced in (B, pT), so each slice needs both shape separation and statistics.
//      This is the most important plot here.
//   3. is the response well behaved in B? (plot 3)
//      Purity, efficiency and the migration matrix on the chosen 4-bin B binning.
//   4. do data and MC look like each other in B at all? (plot 4)
//
// RUN (on LLR, after source setup_roounfold_env.sh):
//
//   root -l -b -q 'plot_z_first_look.C("qcd","pythia","upartv2_B")'
//
// The third argument is the production tag (OUT_TAG in the run scripts). It defaults to the
// "_upartv2_B" production rather than the nominal, so this macro cannot silently read a
// nominal file that has no balance histograms and report empty plots.

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <iomanip>
#include <iostream>
#include <vector>

#include "observables.h"
#include "result_paths.h"
// TFColor / styleTemplate / styleData: the template fit's own palette and stylers. Included
// rather than copied so the 2B/1B/0B colours here are literally the ones the fit plots use.
#include "Help_Functions.h"   // obsProdTag(): which production carries which observable

namespace {

const char *kBase    = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
const char *kBtagTag = "btagWP0712";

// Analysis palette: red = measurement, blue = particle-level MC, green/purple for the rest.
// Data here is not a measurement yet (no fit, no unfolding), but it IS the data, so it keeps
// red; the MC flavour templates take the "rest" colours.
const int kColData = kRed + 1;
const int kCol2b   = kBlue + 1;    // true 2b -- the signal
const int kCol1b   = kGreen + 2;   // single b
const int kCol0b   = kMagenta + 2; // mistagged light/charm

TString sampleSubdir(const TString &sample, const TString &generator)
{
    const bool herwig = (generator == "herwig");
    if (sample == "qcd")  return herwig ? "QCDHerwig"  : "QCD";
    if (sample == "bjet") return herwig ? "bJetHerwig" : "bJet";
    return "";
}

// Sum one named histogram over every per-block file that exists. Same "take what is on
// disk" rule as aggChunkFiles() in apply_unfolding_2d.C -- block counts differ per sample.
template <class TH>
TH *sumOverBlocks(const std::vector<TString> &files, const char *hname, const char *newname)
{
    TH *out = nullptr;
    int used = 0, missing = 0;
    for (const TString &f : files) {
        TFile *fin = TFile::Open(f);
        if (!fin || fin->IsZombie()) { delete fin; continue; }
        TH *h = dynamic_cast<TH *>(fin->Get(hname));
        if (!h) { missing++; fin->Close(); delete fin; continue; }
        if (!out) { out = (TH *) h->Clone(newname); out->SetDirectory(nullptr); }
        else        out->Add(h);
        used++;
        fin->Close();
        delete fin;
    }
    if (!out)
        std::cerr << "ERROR: '" << hname << "' not found in any of the " << files.size()
                  << " input files. If this is a B histogram, the production predates the "
                  << "B observable (or predates the z -> B rename) -- re-run "
                  << "run_agg_ntuple_chunks.sh." << std::endl;
    else if (missing)
        std::cerr << "WARNING: '" << hname << "' missing from " << missing << " of "
                  << (used + missing) << " files" << std::endl;
    return out;
}

// "both" = qcd + bjet, the nominal configuration. The histogram NAMES are identical in the
// two samples' files (only the filename carries the sample tag), so summing across both
// file lists is all "both" needs -- the same operation apply_unfolding_2d.C performs.
//
// ⚠️ Physics caveat, inherited: "both" adds a b-enriched sample onto inclusive QCD, which
// double-counts b jets unless the sample weights already handle the overlap. The bjet
// sample is there for TEMPLATE STATISTICS. Comparing qcd against both shows how much it
// moves.
std::vector<TString> samplesIn(const TString &sample)
{
    if (sample == "both") return {"qcd", "bjet"};
    return {sample};
}

std::vector<TString> mcFiles(const TString &sample, const TString &generator,
                             const TString &tag, bool rmatrix)
{
    std::vector<TString> files;
    for (const TString &s : samplesIn(sample)) {
        const TString subdir = sampleSubdir(s, generator);
        if (subdir.Length() == 0) return {};
        int found = 0;
        for (int b = 0; b < 50; ++b) {
            TString f = TString::Format(
                "%s/%s/agg_ntuple_chunks/block_%04d/%sRun3_%s_template_for_fit_histos_3D_%s_f%s_%s.root",
                kBase, subdir.Data(), b, rmatrix ? "RMatrix_" : "", kBtagTag, s.Data(),
                rmatrix ? "" : "MCGEN", tag.Data());
            if (gSystem->AccessPathName(f) == kFALSE) { files.push_back(f); ++found; }
        }
        if (found == 0)
            std::cerr << "WARNING: no " << (rmatrix ? "RMatrix" : "template")
                      << " files for sample '" << s << "' under " << subdir
                      << " with tag '" << tag << "'" << std::endl;
    }
    return files;
}

// Data is split over 5 primary datasets x 10 blocks.
std::vector<TString> dataFiles(const TString &tag)
{
    std::vector<TString> files;
    for (int pd = 0; pd < 5; ++pd)
        for (int b = 0; b < 10; ++b) {
            TString f = TString::Format(
                "%s/HardProbes/agg_template_chunks/HardProbes%d/block_%04d/"
                "Run3_%s_template_for_fit_histos_3D_data_fMCGEN_%s.root",
                kBase, pd, b, kBtagTag, tag.Data());
            if (gSystem->AccessPathName(f) == kFALSE) files.push_back(f);
        }
    return files;
}

// ---- jet-pT slicing -------------------------------------------------------------------
// jtpt_binsVector = {80, 100, 120} gives two bins, and they are NOT symmetric:
//   bin 1  80 < pT < 100     the lower bin
//   bin 2  pT > 100          OPEN-ENDED -- jtpt_fill() folds every jet above 120 into it,
//                            so it is "100 to infinity", not "100 to 120"
// The EEC results are quoted per pT bin (its plots carry _pt1/_pt2), so these do too.
// ptbin 0 means integrated over the whole range, 80 -> inf.
//
// There is no ROOT underflow/overflow to worry about on this axis: the pT cut in
// passRecoJetKinematics keeps everything below 80 out, and jtpt_fill() folds the high side
// into bin 2, so bins 1..N hold every filled jet (verified: both are exactly 0).
int ptLo(int ptbin) { return ptbin == 0 ? 1 : ptbin; }
int ptHi(int ptbin, TH1 *h) { return ptbin == 0 ? h->GetNbinsZ() : ptbin; }
// TH2D(observable, pT): the pT axis is Y, not Z.
int ptHiY(int ptbin, TH1 *h) { return ptbin == 0 ? h->GetNbinsY() : ptbin; }

TString ptLabel(int ptbin)
{
    if (ptbin == 0) return Form("p_{T} > %g GeV (all)", jtpt_binsVector[0]);
    // Open-ended last bin -> quote it as a threshold, the same rule pt_label() uses.
    if (ptbin == jtpt_bins) return Form("p_{T} > %g GeV", jtpt_binsVector[ptbin - 1]);
    return Form("%g < p_{T} < %g GeV", jtpt_binsVector[ptbin - 1], jtpt_binsVector[ptbin]);
}

TString ptSuffix(int ptbin) { return ptbin == 0 ? "_ptall" : Form("_pt%d", ptbin); }

// Project a TH3D(m_2B, X, pT) onto its X (observable) axis, over all m_2B, for one pT bin.
TH1D *projX(TH3D *h, const char *name, int ptbin)
{
    if (!h) return nullptr;
    TH1D *p = (TH1D *) h->ProjectionY(name, 1, h->GetNbinsX(), ptLo(ptbin), ptHi(ptbin, h));
    p->SetDirectory(nullptr);
    return p;
}

// Project onto m_2B for one bin of the observable axis, for one pT bin.
TH1D *projMB(TH3D *h, int ix, const char *name, int ptbin)
{
    if (!h) return nullptr;
    TH1D *p = (TH1D *) h->ProjectionX(name, ix, ix, ptLo(ptbin), ptHi(ptbin, h));
    p->SetDirectory(nullptr);
    return p;
}

void styleShape(TH1 *h, int col, int style = 1)
{
    if (!h) return;
    h->SetLineColor(col);
    h->SetLineWidth(2);
    h->SetLineStyle(style);
    h->SetMarkerColor(col);
    h->SetStats(0);
}

void unitArea(TH1 *h)
{
    if (h && h->Integral() > 0) h->Scale(1. / h->Integral());
}

} // anonymous namespace


static void z_first_look_one(const TString &sample, const TString &generator,
                             const TString &tag, int ptbin,
                             const TString &observable)
{
    // A "_noeecw" production was filled with the EEC weight OFF, so h3D_* holds yields
    // (dN/dB), not an energy-weighted correlator. The labels have to say which, or the two
    // sets of plots are indistinguishable once they are out of their folders.
    const bool eecOff  = tag.Contains("noeecw");
    const TString wLab = eecOff ? "yields" : "EEC weighted";
    // Reject a bad sample/generator string rather than silently falling back to a default,
    // the same rule apply_unfolding_2d.C follows.
    if (sample != "qcd" && sample != "bjet" && sample != "both") {
        std::cerr << "ERROR: sample must be qcd | bjet | both (got '" << sample << "')" << std::endl;
        return;
    }
    if (generator != "pythia" && generator != "herwig") {
        std::cerr << "ERROR: generator must be pythia | herwig (got '" << generator << "')" << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Any observable, not just the balance: obsByName() supplies the binning, the axis title and the
    // "_<obs>" histogram-name suffix, so this macro reads whichever axis it is asked for.
    // Its own name is historical -- the balance was the first observable that needed a look before the
    // fit existed, and every other one needs the same four questions answered.
    const ObsDef obsd = obsByName(observable);
    if (obsd.nbins == 0) {
        std::cerr << "ERROR: unknown observable '" << observable << "' (use dr | B)"
                  << std::endl;
        return;
    }
    // The production tag and the observable have to agree, or this draws one observable's
    // histograms out of another observable's production -- which is exactly the mistake
    // obsProdTag() exists to prevent. It returns "_upartv2_3obs"; the tag here is an OUT_TAG
    // with no leading underscore and possibly a "noeecw_" in front, so compare on the stem.
    {
        TString stem = obsProdTag(observable);
        if (stem.BeginsWith("_")) stem.Remove(0, 1);
        if (!tag.Contains(stem)) {
            std::cerr << "ERROR: tag '" << tag << "' is not a production of observable '"
                      << observable << "' (expected it to contain '" << stem << "'). "
                      << "Refusing to read one observable's histograms from another's "
                      << "production." << std::endl;
            return;
        }
    }

    std::cout << "\n================ " << observable << " first look: " << sample << " " << generator
              << ", tag '" << tag << "', " << ptLabel(ptbin)
              << " ================" << std::endl;

    std::vector<TString> f_tmpl = mcFiles(sample, generator, tag, false);
    std::vector<TString> f_rm   = mcFiles(sample, generator, tag, true);
    std::vector<TString> f_data = dataFiles(tag);
    std::cout << "  MC template files: " << f_tmpl.size()
              << ", MC RMatrix files: " << f_rm.size()
              << ", data files: " << f_data.size() << std::endl;
    if (f_tmpl.empty()) {
        std::cerr << "ERROR: no MC template files for tag '" << tag << "'. "
                  << "Run:  OUT_TAG=" << tag << " SAMPLE=" << sample
                  << " GENERATOR=" << generator << " ./run_agg_ntuple_chunks.sh" << std::endl;
        return;
    }

    // ---- MC flavour templates, in B and in dR for reference ----
    TH3D *h3_2b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h3D_bb"), "s_2b_B");
    TH3D *h3_1b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h3D_b"),  "s_1b_B");
    TH3D *h3_0b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h3D_0b"), "s_0b_B");
    if (!h3_2b) return;
    TH3D *h3_2b_dr = sumOverBlocks<TH3D>(f_tmpl, "h3D_bb", "s_2b_dr");

    TH3D *h3_data = f_data.empty() ? nullptr
                                   : sumOverBlocks<TH3D>(f_data, obsd.n("h3D_data"), "s_data_B");
    if (!h3_data)
        std::cout << "  NOTE: no data for this tag yet -- MC-only plots will be made."
                  << std::endl;

    // "<obs>_first_look_...": the folder is named from the observable, so the balance now
    // writes into "B_first_look_..." where it used to write "z_first_look_...". The old
    // folders are left where they are rather than moved.
    TString outdir = TString::Format("%s/results/%s_first_look_%s_%s_%s",
                                     kBase, observable.Data(), sample.Data(),
                                     generator.Data(), tag.Data());
    gSystem->mkdir(outdir, true);

    // ---- 1. Do the flavour templates separate in B? ----
    {
        TH1D *p2b = projX(h3_2b, "p2b", ptbin), *p1b = projX(h3_1b, "p1b", ptbin), *p0b = projX(h3_0b, "p0b", ptbin);
        styleShape(p2b, kCol2b); styleShape(p1b, kCol1b); styleShape(p0b, kCol0b);

        std::cout << "\n--- 1. " << (eecOff ? "YIELD" : "EEC-weighted")
                  << " " << observable << " shapes by truth flavour (unit area) ---" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "    " << observable << " bin      |    2b    |    1b    |    0b" << std::endl;
        std::cout << "    -------------+----------+----------+---------" << std::endl;
        TH1D *n2b = (TH1D *) p2b->Clone("n2b"); unitArea(n2b);
        TH1D *n1b = p1b ? (TH1D *) p1b->Clone("n1b") : nullptr; unitArea(n1b);
        TH1D *n0b = p0b ? (TH1D *) p0b->Clone("n0b") : nullptr; unitArea(n0b);
        for (int i = 1; i <= n2b->GetNbinsX(); ++i)
            std::cout << "    " << std::setw(5) << n2b->GetXaxis()->GetBinLowEdge(i)
                      << " - " << std::setw(5) << n2b->GetXaxis()->GetBinUpEdge(i)
                      << " | " << std::setw(8) << n2b->GetBinContent(i)
                      << " | " << std::setw(8) << (n1b ? n1b->GetBinContent(i) : 0.)
                      << " | " << std::setw(8) << (n0b ? n0b->GetBinContent(i) : 0.) << std::endl;

        TCanvas c("c1", "", 900, 700);
        n2b->SetTitle(Form(";%s;normalised", obsd.axis.Data()));
        double mx = n2b->GetMaximum();
        if (n1b) mx = std::max(mx, n1b->GetMaximum());
        if (n0b) mx = std::max(mx, n0b->GetMaximum());
        n2b->SetMaximum(1.45 * mx);
        n2b->SetMinimum(0.);
        n2b->Draw("hist");
        if (n1b) n1b->Draw("hist same");
        if (n0b) n0b->Draw("hist same");
        TLegend leg(0.16, 0.72, 0.60, 0.89);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(n2b, "true 2b (signal)", "l");
        if (n1b) leg.AddEntry(n1b, "1b", "l");
        if (n0b) leg.AddEntry(n0b, "0b (mistag)", "l");
        leg.Draw();
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.032);
        lat.DrawLatex(0.16, 0.93, Form("%s %s, %s, reco level, %s",
                                       sample.Data(), generator.Data(), ptLabel(ptbin).Data(), wLab.Data()));
        c.SaveAs(outdir + "/1_" + observable + "_templates_by_flavour" + ptSuffix(ptbin) + ".pdf");
        c.SaveAs(outdir + "/1_" + observable + "_templates_by_flavour" + ptSuffix(ptbin) + ".png");
    }

    // ---- 2. m_2B inside each B bin: what the template fit will actually fit ----
    {
        const int nz = h3_2b->GetNbinsY();
        TCanvas c("c2", "", 1300, 380);
        c.Divide(nz, 1, 0.001, 0.001);
        std::cout << Form("\n--- 2. m_2B separation inside each %s bin ---", observable.Data()) << std::endl;
        std::cout << "    " << observable << " bin      |   N(2b)   |   N(1b)   |   N(0b)   | S/(S+B)" << std::endl;
        std::cout << "    -------------+-----------+-----------+-----------+--------" << std::endl;
        // Keep the projections alive until the canvas is written.
        std::vector<TH1D *> keep;
        for (int iz = 1; iz <= nz; ++iz) {
            TH1D *m2b = projMB(h3_2b, iz, Form("m2b_%d", iz), ptbin);
            TH1D *m1b = projMB(h3_1b, iz, Form("m1b_%d", iz), ptbin);
            TH1D *m0b = projMB(h3_0b, iz, Form("m0b_%d", iz), ptbin);
            const double n2 = m2b ? m2b->Integral() : 0.;
            const double n1 = m1b ? m1b->Integral() : 0.;
            const double n0 = m0b ? m0b->Integral() : 0.;
            std::cout << "    " << std::setw(5) << h3_2b->GetYaxis()->GetBinLowEdge(iz)
                      << " - " << std::setw(5) << h3_2b->GetYaxis()->GetBinUpEdge(iz)
                      << " | " << std::setw(9) << n2 << " | " << std::setw(9) << n1
                      << " | " << std::setw(9) << n0
                      << " | " << std::setw(7)
                      << ((n2 + n1 + n0) > 0 ? n2 / (n2 + n1 + n0) : 0.) << std::endl;

            styleShape(m2b, kCol2b); styleShape(m1b, kCol1b); styleShape(m0b, kCol0b);
            unitArea(m2b); unitArea(m1b); unitArea(m0b);
            c.cd(iz);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14);
            double mx = m2b ? m2b->GetMaximum() : 1.;
            if (m1b) mx = std::max(mx, m1b->GetMaximum());
            if (m0b) mx = std::max(mx, m0b->GetMaximum());
            if (m2b) {
                m2b->SetTitle(";m_{2B} [GeV];normalised");
                m2b->SetMaximum(1.5 * mx); m2b->SetMinimum(0.);
                m2b->GetXaxis()->SetTitleSize(0.055); m2b->GetYaxis()->SetTitleSize(0.055);
                m2b->Draw("hist");
            }
            if (m1b) m1b->Draw("hist same");
            if (m0b) m0b->Draw("hist same");
            TLatex *lt = new TLatex(); lt->SetNDC(); lt->SetTextSize(0.055);
            // Name the observable, not "B": these panels are the fit's input and end up in
            // talks on their own, where "0.30 < B < 0.50" on a dr plot would be simply wrong.
            lt->DrawLatex(0.20, 0.90, Form("%.2f < %s < %.2f",
                                           h3_2b->GetYaxis()->GetBinLowEdge(iz),
                                           obsSymbol(observable).Data(),
                                           h3_2b->GetYaxis()->GetBinUpEdge(iz)));
            if (iz == 1) {
                TLegend *lg = new TLegend(0.45, 0.62, 0.92, 0.86);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.05);
                lg->AddEntry(m2b, "2b", "l");
                if (m1b) lg->AddEntry(m1b, "1b", "l");
                if (m0b) lg->AddEntry(m0b, "0b", "l");
                lg->Draw();
            }
            keep.push_back(m2b); keep.push_back(m1b); keep.push_back(m0b);
        }
        c.SaveAs(outdir + "/2_mB_in_" + observable + "_bins" + ptSuffix(ptbin) + ".pdf");
        c.SaveAs(outdir + "/2_mB_in_" + observable + "_bins" + ptSuffix(ptbin) + ".png");
    }

    // ---- 3. Response quality in B: purity, efficiency, migration ----
    if (!f_rm.empty()) {
        TH2D *pur_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_purity_numerator_tf"), "pn_B");
        TH2D *pud_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_purity_denominator_tf"), "pd_B");
        TH2D *efn_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_efficiency_numerator_tf"), "en_B");
        TH2D *efd_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_efficiency_denominator_tf"), "ed_B");
        // dR reference, same quantities on its own binning.
        TH2D *pur_d = sumOverBlocks<TH2D>(f_rm, "h_full_purity_numerator_tf", "pn_d");
        TH2D *pud_d = sumOverBlocks<TH2D>(f_rm, "h_full_purity_denominator_tf", "pd_d");

        // ⚠️ Ratios are RECOMPUTED from the summed numerators and denominators, never
        // summed themselves -- adding N stored ratios gives N x the true value. Same rule
        // as apply_unfolding_2d.C's "How 'both' is summed".
        auto ratio1D = [&](TH2D *num, TH2D *den, const char *name) -> TH1D * {
            if (!num || !den) return nullptr;
            TH1D *n = num->ProjectionX(Form("%s_n", name), ptLo(ptbin), ptHiY(ptbin, num));
            TH1D *d = den->ProjectionX(Form("%s_d", name), ptLo(ptbin), ptHiY(ptbin, den));
            TH1D *r = (TH1D *) n->Clone(name);
            r->Divide(n, d, 1., 1., "b");
            r->SetDirectory(nullptr);
            return r;
        };
        TH1D *purity_B = ratio1D(pur_B, pud_B, "purity_B");
        TH1D *eff_B    = ratio1D(efn_B, efd_B, "eff_B");
        TH1D *purity_d = ratio1D(pur_d, pud_d, "purity_d");

        std::cout << Form("\n--- 3. Response quality in %s (dR's own values for scale) ---", observable.Data()) << std::endl;
        if (purity_B) {
            std::cout << "    " << observable << " bin      |  purity  | efficiency" << std::endl;
            std::cout << "    -------------+----------+-----------" << std::endl;
            for (int i = 1; i <= purity_B->GetNbinsX(); ++i)
                std::cout << "    " << std::setw(5) << purity_B->GetXaxis()->GetBinLowEdge(i)
                          << " - " << std::setw(5) << purity_B->GetXaxis()->GetBinUpEdge(i)
                          << " | " << std::setw(8) << purity_B->GetBinContent(i)
                          << " | " << std::setw(10) << (eff_B ? eff_B->GetBinContent(i) : 0.)
                          << std::endl;
        }
        if (purity_d) {
            std::cout << "    dR purity, for scale:";
            for (int i = 1; i <= purity_d->GetNbinsX(); ++i)
                std::cout << " " << std::setprecision(3) << purity_d->GetBinContent(i);
            std::cout << std::endl;
        }

        if (purity_B) {
            TCanvas c("c3", "", 900, 700);
            styleShape(purity_B, kCol2b);
            styleShape(eff_B, kCol1b);
            purity_B->SetTitle(Form(";%s;fraction", obsd.axis.Data()));
            purity_B->SetMinimum(0.); purity_B->SetMaximum(1.35);
            purity_B->Draw("hist");
            if (eff_B) eff_B->Draw("hist same");
            TLegend leg(0.16, 0.74, 0.62, 0.89);
            leg.SetBorderSize(0); leg.SetFillStyle(0);
            leg.AddEntry(purity_B, "purity (reco-binned)", "l");
            if (eff_B) leg.AddEntry(eff_B, "efficiency (gen-binned)", "l");
            leg.Draw();
            TLine l(obsd.bins[0], 1., obsd.bins[obsd.nbins], 1.);
            l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.Draw();
            c.SaveAs(outdir + "/3_" + observable + "_purity_efficiency" + ptSuffix(ptbin) + ".pdf");
            c.SaveAs(outdir + "/3_" + observable + "_purity_efficiency" + ptSuffix(ptbin) + ".png");
        }
    }

    // ---- 4. Data vs MC in B, reco level, shapes only ----
    if (h3_data) {
        TH1D *pdata = projX(h3_data, "pdata", ptbin);
        TH1D *pmc   = projX(h3_2b, "pmc_sum", ptbin);
        if (h3_1b) pmc->Add(projX(h3_1b, "pmc_1b", ptbin));
        if (h3_0b) pmc->Add(projX(h3_0b, "pmc_0b", ptbin));
        styleShape(pdata, kColData); styleShape(pmc, kCol2b);
        pdata->SetMarkerStyle(20);

        std::cout << Form("\n--- 4. Data vs MC in %s (reco level, unit area, NO fit/corrections) ---", observable.Data())
                  << std::endl;
        std::cout << "    " << observable << " bin      |    data  |    MC    |  data/MC" << std::endl;
        std::cout << "    -------------+----------+----------+---------" << std::endl;
        unitArea(pdata); unitArea(pmc);
        for (int i = 1; i <= pdata->GetNbinsX(); ++i) {
            const double d = pdata->GetBinContent(i), m = pmc->GetBinContent(i);
            std::cout << "    " << std::setw(5) << pdata->GetXaxis()->GetBinLowEdge(i)
                      << " - " << std::setw(5) << pdata->GetXaxis()->GetBinUpEdge(i)
                      << " | " << std::setw(8) << d << " | " << std::setw(8) << m
                      << " | " << std::setw(8) << (m > 0 ? d / m : 0.) << std::endl;
        }

        TCanvas c("c4", "", 900, 700);
        pmc->SetTitle(Form(";%s;normalised", obsd.axis.Data()));
        pmc->SetMinimum(0.);
        pmc->SetMaximum(1.45 * std::max(pmc->GetMaximum(), pdata->GetMaximum()));
        pmc->Draw("hist");
        pdata->Draw("e same");
        TLegend leg(0.16, 0.74, 0.65, 0.89);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(pdata, "data (2b-tagged, pre-fit)", "lep");
        leg.AddEntry(pmc, Form("%s MC, all flavours", sample.Data()), "l");
        leg.Draw();
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030);
        lat.DrawLatex(0.16, 0.93, Form("%s -- reco level, %s, no fit/unfolding/UParT SF", ptLabel(ptbin).Data(), wLab.Data()));
        c.SaveAs(outdir + "/4_" + observable + "_data_vs_mc" + ptSuffix(ptbin) + ".pdf");
        c.SaveAs(outdir + "/4_" + observable + "_data_vs_mc" + ptSuffix(ptbin) + ".png");
    }

    // ---- 5. RAW YIELDS per B bin ----------------------------------------------------
    // The h_count_* histograms are the yields: same jets, same selection, same binning as
    // the h3D_* templates, but filled with the tree weight ONLY -- no EEC weight. That is
    // the difference between "how many jets" and "how much EEC", and it is why none of the
    // plots above is a yield (they are all EEC-weighted or unit-area normalised).
    //
    // For DATA the tree weight is 1, so h_count_data IS the raw jet count and its error is
    // sqrt(N). ⚠️ The trigger prescale is NOT in it -- the prescale multiplies `eec`, not
    // the count -- so this is the number of jets actually recorded, not a luminosity-scaled
    // yield. For MC the tree weight is the sample weight, so the count is a weighted
    // prediction; GetEntries() is printed next to it as the raw statistics behind it.
    {
        std::cout << Form("\n--- 5. RAW YIELDS per %s bin (h_count_*, NO EEC weight) ---", observable.Data()) << std::endl;

        TH3D *c_data = f_data.empty() ? nullptr
                                      : sumOverBlocks<TH3D>(f_data, obsd.n("h_count_data"), "c_data_B");
        TH3D *c_2b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h_count_bb"), "c_2b_B");
        TH3D *c_1b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h_count_b"),  "c_1b_B");
        TH3D *c_0b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h_count_0b"), "c_0b_B");

        TH1D *y_data = projX(c_data, "y_data", ptbin);
        TH1D *y_2b   = projX(c_2b,   "y_2b", ptbin);
        TH1D *y_1b   = projX(c_1b,   "y_1b", ptbin);
        TH1D *y_0b   = projX(c_0b,   "y_0b", ptbin);

        std::cout << std::fixed;
        if (y_data) {
            std::cout << "\n  DATA -- raw b-tagged 2-SV jet counts (prescale NOT applied):\n";
            std::cout << "    " << observable << " bin      |      jets  |  sqrt(N)  | rel. stat" << std::endl;
            std::cout << "    -------------+------------+-----------+----------" << std::endl;
            double tot = 0.;
            for (int i = 1; i <= y_data->GetNbinsX(); ++i) {
                const double N = y_data->GetBinContent(i);
                tot += N;
                std::cout << "    " << std::setprecision(3)
                          << std::setw(5) << y_data->GetXaxis()->GetBinLowEdge(i)
                          << " - " << std::setw(5) << y_data->GetXaxis()->GetBinUpEdge(i)
                          << " | " << std::setw(10) << std::setprecision(0) << N
                          << " | " << std::setw(9) << std::setprecision(1) << std::sqrt(N)
                          << " | " << std::setw(8) << std::setprecision(4)
                          << (N > 0 ? std::sqrt(N) / N : 0.) << std::endl;
            }
            std::cout << "    total: " << std::setprecision(0) << tot << " jets" << std::endl;
        }

        if (y_2b) {
            std::cout << "\n  MC (" << sample << " " << generator
                      << ") -- weighted counts, with raw entries behind them:\n";
            std::cout << "    " << observable << " bin      |   2b (wt) |   1b (wt) |   0b (wt) | S/(S+B) | 2b eff. entries"
                      << std::endl;
            std::cout << "    -------------+-----------+-----------+-----------+---------+----------------"
                      << std::endl;
            // Effective entries per bin, (sum w)^2 / sum w^2, from the content and its
            // error. This is the statistical power actually behind a weighted bin --
            // GetEntries() on the projection would give the TOTAL over all bins, the same
            // number four times over, which is why it is not used here.
            auto effEntries = [](TH1D *h, int i) -> double {
                const double c = h->GetBinContent(i), e = h->GetBinError(i);
                return (e > 0.) ? (c * c) / (e * e) : 0.;
            };
            for (int i = 1; i <= y_2b->GetNbinsX(); ++i) {
                const double n2 = y_2b->GetBinContent(i);
                const double n1 = y_1b ? y_1b->GetBinContent(i) : 0.;
                const double n0 = y_0b ? y_0b->GetBinContent(i) : 0.;
                std::cout << "    " << std::setprecision(3)
                          << std::setw(5) << y_2b->GetXaxis()->GetBinLowEdge(i)
                          << " - " << std::setw(5) << y_2b->GetXaxis()->GetBinUpEdge(i)
                          << " | " << std::setw(9) << std::setprecision(4) << n2
                          << " | " << std::setw(9) << n1
                          << " | " << std::setw(9) << n0
                          << " | " << std::setw(7)
                          << ((n2 + n1 + n0) > 0 ? n2 / (n2 + n1 + n0) : 0.)
                          << " | " << std::setw(15) << std::setprecision(0)
                          << effEntries(y_2b, i) << std::endl;
            }
            std::cout << "\n  ⚠️ MC counts are WEIGHTED (sample weight), so they are a prediction,"
                      << " not a number of\n     simulated jets, and they are NOT normalised to the"
                      << " data luminosity -- do not\n     compare the MC and data columns directly."
                      << std::endl;
        }

        // Yields plot: the three MC flavour components stacked, with data on top.
        //
        // Colours and styles come from Help_Functions.h -- TFColor::c2B/c1B/c0B and
        // styleTemplate/styleData -- i.e. exactly what the template fit draws its own plots
        // with: 2B red, 1B blue, 0B green, data black closed circles. Someone flipping
        // between this plot and the fit's m_2B plots is looking at the same categories, so
        // they have to be the same colours.
        //
        // Stack order matches the fit's: 2B at the bottom, then 1B, 0B on top.
        //
        // ⚠️ MC IS SCALED TO THE DATA. The MC counts are weighted predictions and are NOT
        // normalised to the data luminosity, so the raw MC and data numbers are not
        // comparable (the printout above says so too). Scaling the SUM to the data integral,
        // with the three components keeping their relative proportions, is the pre-fit
        // normalisation the template fit itself starts from -- what the plot then shows is
        // whether the SHAPES agree, which is the question worth asking before the fit.
        if (y_data) {
            const Float_t font_scale  = 1200. / 800.;
            const Style_t font_code   = 43;
            const Float_t label_size  = 15. * font_scale;
            const Float_t title_size  = 15. * font_scale;
            const Float_t legend_size = 14. * font_scale;

            const double n_jets = y_data->Integral();

            // Scale the MC sum to the data, preserving the flavour proportions.
            double mc_sum = 0.;
            for (TH1D *h : {y_2b, y_1b, y_0b}) if (h) mc_sum += h->Integral();
            const double mc_scale = (mc_sum > 0.) ? n_jets / mc_sum : 1.;
            for (TH1D *h : {y_2b, y_1b, y_0b}) if (h) h->Scale(mc_scale);

            std::cout << "    MC scaled to data by x" << mc_scale
                      << " (sum of 2b+1b+0b normalised to the data integral)" << std::endl;

            TCanvas c("c5", "", 800, 700);
            c.SetTicks(1, 1);
            c.SetLeftMargin(0.15); c.SetRightMargin(0.05);
            c.SetTopMargin(0.09);  c.SetBottomMargin(0.15);

            if (y_2b) styleTemplate(y_2b, TFColor::c2B());
            if (y_1b) styleTemplate(y_1b, TFColor::c1B());
            if (y_0b) styleTemplate(y_0b, TFColor::c0B());
            styleData(y_data);

            THStack st("st_raw_yields", "");
            if (y_2b) st.Add(y_2b);      // bottom
            if (y_1b) st.Add(y_1b);
            if (y_0b) st.Add(y_0b);      // top

            // Frame from the data histogram, so the axis is the observable's own binning and
            // the range covers whichever of data/stack is taller.
            TH1D *frame = (TH1D *) y_data->Clone("h_raw_frame");
            frame->SetDirectory(nullptr);
            frame->Reset();
            frame->SetTitle("");
            frame->SetMinimum(0.);
            frame->SetMaximum(1.45 * std::max(y_data->GetMaximum(), st.GetMaximum()));
            frame->GetYaxis()->SetTitle("2b jets per bin");
            frame->GetXaxis()->SetTitle(obsd.axis);
            for (TAxis *ax : { frame->GetXaxis(), frame->GetYaxis() }) {
                ax->CenterTitle(true);
                ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
                ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
            }
            frame->GetXaxis()->SetTitleOffset(1.15);
            frame->GetYaxis()->SetTitleOffset(1.9);
            frame->GetYaxis()->SetNdivisions(505);
            // ROOT pads an axis it was not given an explicit range for -- B would be drawn
            // out to 1.1, past where the observable can go. Same fix as the other macros.
            if (observable != "dr") frame->GetXaxis()->SetRange(1, frame->GetNbinsX());
            frame->Draw("AXIS");
            st.Draw("hist same");
            y_data->Draw("PE X0 same");

            TLatex cms; cms.SetNDC();
            cms.SetTextFont(62); cms.SetTextSize(0.045); cms.DrawLatex(0.15, 0.945, "CMS");
            cms.SetTextFont(52); cms.SetTextSize(0.036); cms.DrawLatex(0.255, 0.945, "Internal");
            cms.SetTextFont(42); cms.SetTextSize(0.036); cms.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

            // Two lines only, by request: which pT bin, and what the plot is. The jet count,
            // the MC scale factor and the "no fit, no corrections" caveat were on the plot
            // and are now only in the printout above -- they have not stopped being true.
            // ⚠️ In particular the MC here IS scaled to the data (mc_scale, printed below),
            // so nothing on the canvas says so any more. Say it in the caption.
            TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.032);
            note.SetTextAlign(13);
            double y_note = 0.87;
            note.DrawLatex(0.19, y_note, ptLabel(ptbin));
            note.DrawLatex(0.19, y_note -= 0.048, "Raw yields");

            TLegend leg(0.62, 0.66, 0.94, 0.87);
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetMargin(0.25);
            leg.SetTextFont(font_code); leg.SetTextSize(legend_size);
            leg.AddEntry(y_data, "Data", "pe1");
            if (y_2b) leg.AddEntry(y_2b, "MC 2B (signal)", "f");
            if (y_1b) leg.AddEntry(y_1b, "MC 1B", "f");
            if (y_0b) leg.AddEntry(y_0b, "MC 0B (mistag)", "f");
            leg.Draw();
            c.RedrawAxis();

            c.SaveAs(outdir + "/5_" + observable + "_raw_yields" + ptSuffix(ptbin) + ".pdf");
            c.SaveAs(outdir + "/5_" + observable + "_raw_yields" + ptSuffix(ptbin) + ".png");
        }
    }

    std::cout << "\nWritten to " << outdir << std::endl;
    std::cout << "(/data_CMS is not on the sshfs mount -- use sync_plots_to_eos.sh or scp)"
              << std::endl;
}

// Entry point: one pass per jet-pT slice.
//
// The EEC results are quoted per pT bin, so these are too. Two passes:
//   _pt1    80 < pT < 100   the lower bin
//   _pt2    pT > 100        the main bin, OPEN-ENDED (jtpt_fill folds >120 into it)
//
// The pT-INTEGRATED pass is deliberately NOT produced. The two pT bins have visibly
// different signal fractions (0.33-0.37 against 0.47-0.48 from the B template fit), so a
// pT-integrated shape is a blend of two populations rather than a measurement of either,
// and having it sitting in the same folder as the real ones only invites it being read as
// a result. z_first_look_one() still accepts ptbin = 0 for it if it is ever wanted for a
// statistics check -- it is just never called that way:
//
//   /* ---- disabled (kept for reference): the pT-integrated pass ----
//      z_first_look_one(sample, generator, tag, 0);   // writes *_ptall.*
//      ---- */
void plot_z_first_look(const TString &sample    = "qcd",
                       const TString &generator = "pythia",
                       const TString &tag       = "upartv2_B",
                       const TString &observable = "B")
{
    for (int ptbin = 1; ptbin <= jtpt_bins; ++ptbin)
        z_first_look_one(sample, generator, tag, ptbin, observable);
}
