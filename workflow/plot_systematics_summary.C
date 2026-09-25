// Relative shift of every Herwig variation with respect to the Pythia8 nominal.
//
//   each correction / step swapped to Herwig ON ITS OWN   (thin lines)
//   the whole workflow run with Herwig                    (thick red)
//
// All curves are (variation - nominal) / nominal, so a flat zero means the swap does not
// move the normalised EEC at all.
//
// Inputs are the two final_*_with_systematics.root files written by
// apply_weights_and_systematics.C -- run it for the Pythia and the Herwig nominal first.
//
// usage: root -l -b -q 'plot_systematics_summary.C("both")'

#include "result_paths.h"   // the ONE definition of every result path -- these two plots read
                            // the final files apply_weights_and_systematics.C wrote, so the
                            // paths must come from the same place. They used to be spelled
                            // out by hand here, which is how they would have missed the
                            // matrix_inversion/ subdirectory.
#include <vector>

// The two final files both plots read: the Pythia nominal and the all-Herwig run.
// test_mode 2 (data) -- a generator comparison of anything else is not a measurement.
inline TString finalBandFile(const TString &sample, const TString &generator,
                             const TString &tfGenerator, bool unfoldBayes)
{
    return resultFolder(sample, generator, unfoldBayes, tfGenerator)
         + "final_" + resultLabel(sample, generator, 2, unfoldBayes, tfGenerator)
         + "_with_systematics.root";
}

void plot_systematics_summary(TString SAMPLE = "both", bool unfoldBayes = false)
{
  // Both plots compare two full productions, so they belong with the production they were
  // made from: methodDir keeps the matrix-inversion pair out of the Bayesian one's way
  // instead of overwriting it under the same name.
  const TString R   = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/"
                    + methodDir(unfoldBayes);

  const TString f_pyt = finalBandFile(SAMPLE, "pythia", "pythia", unfoldBayes);
  const TString f_her = finalBandFile(SAMPLE, "herwig", "herwig", unfoldBayes);

  TFile *fp = TFile::Open(f_pyt);
  TFile *fh = TFile::Open(f_her);
  if (!fp || fp->IsZombie()) { std::cerr << "MISSING " << f_pyt << std::endl; return; }
  if (!fh || fh->IsZombie()) { std::cerr << "MISSING " << f_her << std::endl; return; }

  TH1D *nom    = (TH1D *) fp->Get("h_nominal");
  TH1D *nom_hw = (TH1D *) fh->Get("h_nominal");
  if (!nom || !nom_hw) { std::cerr << "ERROR: h_nominal not found" << std::endl; return; }

  // h_delta_* is (nominal - variation); flip the sign so every curve reads as the shift
  // the Herwig swap produces.
  struct Src { const char *hist; const char *label; const char *hex; int style; };
  std::vector<Src> srcs = {
    { "h_delta_template_fit",    "template fit",        "#4C72B0", 2 },
    { "h_delta_unfolding_model", "unfolding MC",        "#4F8F52", 3 },
    { "h_delta_eff_2sv_btag",    "2SV+UParT efficiency","#8C6BB1", 4 },
    { "h_delta_eec_weight",      "EEC weight",          "#C97430", 5 },
  };

  const int nb = nom->GetNbinsX();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c = new TCanvas("c_syst_summary", "", 800, 600);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.04);
  c->SetTopMargin(0.09);  c->SetBottomMargin(0.14);
  c->SetTicks(1, 1);

  TH1D *frame = (TH1D *) nom->Clone("frame_summary");
  frame->SetDirectory(nullptr);
  frame->Reset();

  // Coherent: the whole workflow run with Herwig.
  TH1D *coh = (TH1D *) nom->Clone("h_rel_all_herwig");
  coh->SetDirectory(nullptr);
  double rmax = 0.;
  for (int i = 1; i <= nb; ++i) {
    const double v = nom->GetBinContent(i);
    const double r = (v != 0.) ? (nom_hw->GetBinContent(i) - v) / v : 0.;
    coh->SetBinContent(i, r); coh->SetBinError(i, 0.);
    rmax = std::max(rmax, std::fabs(r));
  }

  std::vector<TH1D *> curves;
  for (const Src &s : srcs) {
    TH1D *d = (TH1D *) fp->Get(s.hist);
    if (!d) { std::cerr << "WARNING: " << s.hist << " missing" << std::endl; continue; }
    TH1D *rel = (TH1D *) nom->Clone(Form("h_rel_%s", s.label));
    rel->SetDirectory(nullptr);
    for (int i = 1; i <= nb; ++i) {
      const double v = nom->GetBinContent(i);
      const double r = (v != 0.) ? -d->GetBinContent(i) / v : 0.; // flip sign
      rel->SetBinContent(i, r); rel->SetBinError(i, 0.);
      rmax = std::max(rmax, std::fabs(r));
    }
    rel->SetLineColor((Color_t) TColor::GetColor(s.hex));
    rel->SetLineWidth(2);
    rel->SetLineStyle(s.style);
    curves.push_back(rel);
  }

  frame->GetYaxis()->SetRangeUser(-1.45 * rmax, 1.75 * rmax);
  frame->GetXaxis()->SetTitle("gen #Delta r_{BB}");
  frame->GetYaxis()->SetTitle("(Herwig #minus Pythia8) / Pythia8");
  frame->GetXaxis()->SetTitleFont(43); frame->GetXaxis()->SetTitleSize(23);
  frame->GetYaxis()->SetTitleFont(43); frame->GetYaxis()->SetTitleSize(23);
  frame->GetXaxis()->SetLabelFont(43); frame->GetXaxis()->SetLabelSize(19);
  frame->GetYaxis()->SetLabelFont(43); frame->GetYaxis()->SetLabelSize(19);
  frame->GetXaxis()->SetTitleOffset(1.3);
  frame->GetYaxis()->SetTitleOffset(1.6);
  frame->GetYaxis()->SetNdivisions(505);
  frame->Draw("AXIS");

  TLine *l0 = new TLine(frame->GetXaxis()->GetXmin(), 0., frame->GetXaxis()->GetXmax(), 0.);
  l0->SetLineColor(kGray + 1); l0->SetLineStyle(2);
  l0->Draw();

  for (TH1D *h : curves) h->Draw("HIST same");

  coh->SetLineColor((Color_t) TColor::GetColor("#C44E52"));
  coh->SetLineWidth(4);
  coh->Draw("HIST same");

  TLegend *leg = new TLegend(0.16, 0.70, 0.62, 0.90);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(43); leg->SetTextSize(19);
  leg->AddEntry(coh, "whole workflow (Herwig)", "l");
  for (size_t k = 0; k < curves.size() && k < srcs.size(); ++k)
    leg->AddEntry(curves[k], srcs[k].label, "l");
  leg->Draw();
  c->RedrawAxis();

  const TString stem = R + "systematics_summary_" + SAMPLE;
  c->Print(stem + ".pdf");
  c->Print(stem + ".png");

  printf("\n %-24s", "bin");
  for (int i = 1; i <= nb; ++i) printf(" %6d", i);
  printf("\n %-24s", "whole workflow [%]");
  for (int i = 1; i <= nb; ++i) printf(" %6.2f", 100 * coh->GetBinContent(i));
  for (size_t k = 0; k < curves.size(); ++k) {
    printf("\n %-24s", srcs[k].label);
    for (int i = 1; i <= nb; ++i) printf(" %6.2f", 100 * curves[k]->GetBinContent(i));
  }
  printf("\n\nWrote %s.{pdf,png}\n", stem.Data());
}


// =====================================================================================
// The final result from the two full analyses -- everything Pythia8, everything Herwig --
// with the generator systematic drawn as a symmetric band, +/- |Herwig - Pythia8|,
// around the Pythia8 nominal.
//
// usage:  root -l -b -q 'plot_systematics_summary.C' -e 'plot_final_generator_band("both")'

Color_t blendW(const char *hex, double a)
{
  TColor *c = gROOT->GetColor(TColor::GetColor(hex));
  if (!c) return kWhite;
  Float_t r, g, b; c->GetRGB(r, g, b);
  return (Color_t) TColor::GetColor(Float_t(1. - a + a * r),
                                    Float_t(1. - a + a * g),
                                    Float_t(1. - a + a * b));
}

void plot_final_generator_band(TString SAMPLE = "both", bool unfoldBayes = false)
{
  // Both plots compare two full productions, so they belong with the production they were
  // made from: methodDir keeps the matrix-inversion pair out of the Bayesian one's way
  // instead of overwriting it under the same name.
  const TString R   = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/"
                    + methodDir(unfoldBayes);

  const TString f_pyt = finalBandFile(SAMPLE, "pythia", "pythia", unfoldBayes);
  const TString f_her = finalBandFile(SAMPLE, "herwig", "herwig", unfoldBayes);

  TFile *fp = TFile::Open(f_pyt);
  TFile *fh = TFile::Open(f_her);
  if (!fp || fp->IsZombie()) { std::cerr << "MISSING " << f_pyt << std::endl; return; }
  if (!fh || fh->IsZombie()) { std::cerr << "MISSING " << f_her << std::endl; return; }

  TH1D *py = (TH1D *) fp->Get("h_nominal");
  TH1D *hw = (TH1D *) fh->Get("h_nominal");
  if (!py || !hw) { std::cerr << "ERROR: h_nominal missing" << std::endl; return; }
  py = (TH1D *) py->Clone("h_full_pythia"); py->SetDirectory(nullptr);
  hw = (TH1D *) hw->Clone("h_full_herwig"); hw->SetDirectory(nullptr);

  const int nb = py->GetNbinsX();

  // ---- identical style to the systematics-curves plot in apply_weights_and_systematics.C
  const Color_t col_nom  = (Color_t) TColor::GetColor("#C44E52"); // measurement red
  const Color_t col_alt  = (Color_t) TColor::GetColor("#4C72B0"); // blue
  const Color_t col_band = blendW("#C44E52", 0.30);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

  const Float_t font_scale  = 1200. / 800.;
  const Style_t font_code   = 43;
  const Float_t label_size  = 15. * font_scale;
  const Float_t title_size  = 15. * font_scale;
  const Float_t legend_size = 15. * font_scale;

  TCanvas *c = new TCanvas("c_gen_band", "", 800, 800);
  TPad *p_top = new TPad("p_top_gb", "", 0., 0.3, 1., 1.);
  TPad *p_bot = new TPad("p_bot_gb", "", 0., 0.,  1., 0.3);
  for (TPad *p : {p_top, p_bot}) {
    p->SetTicks(1, 0);
    p->SetFillColor(0);
  }
  p_top->SetMargin(0.1, 0.1, 0.0,  0.1);   // left, right, bottom, top
  p_bot->SetMargin(0.1, 0.1, 0.23, 0.0);
  c->cd();
  p_top->Draw(); p_bot->Draw();

  // Band: +/- |Herwig - Pythia8| around the Pythia8 result, and its relative version.
  TGraphErrors *g_band     = new TGraphErrors(nb);
  TGraphErrors *g_band_rel = new TGraphErrors(nb);
  g_band->SetName("g_generator_band");
  g_band_rel->SetName("g_generator_band_rel");
  double ymax = 0., rmax = 0.;
  for (int i = 1; i <= nb; ++i) {
    const double v = py->GetBinContent(i);
    const double d = std::fabs(hw->GetBinContent(i) - v);
    g_band->SetPoint(i - 1, py->GetBinCenter(i), v);
    g_band->SetPointError(i - 1, 0.5 * py->GetBinWidth(i), d);
    g_band_rel->SetPoint(i - 1, py->GetBinCenter(i), 0.);
    g_band_rel->SetPointError(i - 1, 0.5 * py->GetBinWidth(i), (v != 0.) ? d / v : 0.);
    ymax = std::max(ymax, v + py->GetBinError(i) + d);
    ymax = std::max(ymax, hw->GetBinContent(i) + hw->GetBinError(i));
    if (v != 0.) rmax = std::max(rmax, d / v + hw->GetBinError(i) / v);
  }
  g_band->SetFillColor(col_band);     g_band->SetLineWidth(0);
  g_band_rel->SetFillColor(col_band); g_band_rel->SetLineWidth(0);

  // ---- top ----
  p_top->cd();
  py->SetTitle("");
  py->GetYaxis()->SetRangeUser(0., ymax * 1.75);
  py->GetYaxis()->SetTitle("EEC(#Delta r)");
  py->GetYaxis()->CenterTitle(true);
  py->GetYaxis()->SetTitleFont(font_code); py->GetYaxis()->SetTitleSize(title_size);
  py->GetYaxis()->SetTitleOffset(1.5);
  py->GetYaxis()->SetLabelFont(font_code); py->GetYaxis()->SetLabelSize(label_size);
  py->GetXaxis()->SetTitleSize(0);
  py->GetXaxis()->SetLabelSize(0);
  py->SetLineColor(col_nom); py->SetMarkerColor(col_nom);
  py->SetMarkerStyle(kFullCircle); py->SetMarkerSize(1); py->SetLineWidth(2);
  // Herwig as a dashed line, like the variation curves on the systematics plots.
  hw->SetLineColor(col_alt); hw->SetMarkerColor(col_alt);
  hw->SetLineStyle(2); hw->SetLineWidth(3);

  py->Draw("AXIS");
  g_band->Draw("2 same");
  hw->Draw("HIST same");
  py->Draw("PE X0 same");

  TLegend *leg = new TLegend(0.15, 0.55, 0.55, 0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.15);
  leg->SetTextFont(font_code);
  leg->SetTextSize(legend_size);
  leg->AddEntry(py,     "full analysis, Pythia8", "pe1");
  leg->AddEntry(hw,     "full analysis, Herwig",  "l");
  leg->AddEntry(g_band, "generator uncertainty",  "f");
  leg->Draw();

  TLatex cms_bl;
  cms_bl.SetNDC();
  cms_bl.SetTextFont(62); cms_bl.SetTextSize(0.045); cms_bl.DrawLatex(0.10, 0.945, "CMS");
  cms_bl.SetTextFont(52); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.205, 0.945, "Internal");
  cms_bl.SetTextFont(42); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

  // ---- bottom: same "difference relative to the nominal" axis as the syst plot ----
  p_bot->cd();
  TH1D *h_frame = (TH1D *) py->Clone("h_gb_frame");
  h_frame->SetDirectory(nullptr);
  h_frame->Reset();
  if (rmax <= 0.) rmax = 0.05;
  h_frame->SetTitle("");
  h_frame->GetYaxis()->SetRangeUser(-1.6 * rmax, 1.6 * rmax);
  h_frame->GetXaxis()->SetTitle("#Delta r");
  h_frame->GetYaxis()->SetTitle("Herwig / Pythia8 - 1");
  h_frame->GetXaxis()->CenterTitle(true);
  h_frame->GetYaxis()->CenterTitle(true);
  h_frame->GetYaxis()->SetNdivisions(505);
  for (TAxis *ax : { h_frame->GetXaxis(), h_frame->GetYaxis() }) {
    ax->SetTitleFont(font_code);
    ax->SetTitleSize(title_size);
    ax->SetLabelFont(font_code);
    ax->SetLabelSize(label_size);
  }
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetTitleOffset(1.5);
  h_frame->Draw("AXIS");
  g_band_rel->Draw("2 same");

  TH1D *h_rel = (TH1D *) hw->Clone("h_rel_hw");
  h_rel->SetDirectory(nullptr);
  for (int i = 1; i <= nb; ++i) {
    const double v = py->GetBinContent(i);
    h_rel->SetBinContent(i, (v != 0.) ? hw->GetBinContent(i) / v - 1. : 0.);
    h_rel->SetBinError(i,   (v != 0.) ? hw->GetBinError(i) / v : 0.);
  }
  h_rel->SetLineColor(col_alt); h_rel->SetMarkerColor(col_alt);
  h_rel->SetLineStyle(2); h_rel->SetLineWidth(3);
  h_rel->Draw("HIST same");

  TLine *l0 = new TLine(h_frame->GetXaxis()->GetXmin(), 0.,
                        h_frame->GetXaxis()->GetXmax(), 0.);
  l0->SetLineStyle(2);
  l0->Draw();
  p_bot->RedrawAxis();

  const TString stem = R + "final_generator_band_" + SAMPLE;
  c->Print(stem + ".pdf");
  c->Print(stem + ".png");

  printf("\n bin | Pythia8 |  Herwig | band (+/-) |  rel\n");
  for (int i = 1; i <= nb; ++i) {
    const double v = py->GetBinContent(i), d = std::fabs(hw->GetBinContent(i) - v);
    printf(" %3d | %7.4f | %7.4f | %10.4f | %5.2f%%\n", i, v, hw->GetBinContent(i), d,
           v != 0. ? 100 * d / v : 0.);
  }
  printf("\nWrote %s.{pdf,png}\n", stem.Data());
}
