#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <iostream>
#include <vector>
#include "TString.h"
#include "tTree.h"

bool skipMC_roc(double pt, double pthat) {
  if (pthat < 0.35 * pt) return true;
  return false;
}

// Build ROC TGraph: x = sig_eff, y = bkg_eff.
// h_sig_num: jets passing the signal selection (numerator)
// h_sig_den: total signal denominator (can differ from h_sig_num, e.g. all 2b vs 2b+2sv)
// h_bkg:     background histogram
TGraph* makeROC(TH1D* h_sig_num, TH1D* h_sig_den, TH1D* h_bkg) {
  int nbins = h_sig_num->GetNbinsX();
  double total_sig = h_sig_den->Integral(0, nbins + 1);
  double total_bkg = h_bkg->Integral(0, nbins + 1);
  if (total_sig == 0 || total_bkg == 0) return new TGraph();

  std::vector<double> sig_eff, bkg_eff;
  for (int ibin = nbins + 1; ibin >= 0; ibin--) {
    double s = h_sig_num->Integral(ibin, nbins + 1);
    double b = h_bkg->Integral(ibin, nbins + 1);
    sig_eff.push_back(s / total_sig);
    bkg_eff.push_back(b / total_bkg);
  }
  TGraph* roc = new TGraph((int)sig_eff.size(), sig_eff.data(), bkg_eff.data());
  return roc;
}

// Convenience overload: numerator and denominator are the same histogram.
TGraph* makeROC(TH1D* h_sig, TH1D* h_bkg) {
  return makeROC(h_sig, h_sig, h_bkg);
}


void calc_btag_roc(TString filename, TString output_folder, TString output_hist,
                   Float_t pT_low, Float_t pT_high, double target_mistag = 1e-3) {

  tTree t;
  t.Init(filename, true);
  t.SetBranchStatus("*", 0);

  std::vector<TString> active_branches = {
    "jtpt", "jteta", "nref", "jtNsvtx", "jtNbHad", "discr_particleNet_BvsAll",
    "weight", "pthat",
    "HLT_HIAK4PFJet40_v1"
  };
  t.SetBranchStatus(active_branches, 1);

  const int   nbins_disc = 1000;
  const double disc_min  = -1.0;
  const double disc_max  =  1.0;

  // discriminator histograms per truth category
  TH1D *h_0b    = new TH1D("h_disc_0b",    "BvsAll score (0b);score;entries",       nbins_disc, disc_min, disc_max);
  TH1D *h_1b    = new TH1D("h_disc_1b",    "BvsAll score (1b);score;entries",       nbins_disc, disc_min, disc_max);
  TH1D *h_2b    = new TH1D("h_disc_2b",    "BvsAll score (2b);score;entries",       nbins_disc, disc_min, disc_max);
  TH1D *h_2b2sv = new TH1D("h_disc_2b2sv", "BvsAll score (2b,#geq2 SV);score;entries", nbins_disc, disc_min, disc_max);
  TH1D *h_all   = new TH1D("h_disc_all",   "BvsAll score (1b+2b);score;entries",    nbins_disc, disc_min, disc_max);

  h_0b->Sumw2(); h_1b->Sumw2(); h_2b->Sumw2(); h_2b2sv->Sumw2(); h_all->Sumw2();

  Long64_t n_events = t.GetEntries();
  for (Long64_t ient = 0; ient < n_events; ient++) {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
    t.GetEntry(ient);

    if (!(t.HLT_HIAK4PFJet40_v1)) continue;

    for (Int_t ijet = 0; ijet < t.nref; ijet++) {
      if (std::abs(t.jteta[ijet]) > 1.9) continue;
      if (skipMC_roc(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < pT_low || t.jtpt[ijet] > pT_high) continue;

      double score = t.discr_particleNet_BvsAll[ijet];
      int    nbHad = t.jtNbHad[ijet];
      int    nSV   = t.jtNsvtx[ijet];
      double w     = t.weight;

      if (nbHad == 0)              h_0b->Fill(score, w);
      if (nbHad >= 1)              h_all->Fill(score, w);
      if (nbHad == 1)              h_1b->Fill(score, w);
      if (nbHad >= 2)              h_2b->Fill(score, w);
      if (nbHad >= 2 && nSV >= 2) h_2b2sv->Fill(score, w);
    }
  }
  std::cout << std::endl;

  // ------- Sanity check: raw 2b+2sv / 2b fraction -------
  double n_2b     = h_2b->Integral(0, h_2b->GetNbinsX() + 1);
  double n_2b2sv  = h_2b2sv->Integral(0, h_2b2sv->GetNbinsX() + 1);
  std::cout << "Weighted:   2b+2sv / 2b = " << n_2b2sv / n_2b
            << "  (" << n_2b2sv << " / " << n_2b << ")" << std::endl;

  // fraction of 2b+2sv / 2b after a given b-tag score cut
  auto svFractionAtCut = [&](double cut) {
    int bin   = h_2b->FindBin(cut);
    int nend  = h_2b->GetNbinsX() + 1;
    double n2b    = h_2b->Integral(bin, nend);
    double n2b2sv = h_2b2sv->Integral(bin, nend);
    double ratio  = (n2b > 0) ? n2b2sv / n2b : 0;
    std::cout << "score > " << cut << "  →  2b+2sv / 2b = " << ratio
              << "  (" << n2b2sv << " / " << n2b << ")" << std::endl;
    return ratio;
  };
  svFractionAtCut(0.99);  // fixed WP for reference

  // ------- Build ROC curves (background = 0b jets) -------
  TGraph* roc_1b    = makeROC(h_1b,    h_0b);  roc_1b->SetName("roc_1b");
  TGraph* roc_2b    = makeROC(h_2b,    h_0b);  roc_2b->SetName("roc_2b");
  TGraph* roc_2b2sv = makeROC(h_2b2sv, h_2b, h_0b);  roc_2b2sv->SetName("roc_2b2sv");
  TGraph* roc_all   = makeROC(h_all,   h_0b);  roc_all->SetName("roc_all");


  //sanity check
  // ------- Working point at score > 0.99 -------
  const double wp_score = 0.99;
  auto wpEff = [](TH1D* h) {
    int bin = h->FindBin(0.99);
    return h->Integral(bin, h->GetNbinsX() + 1) / h->Integral(0, h->GetNbinsX() + 1);
  };
  double wp_mistag    = wpEff(h_0b);
  double wp_sig_1b    = wpEff(h_1b);
  double wp_sig_2b    = wpEff(h_2b);
  double wp_sig_2b2sv = h_2b2sv->Integral(h_2b2sv->FindBin(wp_score), h_2b2sv->GetNbinsX()+1)
                       / h_2b->Integral(0, h_2b->GetNbinsX()+1);
  double wp_sig_all   = wpEff(h_all);

  // ------- Score threshold corresponding to a target mistag rate -------
  auto findThreshold = [&](double target_mistag) {
    double total_bkg = h_0b->Integral(0, h_0b->GetNbinsX() + 1);
    for (int ibin = h_0b->GetNbinsX() + 1; ibin >= 0; ibin--) {
      double mistag = h_0b->Integral(ibin, h_0b->GetNbinsX() + 1) / total_bkg;
      if (mistag >= target_mistag) {
        double score_threshold = h_0b->GetXaxis()->GetBinLowEdge(ibin);
        double eff_1b    = h_1b->Integral(ibin, h_1b->GetNbinsX()+1) / h_1b->Integral(0, h_1b->GetNbinsX()+1);
        double eff_2b    = h_2b->Integral(ibin, h_2b->GetNbinsX()+1) / h_2b->Integral(0, h_2b->GetNbinsX()+1);
        double eff_2b2sv = h_2b2sv->Integral(ibin, h_2b2sv->GetNbinsX()+1) / h_2b->Integral(0, h_2b->GetNbinsX()+1);
        double eff_all   = h_all->Integral(ibin, h_all->GetNbinsX()+1) / h_all->Integral(0, h_all->GetNbinsX()+1);
        std::cout << "\nThreshold for mistag rate " << target_mistag << ":" << std::endl;
        std::cout << "  score > " << score_threshold << "  (actual mistag = " << mistag << ")" << std::endl;
        std::cout << "  sig eff 1b     = " << eff_1b    << std::endl;
        std::cout << "  sig eff 2b     = " << eff_2b    << std::endl;
        std::cout << "  sig eff 2b+2sv = " << eff_2b2sv << std::endl;
        std::cout << "  sig eff 1b+2b  = " << eff_all   << std::endl;
        return score_threshold;
      }
    }
    return -999.0;
  };
  double score_at_target = findThreshold(target_mistag);
  if (score_at_target > -999)
    svFractionAtCut(score_at_target);

  // efficiencies at the target mistag WP (for the second marker on the plot)
  auto effAtScore = [](TH1D* h_num, TH1D* h_den, double score) {
    int bin = h_num->FindBin(score);
    return h_num->Integral(bin, h_num->GetNbinsX()+1) / h_den->Integral(0, h_den->GetNbinsX()+1);
  };
  double mt_mistag    = effAtScore(h_0b,    h_0b,  score_at_target);
  double mt_sig_1b    = effAtScore(h_1b,    h_1b,  score_at_target);
  double mt_sig_2b    = effAtScore(h_2b,    h_2b,  score_at_target);
  double mt_sig_2b2sv = effAtScore(h_2b2sv, h_2b,  score_at_target);
  double mt_sig_all   = effAtScore(h_all,   h_all, score_at_target);

  // pastel colors
  const Color_t col_1b    = TColor::GetColor("#7EB8D4"); // pastel blue
  const Color_t col_2b    = TColor::GetColor("#F4A4A4"); // pastel red
  const Color_t col_2b2sv = TColor::GetColor("#90C98F"); // pastel green
  const Color_t col_all   = TColor::GetColor("#C3A6D8"); // pastel purple

  auto styleGraph = [](TGraph* g, Color_t col, const char* title) {
    g->SetLineColor(col); g->SetLineWidth(3);
    g->SetMarkerColor(col);
    g->SetTitle(title);
  };
  // x axis: sig efficiency = P(pass | gen-Nb)
  // y axis: bkg efficiency = P(pass | gen-0b)
  const char* axes = "b-tag ROC;"
    "signal efficiency;"
    "mistag rate ";
  styleGraph(roc_1b,    col_1b,    axes);
  styleGraph(roc_2b,    col_2b,    axes);
  styleGraph(roc_2b2sv, col_2b2sv, axes);
  styleGraph(roc_all,   col_all,   axes);

  // ------- helper lambda to dress a canvas -------
  auto drawROC = [&](TCanvas* cv, bool logy = false) {
    cv->SetGrid();
    if (logy) cv->SetLogy();

    roc_2b2sv->Draw("AL");
    roc_2b2sv->GetXaxis()->SetLimits(0, 1);
    roc_2b2sv->GetYaxis()->SetRangeUser(logy ? 1e-6 : 0, 1);
    roc_2b->Draw("L same");
    roc_1b->Draw("L same");
    roc_all->Draw("L same");

    // diagonal (random classifier) — start from 1e-6 on log scale to avoid log(0)
    TLine* diag = new TLine(logy ? 1e-6 : 0, logy ? 1e-6 : 0, 1, 1);
    diag->SetLineStyle(2); diag->SetLineColor(kGray+1);
    diag->Draw();

    // target mistag WP markers — diamond (style 33), same curve colors
    if (score_at_target > -999) {
      auto makeMTmarker = [](double x, double y, Color_t col) {
        TMarker* m = new TMarker(x, y, 33);
        m->SetMarkerColor(col); m->SetMarkerSize(2.5);
        m->Draw();
      };
      makeMTmarker(mt_sig_1b,    mt_mistag, col_1b);
      makeMTmarker(mt_sig_2b,    mt_mistag, col_2b);
      makeMTmarker(mt_sig_2b2sv, mt_mistag, col_2b2sv);
      makeMTmarker(mt_sig_all,   mt_mistag, col_all);

      TLine* mt_line = new TLine(0, mt_mistag, 1, mt_mistag);
      mt_line->SetLineStyle(3); mt_line->SetLineColor(kGray+2); mt_line->SetLineWidth(1);
      mt_line->Draw();
    }

    // legend — bottom-right
    TLegend* leg = new TLegend(0.55, 0.12, 0.88, 0.45);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.033);
    leg->AddEntry(roc_1b,    "1b vs 0b",                     "l");
    leg->AddEntry(roc_2b,    "#geq 2b vs 0b",                "l");
    leg->AddEntry(roc_2b2sv, "#geq 2b (#geq 2 SV) vs 0b",  "l");
    leg->AddEntry(roc_all,   "1b + #geq 2b vs 0b",          "l");
    if (score_at_target > -999) {
      TMarker* mt_dummy = new TMarker(0, 0, 33);
      mt_dummy->SetMarkerColor(kGray+2); mt_dummy->SetMarkerSize(2.5);
      leg->AddEntry(mt_dummy, Form("WP: score > %.3f", score_at_target), "p");
    }
    leg->Draw();

    // definitions box — top-left
    TPaveText* def = new TPaveText(0.12, 0.66, 0.62, 0.88, "NDC");
    def->SetBorderSize(0); def->SetFillStyle(0); def->SetTextAlign(12); def->SetTextSize(0.028);
    def->AddText("Signal efficiency: True Positives / All Positives");
    def->AddText("Mistag rate: False Positives / All Negatives");
    def->AddText("score = discr_particleNet_BvsAll,  t scanned in [#minus1, 1]");
    def->AddText(Form("%.0f < p_{T} < %.0f GeV,  |#eta| < 1.9", (double)pT_low, (double)pT_high));
    def->Draw();
  };

  // ------- linear scale canvas -------
  TCanvas* c = new TCanvas("c_roc", "b-tag ROC curves", 800, 700);
  drawROC(c);
  c->SaveAs((output_folder + output_hist + "_roc.pdf").Data());

  // ------- log-y canvas -------
  TCanvas* c_logy = new TCanvas("c_roc_logy", "b-tag ROC curves (log y)", 800, 700);
  drawROC(c_logy, true);
  c_logy->SaveAs((output_folder + output_hist + "_roc_logy.pdf").Data());

  // ------- Write to ROOT file -------
  TFile outFile((output_folder + output_hist + ".root").Data(), "RECREATE");
  h_0b->Write();    h_1b->Write();    h_2b->Write();
  h_2b2sv->Write(); h_all->Write();
  roc_1b->Write();  roc_2b->Write();
  roc_2b2sv->Write(); roc_all->Write();
  c->Write();
  c_logy->Write();
  outFile.Close();

  std::cout << "Output written to: " << output_folder + output_hist + ".root" << std::endl;

  delete c;
  delete c_logy;
}


void btag_roc(Int_t dataType = 1, Float_t pT_low = 80, Float_t pT_high = 140, double target_mistag = 1e-3) {

  TString output_folder = TString(gSystem->DirName(__FILE__)) + "/";
  TString filename, output_hist;

  if (dataType == 1) {
    filename    = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    output_hist = "btag_roc_bjet";
    std::cout << "MC bjet" << std::endl;
  } else if (dataType == 2) {
    filename    = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    output_hist = "btag_roc_qcd";
    std::cout << "MC dijet" << std::endl;
  } else {
    std::cerr << "ROC curve requires MC (dataType 1=bjet, 2=dijet)" << std::endl;
    return;
  }

  calc_btag_roc(filename, output_folder, output_hist, pT_low, pT_high, target_mistag);
  std::cout << "Done :)" << std::endl;
}
