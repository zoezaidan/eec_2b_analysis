#include <TFile.h>
#include <TH2.h>
#include <TH3D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>
#include "binning_histos_small.h"

// pt bin labels matching jtpt_binsVector = {80, 100, 120, 140}
const int N_PT_BINS_STK = 3;
const char* pt_labels_stk[N_PT_BINS_STK] = {"80-100 GeV", "100-120 GeV", "120-140 GeV"};

// Project a 3D histogram onto one axis for a given pt bin (idx=-1 = all integrated)
TH1D* project_stk(TH3D* h3, const char* axis, int idx, const char* name) {
  if (idx >= 0)
    h3->GetZaxis()->SetRange(idx + 1, idx + 1);
  else
    h3->GetZaxis()->SetRange(0, 0);
  TH1D* h = (TH1D*) h3->Project3D(Form("proj_%s_%s", name, axis));
  h->SetName(name);
  h3->GetZaxis()->SetRange(0, 0);
  return h;
}

// Draw one pad: overlapping filled histograms (0B+1B+2B, 1B+2B, 2B) with data on top
// Drawing order: largest first so smaller ones are visible on top
void draw_stack_pad(TH1D* h_0b, TH1D* h_1b, TH1D* h_2b, TH1D* h_data,
                    const char* xtitle, const char* var_name,
                    const char* pt_label, TCanvas* c, int pad) {

  c->cd(pad);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.15);
  if (TString(var_name).Contains("dr")) gPad->SetLogx();

  // Clone so originals are untouched
  TH1D* hs_0b = (TH1D*) h_0b->Clone(Form("stk0b_%s_p%d", var_name, pad));
  TH1D* hs_1b = (TH1D*) h_1b->Clone(Form("stk1b_%s_p%d", var_name, pad));
  TH1D* hs_2b = (TH1D*) h_2b->Clone(Form("stk2b_%s_p%d", var_name, pad));

  // Joint normalisation: scale all three by 1 / (I_0b + I_1b + I_2b)
  // so MC fractions reflect true composition and are on the same scale as data
  double total = hs_0b->Integral() + hs_1b->Integral() + hs_2b->Integral();
  if (total > 0) {
    hs_0b->Scale(1.0 / total);
    hs_1b->Scale(1.0 / total);
    hs_2b->Scale(1.0 / total);
  }

  // Build cumulative sums: total = 0B+1B+2B, h12 = 1B+2B, h2 = 2B only
  TH1D* h_total = (TH1D*) hs_0b->Clone(Form("htotal_%s_p%d", var_name, pad));
  h_total->Add(hs_1b);
  h_total->Add(hs_2b);

  TH1D* h_12 = (TH1D*) hs_1b->Clone(Form("h12_%s_p%d", var_name, pad));
  h_12->Add(hs_2b);

  // Colors matching HEP convention:
  //   h_total (0B+1B+2B) = pale yellow background — visible as thin band = Light Jets
  //   h_12    (1B+2B)    = cyan                   — visible as cyan band  = Single B
  //   hs_2b   (2B only)  = medium purple           — fills the core        = Double B
  h_total->SetFillColor(kYellow-9); h_total->SetLineColor(kYellow-7); h_total->SetLineWidth(1);
  h_12->SetFillColor(kCyan-7);      h_12->SetLineColor(kCyan-5);      h_12->SetLineWidth(1);
  hs_2b->SetFillColor(kBlue-5);     hs_2b->SetLineColor(kBlue-4);     hs_2b->SetLineWidth(1);

  // Draw largest first so smaller ones are visible on top
  h_total->SetTitle(Form(";%s;Entries", xtitle));
  h_total->GetXaxis()->SetTitleSize(0.050);
  h_total->GetYaxis()->SetTitleSize(0.050);
  h_total->GetXaxis()->SetLabelSize(0.035);
  h_total->GetYaxis()->SetLabelSize(0.035);
  h_total->GetXaxis()->SetTitleOffset(1.2);
  h_total->GetYaxis()->SetTitleOffset(1.3);
  h_total->GetXaxis()->SetLabelOffset(0.01);
  h_total->GetYaxis()->SetLabelOffset(0.01);
  if (TString(var_name).Contains("dr")) {
    // Use exact bin edges from binning_histos_small.h
    // dr_binsVector[18] = last edge with meaningful data (~0.537)
    h_total->GetXaxis()->SetRangeUser(dr_binsVector[0], dr_binsVector[bins_dr]);
    h_total->GetXaxis()->SetMoreLogLabels();
    h_total->GetXaxis()->SetNoExponent();
  } else {
    h_total->GetXaxis()->SetRangeUser(0., 7.);
    h_total->GetXaxis()->SetNdivisions(7, 0, 0, true);
  }
  h_total->Draw("hist");
  h_12->Draw("hist same");
  hs_2b->Draw("hist same");

  // Data: draw as filled black triangles (no normalisation)
  TH1D* hd = nullptr;
  if (h_data && h_data->Integral() > 0) {
    hd = (TH1D*) h_data->Clone(Form("stkdata_%s_p%d", var_name, pad));
    // Fold overflow into last bin (mB only — for dr, empty large bins are physics)
    if (!TString(var_name).Contains("dr")) {
      int nb = hd->GetNbinsX();
      double ov     = hd->GetBinContent(nb + 1);
      double ov_err = hd->GetBinError(nb + 1);
      hd->SetBinContent(nb, hd->GetBinContent(nb) + ov);
      hd->SetBinError(nb, std::sqrt(hd->GetBinError(nb)*hd->GetBinError(nb) + ov_err*ov_err));
      hd->SetBinContent(nb + 1, 0);
      hd->SetBinError(nb + 1, 0);
    }
    hd->Scale(1.0 / hd->Integral());
    hd->SetMarkerStyle(22);  // filled triangle up
    hd->SetMarkerSize(1.0);
    hd->SetMarkerColor(kBlack);
    hd->SetLineColor(kBlack);
    hd->SetLineWidth(1);
    hd->Draw("E1 same");
  }

  // Legend (top right) — names match what each visible band represents
  TLegend* leg = new TLegend(0.17, 0.60, 0.52, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite, 0.2);
  leg->SetTextSize(0.030);
  leg->AddEntry(h_total, "All", "f");
  leg->AddEntry(h_12,    "BB+B",   "f");
  leg->AddEntry(hs_2b,   "BB",   "f");
  if (hd) leg->AddEntry(hd, "Data", "lep");
  leg->Draw();

  // pT label: large, top of the pad
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.SetTextFont(62);
  lat.DrawLatex(0.18, 0.88, Form("p_{T}^{jet} #in [%s]", pt_label));
}


void plot_stacked_templates(TString mc_file,
                            TString data_file = "",
                            TString output_pdf = "stacked_templates.pdf") {

  mc_file.ReplaceAll("$mydata", "/data_CMS/cms/zaidan/analysis_lise");
  data_file.ReplaceAll("$mydata", "/data_CMS/cms/zaidan/analysis_lise");

  TFile* fmc = TFile::Open(mc_file);
  if (!fmc || fmc->IsZombie()) {
    std::cerr << "Cannot open MC file: " << mc_file << std::endl;
    return;
  }
  TH3D* h3_0b = (TH3D*) fmc->Get("h3D_0b");
  TH3D* h3_1b = (TH3D*) fmc->Get("h3D_b");
  TH3D* h3_2b = (TH3D*) fmc->Get("h3D_bb");
  if (!h3_0b || !h3_1b || !h3_2b) {
    std::cerr << "Cannot find h3D_0b / h3D_b / h3D_bb in MC file." << std::endl;
    fmc->ls();  // print what's actually in the file
    return;
  }

  TH3D* h3_data = nullptr;
  TFile* fdata  = nullptr;
  if (data_file != "") {
    fdata = TFile::Open(data_file);
    if (!fdata || fdata->IsZombie()) {
      std::cerr << "Cannot open data file: " << data_file << std::endl;
    } else {
      h3_data = (TH3D*) fdata->Get("h_count_data");
      if (!h3_data) {
        std::cerr << "Cannot find h3D_data in data file." << std::endl;
        fdata->ls();  // print what's actually in the file
      } else {
        std::cout << "Data histogram found: " << h3_data->GetEntries() << " entries." << std::endl;
      }
    }
  }

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas("c", "c", 100, 100);
  c->Print((output_pdf + "[").Data());  // open PDF

  // Loop over pt bins + integrated
  for (int ipt = 0; ipt < N_PT_BINS_STK + 1; ipt++) {
    bool integrated = (ipt == N_PT_BINS_STK);
    const char* pt_label = integrated ? "80-140 GeV (integrated)" : pt_labels_stk[ipt];
    int idx = integrated ? -1 : ipt;

    // Project raw (un-normalised) so the stack fractions reflect true MC composition
    TH1D* h_0b_dr = project_stk(h3_0b, "y", idx, Form("h0b_dr_pt%d", ipt));
    TH1D* h_1b_dr = project_stk(h3_1b, "y", idx, Form("h1b_dr_pt%d", ipt));
    TH1D* h_2b_dr = project_stk(h3_2b, "y", idx, Form("h2b_dr_pt%d", ipt));

    TH1D* h_0b_mB = project_stk(h3_0b, "x", idx, Form("h0b_mB_pt%d", ipt));
    TH1D* h_1b_mB = project_stk(h3_1b, "x", idx, Form("h1b_mB_pt%d", ipt));
    TH1D* h_2b_mB = project_stk(h3_2b, "x", idx, Form("h2b_mB_pt%d", ipt));

    TH1D* h_data_dr = nullptr;
    TH1D* h_data_mB = nullptr;
    if (h3_data) {
      h_data_dr = project_stk(h3_data, "y", idx, Form("hdata_dr_pt%d", ipt));
      h_data_mB = project_stk(h3_data, "x", idx, Form("hdata_mB_pt%d", ipt));
    }

    TCanvas* cpt = new TCanvas(Form("cstk_pt%d", ipt),
                               Form("Stacked templates – %s", pt_label), 1200, 600);
    cpt->cd();  // make cpt the active canvas before dividing
    cpt->Divide(2, 1);
    draw_stack_pad(h_0b_dr, h_1b_dr, h_2b_dr, h_data_dr,
                   "#DeltaR", "dr", pt_label, cpt, 1);
    draw_stack_pad(h_0b_mB, h_1b_mB, h_2b_mB, h_data_mB,
                   "m_{B} [GeV]", "mB", pt_label, cpt, 2);
    cpt->Print(output_pdf.Data());
    delete cpt;
  }

  c->Print((output_pdf + "]").Data());  // close PDF
  std::cout << "Saved: " << output_pdf << std::endl;

  fmc->Close();
  if (fdata) fdata->Close();
}
