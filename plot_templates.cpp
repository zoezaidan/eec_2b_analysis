#include <TFile.h>
#include <TH2.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>
#include <iostream>
#include "binning_histos_small.h"

// pt bin labels matching jtpt_binsVector = {80, 100, 120, 140}
const int N_PT_BINS = 3;
const char* pt_labels[N_PT_BINS] = {"80-100 GeV", "100-120 GeV", "120-140 GeV"};

// project a 3D histo onto one axis for a given pt bin (ipt=-1 = all bins integrated)
TH1D* project(TH3D* h3, const char* axis, int ipt, const char* name) {
  if (ipt >= 0)
    h3->GetZaxis()->SetRange(ipt + 1, ipt + 1);
  else
    h3->GetZaxis()->SetRange(0, 0);
  TH1D* h = (TH1D*) h3->Project3D(Form("proj_%s_%s", name, axis));
  h->SetName(name);
  h3->GetZaxis()->SetRange(0, 0);
  return h;
}

void normalize(TH1D* h) {
  if (h && h->Integral() > 0) h->Scale(1.0 / h->Integral());
}

// Stack 1B (bottom) + 2B (top), jointly normalised to 1, overlaid with data shape.
void draw_stack(TH1D* h_b_raw, TH1D* h_bb_raw, TH1D* h_data,
                const char* xtitle, const char* var_name,
                const char* pt_label, TCanvas* c, int pad) {

  c->cd(pad);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  if (TString(var_name).Contains("dr")) gPad->SetLogx();

  // Clone so we don't modify the originals (which are already normalised individually)
  TH1D* hs_b  = (TH1D*) h_b_raw ->Clone(Form("stk_b_%s_p%d",  var_name, pad));
  TH1D* hs_bb = (TH1D*) h_bb_raw->Clone(Form("stk_bb_%s_p%d", var_name, pad));

  // Joint normalisation: scale both by 1 / (integral_1B + integral_2B)
  double total = hs_b->Integral() + hs_bb->Integral();
  if (total > 0) { hs_b->Scale(1.0 / total); hs_bb->Scale(1.0 / total); }

  hs_b ->SetFillColor(kAzure+1);   hs_b ->SetLineColor(kAzure+2);
  hs_bb->SetFillColor(kOrange+7);  hs_bb->SetLineColor(kOrange+8);

  THStack* stack = new THStack(Form("stk_%s_p%d", var_name, pad),
                               Form(";%s;Normalised entries", xtitle));
  stack->Add(hs_b);   // 1B at the bottom
  stack->Add(hs_bb);  // 2B on top

  stack->Draw("hist");
  stack->GetXaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitleSize(0.05);

  TLegend* leg = new TLegend(0.55, 0.62, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hs_b,  "1B (jtNbHad=1)", "f");
  leg->AddEntry(hs_bb, "2B (jtNbHad=2)", "f");

  if (h_data) {
    TH1D* hd = (TH1D*) h_data->Clone(Form("stk_data_%s_p%d", var_name, pad));
    hd->SetMarkerStyle(20);
    hd->SetMarkerSize(0.8);
    hd->SetLineColor(kBlack);
    hd->Draw("E1 same");
    leg->AddEntry(hd, "Data HighEG (HLT Jet80)", "lep");
  }
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.042);
  lat.DrawLatex(0.18, 0.88, Form("p_{T}^{jet} #in [%s]", pt_label));
}


void draw_shapes(TH1D* h_b, TH1D* h_bb, TH1D* h_data,
                 const char* xtitle, const char* var_name,
                 const char* pt_label, TCanvas* c, int pad) {

  c->cd(pad);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  if (TString(var_name).Contains("dr")) gPad->SetLogx();

  h_b->SetLineColor(kAzure+1);
  h_b->SetLineWidth(2);
  h_b->SetFillStyle(0);

  h_bb->SetLineColor(kOrange+7);
  h_bb->SetLineWidth(2);
  h_bb->SetFillStyle(0);

  // find y max across all three
  double ymax = h_b->GetMaximum();
  if (h_bb->GetMaximum() > ymax) ymax = h_bb->GetMaximum();
  if (h_data && h_data->GetMaximum() > ymax) ymax = h_data->GetMaximum();

  h_b->SetMaximum(ymax * 1.3);
  h_b->GetXaxis()->SetTitle(xtitle);
  h_b->GetYaxis()->SetTitle("Normalised entries");
  h_b->GetXaxis()->SetTitleSize(0.05);
  h_b->GetYaxis()->SetTitleSize(0.05);
  h_b->Draw("hist");
  h_bb->Draw("hist same");

  TLegend* leg = new TLegend(0.55, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_b,  "MC b  (jtNbHad=1)", "l");
  leg->AddEntry(h_bb, "MC bb (jtNbHad=2)", "l");

  if (h_data) {
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->SetLineColor(kBlack);
    h_data->Draw("E1 same");
    leg->AddEntry(h_data, "Data HighEG (HLT Jet80)", "lep");
  }
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.042);
  lat.DrawLatex(0.18, 0.88, Form("p_{T}^{jet} #in [%s]", pt_label));
}


void plot_templates(TString mc_file,
                    TString data_file = "",
                    TString output_pdf = "template_shapes.pdf") {

  mc_file.ReplaceAll("$mydata", "/data_CMS/cms/zaidan/analysis_lise");
  data_file.ReplaceAll("$mydata", "/data_CMS/cms/zaidan/analysis_lise");

  TFile* fmc = TFile::Open(mc_file);
  if (!fmc || fmc->IsZombie()) {
    std::cerr << "Cannot open MC file: " << mc_file << std::endl;
    return;
  }
  TH3D* h3_b  = (TH3D*) fmc->Get("h3D_b");
  TH3D* h3_bb = (TH3D*) fmc->Get("h3D_bb");
  if (!h3_b || !h3_bb) {
    std::cerr << "Cannot find h3D_b or h3D_bb in MC file." << std::endl;
    return;
  }

  TH3D* h3_data = nullptr;
  TFile* fdata  = nullptr;
  if (data_file != "") {
    fdata = TFile::Open(data_file);
    if (!fdata || fdata->IsZombie()) {
      std::cerr << "Cannot open data file: " << data_file << std::endl;
    } else {
      h3_data = (TH3D*) fdata->Get("h3D_data");
      if (!h3_data) std::cerr << "Cannot find h3D_data in data file." << std::endl;
    }
  }

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas("c", "c", 100, 100);
  c->Print((output_pdf + "[").Data());  // open pdf

  // loop over pt bins + integrated
  for (int ipt = 0; ipt < N_PT_BINS + 1; ipt++) {
    bool integrated = (ipt == N_PT_BINS);
    const char* pt_label = integrated ? "80-140 GeV (integrated)" : pt_labels[ipt];
    int idx = integrated ? -1 : ipt;

    TH1D* h_b_dr  = project(h3_b,  "y", idx, Form("h_b_dr_pt%d",  ipt));
    TH1D* h_bb_dr = project(h3_bb, "y", idx, Form("h_bb_dr_pt%d", ipt));
    TH1D* h_b_mB  = project(h3_b,  "x", idx, Form("h_b_mB_pt%d",  ipt));
    TH1D* h_bb_mB = project(h3_bb, "x", idx, Form("h_bb_mB_pt%d", ipt));

    normalize(h_b_dr);  normalize(h_bb_dr);
    normalize(h_b_mB);  normalize(h_bb_mB);

    TH1D* h_data_dr = nullptr;
    TH1D* h_data_mB = nullptr;
    if (h3_data) {
      h_data_dr = project(h3_data, "y", idx, Form("h_data_dr_pt%d", ipt));
      h_data_mB = project(h3_data, "x", idx, Form("h_data_mB_pt%d", ipt));
      normalize(h_data_dr);
      normalize(h_data_mB);
    }

    TCanvas* cpt = new TCanvas(Form("cpt%d", ipt), Form("pt bin %d", ipt), 1400, 500);
    cpt->Divide(2, 1);
    draw_shapes(h_b_dr, h_bb_dr, h_data_dr, "#DeltaR",     "dr", pt_label, cpt, 1);
    draw_shapes(h_b_mB, h_bb_mB, h_data_mB, "m_{B} [GeV]", "mB", pt_label, cpt, 2);
    cpt->Print(output_pdf.Data());
    delete cpt;

    // --- Stacked histogram: 1B + 2B jointly normalised, vs data shape ---
    // Project fresh (un-normalised) copies so the stack fractions are meaningful
    TH1D* hraw_b_dr  = project(h3_b,  "y", idx, Form("hraw_b_dr_pt%d",  ipt));
    TH1D* hraw_bb_dr = project(h3_bb, "y", idx, Form("hraw_bb_dr_pt%d", ipt));
    TH1D* hraw_b_mB  = project(h3_b,  "x", idx, Form("hraw_b_mB_pt%d",  ipt));
    TH1D* hraw_bb_mB = project(h3_bb, "x", idx, Form("hraw_bb_mB_pt%d", ipt));

    TCanvas* cstk = new TCanvas(Form("cstk%d", ipt),
                                Form("Stack 1B+2B vs data – pt bin %d", ipt), 1400, 500);
    cstk->Divide(2, 1);
    draw_stack(hraw_b_dr, hraw_bb_dr, h_data_dr, "#DeltaR",     "dr", pt_label, cstk, 1);
    draw_stack(hraw_b_mB, hraw_bb_mB, h_data_mB, "m_{B} [GeV]", "mB", pt_label, cstk, 2);
    cstk->Print(output_pdf.Data());
    delete cstk;
  }

  c->Print((output_pdf + "]").Data());  // close pdf
  std::cout << "Saved: " << output_pdf << std::endl;

  fmc->Close();
  if (fdata) fdata->Close();
}
