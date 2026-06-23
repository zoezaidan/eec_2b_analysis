#include "rootlogon.C"
#include "CMSStyle.C"


void plotNice_bdt_roc(){

  setCMSStyle();

  auto *fin = new TFile("bdt_roc_from_forest.root");

  auto *g1b = (TGraph*) fin->Get("roc1b");
  auto *g2b = (TGraph*) fin->Get("roc2b");
  auto *gAll = (TGraph*) fin->Get("roc_all");
  

  auto hFrame = new TH1F ("hFrame","hFrame",1,0.4,1);
  hFrame->GetXaxis()->SetRangeUser(0.4,1.);
  hFrame->GetYaxis()->SetRangeUser(0.002,1.);

  TCanvas *c=new TCanvas("c","c",600,600);
  c->SetGridx();
  c->SetGridy();
  c->SetLogy();


  g1b->SetLineColor(kblue);
  g2b->SetLineColor(kred);
  gAll->SetLineColor(kgreen);

  g1b->SetLineWidth(3);
  g2b->SetLineWidth(3);
  gAll->SetLineWidth(3);
  
  hFrame->GetXaxis()->SetTitle(g1b->GetXaxis()->GetTitle());
  hFrame->GetYaxis()->SetTitle(g1b->GetYaxis()->GetTitle());
  
  hFrame->Draw();
  g1b->Draw("l");
  g2b->Draw("l");
  gAll->Draw("l");

  TLegend* leg = new TLegend(0.2, 0.55, 0.45, 0.8);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(gAll, "all B (1B + #geq2B)", "l");
  leg->AddEntry(g1b,   "1B", "l");
  leg->AddEntry(g2b,   "#geq2B", "l");
  leg->Draw();


  
  TLatex* btagLbl = new TLatex();
  btagLbl->SetNDC();
  btagLbl->SetTextFont(42);
  btagLbl->SetTextSize(0.04);
  btagLbl->DrawLatex(0.2, 0.85, "b-tag WP > 0.872 (UParT)");
  
  drawCMSLabel(c, "Internal Simulation", "2024 pp (5.36 TeV)");
}
