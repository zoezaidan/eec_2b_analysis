#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "RooUnfold.h"
#include "TLine.h"

void plot_purity_efficiency_response()
{
  gStyle->SetOptStat(1110);
  gStyle->SetPalette(kBird);

  TFile *f = TFile::Open("/data_CMS/cms/zaidan/analysis_lise/Run3/RMatrix_Run3_btagWP872_template_for_fit_histos_3D_qcd_f.root","READ"); //path

  if(!f || f->IsZombie())
    {std::cout << "Cannot open file" << std::endl;  
      return;}

    TH2D *h_purity_num = (TH2D*)f->Get("h_full_purity_numerator_tf");
    TH2D *h_purity_den = (TH2D*)f->Get("h_full_purity_denominator_tf");
    TH2D *h_efficiency_num = (TH2D*)f->Get("h_full_efficiency_numerator_tf");
    TH2D *h_efficiency_den = (TH2D*)f->Get("h_full_efficiency_denominator_tf");
    TH2D *h_purity = (TH2D*)f->Get("h_full_purity_tf");
    TH2D *h_efficiency = (TH2D*)f->Get("h_full_efficiency_tf");

    if(!h_purity_num || !h_purity_den || !h_purity || !h_efficiency_num || !h_efficiency_den || !h_efficiency)
      {std::cout << "Histograms not found" << std::endl;
	    return;}

    double ptBins[4] = {80, 100, 120, 200};

    //TCanvas *cPurity = new TCanvas("cPurity", "Purity vs jtpt and dR", 1800, 600);
    //TCanvas *cEfficiency = new TCanvas("cEfficiency", "Efficiency vs jtpt and dR", 1800, 600);

    TCanvas *cPurity_jtpt_dr = new TCanvas("cPurity_jtpt_dr", "Purity vs jtpt and dR", 600, 600);
    TCanvas *cEfficiency_jtpt_dr = new TCanvas("cEfficiency_jtpt_dr", "Efficiency vs jtpt and dR", 600, 600);
    //TH2D *h2_purity_num = (TH2D*)h_purity_num->Project3D("zy");
    //TH2D *h2_purity_den = (TH2D*)h_purity_den->Project3D("zy");
   
    //TH2D *h2_efficiency_num = (TH2D*)h_efficiency_num->Project3D("zy");
    //TH2D *h2_efficiency_den = (TH2D*)h_efficiency_den->Project3D("zy");
   

    //TH2D *ratio_purity = (TH2D*)h2_purity_num->Clone();
    //ratio_purity->Divide(h2_purity_den);
    //TH2D *ratio_efficiency = (TH2D*)h2_efficiency_num->Clone();
    //ratio_efficiency->Divide(h2_efficiency_den);

    //cPurity->Divide(3,1);
    //cEfficiency->Divide(3,1);

    /* 
    for (int ipt=0; ipt<3; ++ipt){
	    double ptMin = ptBins[ipt];
	    double ptMax = ptBins[ipt+1];

	    int zmin_purity = h_purity->GetZaxis()->FindBin(ptMin + 1e-6);
	    int zmax_purity = h_purity->GetZaxis()->FindBin(ptMax - 1e-6);
      int zmin_efficiency = h_efficiency->GetZaxis()->FindBin(ptMin + 1e-6);
	    int zmax_efficiency = h_efficiency->GetZaxis()->FindBin(ptMax - 1e-6);

	    h_purity->GetZaxis()->SetRange(zmin_purity,zmax_purity);
	    h_efficiency->GetZaxis()->SetRange(zmin_efficiency,zmax_efficiency);

      TH2D *h2_purity = (TH2D*)h_purity->Project3D("xy");
      h2_purity->SetName(Form("purity_pt%d",ipt));
      TH2D *h2_efficiency = (TH2D*)h_efficiency->Project3D("xy");
      h2_efficiency->SetName(Form("efficiency_pt%d",ipt));  

    
      cPurity->cd(ipt+1);
      gPad->SetRightMargin(0.15);
      gPad->SetLogx();

      h2_purity->SetTitle( Form("%.0f < p_{T} < %.0f GeV", ptMin, ptMax));
      h2_purity->GetYaxis()->SetTitle("m_{B}");
      h2_purity->GetXaxis()->SetTitle("#DeltaR");
      h2_purity->GetZaxis()->SetTitle("Purity");
      h2_purity->GetZaxis()->SetRangeUser(0.6,1.0);
      h2_purity->SetMinimum(0.0);
      h2_purity->SetMaximum(1.0);
      h2_purity ->Draw("COLZ");

      cEfficiency->cd(ipt+1);
      gPad->SetRightMargin(0.15);
      gPad->SetLogx();
      h2_efficiency->SetTitle( Form("%.0f < p_{T} < %.0f GeV", ptMin, ptMax));
      h2_efficiency->GetYaxis()->SetTitle("m_{B}");
      h2_efficiency->GetXaxis()->SetTitle("#DeltaR");
      h2_efficiency->GetZaxis()->SetTitle("Efficiency");
      h2_efficiency->GetZaxis()->SetRangeUser(0.6,1.0);
      h2_efficiency->SetMinimum(0.0);
      h2_efficiency->SetMaximum(1.0); 
      h2_efficiency ->Draw("COLZ");
  
    }

    
    cPurity_jtpt_dr->cd(0);
    gPad->SetRightMargin(0.15);
    gPad->SetLogx();

    ratio_purity->GetYaxis()->SetTitle("jtpt");
    ratio_purity->GetXaxis()->SetTitle("#DeltaR");
    ratio_purity->GetZaxis()->SetTitle("Purity");
    ratio_purity->GetZaxis()->SetRangeUser(0.6,1.0);
    ratio_purity->SetMinimum(0.0);
    ratio_purity->SetMaximum(1.0);
    ratio_purity ->Draw("COLX");

    cEfficiency_jtpt_dr->cd(0);
    gPad->SetRightMargin(0.15);
    gPad->SetLogx();
    ratio_efficiency->GetYaxis()->SetTitle("jtpt");
    ratio_efficiency->GetXaxis()->SetTitle("#DeltaR");
    ratio_efficiency->GetZaxis()->SetTitle("Efficiency");
    ratio_efficiency->GetZaxis()->SetRangeUser(0.6,1.0);
    ratio_efficiency->SetMinimum(0.0);
    ratio_efficiency->SetMaximum(1.0);
    ratio_efficiency ->Draw("COLX");
*/

    cPurity_jtpt_dr->cd(0);
    gPad->SetRightMargin(0.15);
    gPad->SetLogx();
    h_purity->SetStats(0);
    h_purity->GetYaxis()->SetTitle("jtpt");
    h_purity->GetXaxis()->SetTitle("#DeltaR");
    h_purity->GetZaxis()->SetTitle("Purity");
    h_purity->GetZaxis()->SetRangeUser(0.6,1.0);
    h_purity->SetMinimum(0.0);
    h_purity->SetMaximum(1.0);
    h_purity ->Draw("COLZ");

    cEfficiency_jtpt_dr->cd(0);
    gPad->SetRightMargin(0.15);
    gPad->SetLogx();
    h_efficiency->SetStats(0);
    h_efficiency->GetYaxis()->SetTitle("jtpt");
    h_efficiency->GetXaxis()->SetTitle("#DeltaR");
    h_efficiency->GetZaxis()->SetTitle("Efficiency");
    h_efficiency->GetZaxis()->SetRangeUser(0.6,1.0);
    h_efficiency->SetMinimum(0.0);
    h_efficiency->SetMaximum(1.0);
    h_efficiency ->Draw("COLZ");

    RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("response_tf_full");

    if (!resp) {
      std::cout << "Response not found!" << std::endl;
      return;}

    TH2 *hResponse = resp->Hresponse();
    hResponse->SetTitle("Response matrix;Detector-level bin;Particle-level bin");

    TCanvas *c = new TCanvas("c","response",800,700);

    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetFixedAspectRatio();
    hResponse ->SetStats(0);
    //hResponse->SetMinimum(0.0);

    // optional normalization style:
    hResponse->Scale(1.0 / hResponse->Integral());

    hResponse->Draw("COLZ");
    double n = hResponse->GetNbinsX();


    c->SaveAs("response_matrix.pdf");
    c->SaveAs("response_matrix.png");






    //cPurity->SaveAs("purity_allPtBins.pdf");
    //cPurity->SaveAs("purity_allPtBins.png");
    //cEfficiency->SaveAs("efficiency_allPtBins.pdf");
    //cEfficiency->SaveAs("efficiency_allPtBins.png");

    cPurity_jtpt_dr->SaveAs("purity_jtpt_dr.pdf");
    cPurity_jtpt_dr->SaveAs("purity_jtpt_dr.png");
    cEfficiency_jtpt_dr->SaveAs("efficiency_jtpt_dr.pdf");
    cEfficiency_jtpt_dr->SaveAs("efficiency_jtpt_dr.png");
}
