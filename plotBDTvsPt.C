#include "CMSStyle.C"

void plotBDTvsPt(){


  setCMSStyle();

  
  TFile *fin = new TFile("BDTvsPt.root");
  TH3F *h1B_sig_bdt_jtpt_bpt = (TH3F*)fin->Get("h1B_sig_bdt_jtpt_bpt");	
  TH3F *h2B_sig_bdt_jtpt_bpt = (TH3F*)fin->Get("h2B_sig_bdt_jtpt_bpt");	
  TH2F *h1B_bkd_bdt_jtpt = (TH2F*)fin->Get("h1B_bkd_bdt_jtpt");	
  TH2F *h2B_bkd_bdt_jtpt = (TH2F*)fin->Get("h2B_bkd_bdt_jtpt");	

  // 0-100 for b-hadrons
  // 80-140 for jets
  const int nPtBins = 6;
  const int nBBins = 25;

  TH1F *h1B_jtpt_sig[nPtBins], *h2B_jtpt_sig[nPtBins]; 
  TH1F *h1B_jtpt_bkd[nPtBins], *h2B_jtpt_bkd[nPtBins]; 
  TH1F *h1B_bpt_sig[nBBins], *h2B_bpt_sig[nBBins]; 


  gStyle->SetOptStat(0);

  auto styleHist = [](TH1F *h){
    h->SetTitle("");
    h->GetXaxis()->SetTitle("BDT Score");
    h->GetYaxis()->SetTitle("Unit Norm.");
    h->GetXaxis()->SetTitleSize(0.09);
    h->GetYaxis()->SetTitleSize(0.10);
    h->GetXaxis()->SetLabelSize(0.075);
    h->GetYaxis()->SetLabelSize(0.075);
    h->GetXaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetTitleOffset(0.95);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
  };

  auto drawPadLabel = [](const char *txt, bool bottomRow=false, bool leftCol=false, bool firstCanvas=false){
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    double baseSize = firstCanvas ? 0.078 : 0.075;
    double bottomSize = firstCanvas ? 0.072 : 0.065;
    double size = bottomRow ? bottomSize : baseSize;
    if(leftCol) size *= 0.92;
    lat.SetTextSize(size);
    lat.DrawLatex(leftCol ? 0.25 : 0.12,0.84,txt);
  };

  auto drawRightSideInfo = [](TCanvas *canv, bool firstCanvas=false){
    canv->cd();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(firstCanvas ? 0.030 : 0.022);
    lat.SetTextAlign(11);

    if(firstCanvas){
      lat.DrawLatex(0.835,0.79,"|#eta_{jet}| < 2");
    } else {
      lat.DrawLatex(0.835,0.79,"100 < p_{T,jet} < 120 GeV");
      lat.DrawLatex(0.835,0.755,"|#eta_{jet}| < 2");
    }

    TLegend *legColor = new TLegend(0.835, firstCanvas ? 0.58 : 0.58, 0.975, firstCanvas ? 0.73 : 0.70);
    legColor->SetBorderSize(0);
    legColor->SetFillStyle(0);
    legColor->SetTextFont(42);
    legColor->SetTextSize(firstCanvas ? 0.034 : 0.030);

    TLine *l1 = new TLine();
    TLine *l2 = new TLine();
    l1->SetLineColor(kblue);
    l2->SetLineColor(kred);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    legColor->AddEntry(l1,"1B","l");
    legColor->AddEntry(l2,"2B","l");
    legColor->Draw();

    if(firstCanvas){
      TLegend *legStyle = new TLegend(0.835,0.38,0.975,0.53);
      legStyle->SetBorderSize(0);
      legStyle->SetFillStyle(0);
      legStyle->SetTextFont(42);
      legStyle->SetTextSize(0.034);

      TLine *lSolid = new TLine();
      TLine *lDash = new TLine();
      lSolid->SetLineColor(kBlack);
      lDash->SetLineColor(kBlack);
      lSolid->SetLineWidth(2);
      lDash->SetLineWidth(2);
      lDash->SetLineStyle(7);
      legStyle->AddEntry(lSolid,"b-decay","l");
      legStyle->AddEntry(lDash,"non-b","l");
      legStyle->Draw();
    }
  };

  auto drawCMSHeader = [](TCanvas *canv, bool firstCanvas=false){
    canv->cd();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextColor(kBlack);

    double xCMS = firstCanvas ? 0.12 : 0.085;
    double xExtra = firstCanvas ? 0.175 : 0.16;

    lat.SetTextAlign(11);
    lat.SetTextFont(61);
    lat.SetTextSize(0.035);
    lat.DrawLatex(xCMS,0.945,"CMS");

    lat.SetTextFont(52);
    lat.SetTextSize(0.028);
    lat.DrawLatex(xExtra,0.945,"Internal Simulation");

    double xLumi = 0.82;
    lat.SetTextAlign(31);
    lat.SetTextFont(42);
    lat.SetTextSize(0.028);
    lat.DrawLatex(xLumi,0.945,"2024 pp (5.36 TeV)");
  };

  auto c=new TCanvas("c","c",900,600);
  auto c2=new TCanvas("c2","c2",1000,1000);
  
  // Manually build equal-frame multi-pad canvases.
  // The outer pads are made larger so that, after their larger margins,
  // the actual histogram frames remain the same size in every panel.
  auto makeEqualFramePads = [](TCanvas *canv, const char *prefix, int ncol, int nrow){
    const double inner = 0.025;
    const double leftOuter = 0.18;
    const double rightOuter = 0.03;
    const double topOuter = 0.03;
    const double bottomOuter = 0.15;

    const double canvasLeft = 0.05;
    const double canvasRight = 0.18;
    const double canvasTop = 0.06;
    const double canvasBottom = 0.03;

    canv->cd();
    TPad *mother = new TPad(Form("%s_mother",prefix),"",canvasLeft,canvasBottom,1.-canvasRight,1.-canvasTop);
    mother->SetMargin(0,0,0,0);
    mother->SetBorderMode(0);
    mother->Draw();

    double colWidth[10];
    double rowHeight[10];

    double sumColWeights = 0.;
    for(int col=0; col<ncol; col++){
      double lm = (col==0)      ? leftOuter  : inner;
      double rm = (col==ncol-1) ? rightOuter : inner;
      colWidth[col] = 1./(1. - lm - rm);
      sumColWeights += colWidth[col];
    }
    for(int col=0; col<ncol; col++) colWidth[col] /= sumColWeights;

    double sumRowWeights = 0.;
    for(int row=0; row<nrow; row++){
      double tm = (row==0)      ? topOuter    : inner;
      double bm = (row==nrow-1) ? bottomOuter : inner;
      rowHeight[row] = 1./(1. - tm - bm);
      sumRowWeights += rowHeight[row];
    }
    for(int row=0; row<nrow; row++) rowHeight[row] /= sumRowWeights;

    double yHigh = 1.;
    for(int row=0; row<nrow; row++){
      double yLow = (row==nrow-1) ? 0. : yHigh - rowHeight[row];
      double xLow = 0.;
      for(int col=0; col<ncol; col++){
        double xHigh = (col==ncol-1) ? 1. : xLow + colWidth[col];
        int ipad = row*ncol + col + 1;

        mother->cd();
        TPad *pad = new TPad(Form("%s_%d",prefix,ipad), "", xLow, yLow, xHigh, yHigh);
        pad->SetNumber(ipad);
        pad->SetLeftMargin(col==0 ? leftOuter : inner);
        pad->SetRightMargin(col==ncol-1 ? rightOuter : inner);
        pad->SetTopMargin(row==0 ? topOuter : inner);
        pad->SetBottomMargin(row==nrow-1 ? bottomOuter : inner);
        pad->SetBorderMode(0);
        pad->Draw();

        xLow = xHigh;
      }
      yHigh = yLow;
    }
  };

  makeEqualFramePads(c, "c_pad", 3, 2);
  makeEqualFramePads(c2, "c2_pad", 5, 5);

  
  for(int ipt=0;ipt<nPtBins;ipt++){
    
    h1B_jtpt_sig[ipt] = (TH1F*)h1B_sig_bdt_jtpt_bpt->ProjectionX(Form("h1B_jtpt_sig%d",ipt),ipt+1,ipt+1,1,100);  
    h2B_jtpt_sig[ipt] = (TH1F*)h2B_sig_bdt_jtpt_bpt->ProjectionX(Form("h2B_jtpt_sig%d",ipt),ipt+1,ipt+1,1,100);  
    h1B_jtpt_bkd[ipt] = (TH1F*)h1B_bkd_bdt_jtpt->ProjectionX(Form("h1B_jtpt_bkd%d",ipt),ipt+1,ipt+1);  
    h2B_jtpt_bkd[ipt] = (TH1F*)h2B_bkd_bdt_jtpt->ProjectionX(Form("h2B_jtpt_bkd%d",ipt),ipt+1,ipt+1);  

    h1B_jtpt_sig[ipt]->Scale(1./h1B_jtpt_sig[ipt]->Integral());
    h2B_jtpt_sig[ipt]->Scale(1./h2B_jtpt_sig[ipt]->Integral());
    h1B_jtpt_bkd[ipt]->Scale(1./h1B_jtpt_bkd[ipt]->Integral());
    h2B_jtpt_bkd[ipt]->Scale(1./h2B_jtpt_bkd[ipt]->Integral());

    TPad *mother = (TPad*)c->FindObject("c_pad_mother");
    TPad *pad = mother ? (TPad*)mother->FindObject(Form("c_pad_%d",ipt+1)) : 0;
    if(!pad){
      std::cout << "Missing pad c_pad_" << ipt+1 << std::endl;
      continue;
    }
    pad->cd();
    pad->SetLogy();
    h1B_jtpt_sig[ipt]->SetMinimum(0.004);
    h1B_jtpt_sig[ipt]->SetMaximum(1.);
    styleHist(h1B_jtpt_sig[ipt]);

    if(ipt%3!=0){
      h1B_jtpt_sig[ipt]->GetYaxis()->SetLabelSize(0);
      h1B_jtpt_sig[ipt]->GetYaxis()->SetTitleSize(0);
    }
    if(ipt<3){
      h1B_jtpt_sig[ipt]->GetXaxis()->SetLabelSize(0);
      h1B_jtpt_sig[ipt]->GetXaxis()->SetTitleSize(0);
    }
    if(ipt>=3 && ipt%3==0){
      h1B_jtpt_sig[ipt]->GetYaxis()->SetTitleSize(0.085);
      h1B_jtpt_sig[ipt]->GetYaxis()->SetTitleOffset(1.2);
    }

    h1B_jtpt_sig[ipt]->SetLineColor(kblue);
    h1B_jtpt_sig[ipt]->Draw("hist,E0");   
    h2B_jtpt_sig[ipt]->SetLineColor(kred);
    h2B_jtpt_sig[ipt]->Draw("hist,E0,same");

    h1B_jtpt_bkd[ipt]->SetLineStyle(7);
    h2B_jtpt_bkd[ipt]->SetLineStyle(7);

    h1B_jtpt_bkd[ipt]->SetLineColor(kblue);
    h1B_jtpt_bkd[ipt]->Draw("hist,E0,same");   
    h2B_jtpt_bkd[ipt]->SetLineColor(kred);
    h2B_jtpt_bkd[ipt]->Draw("hist,E0,same");
    drawPadLabel(Form("%d < jet p_{T} < %d GeV",80+20*ipt,100+20*ipt), ipt>=3, ipt%3==0, true);
  }
  
 
  for(int iB=0;iB<nBBins;iB++){

      
    h1B_bpt_sig[iB] = (TH1F*)h1B_sig_bdt_jtpt_bpt->ProjectionX(Form("h1B_bpt_sig%d",iB),2,2,2*iB+1,2*(iB+1));  
    h2B_bpt_sig[iB] = (TH1F*)h2B_sig_bdt_jtpt_bpt->ProjectionX(Form("h2B_bpt_sig%d",iB),2,2,2*iB+1,2*(iB+1));  
    
    h1B_bpt_sig[iB]->Scale(1./h1B_bpt_sig[iB]->Integral());
    h2B_bpt_sig[iB]->Scale(1./h2B_bpt_sig[iB]->Integral());

    TPad *mother = (TPad*)c2->FindObject("c2_pad_mother");
    TPad *pad = mother ? (TPad*)mother->FindObject(Form("c2_pad_%d",iB+1)) : 0;
    if(!pad){
      std::cout << "Missing pad c2_pad_" << iB+1 << std::endl;
      continue;
    }
    pad->cd();
    pad->SetLogy();
    h1B_bpt_sig[iB]->SetMinimum(0.004);
    h1B_bpt_sig[iB]->SetMaximum(1.);
    styleHist(h1B_bpt_sig[iB]);
    if(iB%5!=0){
      h1B_bpt_sig[iB]->GetYaxis()->SetLabelSize(0);
      h1B_bpt_sig[iB]->GetYaxis()->SetTitleSize(0);
    }
    if(iB<20){
      h1B_bpt_sig[iB]->GetXaxis()->SetLabelSize(0);
      h1B_bpt_sig[iB]->GetXaxis()->SetTitleSize(0);
    }
    if(iB>=20 && iB%5==0){
      h1B_bpt_sig[iB]->GetYaxis()->SetTitleSize(0.085);
      h1B_bpt_sig[iB]->GetYaxis()->SetTitleOffset(1.1);
    }
    h1B_bpt_sig[iB]->SetLineColor(kblue);
    h1B_bpt_sig[iB]->Draw("hist,E0");
    h2B_bpt_sig[iB]->SetLineColor(kred);
    h2B_bpt_sig[iB]->Draw("hist,E0,same");


    drawPadLabel(Form("%d < B p_{T} < %d GeV",4*iB,4*(iB+1)), iB>=20, iB%5==0, false);
  }

  
  
  drawCMSHeader(c,true);
  drawCMSHeader(c2,false);

  drawRightSideInfo(c,true);
  drawRightSideInfo(c2,false);
  
  //drawCMSLabel(c, "Internal Simulation", "2024 pp (5.36 TeV)");

  
}
