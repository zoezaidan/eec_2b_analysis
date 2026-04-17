
// ROOT includes 
#include <TString.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TError.h>
#include <TObject.h>
#include <TString.h>
#include <TFile.h>
#include <TStyle.h>
#include <THStack.h>
#include <TMath.h>
#include <TSystem.h>

// RooFit related 
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooFit.h>


// name for outpur directory (in path) 
TString sDirname = "TemplateFitOutput"; // has to be before Help_Functions.h

std::unique_ptr<TCanvas> draw_template_fit_result(
    TString fout_name,
    TFile* foutputPlots,
    TString &dataset,
    TString &folder,
    TString &pT_selection, Int_t pt_bin = 0, bool isIntegDeltaR = false);
    
void DrawCommonTextTopRight(TPad*pad,  int ibin_dr, int ibin_pt, const char* extra="") {
    /*Draw dr range, pt range + default Build Legend*/    
    // you can pass dirctly a canvas pointer 
    pad->cd(); 
    // draw auto legend first
    auto leg = pad->BuildLegend(0.6, 0.5, 0.88, 0.75); // 0.6, 0.5, 0.88, 0.7
        // optional styling
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

    // -- text for other information
    double right = 1.0 - gPad->GetRightMargin();
    double top   = 1.0 - gPad->GetTopMargin();

    TLatex latex;
    latex.SetNDC(); // relative coorinates
    latex.SetTextSize(0.028); // 0.035
    latex.SetTextAlign(13); // left top alignment
    latex.SetTextFont(42); // font helvatic normal (42). bold(62)

    double x = right - 0.35; //  0.3 shift left from right edge
    double y = top - 0.03;

    // legend text setting: read bins to numbers 
    double pt_first = 0;
    double pt_last = 0;
    double dr_first = 0;
    double dr_last = 0;
    if(!ibin_dr){dr_first = 0; dr_last = 1;} // converted to infintiy 
        else { dr_first =  dr_binsVector[ibin_dr -1]; dr_last  = dr_binsVector[ibin_dr];}

    if (!ibin_pt){pt_first = jtpt_binsVector[0];  pt_last = jtpt_binsVector[jtpt_bins];} 
        else{pt_first = jtpt_binsVector[ibin_pt-1]; pt_last = jtpt_binsVector[ibin_pt]; }

        cout << "pt bin #"<< ibin_pt << "pt first and last are: " << pt_first << ", " << pt_last << endl;

    if (dr_last == 1.0) latex.DrawLatex(x, y, Form("%g < #DeltaR < #infty", dr_first));
    else latex.DrawLatex(x, y, Form("%g < #DeltaR < %g", dr_first, dr_last));
    latex.DrawLatex(x, y - 0.04, Form("%g < p_{T} < %g GeV", pt_first, pt_last));
    latex.DrawLatex(x, y - 0.08, extra);
}


void AddRatioPlot(TH1* h1, TH1* h2) {
    // draw ratio of two hists
    TH1* ratio = (TH1*)h1->Clone("ratio");
    ratio->Reset();
    ratio->Divide(h1, h2);

    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("Data/Fit");
    ratio->GetYaxis()->CenterTitle(true);

    // cout << "ratio max and min "<< ratio->GetMaximum() << ", " << ratio->GetMinimum() <<endl;
    double max = 0; 
    if (ratio->GetMaximum() > 10) max = 5; else max = ratio->GetMaximum() + 0.2; // avoid very very large values when fit is almost zero.
    ratio->SetMaximum( max ); // max
    ratio->SetMinimum( ratio->GetMinimum() - 0.2 ); // min

        // larger label
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetYaxis()->SetTitleOffset(0.5);

        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->GetXaxis()->SetTitleOffset(0.8);

    ratio->SetMarkerStyle(20);
    ratio->SetLineColor(kBlack);
    ratio->Draw("EP");

    // unity line
    TLine* line = new TLine(
        ratio->GetXaxis()->GetXmin(), 1.0,
        ratio->GetXaxis()->GetXmax(), 1.0
    );
    line->SetLineStyle(2);
    line->Draw();
}

//Draws the result of the template fit
std::unique_ptr<TCanvas> draw_template_fit_result(
    TString fout_name,
    TFile* foutputPlots,
    TString &dataset,
    TString &folder,
    TString &pT_selection, Int_t pt_bin)
{
    
    // ----------------------------------
    // --- WORK IN PROGRESS -------------
    // ----------------------------------


    // pt_bin is kept with convenstion starting from 0 (for integarted interval, then then the rest) 

    // For trvial tests 
    TString trivialMC_label = fout_name.Contains("trivialMC") ? "_trivialMC": "";

    //Get fractions for the jtpt bin pt_bin
    TFile *file = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "read");
    if (!file) {Error("Input File:", "File does not exist'%s'", file->GetName());return nullptr;}
    cout << "input file name " << file->GetName() << endl;

    // -- Read output of template fit files 
    TH1D *h;
    TH1D *htrue;
    TH1D *hbkg;
    TH1D *hbkg_true;
    TString sname_canvas;
        sname_canvas = Form("c_%s_ptbin_%d%s", dataset.Data(), pt_bin,  trivialMC_label.Data());

    // -- Read input TH2D  summary of signa and background fractions 
    //signal fraction
    // Fit result (uncertaintiy already at the point)
    TH2D *h_2D = (TH2D*)file->Get("h_sig_fraction"); h_2D->SetDirectory(nullptr);
        h = (TH1D*) h_2D->ProjectionX("h", pt_bin+1 , pt_bin+1);
        h->GetXaxis()->SetTitle("Bin (\\#Delta\\ r)");
    // True fractions (uncertaintiy already at the point)
    TH2D *htrue_2D = (TH2D*)file->Get("h_sig_frac_true"); htrue_2D->SetDirectory(nullptr);
        htrue = (TH1D*)htrue_2D->ProjectionX("htrue", pt_bin+1 , pt_bin+1);
        
        // and its uncertainity 
        // TH2D *htrue_err_2D = (TH2D*)file->Get("h_sig_frac_true_error");
        // TH1D* htrue_err = (TH1D*) htrue_err_2D->ProjectionX("htrue_err", pt_bin+1 , pt_bin+1);
        // for (int ibin = 1; ibin <= htrue_err->GetNbinsX(); ibin++)
        // {
        //     htrue->SetBinError(ibin, htrue_err->GetBinContent(ibin));
        // }

    //background fraction
    TH2D *hbkg_2D = (TH2D*)file->Get("h_bkg_fraction"); hbkg_2D->SetDirectory(nullptr);
        hbkg = (TH1D*)hbkg_2D->ProjectionX("hbkg", pt_bin+1, pt_bin+1);
    TH2D *hbkg_true_2D = (TH2D*)file->Get("h_bkg_frac_true");   hbkg_true_2D->SetDirectory(nullptr);
        hbkg_true = (TH1D*)hbkg_true_2D->ProjectionX("hbkg_true",pt_bin+1, pt_bin+1);
        // // and its uncertainity 
        // TH2D *hbkg_true_err_2D = (TH2D*)file->Get("h_bkg_frac_true_error");
        // TH1D* hbkg_true_err = (TH1D*) hbkg_true_err_2D->ProjectionX("hbkg_true_err", pt_bin,pt_bin);
        // for (int ibin = 1; ibin <= hbkg_true_err->GetNbinsX(); ibin++)
        // {
        //     /// Assign the uncertainity to the iD hist of true signal Vs. DeltaR 
        //     hbkg_true->SetBinError(ibin, hbkg_true_err->GetBinContent(ibin));
        // }

    // -- deattach hists from input root file: for drawing 
    h->SetDirectory(nullptr);
    htrue->SetDirectory(nullptr);
    hbkg->SetDirectory(nullptr);
    hbkg_true->SetDirectory(nullptr);

    // -- Output summary canvas
    auto c = std::make_unique<TCanvas>(sname_canvas,"Template fit result", 900, 1100);// 800 x 800
        gROOT->GetListOfCanvases()->Remove(c.get()); // to save it later 
        c->SetTitle("");

    // --- Pads
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);
    // -- dettach pads from root ownership
    pad1->SetBit(kCanDelete, false);
    pad2->SetBit(kCanDelete, false);

    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.08);
    pad1->SetLogx();
    
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.12);
    pad2->SetLogx();
    pad1->Draw();
    pad2->Draw();
    
    //--- pad 1
    pad1->cd();

    //Draw results
    h->SetStats(0);
    h->SetTitle("");
    h->GetYaxis()->SetRangeUser(0,1);
    // h->GetXaxis()->CenterTitle(true);
    h->GetXaxis()->SetTitle("");
    h->GetYaxis()->SetTitle("Signal or Background fraction");
    h->GetYaxis()->CenterTitle(true);
    h->SetMarkerColor(kBlue);
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    
    h->SetDirectory(nullptr);
    h->Draw("EH");
    
    htrue->SetDirectory(nullptr);
    htrue->SetMarkerColor(kRed);
    htrue->SetLineColor(kRed);
    htrue->SetLineStyle(9);
    htrue->SetLineWidth(2);
    htrue->Draw("same");

    hbkg->SetDirectory(nullptr);
    hbkg->SetStats(0);
    hbkg->SetMarkerColor(kCyan+2);
    hbkg->SetLineColor(kCyan+2);
    hbkg->Draw("EH SAME");
    hbkg_true->SetDirectory(nullptr);
    hbkg_true->SetMarkerColor(kMagenta);
    hbkg_true->SetLineColor(kMagenta);
    hbkg_true->SetLineStyle(9);
    hbkg_true->Draw("same");

    

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    // if(!isIntegDeltaR) test_info_text->DrawLatex(0.15, 0.5, Form("%g < p_{T} < %g GeV", jtpt_binsVector[pt_bin-1], jtpt_binsVector[pt_bin]));
    // else{ 
    //     TString srangeDeltaR = Form("#DeltaR [Full range]");
    //     test_info_text->DrawLatex(0.15, 0.5, srangeDeltaR);
    // }
    // test_info_text->DrawLatex(0.15, 0.45, dataset);
    // test_info_text->DrawLatex(0.15, 0.4, "Template from MC bjet");
    test_info_text->Draw("SAME");


    TLegend *leg = new TLegend(0.406,0.39,0.7,0.55, ""); //get position from legend.C file (lower left corner is 0,0)
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h, "fitted signal");
    leg->AddEntry(htrue, "MC signal");
    leg->AddEntry(hbkg, "fitted B background");
    leg->AddEntry(hbkg_true, "MC B background");
    leg->Draw("same");

    pad2->cd();

  // Create ratio histogram
    TH1D* h_ratio = (TH1D*) h->Clone("h_ratio");
    h_ratio->SetDirectory(nullptr);
    h_ratio->Divide(htrue);
    // for (int ibin = 1; ibin <= h_ratio->GetNbinsX() ; ibin++)
    // {
    //     h_ratio->SetBinError(ibin, h->GetBinError(ibin)/htrue->GetBinContent(ibin));
    //     // cout << "bin error = "<<  h->GetBinError(ibin)/htrue->GetBinContent(ibin)  << endl;
    // }

    h_ratio->SetStats(0);
    h_ratio->SetLineColor(kBlue);
    h_ratio->SetLineWidth(2);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlue);
    // h_ratio->SetFillColor(kBlue-6);
    // h_ratio->SetFillStyle(3018);


    TH1D* h_ratio2 = (TH1D*) hbkg->Clone("h_ratio2");
    h_ratio2->SetDirectory(nullptr);
    h_ratio2->Divide(hbkg_true);
    // for (int ibin = 1; ibin <= h_ratio2->GetNbinsX() ; ibin++)
    // {
    //     h_ratio2->SetBinError(ibin, hbkg->GetBinError(ibin)/hbkg_true->GetBinContent(ibin));
    //     cout << "bkg hist bin error = "<<  hbkg->GetBinError(ibin)/hbkg_true->GetBinContent(ibin)  << endl;
    // }

    h_ratio2->SetStats(0);
    h_ratio2->SetLineColor(kCyan+2);
    h_ratio2->SetLineWidth(2);
    h_ratio2->SetMarkerStyle(4);
    h_ratio2->SetMarkerColor(kCyan+2);
    // h_ratio->SetFillColor(kCyan-6);
    // h_ratio->SetFillStyle(3018);

    // -- styles for the ratio plot 
    h_ratio->GetYaxis()->SetTitle("ratio");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.10);
    h_ratio->GetYaxis()->SetLabelSize(0.09);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    
       h_ratio->GetXaxis()->SetTitle("#DeltaR");
     // else { h_ratio->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");}
    h_ratio->GetXaxis()->SetTitleSize(0.12);
    h_ratio->GetXaxis()->SetLabelSize(0.10);
    // h_ratio->GetXaxis()->SetRangeUser(0,12);
    
    h_ratio->SetMinimum(0.1); //-5
    h_ratio->SetMaximum(1.9);// 10
    
    h_ratio->Draw("HIST PE1 L  ");
    h_ratio2->Draw("HIST PE1 L same");

    // Reference line at 1
    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw();
 

    TLegend *leg2 = new TLegend(0.406,0.7,0.7,0.9, ""); //get position from legend.C file (lower left corner is 0,0) // 0.406,0.39,0.7,0.55, ""
    leg2->SetTextSize(0.03);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.50);
    leg2->AddEntry(h_ratio, "sig ratio");
    leg2->AddEntry(h_ratio2, "bkg ratio");
    leg2->Draw("same");
    
    TString ptbin_name;
    // if(!isIntegDeltaR) 
    // ------------------------------- range has to be modified  -------------------------------
        {ptbin_name= Form("%g_%g", jtpt_binsVector[pt_bin-1], jtpt_binsVector[pt_bin]);}
    // else if (isIntegDeltaR) { cout << "HELLO   --- This is integarted bin "<< endl; ptbin_name= Form("IntegDeltaR_vspTbins");}
    c->Print( folder + sDirname + "/" + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".pdf");
    c->SaveAs( folder + sDirname + "/"  + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".png");




    // -- Write plotted canvas to output file 
     if (!foutputPlots) {
        std::cerr << "Invalid output file pointer!" << std::endl;
        return nullptr;
    }
    foutputPlots->cd();
    c->Write();
    foutputPlots->Write();


    file->Close(); delete file;
    return c;
}


