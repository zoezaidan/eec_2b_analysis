
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

// C++ includes 
#include <vector>
    #include <TLegend.h>
#include <TObject.h>
#include <vector>
#include <string>


// name for outpur directory (in path) 
TString sDirname = "TemplateFitOutput"; // has to be before Help_Functions.h

std::unique_ptr<TCanvas> draw_template_fit_result(
    TString fout_name,
    TFile* foutputPlots,
    TString &dataset,
    TString &folder,
    TString &pT_selection,
     Int_t pt_bin = 0);

TLegend* CreateLegend(
    double x1, double y1, double x2, double y2,
    const std::vector<TObject*>& objects,
    const std::vector<std::string>& options,
    const std::vector<std::string>& labels = {})
{
    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (size_t i = 0; i < objects.size(); ++i) {
        TObject* obj = objects[i];

        // decide label
        const char* label = nullptr;
        if (!labels.empty() && i < labels.size() && labels[i] != "") {
            label = labels[i].c_str();  // custom label
        } else {
            label = obj->GetTitle();   // fallback to object title
        }

        const char* opt = (i < options.size()) ? options[i].c_str() : "l";

        leg->AddEntry(obj, label, opt);
    }

    return leg;
}


void DrawCommonTextTopRight(TPad*pad,  int ibin_dr, int ibin_pt, bool useDeaultLegend =true,const char* extra="") {
    /*Draw dr range, pt range + default Build Legend*/    
    // you can pass dirctly a canvas pointer 
    pad->cd(); 
    // draw auto legend first
    TLegend* leg = nullptr;
    if(useDeaultLegend){
        leg = pad->BuildLegend(0.59, 0.5, 0.87, 0.75); // 0.6, 0.5, 0.88, 0.7
        // optional styling
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
    }

    // -------------- pt dr legend ----------------------------
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

        // cout << "pt bin #"<< ibin_pt << "pt first and last are : " << pt_first << ", " << pt_last << endl;

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
    TString &pT_selection, Int_t ibin_pt)
{

    // ----------------------------------
    // --- WORK IN PROGRESS -------------
    // ----------------------------------


    // ibin_pt is kept with convenstion starting from 0 (for integarted interval, then then the rest) 

    // For trvial tests 
    TString trivialMC_label = fout_name.Contains("trivialMC") ? "_trivialMC": "";

    // -- general style 
    gStyle->SetOptStat(0);

    // -- output directory 
    TString sresultDir = Form("%s/FitResult_Summary_S_B_fractions", sDirname.Data());// aprent already exist
    gSystem->mkdir(sresultDir, kTRUE); 

    //Get fractions for the jtpt bin ibin_pt
    TFile *file = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "read");
    if (!file) {Error("Input File:", "File does not exist'%s'", file->GetName());return nullptr;}
    cout << "input file name " << file->GetName() << endl;

    // -- Read output of template fit files 
    TH1D *h;
    TH1D *htrue;
    TH1D *hbkg;
    TH1D *hbkg_true;
    TString sname_canvas;
        sname_canvas = Form("c_%s_ptbin_%d%s", dataset.Data(), ibin_pt,  trivialMC_label.Data());

    // -- Read input TH2D  summary of signa and background fractions 
    //signal fraction
    // Fit result (uncertaintiy already at the point)
    TH2D *h_2D = (TH2D*)file->Get("h_sig_fraction"); h_2D->SetDirectory(nullptr);
        h = (TH1D*) h_2D->ProjectionX("h", ibin_pt+1 , ibin_pt+1);
        h->GetXaxis()->SetTitle("Bin (\\#Delta\\ r)");
    // True fractions (uncertaintiy already at the point)
    TH2D *htrue_2D = (TH2D*)file->Get("h_sig_frac_true"); htrue_2D->SetDirectory(nullptr);
        htrue = (TH1D*)htrue_2D->ProjectionX("htrue", ibin_pt+1 , ibin_pt+1);
        
        // and its uncertainity 
        // TH2D *htrue_err_2D = (TH2D*)file->Get("h_sig_frac_true_error");
        // TH1D* htrue_err = (TH1D*) htrue_err_2D->ProjectionX("htrue_err", ibin_pt+1 , ibin_pt+1);
        // for (int ibin = 1; ibin <= htrue_err->GetNbinsX(); ibin++)
        // {
        //     htrue->SetBinError(ibin, htrue_err->GetBinContent(ibin));
        // }

    //background fraction
    TH2D *hbkg_2D = (TH2D*)file->Get("h_bkg_fraction"); hbkg_2D->SetDirectory(nullptr);
        hbkg = (TH1D*)hbkg_2D->ProjectionX("hbkg", ibin_pt+1, ibin_pt+1);
    TH2D *hbkg_true_2D = (TH2D*)file->Get("h_bkg_frac_true");   hbkg_true_2D->SetDirectory(nullptr);
        hbkg_true = (TH1D*)hbkg_true_2D->ProjectionX("hbkg_true",ibin_pt+1, ibin_pt+1);
        // // and its uncertainity 
        // TH2D *hbkg_true_err_2D = (TH2D*)file->Get("h_bkg_frac_true_error");
        // TH1D* hbkg_true_err = (TH1D*) hbkg_true_err_2D->ProjectionX("hbkg_true_err", ibin_pt,ibin_pt);
        // for (int ibin = 1; ibin <= hbkg_true_err->GetNbinsX(); ibin++)
        // {
        //     /// Assign the uncertainity to the iD hist of true signal Vs. DeltaR 
        //     hbkg_true->SetBinError(ibin, hbkg_true_err->GetBinContent(ibin));
        // }

    // -- Label vertically dr bins axis 


    // -- deattach hists from input root file: for drawing 
    h->SetDirectory(nullptr);
    htrue->SetDirectory(nullptr);
    hbkg->SetDirectory(nullptr);
    hbkg_true->SetDirectory(nullptr);

    // -- Set styles 
    h ->SetLineColor(kRed +2);
    htrue ->SetLineColor(kBlue +1);
    hbkg ->SetLineColor(kGreen +2);
    hbkg_true ->SetLineColor(kMagenta +2);

    h ->SetLineWidth(2);
    htrue ->SetLineWidth(2);
    hbkg ->SetLineWidth(2);
    hbkg_true ->SetLineWidth(2);

    // -- Output summary canvas
    auto c = std::make_unique<TCanvas>(sname_canvas,"Template fit result", 900, 1100);// 800 x 800
        gROOT->GetListOfCanvases()->Remove(c.get()); // to save it later 
        c->SetTitle("");

    // --- Pads
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);
    // -- dettach pads from root ownership
    pad1->SetBit(kCanDelete, false); // not sure if neeed 
    pad2->SetBit(kCanDelete, false); // not sure if neeed 

    // --  set pads margins 
    pad1->SetBottomMargin(0.01);
    pad1->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.08);
    
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.12);
    pad1->Draw();
    pad2->Draw();
    
    //--- Draw on pad 1 
    pad1->cd();
    //Draw results
        h->SetStats(0);
        h->SetTitle("");
        h->GetYaxis()->SetRangeUser(0,1);
        // h->GetXaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("Signal or Background fraction");
        h->GetYaxis()->CenterTitle(true);
        h->Draw("Hist E");
        htrue->Draw("Hist E same");
        hbkg->Draw("Hist E same");
        hbkg_true->Draw("Hist E same");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    // if(!isIntegDeltaR) test_info_text->DrawLatex(0.15, 0.5, Form("%g < p_{T} < %g GeV", jtibin_ptsVector[ibin_pt-1], jtibin_ptsVector[ibin_pt]));
    // else{ 
    //     TString srangeDeltaR = Form("#DeltaR [Full range]");
    //     test_info_text->DrawLatex(0.15, 0.5, srangeDeltaR);
    // }
    // test_info_text->DrawLatex(0.15, 0.45, dataset);
    // test_info_text->DrawLatex(0.15, 0.4, "Template from MC bjet");
    test_info_text->Draw("SAME");


    TLegend *leg = new TLegend(0.32,0.39, 0.61,0.55, ""); //get position from legend.C file (lower left corner is 0,0)
    // 0.406 ,0.39,0.7,0.55, ""
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(hbkg_true, "Bkg. frac. (True)");
    leg->AddEntry(hbkg, "Bkg. frac. (Fit)");
    leg->AddEntry(htrue, "Sig. frac. (True)");
    leg->AddEntry(h, "Sig. frac. (Fit)");
    // Add pt range  in leg title 
    double pt_first = 0, pt_last = 0;
        if (!ibin_pt){pt_first = jtpt_binsVector[0];  pt_last = jtpt_binsVector[jtpt_bins];} 
        else{pt_first = jtpt_binsVector[ibin_pt-1]; pt_last = jtpt_binsVector[ibin_pt]; }
        leg->SetHeader(Form("%g < p_{T} < %g GeV", pt_first, pt_last), "C"); //centered 
    leg->Draw("same");

    // ----------------------------------------------------------

    pad2->cd();
    // Create ratio histogram
    TH1D* h_ratio = (TH1D*) h->Clone("h_ratio"); h_ratio->Reset();
    h_ratio->Sumw2();
    h_ratio->Divide(h, htrue);
    h_ratio->SetDirectory(nullptr);
    // ratio color is same as the true one 
    h_ratio->SetLineColor(kBlue + 1);
    h_ratio->SetLineWidth(2);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlue + 1);

    TH1D* h_ratio2 = (TH1D*) hbkg->Clone("h_ratio2"); h_ratio2->Reset();
    h_ratio2->Sumw2();
    h_ratio2->Divide(hbkg, hbkg_true);
    h_ratio2->SetDirectory(nullptr);
    h_ratio2->SetLineColor(kMagenta+1);
    h_ratio2->SetLineWidth(2);
    h_ratio2->SetMarkerStyle(4);
    h_ratio2->SetMarkerColor(kMagenta+1);

    // -- styles for the ratio plot 
    h_ratio->GetYaxis()->SetTitle("Fit/True");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.10);
    h_ratio->GetYaxis()->SetLabelSize(0.09);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->GetXaxis()->SetTitle("Bin(#DeltaR)");
    h_ratio->GetXaxis()->SetTitleSize(0.12);
    h_ratio->GetXaxis()->SetLabelSize(0.10);
    
    // -- Set pad 2 y axis limit 
    double max = 0, min = 0; 
    if (h_ratio->GetMaximum() > 10) max = 5; else max = h_ratio->GetMaximum() + 0.2; // avoid very very large values when fit is almost zero.
    if (h_ratio->GetMinimum() < 0.01) min = 0; else min = h_ratio->GetMinimum()- 0.2;
    h_ratio->SetMaximum( max ); // max
    h_ratio->SetMinimum( min); // min

    // Draw lower Pad 
    h_ratio->Draw("HIST E  ");
    h_ratio2->Draw("HIST E  same");

    // Reference line at 1
    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw();
 

    TLegend *leg2 = new TLegend(0.7, 0.7, 0.85, 0.9, "");
     //get position from legend.C file (lower left corner is 0,0) // 0.406,0.39,0.7,0.55, ""
    leg2->SetTextSize(0.05);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.50);
    leg2->AddEntry(h_ratio, "Sig.");
    leg2->AddEntry(h_ratio2, "Bkg.");
    leg2->Draw("same");
    
    TString ptbin_name;
    // if(!isIntegDeltaR) 
    // ------------------------------- range has to be modified  -------------------------------
    if (ibin_pt > 0){ ptbin_name = Form("%g_%g", jtpt_binsVector[ibin_pt-1], jtpt_binsVector[ibin_pt]);}
    else ptbin_name = Form("%g_%g", jtpt_binsVector[0], jtpt_binsVector[bins_pt]);

    // else if (isIntegDeltaR) { cout << "HELLO   --- This is integarted bin "<< endl; ptbin_name= Form("IntegDeltaR_vspTbins");}
    c->Print( folder + sresultDir + "/" + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".pdf");
    c->SaveAs( folder + sresultDir + "/"  + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".png");


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


