
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
TString sDirname = "TemplateFitOutput_Rebinned"; // Must use before Help_Functions.h
// for test rebinning: add _Rebinned


std::unique_ptr<TCanvas> draw_template_fit_result(
    TString fout_name,
    TFile* foutputPlots,
    TString &dataset,
    TString &folder,
    TString &pT_selection,
     Int_t pt_bin = 0);

// -----------------------
TH3D* CutAndRebinY(TH3D* h3, double yMin, int newNyBins, const double* newyBins, const double* xBins, const double* zBins)
{
    /* This fucntion takes TH3D input hitograms with axis y: dr that has binning of dr from ~0-1, and 
        convert it into TH3D that start at lower edge: yMin (without doing underflow): it is cut
        then rebin the other bins based on the given array yBins.
        Assuming the array is already at the edges of the old binnings.
    */

    if (!h3 || !newyBins || newNyBins < 2) return nullptr;

    // --- Y-axis cut of the old hist
    TAxis* yAxis = h3->GetYaxis();
    int yMinBin = yAxis->FindBin(yMin); // Find bin take it even if it is  the lower edge of next bin! 
    int yMaxBin = yAxis->GetNbins();

    cout << "yMinBin = " << yMinBin << endl;

    // --- Build new histogram with variable Y binning
    TH3D* hNew = new TH3D(
        Form("%s_rebinnedY", h3->GetName()),
        h3->GetTitle(),
        // X axis (kept unchanged)
        h3->GetXaxis()->GetNbins(),
        xBins,  
        // Y axis (variable binning)
        newNyBins,
        newyBins,
        // Z axis (kept unchanged)
        h3->GetZaxis()->GetNbins(),
        zBins
    );
    // hNew->Reset();
    hNew->Sumw2();

    // --- Fill new histogram
    for (int ix = 1; ix <= h3->GetXaxis()->GetNbins(); ix++) {
        for (int iy = yMinBin; iy <= yMaxBin; iy++) {
            for (int iz = 1; iz <= h3->GetZaxis()->GetNbins(); iz++) {

                double content = h3->GetBinContent(ix, iy, iz);
                double error   = h3->GetBinError(ix, iy, iz);

                if (content == 0 && error == 0) continue;
                // test using lowedge intead of bin center for x and z 
                double x = h3->GetXaxis()->GetBinLowEdge(ix) + 1e-03; //  + 1e-12
                double z = h3->GetZaxis()->GetBinLowEdge(iz) + 1e-03;
                double y = h3->GetYaxis()->GetBinLowEdge(iy) + 1e-03; // use low edge instead of center, and add this epsilon to avoid the abiguity of bin edge (we use up to 9 digits)


                // The new bins 
                int xbin =  hNew->GetXaxis()->FindBin(x);
                int zbin =  hNew->GetZaxis()->FindBin(z);
                int ybin = -1;

                // USe global bin number for add content 
                int newBin;

                // add overflow old bins to the last new bin
                if (y >= newyBins[newNyBins]){
                    ybin = newNyBins;
                }
                else {
                    ybin =  hNew->GetYaxis()->FindBin(y);  
                }
                // global bin 
                newBin =  hNew->GetBin(xbin, ybin, zbin);
                // accumulate content
                hNew->AddBinContent(newBin, content);

                // combine errors in quadrature
                double oldErr = hNew->GetBinError(newBin);
                hNew->SetBinError(newBin,
                std::sqrt(oldErr * oldErr + error * error));
            
                // -- debug 
                // if (iy == SOME_BIN && ix == SOME_BIN && iz == SOME_BIN) {
                //         cout << "EVENT: x=" << x << " y=" << y << " z=" << z << endl;
                //         cout << " -> xbin=" << xbin
                //              << " ybin=" << ybin
                //              << " zbin=" << zbin << endl;
                //     }

            }
        }
    }

    return hNew;
}


// TH3D* CutMergeAndRebinY(TH3D* h3,
//                         const double* xBins,
//                         const double* newyBins,
//                         const double* zBins,
//                         int newNyBins)
// {
//     if (!h3 || !newyBins || newNyBins < 2) return nullptr;

//     TAxis* yA = h3->GetYaxis();

//     // =====================================================
//     // STEP 1: CUT
//     // =====================================================
//     double yCut = 0.004079;

//     int startBin = 1;
//     for (int i = 1; i <= yA->GetNbins(); i++) {
//         if (yA->GetBinUpEdge(i) > yCut) {
//             startBin = i;
//             break;
//         }
//     }

//     int endBin = yA->GetNbins();

//     // =====================================================
//     // STEP 2: CREATE FINAL HIST DIRECTLY (NO ROOT REBIN)
//     // =====================================================
//     TH3D* hFinal = new TH3D(
//         "hFinal",
//         h3->GetTitle(),

//         h3->GetNbinsX(),
//         xBins,

//         newNyBins,
//         newyBins,

//         h3->GetNbinsZ(),
//         zBins
//     );

//     hFinal->Sumw2();

//     TAxis* yNew = hFinal->GetYaxis();

//     // =====================================================
//     // STEP 3: FILL DIRECTLY INTO FINAL HIST
//     // =====================================================
//     for (int ix = 1; ix <= h3->GetNbinsX(); ix++) {
//         for (int iz = 1; iz <= h3->GetNbinsZ(); iz++) {

//             for (int iy = startBin; iy <= endBin; iy++) {

//                 double c = h3->GetBinContent(ix, iy, iz);
//                 double e = h3->GetBinError(ix, iy, iz);

//                 if (c == 0 && e == 0) continue;

//                 double yLow  = yA->GetBinLowEdge(iy);
//                 double yHigh = yA->GetBinUpEdge(iy);

//                 // distribute into variable bins
//                 for (int j = 1; j <= newNyBins; j++) {

//                     double nLow  = yNew->GetBinLowEdge(j);
//                     double nHigh = yNew->GetBinUpEdge(j);

//                     double overlap = std::min(yHigh, nHigh)
//                                    - std::max(yLow, nLow);

//                     if (overlap <= 0) continue;

//                     double frac = overlap / (yHigh - yLow);

//                     int bin = hFinal->GetBin(ix, j, iz);

//                     hFinal->AddBinContent(bin, c * frac);

//                     hFinal->SetBinError(bin,
//                         std::sqrt(
//                             std::pow(hFinal->GetBinError(bin),2) +
//                             std::pow(e * frac,2)
//                         )
//                     );
//                 }
//             }
//         }
//     }

//     return hFinal;
// }


// 
// //Rebin the dr axis to make the bins start later
// void FuncRebin_dr(TH3D* &h){

//     Int_t first_bin = 3;
//     Int_t bins_dr = h->GetNbinsY();
//     Int_t bins_mb = h->GetNbinsX();
//     Int_t bins_pt = h->GetNbinsZ();

//     for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
//         for(Int_t imb_bin = 1; imb_bin <= bins_mb; imb_bin++){
//             Float_t first_bin_content = 0;
//             // for(Int_t ibin_dr = 1; ibin_dr <= first_bin; ibin_dr++){
//             //     first_bin_content += h->GetBinContent(imb_bin, ibin_dr, ibin_pt);
//             // }

//             first_bin_content  = h->GetBinContent(imb_bin, first_bin, ibin_pt); // cut 
//             h->GetYaxis()->SetRange(first_bin, bins_dr);
//             h->SetBinContent(imb_bin, first_bin, ibin_pt, first_bin_content);
//         }
//     }
// }

// -----------------------
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

// -----------------------
void DrawCommonTextTopRight(TPad*pad,  int ibin_dr, int ibin_pt, const double* newyBins, bool useDeaultLegend =true,const char* extra="") {
    /*Draw dr range, pt range + default Build Legend*/    
    // you can pass dirctly a canvas pointer 
    // Pass ybins array: to include changes when rebinning dr axis

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
        else { dr_first =  newyBins[ibin_dr -1]; dr_last  = newyBins[ibin_dr];}

    if (!ibin_pt){pt_first = jtpt_binsVector[0];  pt_last = jtpt_binsVector[jtpt_bins];} 
        else{pt_first = jtpt_binsVector[ibin_pt-1]; pt_last = jtpt_binsVector[ibin_pt]; }

        // cout << "pt bin #"<< ibin_pt << "pt first and last are : " << pt_first << ", " << pt_last << endl;

    if (dr_last == 1.0) latex.DrawLatex(x, y, Form("%g < #DeltaR < #infty", dr_first));
    else latex.DrawLatex(x, y, Form("%g < #DeltaR < %g", dr_first, dr_last));
    latex.DrawLatex(x, y - 0.04, Form("%g < p_{T} < %g GeV", pt_first, pt_last));
    latex.DrawLatex(x, y - 0.08, extra);
}

// -----------------------
void AddRatioPlot(TH1* h1, TH1* h2, const char* drawoption="EP", int linecolor = kBlack) {
    // draw ratio of two hists
    TH1* ratio = (TH1*)h1->Clone(Form("ratio_%s_%s", h1->GetName(), h2->GetName() ));
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
    ratio->SetLineColor(linecolor);
    ratio->Draw(drawoption);

    // unity line
    TLine* line = new TLine(
        ratio->GetXaxis()->GetXmin(), 1.0,
        ratio->GetXaxis()->GetXmax(), 1.0
    );
    line->SetLineStyle(2);
    line->Draw();
}

// -----------------------
//Draws the result of the template fit
std::unique_ptr<TCanvas> draw_template_fit_result(
    TString fout_name,
    TFile* foutputPlots,
    TString &dataset,
    TString &folder,
    TString &pT_selection,
    Int_t ibin_pt,
    bool rebin_dr = true)
{
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

    //background fraction
    TH2D *hbkg_2D = (TH2D*)file->Get("h_bkg_fraction"); hbkg_2D->SetDirectory(nullptr);
        hbkg = (TH1D*)hbkg_2D->ProjectionX("hbkg", ibin_pt+1, ibin_pt+1);
    TH2D *hbkg_true_2D = (TH2D*)file->Get("h_bkg_frac_true");   hbkg_true_2D->SetDirectory(nullptr);
        hbkg_true = (TH1D*)hbkg_true_2D->ProjectionX("hbkg_true",ibin_pt+1, ibin_pt+1);

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
        c->SetTitle("Bin1 for full dr range, then differential bins");

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

    TLegend *leg = new TLegend(0.32,0.39, 0.61,0.55, ""); //get position from legend.C file (lower left corner is 0,0)
    // 0.406 ,0.39,0.7,0.55, ""
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(hbkg_true, "Bkg. frac. (MC)");
    leg->AddEntry(hbkg, "Bkg. frac. (Data)");
    leg->AddEntry(htrue, "Sig. frac. (MC)");
    leg->AddEntry(h, "Sig. frac. (Data)");
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
    h_ratio->GetYaxis()->SetTitle("Data/MC");
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
    // ------------------------------- range has to be modified  -------------------------------
    if (ibin_pt > 0){ ptbin_name = Form("%g_%g", jtpt_binsVector[ibin_pt-1], jtpt_binsVector[ibin_pt]);}
    else ptbin_name = Form("%g_%g", jtpt_binsVector[0], jtpt_binsVector[bins_pt]);
    c->SaveAs( folder + sresultDir + "/" + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".pdf");
    c->SaveAs( folder + sresultDir + "/"  + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".png");


    // ---------- 
    // -- Another canvas: for each pt, draw only differntial bins, without integarted dr result + axis is absolute Dr values
    auto c_drval = std::make_unique<TCanvas>(Form("drval_%s", sname_canvas.Data()),"Template fit result vs. dr intervals", 900, 1100);// 800 x 800
        gROOT->GetListOfCanvases()->Remove(c_drval.get()); // to save it later 
        c_drval->SetTitle("");
        // c_drval->SetLogx();

        // --- Pads
        TPad *pad1_val = new TPad("pad1_val", "pad1_val", 0, 0.30, 1, 1.0);
        TPad *pad2_val = new TPad("pad2_val", "pad2_val", 0, 0.00, 1, 0.30);

        // -- convert bin number to absolute dr values 
        // -- signal
        Int_t N_dr_bins; // 
        const double* binsvector = nullptr;
        if (rebin_dr) { binsvector = dr_binsVector_wider; N_dr_bins = bins_dr_wider;}
        else { binsvector = dr_binsVector; N_dr_bins = bins_dr;}


        TH1D* h_dr = new TH1D("h_dr", "h_dr", N_dr_bins, binsvector); 
        h_dr->GetXaxis() ->SetTitle("#DeltaR"); 
        h_dr->Reset();
        cout << "\n old histogram binning with integrated dr:  #bins = " << h->GetNbinsX()<< endl;
        cout << " and the new dr axis without integarted dr: #bins = " << h_dr->GetNbinsX()<< endl;       
        
        // - and the rest of histograms 
        TH1D* htrue_dr = (TH1D*) h_dr->Clone("htrue_dr");htrue_dr->Reset();
        TH1D* hbkg_dr = (TH1D*) h_dr->Clone("hbkg_dr");hbkg_dr->Reset();        
        TH1D* hbkg_true_dr = (TH1D*) h_dr->Clone("hbkg_true_dr");hbkg_true_dr->Reset();

        // change the axes here -----------------------
        for (int i = 1; i <= N_dr_bins; i++)
        {
            h_dr->SetBinContent(i, h->GetBinContent(i+1));
            h_dr->SetBinError(i, h->GetBinError(i+1));

            htrue_dr->SetBinContent(i, htrue ->GetBinContent(i+1));
            htrue_dr->SetBinError(i, htrue ->GetBinError(i+1));

            hbkg_dr->SetBinContent(i, hbkg ->GetBinContent(i+1));
            hbkg_dr->SetBinError(i, hbkg ->GetBinError(i+1));

            hbkg_true_dr ->SetBinContent(i, hbkg_true ->GetBinContent(i+1));
            hbkg_true_dr ->SetBinError(i, hbkg_true ->GetBinError(i+1));
        }

        h_dr->SetDirectory(nullptr);
        htrue_dr->SetDirectory(nullptr);
        hbkg_dr->SetDirectory(nullptr);
        hbkg_true_dr->SetDirectory(nullptr);

    // -- Set styles 
    h_dr ->SetLineColor(kRed +2);
    htrue_dr ->SetLineColor(kBlue +1);
    hbkg_dr ->SetLineColor(kGreen +2);
    hbkg_true_dr ->SetLineColor(kMagenta +2);

    h_dr ->SetLineWidth(2);
    htrue_dr ->SetLineWidth(2);
    hbkg_dr ->SetLineWidth(2);
    hbkg_true_dr ->SetLineWidth(2);

        // -- Draw the new histograms of the new axis 
         // --  set pads margins 
            pad1_val->SetBottomMargin(0.01);
            pad1_val->SetLeftMargin(0.12);
            pad1_val->SetTopMargin(0.08);
            
            pad2_val->SetTopMargin(0.05);
            pad2_val->SetBottomMargin(0.35);
            pad2_val->SetLeftMargin(0.12);
            pad1_val->Draw();
            pad2_val->Draw();
            // 
            pad1_val ->cd();
            pad1_val->SetLogx();

                h_dr->SetStats(0);
                h_dr->SetTitle("");
                h_dr->GetYaxis()->SetRangeUser(0,1);
                h_dr->GetXaxis()->SetTitle("");
                h_dr->GetYaxis()->SetTitle("Signal or Background fraction");
                h_dr->GetYaxis()->CenterTitle(true);
                h_dr->Draw("Hist E");
                htrue_dr->Draw("Hist E same");
                hbkg_dr->Draw("Hist E same");
                hbkg_true_dr->Draw("Hist E same");
                TLegend *leg_val = new TLegend(0.32,0.39, 0.61,0.55, "");
                    leg_val->SetTextSize(0.03);
                    leg_val->SetFillStyle(0);
                    leg_val->SetBorderSize(0);
                    leg_val->SetMargin(0.50);
                    leg_val->AddEntry(hbkg_true_dr, "Bkg. frac. (MC)");
                    leg_val->AddEntry(hbkg_dr, "Bkg. frac. (Data)");
                    leg_val->AddEntry(htrue_dr, "Sig. frac. (MC)");
                    leg_val->AddEntry(h_dr, "Sig. frac. (Data)");
                    leg_val->SetHeader(Form("%g < p_{T} < %g GeV", pt_first, pt_last), "C"); //defined before
                    leg_val->Draw("same");

            pad2_val->cd();
            pad2_val->SetLogx();

            TH1D* h_ratio_dr = (TH1D*) h_dr->Clone("h_ratio_dr"); h_ratio_dr->Reset();
            h_ratio_dr->Divide(h_dr, htrue_dr);
            h_ratio_dr->SetDirectory(nullptr);
            h_ratio_dr->SetLineColor(kBlue + 1);
            h_ratio_dr->SetLineWidth(2);
            h_ratio_dr->SetMarkerStyle(20);
            h_ratio_dr->SetMarkerColor(kBlue + 1);

            TH1D* h_ratio2_dr = (TH1D*) hbkg_dr->Clone("h_ratio2_dr"); h_ratio2_dr->Reset();
            h_ratio2_dr->Divide(hbkg_dr, hbkg_true_dr);
            h_ratio2_dr->SetDirectory(nullptr);
            h_ratio2_dr->SetLineColor(kMagenta+1);
            h_ratio2_dr->SetLineWidth(2);
            h_ratio2_dr->SetMarkerStyle(4);
            h_ratio2_dr->SetMarkerColor(kMagenta+1);

            // -- styles for the ratio plot 
            h_ratio_dr->GetYaxis()->SetTitle("Data/MC");
            h_ratio_dr->GetYaxis()->SetNdivisions(505);
            h_ratio_dr->GetYaxis()->SetTitleSize(0.10);
            h_ratio_dr->GetYaxis()->SetLabelSize(0.09);
            h_ratio_dr->GetYaxis()->SetTitleOffset(0.5);
            h_ratio_dr->GetXaxis()->SetTitle("#DeltaR");
            h_ratio_dr->GetXaxis()->SetTitleSize(0.12);
            h_ratio_dr->GetXaxis()->SetLabelSize(0.10);

            if (h_ratio_dr->GetMaximum() > 10) max = 5; else max = h_ratio_dr->GetMaximum() + 0.2; // avoid very very large values when fit is almost zero.
            if (h_ratio_dr->GetMinimum() < 0.01) min = 0; else min = h_ratio_dr->GetMinimum()- 0.2;
            h_ratio_dr->SetMaximum( max ); // max
            h_ratio_dr->SetMinimum( min); // min
            // Draw lower Pad 
            h_ratio_dr->Draw("HIST E  ");
            h_ratio2_dr->Draw("HIST E  same");

            // Reference line at 1
            TLine *line_val = new TLine(h_ratio_dr->GetXaxis()->GetXmin(), 1.0, h_ratio_dr->GetXaxis()->GetXmax(), 1.0);
            line_val->SetLineStyle(2);
            line_val->Draw();
         

            TLegend *leg2_val = new TLegend(0.7, 0.7, 0.85, 0.9, "");
            leg2_val->SetTextSize(0.05);
            leg2_val->SetFillStyle(0);
            leg2_val->SetBorderSize(0);
            leg2_val->SetMargin(0.50);
            leg2_val->AddEntry(h_ratio_dr, "Sig.");
            leg2_val->AddEntry(h_ratio2_dr, "Bkg.");
            leg2_val->Draw("same");

            c_drval->SaveAs( folder + sresultDir + "/" + trivialMC_label + "sign_frac_result_differentialDronly_" + dataset + "_" + ptbin_name + ".png");
                c_drval->SaveAs( folder + sresultDir + "/" + trivialMC_label + "sign_frac_result_differentialDronly_" + dataset + "_" + ptbin_name + ".pdf");


    // -- Write plotted canvas to output file 
     if (!foutputPlots) {
        std::cerr << "Invalid output file pointer!" << std::endl;
        return nullptr;
    }
    foutputPlots->cd();
    c->Write();
    c_drval->Write();
    foutputPlots->Write();


    file->Close(); delete file;
    return c;
}


