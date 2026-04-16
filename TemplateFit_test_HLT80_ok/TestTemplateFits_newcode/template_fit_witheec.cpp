#include "tTree.h"
#include "TFile.h"
#include "binning_histos_all.h"

// For the output to be saved: 
TString sDirname = "TemplateFitOutput_HLT80_HGdata_bjet_witheec";

//Rebin the mB axis to have a bigger overflow bin
void rebin_mass(TH3D* &h){
    Int_t bins_mb = h->GetNbinsX();
    Int_t bins_pt = h->GetNbinsZ();
    Int_t bins_dr = h->GetNbinsY();

    Int_t last_mb_bin = 12;

    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){
            Float_t last_bin = 0;
            for(Int_t imb_bin = last_mb_bin; imb_bin <= bins_mb; imb_bin++){
                last_bin += h->GetBinContent(imb_bin, ibin_dr, ibin_pt);
            }
            h->GetXaxis()->SetRange(1, last_mb_bin);
            h->SetBinContent(last_mb_bin, ibin_dr, ibin_pt);
        }
    } 
}

//Rebin the dr axis to make the bins start later
void rebin_dr(TH3D* &h){

    Int_t first_bin = 3;
    Int_t bins_dr = h->GetNbinsY();
    Int_t bins_mb = h->GetNbinsX();
    Int_t bins_pt = h->GetNbinsZ();

    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t imb_bin = 1; imb_bin <= bins_mb; imb_bin++){
            Float_t first_bin_content = 0;
            for(Int_t ibin_dr = 1; ibin_dr <= first_bin; ibin_dr++){
                first_bin_content += h->GetBinContent(imb_bin, ibin_dr, ibin_pt);
            }
            h->GetYaxis()->SetRange(first_bin, bins_dr);
            h->SetBinContent(imb_bin, first_bin, ibin_pt, first_bin_content);
            
            
        }
    }
}

void do_template_fit (TString &dataset_HighEG, TString &dataset_LowEG, bool& alsoLowEG, TString &templates_qcd, 
                      TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name,  
                      TString dataLabel = "Data" ){
    
// Dataset: for now HGtrigered data. 
//-- Samples: qcd + bjet 
// -- Template fit pdfs : includes (0B, 1B) as bkg (added effectively), 2B as sig  

    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");
    TString pdfBase = fout_name; pdfBase.ReplaceAll(".root", "");
    TString pdfName = Form("%s/%s.pdf", sDirname.Data(), pdfBase.Data());

    // -- Using new naming convention
    TString namehData = "h3D_data"; // h3D_data with eec on the weight 

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B  = "h3D_bb";
    TString nameh1B = "h3D_b";
    TString nameh0B = "h3D_0b";
   
    //-- QCD sample:
    TFile *file_qcd= new TFile(folder + templates_qcd, "read");
    if (!file_qcd || file_qcd->IsZombie()) {Error("Input File:", "File does not exist: %s", (folder + templates_qcd).Data()); return;}

    TH3D *h3D_b = (TH3D*)file_qcd->Get(nameh1B)->Clone("h3D_b");
    if(!h3D_b) {Error("Get:", "histogram does not exist: h3D_b"); return;}
    
    TH3D *h3D_bb = (TH3D*)file_qcd->Get(nameh2B)->Clone("h3D_bb");
    if(!h3D_bb){Error("Get:", "histogram does not exist: h3D_bb"); return;}
    
    TH3D *h3D_nob = (TH3D*)file_qcd->Get(nameh0B)->Clone(nameh0B);
    if(!h3D_nob){Error("Get:", "histogram does not exist: h3D_nob"); return;}

    //-- Bjet sample: Templates: Open file and get histograms from bjet: templates 
    TFile *file_bjet = new TFile(templates_bjet, "read");
    if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}

    TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
    if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
    
    TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
    if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}
    
   
    //Open dataset:
    TFile *file_data_HighEG = new TFile(dataset_HighEG, "read");
    if (!file_data_HighEG) {Error("Input File:", "File does not exist'%s'", file_data_HighEG->GetName());return;}
    cout << "file HighEG name " << file_data_HighEG->GetName() << endl;
    
    TH3D *h3D_data;
    h3D_data= (TH3D*)file_data_HighEG->Get(namehData)->Clone("h3D_data");
    if(!h3D_data){Error("Get:", "histogram does not exist '%s' ", h3D_data->GetName()); return;}

    // -- Add also the lower triggerd data //ZOE CHECK IF THIS IS CORRECT
    if(alsoLowEG){
        TFile *file_data_LowEG = new TFile(dataset_LowEG, "read");
        if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
        cout << "file LowEG name " << file_data_LowEG->GetName() << endl;
        
        TH3D *h3D_dataLowEG;
        h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
        if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ", h3D_dataLowEG->GetName()); return;}
        
        h3D_data->Add(h3D_dataLowEG); 
        }

    h3D_data->SetMarkerStyle(20);
    h3D_data->SetMarkerColor(0);
    h3D_data->SetLineColor(0);

    //Get number of bin entries
    Int_t bins_pt = h3D_data->GetNbinsZ();
    Int_t bins_dr = h3D_data->GetNbinsY(); 
    Int_t mb_bins = h3D_data->GetNbinsX();
    cout << "-- Data hist binning:" << endl;
    cout << "pt bins = "<< bins_pt << endl;
    cout << "dr bins = "<< bins_dr << endl;
    cout << "mb bins = "<< mb_bins << endl;

    // -- Draw THstach for mb projection 
    // Get projections X 
    // hists names as for slices in bins (but outside {} so it is ok)
    TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX("h_data_mb", 1,bins_dr, 1, bins_pt ); // all bins
    h_data_mb->SetLineColor(kBlack);
    h_data_mb->SetLineWidth(2);
    h_data_mb->SetMarkerColor(kBlack);
    h_data_mb->SetTitle(dataLabel);
    h_data_mb->SetMarkerStyle(24); // open circle, matching slides
    //TH1D *h_b = (TH1D *) h3D_b->ProjectionX("h_b");
    TH1D *h_b_bjet = (TH1D *) h3D_b_bjet->ProjectionX("h_b_bjet");
    //TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX("h_bb");
    TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet->ProjectionX("h_bb_bjet");
    //TH1D *h_nob = (TH1D *) h3D_nob->ProjectionX("h_nob");
  
    // -- Set titles for the hisograms  to be used in the default legends 
    //h_b->SetTitle("1B (QCD)");
    
    h_b_bjet->SetTitle("1B (bjet)");
    //h_bb->SetTitle("2B (QCD)");
    h_bb_bjet->SetTitle("2B (bjet)");
    //h_nob->SetTitle("0B");
    
    

    
    // -- Deattach hists from root file 
    h_data_mb->SetDirectory(0);
    //h_b->SetDirectory(0);
    // h_bb->SetDirectory(0);
    h_b_bjet->SetDirectory(0);
    h_bb_bjet->SetDirectory(0);
    
    //h_nob->SetDirectory(0);
    // Set styles of hist to be drawn:
    //h_b->SetFillColor(kRed-10); h_b->SetLineColor(kRed-10);        
    //h_bb->SetFillColor(kBlue-10); h_bb->SetLineColor(kBlue-10);
    h_b_bjet->SetFillColor(kTeal-8);  h_b_bjet->SetLineColor(kTeal-6);   // 1B
    h_bb_bjet->SetFillColor(kRed-9);  h_bb_bjet->SetLineColor(kRed-3);   // 2B
   //  h_nob->SetFillColor(kRed); h_nob->SetLineColor(kRed);
    
    
    // =========================================================
    // PLOT 1 (pre-fit, non-normalized, integrated over all pT and deltaR)
    // Shows the raw m_2B distribution: stacked 1B + 2B templates from bjet MC
    // scaled to match the total data yield, overlaid with the data histogram.
    // Purpose: quick check that the template shapes and overall normalisation
    // are reasonable before running the fit.
    // =========================================================
    // Scale templates to match the total data integral
    double scale_to_data = (h_b_bjet->Integral(1, mb_bins) + h_bb_bjet->Integral(1, mb_bins)) > 0
    ? h_data_mb->Integral(1, mb_bins) / (h_b_bjet->Integral(1, mb_bins) + h_bb_bjet->Integral(1, mb_bins))
    : 1.0;
    TH1D* h_b_bjet_scaled  = (TH1D*) h_b_bjet->Clone("h_b_bjet_scaled");   h_b_bjet_scaled ->Scale(scale_to_data);
    TH1D* h_bb_bjet_scaled = (TH1D*) h_bb_bjet->Clone("h_bb_bjet_scaled");  h_bb_bjet_scaled->Scale(scale_to_data);
    
    THStack hstack_All("hstack_All","Mass stacked histogram");
    hstack_All.SetTitle("Before fit;m_{2B} [GeV];");
    hstack_All.Add(h_bb_bjet_scaled); // 2B first
    hstack_All.Add(h_b_bjet_scaled);  // 1B second
   
    
    auto canva = new TCanvas("All_templates_Data", "m_{2B}: data vs templates (non-normalized, scaled to data)", 800, 800);
    canva->cd();
    hstack_All.Draw("hist E");
    h_data_mb->Draw("HIST E same");
    gPad->Modified(); gPad->Update();
    canva->Modified(); canva->Update();
    fout->cd();
    canva->BuildLegend(0.6, 0.65, 0.9, 0.9);
    canva->Write();
    canva->Print(pdfName + "("); // open PDF — this is page 1

    
    // =========================================================
    // PLOT 2 (pre-fit, normalized, integrated over all pT and deltaR)
    
    // Same templates and data but each divided by their own total integral
    // so that all distributions have unit area.
    // 0B from QCD is also included in the stack (normalized to the same total).
    // Purpose: compare the shapes of 0B, 1B, 2B templates against data
    // independent of normalisation, to validate the template shapes.
    // =========================================================
    THStack hstack_All_norm("hstack_All_norm","Mass stacked histogram");
    hstack_All_norm.SetTitle("Before fit;m_{2B} [GeV];");
    // normalize all template components to the same denominator (1B + 2B integral)
    double total_int_1b_2b_bjet = h_b_bjet->Integral(1, mb_bins) + h_bb_bjet->Integral(1, mb_bins);
    TH1D* h_b_bjet_normstack  = (TH1D*) h_b_bjet->Clone("h_b_bjet_normstack");   h_b_bjet_normstack ->Scale(1./total_int_1b_2b_bjet);
    TH1D* h_bb_bjet_normstack = (TH1D*) h_bb_bjet->Clone("h_bb_bjet_normstack");  h_bb_bjet_normstack->Scale(1./total_int_1b_2b_bjet);
    // 0B from QCD, projected over all pT and deltaR bins, normalized to same denominator
    TH1D* h_nob_proj = (TH1D*) h3D_nob->ProjectionX("h_nob_proj", 1, bins_dr, 1, bins_pt);
    h_nob_proj->SetTitle("0B (QCD)");
    h_nob_proj->SetFillColor(kBlue-5);
    h_nob_proj->SetLineColor(kBlue-5); // 0B
    TH1D* h_nob_normstack = (TH1D*) h_nob_proj->Clone("h_nob_normstack");
    h_nob_normstack->Scale(1./total_int_1b_2b_bjet);
    // stack order: 0B at back, then 1B, then 2B on top
    hstack_All_norm.Add(h_bb_bjet_normstack);
    hstack_All_norm.Add(h_b_bjet_normstack);
    hstack_All_norm.Add(h_nob_normstack);
    // data self-normalized: divided by its own integral
    TH1D* hdata_selfnorm = (TH1D*) h_data_mb->Clone("hdata_selfnorm");
    hdata_selfnorm->Scale(1./h_data_mb->Integral(1, mb_bins));
    hdata_selfnorm->SetMarkerStyle(24); // open circle
    hdata_selfnorm->SetMarkerColor(kBlack); hdata_selfnorm->SetLineColor(kBlack); hdata_selfnorm->SetLineWidth(2);
    
    auto canva_stack_norm = new TCanvas("All_templates_Data_stacked_normalized", "m_{2B}: data vs templates (normalized)", 800, 800);
    canva_stack_norm->cd();
    hstack_All_norm.Draw("hist E");
    hdata_selfnorm->Draw("PE same");
    gPad->Modified(); gPad->Update();
    canva_stack_norm->Modified(); canva_stack_norm->Update();
    fout->cd();
    canva_stack_norm->BuildLegend(0.6, 0.65, 0.9, 0.9);
    canva_stack_norm->Write();
    canva_stack_norm->Print(pdfName); // page 2
    
    
    // return; // removed: allow fit loop to run
    /*
    // -- Combined QCD + bjet
    // Another plot for the distributions  before slicing, but combined samples 
    auto canvatotal = new TCanvas("All_templates_combined",Form(""), 800, 800 );
    canvatotal->cd();
    // canvatotal->SetLogy();
    // Get combined versions + their normalization
    TH1D* normhc_nob = (TH1D*) h_nob->Clone("normhc_nob"); normhc_nob->Scale(1./normhc_nob->Integral());
    TH1D* hc_1b = (TH1D*) h_b->Clone("hc_1b"); hc_1b->Add(h_b_bjet); hc_1b->SetTitle("1B");
    TH1D* normhc_1b = (TH1D*) hc_1b->Clone("normhc_1b");   normhc_1b->Scale(1./normhc_1b->Integral());
    TH1D* hc_2b = (TH1D*) h_bb ->Clone("hc_2b"); hc_2b->Add(h_bb_bjet); hc_2b->SetTitle("2B");
    TH1D* normhc_2b = (TH1D*) hc_2b->Clone("normhc_2b");   normhc_2b->Scale(1./normhc_2b->Integral());
    
    // - draw stack 
    THStack hstack_All_combined("hstack_All_combined","Mass stacked histogram");
    hstack_All_combined.SetTitle(";m_{2B} [GeV];");
    hstack_All_combined.Add(normhc_nob);
    hstack_All_combined.Add(normhc_2b);
    hstack_All_combined.Add(normhc_1b);
    hstack_All_combined.Draw("hist E"); // templates stack
    // -- normalizaed data 
    TH1D* normhc_data = (TH1D*) h_data_mb->Clone("hd_norm");
    normhc_data->Scale(1./normhc_data->Integral());
    normhc_data->Draw("HIST E same"); 
    normhc_data->SetTitle("data");
    canvatotal->SetTitle("Normalized distributions");         
    gPad->Modified();   
    gPad->Update();
    canvatotal->Modified();
    canvatotal->Update();
    fout->cd();
    canvatotal->BuildLegend(0.6, 0.65, 0.9, 0.9);
    canvatotal->Write();
    canvatotal->Print(Form("%s/%s.png", sDirname.Data(), canvatotal->GetName()));
    */
   
   // return;
   
   // -- Output hist declarations: empty histograms 
   // -- For the fit result: per (deltaR, jtpt) bin
   //Create the histograms storing the fit parameters in TH2D (the dr are stored in the X, the jtpt in the Y)
   TH2D *h_sig_fraction = (TH2D*) h3D_data->Project3D("zy")->Clone("h_sig_fraction");// just to get the axes from hist
   h_sig_fraction->SetTitle("x=#DeltaR, y=jtpt, fraction"); 
   h_sig_fraction->Reset();
   TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
   TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
   TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error");
   //Store the true parameters and errors
   TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
   TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");
   TH2D *h_bkg_frac_true = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true");
   TH2D *h_bkg_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true_error");
   
   // -- And for the integarted bin in deltaR, Fitparaemter Vs. ptbin 
   TH1D *h_sig_fraction_DeltaRIntBin = (TH1D*) h3D_data->ProjectionZ("h_sig_fraction_DeltaRIntBin"); h_sig_fraction_DeltaRIntBin->Reset();
   TH1D *h_bkg_fraction_DeltaRIntBin = (TH1D *) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_fraction_DeltaRIntBin");
   TH1D *h_sig_true_fraction_DeltaRIntBin = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_sig_true_fraction_DeltaRIntBin");
   TH1D *h_bkg_true_fraction_DeltaRIntBin = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_true_fraction_DeltaRIntBin");
   TH1D *h_sig_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_sig_fraction_DeltaRIntBin_error");
   TH1D *h_bkg_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_fraction_DeltaRIntBin_error");
   TH1D *h_sig_true_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_sig_true_fraction_DeltaRIntBin_error");
   TH1D *h_bkg_true_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_true_fraction_DeltaRIntBin_error");
   
   // --- Vector to test the convergence
   std::vector <std::pair<int, int>> non_converge_bins;
   
   
   // Fitting - loop over dr and jtpt entries
   for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
       // for(Int_t ibin_pt = 1; ibin_pt <= 1; ibin_pt++){
        // for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){
            
        // for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){
            for(Int_t ibin_dr = 0; ibin_dr <= 0; ibin_dr++){
                // Bin = 0: integaretd deltaR range, other bins are the DeltaR bins 
                
                // The slice range for projection is different for integarted Vs. one bin slice
                Int_t SliceFirstbin_dr = ibin_dr;
                Int_t SliceLastbin_dr =  ibin_dr;
                
                
                // For integarted bin only: bin = 0
                if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = bins_dr;}
                // -- test rebinned dr in 2 
                /*else{
                    
                SliceLastbin_dr =  ibin_dr + 1;
                ibin_dr++; // to make step of 2  
                }*/
               
               // Make slices for data
               TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
               h_data_mb->SetMarkerColor(kBlack);
               h_data_mb->SetMarkerStyle(20);
               h_data_mb->SetLineColor(kBlack);
               h_data_mb->GetXaxis()->SetTitle("m_{B} [GeV]");
               // Make slices for QCD
               TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
               TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
               TH1D *h_nob = (TH1D *) h3D_nob->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
               // TH1D *h_nob = (TH1D*) h_data_mb->Clone("h_nob_zeroinputfortest"); h_nob->Reset(); // zero test (disabled)
               
               //  Make slices for Bjet 
               TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
               TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
               
               
               
               // -- Deattach hists from root file
               h_data_mb->SetDirectory(0);
               h_b->SetDirectory(0);
               h_bb->SetDirectory(0);
               h_b_bjet->SetDirectory(0);
               h_bb_bjet->SetDirectory(0);
               // h_nob->SetDirectory(0);
               
               
               TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);
               TString srangeDeltaR = Form("#DeltaR [Full range]");
               // pT bin label for plot titles and TLatex annotations
               TString srangePt = Form("p_{T} [%.0f-%.0f] GeV", jtpt_binsVector[ibin_pt-1], jtpt_binsVector[ibin_pt]);
               
               // Calculate true fractions to be used as initial values for the fit (the true fractions are the QCD qcd ones)            
               // double int0 = h_bb->Integral(1, mb_bins, "width");// sig
               // double int1 = h_b->Integral(1, mb_bins, "width");// bkg1
               // double int3 = h_nob->Integral(1, mb_bins, "width"); // bkg 3
               
               // From bjet sample 
               double int0_bjet = h_bb_bjet->Integral(1, mb_bins, "width");// sig
               double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");// bkg1
               
               // For MC test case only
               // double integral_inputdata1 = h_data_mb->Integral(1, mb_bins,"width");
               // double integral_inputdata2 = h_data_bjet_mb->Integral(1, mb_bins,"width");
               // cout << "Before Samples additions:\n QCD data integral = "<< integral_inputdata1 << endl;
               // cout << "Bjet data integral = "<< integral_inputdata2 << endl;
               
               // Add QCD and bjet sample for the pseudodata 
               // h_data_mb->Add(h_data_bjet_mb);
               
               // For data: After the sum of QCD + bjet 
               double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
               
               cout << "Total input data integral = "<< integral_inputdata << endl;
               // std::cout << "int0 h_bb=" << int0 << std::endl;
               // std::cout << "int1 h_b =" << int1 << std::endl;
               // cout << "2B/1B in QCD sample = " << int0/int1 << endl;
               // cout << "QCD: 2B/(1B + 2B) = " << int0/(int0 + int1) << endl;
               
               std::cout << "int0 h_bb_bjet=" << int0_bjet << std::endl;
               std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
               cout <<  "2B/1B in  Bjet sample = " << int0_bjet/int1_bjet << endl;
               // cout << "2B/(2B + 1B) in qcd sample = " << int0/(int0 + int1) << endl;
               
               // std::cout << "int3 h_nob =" << int3 << std::endl;
               
               
               // --- For now: the TRUE signal (2B) fraction is used from the QCD sample.
               // the fraction of 2B/1B in Bjet and QCD sample is exactly the same
               // double sig_fraction_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int2 + int3 )); 
               // double bkg_fraction_b_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int2 + int3 ));
               
               
               // No contiubtion from hmore2B 
               /*
               double sig_fraction_true = (int0  + int1 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int3 )); 
               double bkg_fraction_b_true = (int0  + int1 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int3 ));
               std::cout << "QCD : 2B sig_fraction_true =" << sig_fraction_true << std::endl;
               std::cout << "QCD: 1B bkg_fraction_b_true =" << bkg_fraction_b_true << std::endl;
               // -- Compute the errors of true fractions for later fill
               double True_bkg_b_err = 0.0;
               double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");
               double True_sig_err = 0.0;
               double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");
               */
              
              // -- test bjet only as true signal 
              double sig_fraction_true = (int0_bjet  + int1_bjet ) == 0 ? 0 : (int0_bjet / (int0_bjet + int1_bjet)); 
              double bkg_fraction_b_true = (int0_bjet  + int1_bjet ) == 0 ? 0 : (int1_bjet / (int0_bjet + int1_bjet));
              std::cout << "bjet : 2B sig_fraction_true =" << sig_fraction_true << std::endl;
              std::cout << "bjet: 1B bkg_fraction_b_true =" << bkg_fraction_b_true << std::endl;  
              double True_bkg_b_err = 0.0;
              double True_bkg_b_integral = h_b_bjet->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");
              double True_sig_err = 0.0;
              double True_sig_integral = h_bb_bjet->IntegralAndError(1, mb_bins, True_sig_err, "width");
              
              // -- bjet-only templates (disabled in favour of QCD+bjet combined below)
              // TH1D* h_sig = (TH1D*) h_bb_bjet->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt));
              // TH1D* h_bkg = (TH1D*) h_b_bjet->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt));
              //     Int_t h_sig_bins = h_sig->GetNbinsX();
              //     Int_t h_bkg_bins = h_bkg->GetNbinsX();
              
              // -- Get Normalized templates as input PDF to the fit (shape matters only)
              // signal: 2B from QCD + Bjet sample combined
              TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
              h_sumsig->Add(h_bb_bjet);
              Int_t h_sig_bins = h_sumsig->GetNbinsX();
              TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized
              // background: 1B from QCD + Bjet sample combined
              TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
              h_sumbkg->Add(h_b_bjet);
              Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
              TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized
              // Total sum bkg = 0B + 1B
              TH1D* h_sumbkg_0b_1b = (TH1D*)  h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
              h_sumbkg_0b_1b->Add(h_nob);
              // write combined signal/bkg hists
              fout->cd();
              h_sumsig->Write();
              h_sumbkg->Write();
              h_sumbkg_0b_1b->Write();
              h_data_mb->Write();
              
              // set titles for fit components (used in legends below)
              h_sumsig->SetTitle("2B: QCD + bjet");
              h_sumbkg->SetTitle("1B: QCD + bjet");
              h_nob->SetTitle("0B");
              h_nob->SetFillColor(kRed);
              h_sumsig->SetFillColor(kBlue-10);
              h_sumbkg->SetFillColor(kRed-10);
              h_sumbkg->SetFillStyle(3244);
              
              // ----- To avoid empty bins if exist!, set them to eps value
              const double eps = 1e-6; 
              for (int i = 1; i <= h_sig_bins; i++){
                  if (h_sig->GetBinContent(i) <= 0){
                      h_sig->SetBinContent(i, eps);
                      cout << 
                      Form("INFO: (ptbin %d, deltaRbin %d) has empty  signal Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i) 
                      << endl; 
                    }
                }
                
                for (int i = 1; i <= h_bkg_bins; i++){
                    if (h_bkg->GetBinContent(i) <= 0){
                        h_bkg->SetBinContent(i, eps);
                        cout <<
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty bkg  Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i)  
                        << endl; 
                    }
                }
                /// Normalize safely for the Roofit
                if(h_sig->Integral(1, h_sig_bins, "width") == 0.0 || h_bkg->Integral(1, h_bkg_bins, "width") == 0.0)
                {
                    cerr<< "Signal or bkg templates has zero integral"<< endl;
                    continue;
                }
                
                h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));
                h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));
                cout << "After Normalization \n True 2B integral = " << h_sig->Integral(1, h_sig_bins, "width") << endl;
                cout << "True 1B integral with width option =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;
                cout << "True 1B integral without width option =  "<<  h_bkg->Integral(1, h_bkg_bins) << endl;
                
                /// Hay que hacer esto ???
                // ----------- Franc. added the hrest to the hsig template ----------- 
                //Add signal and rest with the correct scaling (from Lida)
                // add together light+c to sig with nominal MC ratio
                // a sig + b bb + c light = n, a+b+c=1, a+c=1-b
                // c'=c/(a+c)=(1-a-b)/(1-b)
                // a'=a/(a+c)=a/(1-b)=1-c'
                // -- Afnan: what are these param?
                // -- Not used in my case: 
                // double cl_frac = (1 - sig_fraction_true - bkg_fraction_b_true) / (1 - bkg_fraction_b_true); // nominal
                // cl_frac = 2*cl_frac; // systematic
                ////// h_sig->Add(h_sig, h_rest, sig_frac, 1-sig_frac);
                
                
                /*
                // If no other than 2B and 1B  are used: 
                /////// For now: since i have pure signal, use charm-light contribution = 0
                double cl_frac = 0.; // systematic
                double sig_frac = 1 - cl_frac;
                */
               
               // -- Effective bkg PDF 
               // Assume : a + b + c = 1; where a, b, c, are fractions in in nominal MC like 
               // a: 2B fraction, b: 1B fraction, c: 0B fraction, in 
               // and For Effective bkg: b`+ c` should = 1, so relate b`to b .. and simialr for c`
               // This result in: b`= b/(b+c), and c`= c/(b+c); one can use b+c = 1-a;  
               // I will write it in terms of the 2B sgnal, 1B bkg true fractions, before use effecive ones 
               double eff_bkg0B = (1- sig_fraction_true - bkg_fraction_b_true)/(1- sig_fraction_true);// c`
               double eff_bkg1B = 1. - eff_bkg0B;// b`
               
               // -- test bjet, no 0B //Zoe
               eff_bkg0B = 0;
               eff_bkg1B = 1; 
               
               // Build effective signal hist: with new relaitve normalization, the integral should = 1
               h_bkg->Add(h_bkg, h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist
               
               double err_int;
               double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins, err_int ,"width");
               cout << "Effective Bkg normalized hist integral +/- uncertainity = "<< int_val << "+/-" << err_int << endl;
               cout << "effective 1B fraction and 0B = "<< eff_bkg1B << ", " << eff_bkg0B << endl;
               // -- Note: the effective bkg integral is Effective Bkg normalized hist integral = 0.997777
               // if (h_bkg->Integral(1,  h_bkg_bins, "width") != 1.) {cout << "Effective BKG PDF is not normalized!"<< endl; return;}
               
               // =========================================================
               // PLOT 3 (pre-fit, non-normalized, per pT bin / deltaR slice)
               // Stacked 1B (green) + 2B (blue) templates from the combined
               // QCD+bjet MC, scaled so that (1B+2B) total matches the data
               // integral in this bin. Data shown as markers.
               // Purpose: check template shapes and scale in each kinematic bin
               // before the fit is run.
               // =========================================================
               ///// Pre-fit stacked canvas: templates scaled to data (before fit)
               {
                   TH1D* h_sig_prefit = (TH1D*) h_sig->Clone(Form("h_sig_prefit_%d_%d", ibin_dr, ibin_pt));
                   h_sig_prefit->Scale(sig_fraction_true * integral_inputdata);
                   h_sig_prefit->SetTitle("2B (bjet)");
                   h_sig_prefit->SetFillColor(kRed-9); h_sig_prefit->SetLineColor(kRed-3);  // 2B
                   TH1D* h_bkg_prefit = (TH1D*) h_bkg->Clone(Form("h_bkg_prefit_%d_%d", ibin_dr, ibin_pt));
                   h_bkg_prefit->Scale(bkg_fraction_b_true * integral_inputdata);
                   h_bkg_prefit->SetTitle("1B (bjet)");
                   h_bkg_prefit->SetFillColor(kTeal-8); h_bkg_prefit->SetLineColor(kTeal-6); // 1B
                   // data: open circle markers, black
                   h_data_mb->SetTitle(dataLabel);
                   h_data_mb->SetMarkerStyle(24); // open circle
                   h_data_mb->SetMarkerColor(kBlack); h_data_mb->SetLineColor(kBlack); h_data_mb->SetLineWidth(2);
                   
                   THStack hstack_prefit(Form("hstack_prefit_%d_%d", ibin_dr, ibin_pt), "");
                   hstack_prefit.SetTitle(";m_{2B} [GeV];Counts");
                   hstack_prefit.Add(h_sig_prefit);  // 2B first
                   hstack_prefit.Add(h_bkg_prefit);  // 1B second
                   
                   hstack_prefit.SetMaximum(1.2 * h_data_mb->GetMaximum());
                   
                   auto canva_prefit = new TCanvas(Form("prefit_stacked_%d_%d", ibin_dr, ibin_pt),
                   Form("Pre-fit | %s | %s", srangePt.Data(), srangeDeltaR.Data()), 800, 800);
                   canva_prefit->cd();
                   hstack_prefit.Draw("hist E");
                   h_data_mb->Draw("hist E same");
                   // pT bin label drawn on the plot (top-left, outside legend area)
                   TLatex latex_prefit;
                   latex_prefit.SetNDC(); latex_prefit.SetTextSize(0.04); latex_prefit.SetTextFont(62);
                   latex_prefit.DrawLatex(0.15, 0.88, srangePt);
                   gPad->Modified(); gPad->Update();
                   canva_prefit->Modified(); canva_prefit->Update();
                   fout->cd();
                   canva_prefit->BuildLegend(0.6, 0.65, 0.9, 0.9, srangeDeltaR);
                   canva_prefit->Write();
                   canva_prefit->Print(pdfName);
                }
                
                ///// Fitting
                // Create the observable
                Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
                Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
                RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
                mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
                
                // Inputs 
                // Create the RooDataHist object for the observed data + templates
                RooDataHist *dh_data_mb = new RooDataHist(Form("dh_data_mb_%d_%d", ibin_dr, ibin_pt), "dh_data_mb", mb, RooFit::Import(*h_data_mb));
                RooDataHist *dh_sig_mb = new RooDataHist(Form("dh_sig_mb_%d_%d", ibin_dr, ibin_pt), "dh_sig_mb", mb, RooFit::Import(*h_sig));
                RooDataHist *dh_bkg_mb = new RooDataHist(Form("h_bkg_mb_%d_%d", ibin_dr, ibin_pt), "dh_bkg_mb", mb, RooFit::Import(*h_bkg));
                
                // Create the RooHistPdf objects for the template PDFs
                RooHistPdf sig_template(Form("sig_template_%d_%d", ibin_dr, ibin_pt), "sig_template", mb, *dh_sig_mb);
                RooHistPdf bkg_template(Form("bkg_template_%d_%d", ibin_dr, ibin_pt), "bkg_template", mb, *dh_bkg_mb);
                
                // Create list of templates
                RooArgList template_list(sig_template, bkg_template, "template_list");
                
                // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
                // debug 
                // sig_fraction_true = 0.5;
                RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true , 0., 1);// , // 1-bkg_fraction_b_true
                
                /*
                // test fix fit par:
                sig_fraction_val.setVal(sig_fraction_true);     // set it to zero
                sig_fraction_val.setConstant(true); // fix it so it won't float in the fit
                */
               
               // Create the composite PDF using a linear combination of the template PDFs
               // RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, true); 
               //fit the scaling when adding the single b with the more b distribution
               
               RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, false); // fasle: dont use extended fit: assume fraction not absolute yields  
               RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name            
               // -- test fit in a smaller  ragne? then fit 
               /*
               mb.setRange("fitRange", 4,  9); // limit the fit range 
               RooFitResult* result = model0.fitTo(*dh_data_mb,
               RooFit::Range("fitRange"),
               RooFit::NormRange("fitRange"),
               RooFit::SumW2Error(true),
               RooFit::Save(),
               RooFit::CloneData(true),
               RooFit::PrintLevel(2),
               RooFit::Strategy(1),
               RooFit::Minos(false));
               
               */
              
              //instead of sign histogram here we would put the data histogram
              Int_t status = result->status();
              result->Print();
              
              std::cout << "covariance matrix:" << std::endl;
              (result->covarianceMatrix().Print());
              
              //Check if it converged for a dr and jtpt bin
              if (status != 0) {
                  std::cout << "\n\n\n\n!!!Fitting for ipt = " << ibin_pt 
                  << ", ix = " << ibin_dr 
                  << " did not converge\n\n\n\n" << std::endl;
                  non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
                  continue;
                }
                
                // Get the fitted parameter values
                double a = sig_fraction_val.getValV();
                double da = sig_fraction_val.getError();
                
                //Print some check
                std::cout << "RooFit result: \n Signal PDF fraction: a = " << a << " \n its uncertainity: da = " << da << std::endl;
                /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
                
                Double_t p0, p1, p2, errP0, errP1, errP2;
                p0 = a; // signal 2B
                // Total bkg: (1-a)
                p1 = (1-a)*eff_bkg1B; // 1B bkg 
                p2 = (1-a)*eff_bkg0B; // 0B bkg 
                errP0 = da;
                errP1 = da *eff_bkg1B;
                errP2 = da *eff_bkg0B;
                
                std::cout << "a: 2B =" << sig_fraction_true << ", a': after fit=" << p0 << std::endl;            
                std::cout << "a'/a  for 2B = " << p0/sig_fraction_true << std::endl;
                
                
                std::cout << "b: 1B =" << bkg_fraction_b_true << ", b': 1B =" << p1 << std::endl;
                std::cout << "b'/b for 1B = " << p1/bkg_fraction_b_true << std::endl;
                
                
                std::cout << "c: 0B =" << (1 - sig_fraction_true - bkg_fraction_b_true) << ", c': 0B after fit =" << p2 << std::endl;
                std::cout << "c'/c = " << p2/(1-sig_fraction_true - bkg_fraction_b_true) << std::endl;
                
                
                
                
                // save the fit
                // rescaling for the signal only (not rest background)
                h_sig_fraction->SetBinContent(ibin_dr, ibin_pt, p0);
                h_sig_fraction->SetBinError(ibin_dr, ibin_pt, errP0);
                h_sig_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP0);
                
                h_bkg_fraction->SetBinContent(ibin_dr, ibin_pt, p1);
                h_bkg_fraction->SetBinError(ibin_dr, ibin_pt, errP1);
                h_bkg_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP1);
                
                //save the true fraction
                h_sig_frac_true->SetBinContent(ibin_dr, ibin_pt, sig_fraction_true);
                h_bkg_frac_true->SetBinContent(ibin_dr, ibin_pt, bkg_fraction_b_true);
                
                // -- Compute siganl and bkg fractions uncertainity 
                // -- sigma True for S or B  = (B x sigmaS + SxsigmaB)/(S+B)²
                // double err_true_frac = (int1 * True_sig_err + int0 * True_bkg_b_err)/(int0+ int1)/(int0+ int1);// int0 for signal, int1 for bkg 
                // when bjet oly 
                double err_true_frac = (int1_bjet * True_sig_err + int0_bjet * True_bkg_b_err)/(int0_bjet+ int1_bjet)/(int0_bjet+ int1_bjet);
                
                
                // -- Save true fraction uncertainity
                h_sig_frac_true_error->SetBinContent(ibin_dr, ibin_pt, err_true_frac);
                h_bkg_frac_true_error->SetBinContent(ibin_dr, ibin_pt, err_true_frac);
                
                // -- And for the integrated bin 
                if(!ibin_dr){
                    h_sig_fraction_DeltaRIntBin->SetBinContent(ibin_pt, p0);
                    h_sig_fraction_DeltaRIntBin->SetBinError(ibin_pt, errP0);
                    h_sig_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, errP0);
                    
                    h_bkg_fraction_DeltaRIntBin->SetBinContent(ibin_pt, p1);
                    h_bkg_fraction_DeltaRIntBin->SetBinError(ibin_pt, errP1);
                    h_bkg_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, errP1);
                    
                    //save the true fraction
                    h_sig_true_fraction_DeltaRIntBin->SetBinContent(ibin_pt, sig_fraction_true);
                    h_bkg_true_fraction_DeltaRIntBin->SetBinContent(ibin_pt, bkg_fraction_b_true);
                    
                    // -- Save true fraction uncertainity
                    h_sig_true_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, err_true_frac);
                    h_bkg_true_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, err_true_frac);
                }
                
                
                
                // -- After fits: save mass distribution re and post-fit 
                TH1D *h_sig_fit = (TH1D*) h_sig->Clone(Form("h_sig_fit_%d_%d", ibin_dr, ibin_pt));
                h_sig_fit->Scale(p0 * integral_inputdata);
                h_sig_fit->SetTitle("After fit: Sig. (2B)");
                h_sig_fit->SetMarkerColor(kRed-3);
                h_sig_fit->SetLineColor(kRed-3);
                h_sig_fit->SetFillColor(kRed-9);   // 2B: kRed-9 fill, kRed-3 border
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

                TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("After fit: bkg (1B+ 0B)");
                h_bkg_fit->SetMarkerColor(kBlue-7);
                h_bkg_fit->SetLineColor(kBlue-7);

                TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("After fit: 1B");
                h_bkg_fit_1b->SetMarkerColor(kTeal-6);
                h_bkg_fit_1b->SetLineColor(kTeal-6);
                h_bkg_fit_1b->SetLineWidth(2);
                h_bkg_fit_1b->SetFillColor(kTeal-8);  // 1B: kTeal-8 fill, kTeal-6 border
                
                TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("After fit: 0B");
                h_bkg_fit_nob->SetFillColor(kBlue-5);
                h_bkg_fit_nob->SetLineColor(kBlue-5);  // 0B: kBlue-5
                h_bkg_fit_nob->SetLineWidth(2);

                // And the post-fit template (fitted sig + bkg) — used as ratio denominator
                TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                h_total_fit->Add(h_bkg_fit);
                h_total_fit->SetTitle("Total fit");
                h_total_fit->SetLineColor(kBlue-7);
                h_total_fit->SetLineWidth(2);
                h_total_fit->SetFillStyle(0);
                
                
                // Save Signal and bkg mass distributions after the fit 
                fout->cd();
                h_sig_fit->Write();
                h_bkg_fit->Write();
                h_total_fit->Write();
                
                // =========================================================
                // PLOT 4 (after-fit, non-normalized, per pT bin / deltaR slice)
                // Stacked 1B (orange) + 2B (red) templates rescaled by the fitted
                // fractions (p1 and p0) times the total data integral.
                // Data shown as markers.
                // Purpose: show how well the fit decomposition reproduces the
                // observed m_2B distribution in absolute counts.
                // =========================================================
                ///// After-fit stacked canvas (non-normalized): 1B + 2B + data + ratio pad
                {
                    h_data_mb->SetTitle(dataLabel);
                    h_data_mb->SetMarkerStyle(24); // open circle
                    h_data_mb->SetMarkerColor(kBlack); h_data_mb->SetLineColor(kBlack); h_data_mb->SetLineWidth(2);

                    // total fit = 1B + 2B (for ratio denominator and line on plot)
                    TH1D* h_total_fit_af = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_af_%d_%d", ibin_dr, ibin_pt));
                    h_total_fit_af->Add(h_bkg_fit_1b);
                    h_total_fit_af->SetTitle("Total fit");
                    h_total_fit_af->SetLineColor(kBlue-7); h_total_fit_af->SetLineWidth(2); h_total_fit_af->SetFillStyle(0);

                    THStack* hstack_afterfit = new THStack(Form("hstack_afterfit_%d_%d", ibin_dr, ibin_pt), "");
                    hstack_afterfit->SetTitle("After fit;m_{2B} [GeV];Counts");
                    hstack_afterfit->Add(h_sig_fit);    // 2B on top
                    hstack_afterfit->Add(h_bkg_fit_1b); // 1B at bottom
                    
                    double ymax_af = TMath::Max(h_data_mb->GetMaximum(), h_sig_fit->GetMaximum() + h_bkg_fit_1b->GetMaximum());
                    hstack_afterfit->SetMaximum(1.2 * ymax_af);

                    // ratio: data / total fit
                    TH1D* h_ratio_af = (TH1D*) h_data_mb->Clone(Form("h_ratio_af_%d_%d", ibin_dr, ibin_pt));
                    h_ratio_af->Divide(h_total_fit_af);
                    h_ratio_af->SetTitle(";m_{2B} [GeV];Data / Fit");
                    h_ratio_af->SetMarkerStyle(24); h_ratio_af->SetMarkerColor(kBlue-7);
                    h_ratio_af->SetLineColor(kBlue-7); h_ratio_af->SetLineWidth(2);
                    h_ratio_af->GetYaxis()->SetRangeUser(0.5, 1.5);
                    h_ratio_af->GetYaxis()->SetNdivisions(504);
                    h_ratio_af->GetYaxis()->SetTitleSize(0.10); h_ratio_af->GetYaxis()->SetTitleOffset(0.45);
                    h_ratio_af->GetXaxis()->SetTitleSize(0.10); h_ratio_af->GetXaxis()->SetLabelSize(0.09);
                    h_ratio_af->GetYaxis()->SetLabelSize(0.09);

                    TString sname_canvas_afterfit = sname_canvas + "_afterfit";
                    auto canva_afterfit = new TCanvas(sname_canvas_afterfit,
                        Form("After fit (non-normalized) | %s | %s", srangePt.Data(), srangeDeltaR.Data()), 800, 900);
                    // upper pad: main plot (~70%)
                    TPad* pad1_af = new TPad("pad1_af", "", 0, 0.28, 1, 1);
                    pad1_af->SetBottomMargin(0.03); pad1_af->SetTopMargin(0.08);
                    pad1_af->Draw(); pad1_af->cd();
                    hstack_afterfit->Draw("hist E");
                    // remove x-axis label from upper pad
                    hstack_afterfit->GetXaxis()->SetLabelSize(0);
                    hstack_afterfit->GetXaxis()->SetTitleSize(0);
                    h_data_mb->Draw("hist E same");
                    h_total_fit_af->Draw("HIST E same"); // total fit line on top
                    TLatex latex_af;
                    latex_af.SetNDC(); latex_af.SetTextSize(0.05); latex_af.SetTextFont(62);
                    latex_af.DrawLatex(0.15, 0.88, srangePt);
                    gPad->Modified(); gPad->Update();
                    pad1_af->BuildLegend(0.6, 0.60, 0.92, 0.88, srangeDeltaR);
                    // lower pad: ratio (~28%)
                    canva_afterfit->cd();
                    TPad* pad2_af = new TPad("pad2_af", "", 0, 0, 1, 0.28);
                    pad2_af->SetTopMargin(0.03); pad2_af->SetBottomMargin(0.38);
                    pad2_af->Draw(); pad2_af->cd();
                    h_ratio_af->Draw("PE");
                    // reference line at 1
                    TLine* line_af = new TLine(h_ratio_af->GetXaxis()->GetXmin(), 1,
                                               h_ratio_af->GetXaxis()->GetXmax(), 1);
                    line_af->SetLineColor(kBlack); line_af->SetLineStyle(2); line_af->Draw();
                    gPad->Modified(); gPad->Update();
                    canva_afterfit->Modified(); canva_afterfit->Update();
                    fout->cd();
                    canva_afterfit->Write();
                    canva_afterfit->Print(pdfName);
                }

                // =========================================================
                // PLOT 5 (after-fit, normalized, per pT bin / deltaR slice)
                // Same post-fit components as plot 4, but each divided by the
                // total data integral so all distributions have the same area.
                // Data also self-normalized. Ratio pad: data_norm / fit_norm.
                // =========================================================
                ///// After-fit stacked canvas (normalized): 1B + 2B + data + ratio pad
                {
                    TH1D* h_sig_fit_normstack    = (TH1D*) h_sig_fit->Clone(Form("h_sig_fit_normstack_%d_%d", ibin_dr, ibin_pt));
                    h_sig_fit_normstack->Scale(1./integral_inputdata);
                    h_sig_fit_normstack->SetTitle("After fit: Sig. (2B)");
                    h_sig_fit_normstack->SetFillColor(kRed-9); h_sig_fit_normstack->SetLineColor(kRed-3);
                    TH1D* h_bkg_fit_1b_normstack = (TH1D*) h_bkg_fit_1b->Clone(Form("h_bkg_fit_1b_normstack_%d_%d", ibin_dr, ibin_pt));
                    h_bkg_fit_1b_normstack->Scale(1./integral_inputdata);
                    h_bkg_fit_1b_normstack->SetTitle("After fit: 1B");
                    h_bkg_fit_1b_normstack->SetFillColor(kTeal-8); h_bkg_fit_1b_normstack->SetLineColor(kTeal-6);
                    TH1D* hdata_norm_af = (TH1D*) h_data_mb->Clone(Form("hdata_norm_af_%d_%d", ibin_dr, ibin_pt));
                    hdata_norm_af->Scale(1./integral_inputdata);
                    hdata_norm_af->SetTitle(dataLabel);
                    hdata_norm_af->SetMarkerStyle(24);
                    hdata_norm_af->SetLineColor(kBlack); hdata_norm_af->SetMarkerColor(kBlack); hdata_norm_af->SetLineWidth(2);

                    // total fit normalized
                    TH1D* h_total_fit_norm = (TH1D*) h_sig_fit_normstack->Clone(Form("h_total_fit_norm_%d_%d", ibin_dr, ibin_pt));
                    h_total_fit_norm->Add(h_bkg_fit_1b_normstack);
                    h_total_fit_norm->SetTitle("Total fit");
                    h_total_fit_norm->SetLineColor(kBlue-7); h_total_fit_norm->SetLineWidth(2); h_total_fit_norm->SetFillStyle(0);

                    THStack* hstack_norm_afterfit = new THStack(Form("hstack_norm_afterfit_%d_%d", ibin_dr, ibin_pt), "");
                    hstack_norm_afterfit->SetTitle(";m_{2B} [GeV];");
                    hstack_norm_afterfit->Add(h_sig_fit_normstack);    // 2B on top
                    hstack_norm_afterfit->Add(h_bkg_fit_1b_normstack); // 1B at bottom
                    double ymax_afn = TMath::Max(hdata_norm_af->GetMaximum(),
                        h_sig_fit_normstack->GetMaximum() + h_bkg_fit_1b_normstack->GetMaximum());
                    hstack_norm_afterfit->SetMaximum(1.2 * ymax_afn);

                    // ratio: data_norm / total_fit_norm
                    TH1D* h_ratio_afn = (TH1D*) hdata_norm_af->Clone(Form("h_ratio_afn_%d_%d", ibin_dr, ibin_pt));
                    h_ratio_afn->Divide(h_total_fit_norm);
                    h_ratio_afn->SetTitle(";m_{2B} [GeV];Data / Fit");
                    h_ratio_afn->SetMarkerStyle(24); h_ratio_afn->SetMarkerColor(kBlue-7);
                    h_ratio_afn->SetLineColor(kBlue-7); h_ratio_afn->SetLineWidth(2);
                    h_ratio_afn->GetYaxis()->SetRangeUser(0.5, 1.5);
                    h_ratio_afn->GetYaxis()->SetNdivisions(504);
                    h_ratio_afn->GetYaxis()->SetTitleSize(0.10); h_ratio_afn->GetYaxis()->SetTitleOffset(0.45);
                    h_ratio_afn->GetXaxis()->SetTitleSize(0.10); h_ratio_afn->GetXaxis()->SetLabelSize(0.09);
                    h_ratio_afn->GetYaxis()->SetLabelSize(0.09);

                    auto canva_stack_norm_afterfit = new TCanvas(Form("stacked_norm_afterfit_%d_%d", ibin_dr, ibin_pt),
                        Form("After fit (normalized) | %s | %s", srangePt.Data(), srangeDeltaR.Data()), 800, 900);
                    TPad* pad1_afn = new TPad("pad1_afn", "", 0, 0.28, 1, 1);
                    pad1_afn->SetBottomMargin(0.03); pad1_afn->SetTopMargin(0.08);
                    pad1_afn->Draw(); pad1_afn->cd();
                    hstack_norm_afterfit->Draw("hist E");
                    hstack_norm_afterfit->GetXaxis()->SetLabelSize(0);
                    hstack_norm_afterfit->GetXaxis()->SetTitleSize(0);
                    hdata_norm_af->Draw("PE same");
                    h_total_fit_norm->Draw("HIST same");
                    TLatex latex_afn;
                    latex_afn.SetNDC(); latex_afn.SetTextSize(0.05); latex_afn.SetTextFont(62);
                    latex_afn.DrawLatex(0.15, 0.88, srangePt);
                    gPad->Modified(); gPad->Update();
                    pad1_afn->BuildLegend(0.6, 0.60, 0.92, 0.88, srangeDeltaR);
                    canva_stack_norm_afterfit->cd();
                    TPad* pad2_afn = new TPad("pad2_afn", "", 0, 0, 1, 0.28);
                    pad2_afn->SetTopMargin(0.03); pad2_afn->SetBottomMargin(0.38);
                    pad2_afn->Draw(); pad2_afn->cd();
                    h_ratio_afn->Draw("PE");
                    TLine* line_afn = new TLine(h_ratio_afn->GetXaxis()->GetXmin(), 1,
                                                h_ratio_afn->GetXaxis()->GetXmax(), 1);
                    line_afn->SetLineColor(kBlack); line_afn->SetLineStyle(2); line_afn->Draw();
                    gPad->Modified(); gPad->Update();
                    canva_stack_norm_afterfit->Modified(); canva_stack_norm_afterfit->Update();
                    fout->cd();
                    canva_stack_norm_afterfit->Write();
                    canva_stack_norm_afterfit->Print(pdfName);
                }
                    
                    cout << "---------------------\n\n\n" << endl; 
                } // loop over deltaR bins 
            }
            
            // Save histograms
            // TH3D 
            // for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
            for (auto h : {h3D_data, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
            
            
            
            // TH2D     
            for (auto h : {h_sig_fraction, h_sig_fraction_error,
                h_bkg_fraction, h_bkg_fraction_error,
                h_sig_frac_true, h_sig_frac_true_error,
                h_bkg_frac_true, h_bkg_frac_true_error 
            }) {
                h->Write();
            }
            
            // TH1D              
            for (auto h :{ 
                h_sig_fraction_DeltaRIntBin,
                h_sig_fraction_DeltaRIntBin_error,
                h_sig_true_fraction_DeltaRIntBin,
                h_sig_true_fraction_DeltaRIntBin_error,
                h_bkg_fraction_DeltaRIntBin,
                h_bkg_fraction_DeltaRIntBin_error,
                h_bkg_true_fraction_DeltaRIntBin,
                h_bkg_true_fraction_DeltaRIntBin_error
            }) {
                h->Write();
            }
            
            
            // //See if some bins did not converge
            for (auto p : non_converge_bins) {
                std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
            }
            
            // Close the PDF ("]" closes without adding a blank page)
            TCanvas("cclose","cclose").Print(pdfName + "]");
            
            
            fout->Close();
                            
            
        }
        //Draws the result of the template fit
        
        std::unique_ptr<TCanvas> draw_template_fit_result(
            TString fout_name,
            TFile* foutputPlots,
            TString &dataset,
            TString &folder,
            TString &pT_selection, Int_t pt_bin = 0, bool isIntegDeltaR = false){
                
                // pt bin is only needed for not integarted deltaR case 
                
                // -- the path to the outpt file 
                
                // For trvial tests 
                TString trivialMC_label = fout_name.Contains("trivialMC") ? "trivialMC": "";
                cout << "For drawing tempaltes, the MC label is " << trivialMC_label << endl;
                
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
                
                if(!isIntegDeltaR){
                    //signal fraction
                    TH2D *h_2D = (TH2D*)file->Get("h_sig_fraction");
                    h_2D->SetDirectory(nullptr);
                    h = (TH1D*)h_2D->ProjectionX("h", pt_bin,pt_bin);
                    h->GetXaxis()->SetTitle("\\#Delta\\ r");
                    TH2D *htrue_2D = (TH2D*)file->Get("h_sig_frac_true");
                    htrue_2D->SetDirectory(nullptr);
                    htrue = (TH1D*)htrue_2D->ProjectionX("htrue", pt_bin,pt_bin);
                    // and its uncertainity 
                    TH2D *htrue_err_2D = (TH2D*)file->Get("h_sig_frac_true_error");
                    TH1D* htrue_err = (TH1D*) htrue_err_2D->ProjectionX("htrue_err", pt_bin,pt_bin);
                    for (int ibin = 1; ibin <= htrue_err->GetNbinsX(); ibin++)
                    {
                        /// Assign the uncertainity to the iD hist of true signal Vs. DeltaR 
                        htrue->SetBinError(ibin, htrue_err->GetBinContent(ibin));
                    }
                    
                    //background fraction
                    TH2D *hbkg_2D = (TH2D*)file->Get("h_bkg_fraction");
                    hbkg_2D->SetDirectory(nullptr);
                    hbkg = (TH1D*)hbkg_2D->ProjectionX("hbkg", pt_bin,pt_bin);
                    
                    TH2D *hbkg_true_2D = (TH2D*)file->Get("h_bkg_frac_true");
                    hbkg_true_2D->SetDirectory(nullptr);
                    hbkg_true = (TH1D*)hbkg_true_2D->ProjectionX("hbkg_true",pt_bin,pt_bin);
                    // and its uncertainity 
                    TH2D *hbkg_true_err_2D = (TH2D*)file->Get("h_bkg_frac_true_error");
                    TH1D* hbkg_true_err = (TH1D*) hbkg_true_err_2D->ProjectionX("hbkg_true_err", pt_bin,pt_bin);
                    for (int ibin = 1; ibin <= hbkg_true_err->GetNbinsX(); ibin++)
                    {
                        /// Assign the uncertainity to the iD hist of true signal Vs. DeltaR 
                        hbkg_true->SetBinError(ibin, hbkg_true_err->GetBinContent(ibin));
                    }
                    
                    // For rebinnng the deltaR nto 2
                    htrue->Rebin(2);
                    hbkg_true->Rebin(2);
                    h->Rebin(2);
                    hbkg->Rebin();
                    
                    sname_canvas = Form("%s_c_%s_ptbin_%d", trivialMC_label.Data() ,dataset.Data(), pt_bin);
                    
                }
                else if (isIntegDeltaR)
                {
                    cout << "Hello: this is integarted deltaR bin" << endl;
                    // -- DeltaR bin inetgarted (bin name in histograms = 0), plot Sig and Bkg fractions Vs. ptbins axis once 
                    h = (TH1D*) file->Get("h_sig_fraction_DeltaRIntBin"); // Signal 
                    h->GetXaxis()->SetTitle("\\#it{p}_{T}\\ [GeV]");
                    htrue = (TH1D*) file->Get("h_sig_true_fraction_DeltaRIntBin");
                    
                    TH1D* htrue_err = (TH1D*) file->Get("h_sig_true_fraction_DeltaRIntBin_error");
                    for (int ibin = 1; ibin <= htrue_err->GetNbinsX(); ibin++)
                    {
                        /// Assign the uncertainity to the iD hist of true signal Vs. DeltaR 
                        htrue->SetBinError(ibin, htrue_err->GetBinContent(ibin));
                    }
                    
                    hbkg = (TH1D*) file->Get("h_bkg_fraction_DeltaRIntBin");   
                    hbkg_true =  (TH1D*) file->Get("h_bkg_true_fraction_DeltaRIntBin");
                    
                    TH1D* hbkg_true_err = (TH1D*) file->Get("h_bkg_true_fraction_DeltaRIntBin_error");
                    for (int ibin = 1; ibin <= hbkg_true_err->GetNbinsX(); ibin++)
                    {
                        /// Assign the uncertainity to the iD hist of true signal Vs. DeltaR 
                        hbkg_true->SetBinError(ibin, hbkg_true_err->GetBinContent(ibin));
                    } 
                    
                    sname_canvas = Form("%s_c_%s_DeltaRIntBin_Allptbins", trivialMC_label.Data() ,dataset.Data());
                    
                }  
                // -- deattach hists from input root file: for drawing 
                h->SetDirectory(nullptr);
                htrue->SetDirectory(nullptr);
                hbkg->SetDirectory(nullptr);
                hbkg_true->SetDirectory(nullptr);
                
                auto c = std::make_unique<TCanvas>(sname_canvas,"Template Fit", 800, 800);
                gROOT->GetListOfCanvases()->Remove(c.get()); // to save it later 
                c->SetTitle(dataset);
                
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
                h->GetYaxis()->SetTitle("fitted parameter");
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
                if(!isIntegDeltaR) test_info_text->DrawLatex(0.15, 0.5, Form("%g < p_{T} < %g GeV", jtpt_binsVector[pt_bin-1], jtpt_binsVector[pt_bin]));
                else{ 
                    TString srangeDeltaR = Form("#DeltaR [Full range]");
                    test_info_text->DrawLatex(0.15, 0.5, srangeDeltaR);
                }
                // test_info_text->DrawLatex(0.15, 0.45, dataset);
                test_info_text->DrawLatex(0.15, 0.4, "Template from MC bjet");
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
                        
                        if(!isIntegDeltaR)  h_ratio->GetXaxis()->SetTitle("#DeltaR");
                        else { h_ratio->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");}
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
                        if(!isIntegDeltaR) {ptbin_name= Form("%g_%g", jtpt_binsVector[pt_bin-1], jtpt_binsVector[pt_bin]);}
                        else if (isIntegDeltaR) { cout << "HELLO   --- This is integarted bin "<< endl; ptbin_name= Form("IntegDeltaR_vspTbins");}
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
        
void template_fit_witheec(){
                        
    TString pT_selection = "80_140";// full range 
    TString folder = "";
    //Get data and mc labels 
    bool norm = true; //Normalise the eec plot
    bool all = true; //Plot also the mc eec
    bool eff_corr = false; //Correct for tagging efficiency
    bool alsoLowEG = false; // Add LowEG data 
    
    // -- Output folder to save the result of the tests 
    gSystem->mkdir(sDirname, kTRUE); // make directory if does not exist for the output files to be all in one place 
    
    TString sfoutputPlots_qcd = "Summary_histo_templatefit_test_14avril.root";
    TFile *foutputPlots_qcd = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_qcd.Data()), "RECREATE"); // For canvas output
    if (!foutputPlots_qcd || foutputPlots_qcd->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }
    
    // --  HG data and bjet with HLT80: using result of Zoe code and condor result 
    TString dataset_HighEG = "/data_CMS/cms/zaidan/analysis_lise/template_for_fit_histos_3D_HighEG_btag.root"; 
    TString dataset_LowEG = ""; 
    if (alsoLowEG) dataset_LowEG = "/data_CMS/cms/zaidan/analysis_lise/template_for_fit_histos_3D_LowEG_btag.root"; // with both LowEG

    TString dataname = "HighEG_allHLT";

    TString templates_qcd = "/data_CMS/cms/zaidan/analysis_lise/template_for_fit_histos_3D_qcd_btag.root"; // not used now 
    TString templates_bjet = "/data_CMS/cms/zaidan/analysis_lise/template_for_fit_histos_3D_bjet_test_btag.root";
    TString fout_name = "TemplateFits_histos_3d_allHLT_HighEGdata_bjet_" + pT_selection +  ".root";
    do_template_fit(dataset_HighEG, dataset_LowEG, alsoLowEG, templates_qcd, templates_bjet,  pT_selection, folder, fout_name,  "HighEG data");
    
    foutputPlots_qcd->Print();
    foutputPlots_qcd->Close();
    delete foutputPlots_qcd;
} 

