#include "tTree.h"
#include "TFile.h"
#include "binning_histos_all.h"


// For the output to be saved: 
TString sDirname = "TemplateFitOutput_HLT80_HGdata_bjet_noeec";


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
void do_template_fit_bjetonly(TString &dataset, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name, bool& alsoLowEG){
    /*
    This fit use sig and bkg templates from the input files.
    */
    // Dataset: for now HGtrigered data. 
    //-- Samples: Dijet + bjet 
    // -- Template fit pdfs : includes (0B, 1B) as bkg (added effectively), 2B as sig  

    bool isSampled = false; // for templates 
    
    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");

    // -- Using new naming convention
    TString namehData = "h_count_data";
    TString namehMC = "";

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h_count_bb";
    TString nameh1B = "h_count_b";

    TString nameh0B = "";
    // -- Additional histograms
    TString namehmore2B = "";

 
    //-- Dijet sample:
    TFile *file_dijet = new TFile(folder + templates, "read");
        if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());}
/*   
    TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b"); // 
    TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
        if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
        if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
    // -- 0B and more B contribution for drawing only
    TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone(nameh0B);
       if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}
*/

    //-- Bjet sample: Templates: Open file and get histograms from bjet: templates 
    TFile *file_bjet = new TFile(templates_bjet, "read");
         if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}
        TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
        TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
        if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
        if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}

/*
    // -- Additional contrinution to draw, not to fit 
        TH3D *h3D_more2b = (TH3D*)file_dijet->Get(namehmore2B)->Clone(namehmore2B);
            if(!h3D_more2b){Error("Get:", "histogram does not exist '%s' ",h3D_more2b->GetName());}
        TH3D *h3D_more2b_bjet = (TH3D*) file_bjet->Get(namehmore2B)->Clone(Form("%s_bjet", namehmore2B.Data() ) );
            if(!h3D_more2b_bjet){Error("Get:", "histogram does not exist '%s' ", h3D_more2b_bjet->GetName());}

*/
    //Open dataset:
    TFile *file_data = new TFile(dataset, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
        cout << "file data name " << file_data->GetName() << endl;
        TH3D *h3D_data;
        h3D_data = (TH3D*)file_data->Get(namehData)->Clone("h3D_data");
        if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}

        // -- Add also the lower triggerd data
        if(alsoLowEG){
            TString LowEGdata_name = "Prepare_Data/Data_Hist/withunder_overFlow_includedibins/hist_3d_4d_aggr_BDT_n1_data_LowEGJet_80_140.root";
            TFile *file_data_LowEG = new TFile(LowEGdata_name, "read");
            if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
            cout << "file data name " << file_data_LowEG->GetName() << endl;
            TH3D *h3D_dataLowEG;
            h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
            if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ",h3D_dataLowEG->GetName()); return;}
            h3D_data->Add(h3D_dataLowEG);
        }

        h3D_data->SetMarkerStyle(20);
        h3D_data->SetMarkerColor(0);
        h3D_data->SetLineColor(0);

        // -- For MC trivial test 
        TH3D *h3D_MC_pseudo;
        TH3D *h3D_MC_pseudo_bjet;
        /*
        if(isSampled){
            h3D_MC_pseudo = (TH3D*)file_dijet->Get(namehMC)->Clone(namehMC);
            h3D_MC_pseudo_bjet = (TH3D*)file_bjet->Get(namehMC)->Clone("h3D_data_bjet");
        }
        else {
            // Trivial closure test
            
            h3D_MC_pseudo = (TH3D*) h3D_b->Clone("h3D_pseudodata");
            h3D_MC_pseudo->Add(h3D_bb);
            h3D_MC_pseudo->Add(h3D_nob);
            
            // then add bjet sample also
            h3D_MC_pseudo_bjet = (TH3D*) h3D_b_bjet->Clone("h3D_data_bjet");
            h3D_MC_pseudo_bjet->Add(h3D_bb_bjet);
        }
        if(!h3D_MC_pseudo){Error("Get:", "histogram does not exist '%s' ",h3D_MC_pseudo->GetName());}
        if(!h3D_MC_pseudo_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_MC_pseudo_bjet->GetName());}
        */

    //Get number of bin entries
    Int_t bins_pt = h3D_data->GetNbinsZ();
    Int_t bins_dr = h3D_data->GetNbinsY(); 
    Int_t mb_bins = h3D_data->GetNbinsX();
    cout << "-- Data hist binning" << endl;
    cout << "pt bins = "<< bins_pt << endl;
    cout << "dr bins = "<< bins_dr << endl;
    cout << "mb bins = "<< mb_bins << endl;

    // -- Draw THstach for mb projection 
    // Get projections X 
    // hists names as for slices in bins (but outside {} so it is ok)
        TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX("h_data_mb", 1,bins_dr, 1, bins_pt ); // all bins 
            h_data_mb->SetTitle("Data");
            h_data_mb->SetLineColor(kBlack);
            h_data_mb->SetLineWidth(3);
            h_data_mb->SetMarkerColor(kBlack);
            h_data_mb->SetMarkerStyle(4);
        // TH1D *h_b = (TH1D *) h3D_b->ProjectionX("h_b");
        TH1D *h_b_bjet = (TH1D *) h3D_b_bjet->ProjectionX("h_b_bjet");
        // TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX("h_bb");
        TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet->ProjectionX("h_bb_bjet");
        // TH1D *h_nob = (TH1D *) h3D_nob->ProjectionX("h_nob");
        // TH1D *h_more2b = (TH1D *) h3D_more2b->ProjectionX("h_more2b");
        // TH1D *h_more2b_bjet = (TH1D *) h3D_more2b_bjet ->ProjectionX("h_more2b_bjet");
            /*h_more2b->SetDirectory(0);
            h_more2b_bjet->SetDirectory(0);
            h_more2b->SetFillColor(kOrange -2 ); 
            h_more2b->SetLineColor(kOrange-2);
            h_more2b_bjet->SetFillColor(kOrange -2 );
            h_more2b_bjet->SetLineColor(kOrange-2);
            h_more2b->SetTitle(">2B (Dijet)");
            h_more2b_bjet->SetTitle(">2B (bjet)");
            */

        // -- Set titles for the hisograms  to be used in the default legends 
        // h_b->SetTitle("1B (dijet)");
        h_b_bjet->SetTitle("1B (bjet)");
        // h_bb->SetTitle("2B (dijet)");
        h_bb_bjet->SetTitle("2B (bjet)");
        // h_nob->SetTitle("0B");

        /*
        // The trivial sum of sig+bkg t compare to full MC 
        TH1D *hsum = (TH1D*) h_b->Clone("hsum");
        hsum->SetTitle("BB signal + 1B bkg + 0B bkg");
        hsum->Add(h_b_bjet);
        hsum->Add(h_bb);
        hsum->Add(h_bb_bjet);
        hsum->Add(h_nob);
        hsum->SetLineWidth(3);
        hsum->SetLineColor(kMagenta+2);
        hsum->SetDirectory(0);
        fout->cd();
        hsum->Write();
        */
       

       // -- Deattach hists from root file 
        h_data_mb->SetDirectory(0);
        // h_b->SetDirectory(0);
        // h_bb->SetDirectory(0);
        h_b_bjet->SetDirectory(0);
        h_bb_bjet->SetDirectory(0);
        // h_more2b->SetDirectory(0);
        // h_nob->SetDirectory(0);
        // Set styles of hist to be drawn:
        // h_b->SetFillColor(kRed-10); h_b->SetLineColor(kRed-10);        
        // h_bb->SetFillColor(kBlue-10); h_bb->SetLineColor(kBlue-10);
        h_b_bjet->SetFillColor(kGreen-2); h_b_bjet->SetLineColor(kGreen-2);        
        h_bb_bjet->SetFillColor(kOrange-6); h_bb_bjet->SetLineColor(kOrange-6); // brown 
        // h_more2b->SetFillColor(kOrange); h_more2b->SetLineColor(kOrange);
        // h_nob->SetFillColor(kRed); h_nob->SetLineColor(kRed);
       

    THStack hstack_All("hstack_All","Mass stacked histogram");
        hstack_All.SetTitle(";m_{2B} [GeV];");
        // hstack_All.Add(h_more2b);
        // hstack_All.Add(h_more2b_bjet);
        // hstack_All.Add(h_nob);
        // hstack_All.Add(h_b);
        // hstack_All.Add(h_bb);
        hstack_All.Add(h_b_bjet);
        hstack_All.Add(h_bb_bjet);
        hstack_All.SetMaximum(20); /// for eec wights 


    auto canva = new TCanvas("All_templates_Data",Form("Templates(mb) and Full MC"), 800, 800 );
        canva->cd();
        // canva->SetLogy();
        hstack_All.Draw("hist E"); // templates stack
        // -- try draw scaled down data histogram?
        TH1D* hd_norm = (TH1D*) h_data_mb->Clone("hd_norm");
            hd_norm->Scale(1./10000);
            hd_norm->Draw("HIST E same"); 
            hd_norm->SetTitle("data/1e+04");
        // h_data_mb->Draw("HIST E same"); // data points 
        // hsum->Draw(" Hist same"); // templates sum          
        gPad->Modified();   
        gPad->Update();
        canva->Modified();
        canva->Update();
        fout->cd();
        canva->BuildLegend(0.8, 0.9, 0.8, 0.9);
        canva->Write();
        canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));

// return;
    // Another canvas for the normalized version stacked
        // stacked histograms normalization is to the total (not self normalized)
        THStack hstack_All_norm("hstack_All_norm","Mass stacked histogram");
            hstack_All_norm.SetTitle(";m_{2B} [GeV];");
            // normalize
            double total_int_1b_2b_bjet = h_b_bjet->Integral(1, mb_bins) +  h_bb_bjet->Integral(1, mb_bins); 
                TH1D* h_b_bjet_normstack = (TH1D*)  h_b_bjet->Clone("h_b_bjet_normstack");  h_b_bjet_normstack ->Scale(1./total_int_1b_2b_bjet);
                TH1D* h_bb_bjet_normstack = (TH1D*)  h_bb_bjet->Clone("h_bb_bjet_normstack");  h_bb_bjet_normstack ->Scale(1./total_int_1b_2b_bjet);
            hstack_All_norm.Add(h_b_bjet_normstack);
            hstack_All_norm.Add(h_bb_bjet_normstack);
            
            TH1D* hdata_selfnorm = (TH1D*) h_data_mb->Clone("hdata_selfnorm");
                hdata_selfnorm->Scale(1./h_data_mb->Integral(1, mb_bins));
            
            auto canva_stack_norm = new TCanvas("All_templates_Data_stacked_normalized",Form(""), 800, 800 );
                canva_stack_norm->cd();
                hstack_All_norm.Draw("hist E"); 
                hdata_selfnorm->Draw("HIST E same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_norm->Modified();
                canva_stack_norm->Update();
                fout->cd();
                canva_stack_norm->BuildLegend(0.8, 0.9, 0.8, 0.9);
                canva_stack_norm->Write();
                canva_stack_norm->Print(Form("%s/%s.png", sDirname.Data(), canva_stack_norm->GetName()));


// return;
/*  
// -- Combined Dijet + bjet
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
        canvatotal->BuildLegend(0.7, 0.9, 0.7, 0.9);
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
            /*
            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
           */
            // set to zero for test 
              TH1D *h_nob = (TH1D*) h_data_mb->Clone("h_nob_zeroinputfortest"); h_nob->Reset();

           //  Make slices for Bjet 
            TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            

      
            // -- Deattach hists from root file 
            h_data_mb->SetDirectory(0);
            // h_b->SetDirectory(0);
            // h_bb->SetDirectory(0);
            h_b_bjet->SetDirectory(0);
            h_bb_bjet->SetDirectory(0);
            // h_nob->SetDirectory(0);
            // h_more2b->SetDirectory(0);

            TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);

            // Calculate true fractions to be used as initial values for the fit (the true fractions are the dijet qcd ones)            
            // double int0 = h_bb->Integral(1, mb_bins, "width");// sig
            // double int1 = h_b->Integral(1, mb_bins, "width");// bkg1
            // double int2 = h_more2b->Integral(1, mb_bins, "width"); // bkg2
            // double int3 = h_nob->Integral(1, mb_bins, "width"); // bkg 3

            // From bjet sample 
            double int0_bjet = h_bb_bjet->Integral(1, mb_bins, "width");// sig
            double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");// bkg1

            // For MC test case only
            // double integral_inputdata1 = h_data_mb->Integral(1, mb_bins,"width");
            // double integral_inputdata2 = h_data_bjet_mb->Integral(1, mb_bins,"width");
                // cout << "Before Samples additions:\n Dijet data integral = "<< integral_inputdata1 << endl;
                // cout << "Bjet data integral = "<< integral_inputdata2 << endl;

            // Add Dijet and bjet sample for the pseudodata 
                // h_data_mb->Add(h_data_bjet_mb);

            // For data: After the sum of dijet + bjet 
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
                
                cout << "Total input data integral = "<< integral_inputdata << endl;
                // std::cout << "int0 h_bb=" << int0 << std::endl;
                // std::cout << "int1 h_b =" << int1 << std::endl;
                // cout << "2B/1B in Dijet sample = " << int0/int1 << endl;
                // cout << "Dijet: 2B/(1B + 2B) = " << int0/(int0 + int1) << endl;

                std::cout << "int0 h_bb_bjet=" << int0_bjet << std::endl;
                std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                cout <<  "2B/1B in  Bjet sample = " << int0_bjet/int1_bjet << endl;
                // cout << "2B/(2B + 1B) in qcd sample = " << int0/(int0 + int1) << endl;
                // std::cout << "int2 h_more2b =" << int2 << std::endl;
                // std::cout << "int3 h_nob =" << int3 << std::endl;


            // --- For now: the TRUE signal (2B) fraction is used from the Dijet sample.
                // the fraction of 2B/1B in Bjet and Dijet sample is exactly the same
            // double sig_fraction_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int2 + int3 )); 
            // double bkg_fraction_b_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int2 + int3 ));
            

            // No contiubtion from hmore2B 
            /*
            double sig_fraction_true = (int0  + int1 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int3 )); 
            double bkg_fraction_b_true = (int0  + int1 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int3 ));
                std::cout << "Dijet : 2B sig_fraction_true =" << sig_fraction_true << std::endl;
                std::cout << "Dijet: 1B bkg_fraction_b_true =" << bkg_fraction_b_true << std::endl;
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

                // -- and their histograms for next normalization 
                TH1D* h_sig = (TH1D*) h_bb_bjet ->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); //  
                TH1D* h_bkg = (TH1D*) h_b_bjet ->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); //
                    Int_t h_sig_bins = h_sig->GetNbinsX();
                    Int_t h_bkg_bins = h_bkg->GetNbinsX();




            /*
            // -- Get Normalized templates as input  PDF to the fit(shape matters only): where PDF integral is 1 (Probability density function)
            //signal: 2B from Dijet + Bjet sample 
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
            h_sumsig->Add(h_bb_bjet);
            Int_t h_sig_bins = h_sumsig->GetNbinsX();
            TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized
            //Add background: here only from 1B 
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
            h_sumbkg->Add(h_b_bjet);
            Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
            TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized
            // Total sum bkg = 0B + 1B  
            TH1D* h_sumbkg_0b_1b = (TH1D*)  h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
            h_sumbkg_0b_1b->Add(h_nob);
            // write signal hists after combine Dijet and bjet sample 
            fout->cd();
            h_sumsig->Write();
            h_sumbkg->Write();
            h_sumbkg_0b_1b->Write();
            // -- Write data hist 
            h_data_mb->Write();
           
            // -- Draw Stack for combined: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before ormalization 
              THStack hstack_templatesforfit(Form("hstack_templatesforfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesforfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesforfit.Add(h_nob);
                    hstack_templatesforfit.Add(h_sumbkg);
                    hstack_templatesforfit.Add(h_sumsig);
                    h_sumsig->SetTitle("2B: Dijet + bjet");
                    h_sumbkg->SetTitle("1B: Dijet + bjet");
                    h_nob->SetTitle("0B");
                    // -- Set fill colors 
                    h_nob->SetFillColor(kRed);
                    h_sumsig->SetFillColor(kBlue -10);
                    h_sumbkg->SetFillColor(kRed-10);
                    h_sumbkg->SetFillStyle(3244);

              auto canva = new TCanvas(Form("All_templates_%d_%d", ibin_dr, ibin_pt),Form("Templates(mb) and Full MC"), 800, 800 );
                    canva->cd();
                    canva->SetLogy();
                    // hstack_templatesforfit.GetYaxis()->SetRange(0, 2e+08);
                    hstack_templatesforfit.Draw("hist E"); 
                    // h_data_mb->SetTitle("Data");
                    // h_data_mb->Draw("PE same");
                    gPad->Modified();   
                    gPad->Update();
                    canva->Modified();
                    canva->Update();
                    fout->cd();
                    canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
                    canva->Write();
                    canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));
*/

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
            
            // -- test bjet, no 0B
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
                h_sig_fit->SetMarkerColor(kRed+2);
                h_sig_fit->SetLineColor(kRed+2);
                h_sig_fit->SetFillColor(kRed+2);
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg 
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("After fit: bkg (1B+ 0B)");
                h_bkg_fit->SetMarkerColor(kCyan+2);
                h_bkg_fit->SetLineColor(kCyan+2);

            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only 
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("After fit: 1B");
                h_bkg_fit_1b->SetMarkerColor(kOrange+2);
                h_bkg_fit_1b->SetLineColor(kOrange+2);
                h_bkg_fit_1b->SetLineWidth(2);
                h_bkg_fit_1b->SetFillColor(kOrange+2);

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only 
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("After fit: 0B");
                h_bkg_fit_nob->SetMarkerColor(kYellow+2);
                h_bkg_fit_nob->SetLineColor(kYellow+2);
                h_bkg_fit_nob->SetLineWidth(2);

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("After fit: sig + bkg");
                    h_total_fit->SetMarkerColor(kGreen+2);
                    h_total_fit->SetLineColor(kGreen+2);
                    h_total_fit->SetMarkerStyle(22);
                    h_total_fit->SetLineWidth(2);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();
            h_bkg_fit->Write();
            h_total_fit->Write();

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
            auto canva_afterfit = new TCanvas(sname_canvas_afterfit,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetMaximum(1.2 * h_data_mb->GetMaximum());
                h_data_mb->SetTitle("Data");
                // h_sumbkg->SetTitle("1B (bkg) template (combined)");
                // h_sumsig->SetTitle("2B (sig) template (combined)");
                // h_sumbkg_0b_1b->SetTitle("Bkg template (0B + 1B) (combined)");
                // h_nob->SetTitle("0B template (Dijet only)");

                // h_sumbkg->SetLineColor(kOrange);
                // h_nob->SetLineColor(kMagenta);
                h_data_mb->SetLineWidth(2);


                // -- Make better visual
                // h_bb->SetLineWidth(3);
                h_bb_bjet->SetLineWidth(3);
                // h_sumsig->SetLineWidth(3);
                h_sig_fit->SetLineWidth(3);

                // h_b->SetLineStyle(9); // close dashes
                h_b_bjet->SetLineStyle(9);
                // h_sumbkg_0b_1b->SetLineStyle(7);
                // h_sumbkg_0b_1b->SetLineWidth(3);
                h_bkg_fit->SetLineWidth(2);
                h_bkg_fit->SetLineStyle(10);  
                
                // Add after fit templates
                h_data_mb->Draw("P E");
                h_total_fit->Draw("P E SAME");

                // h_bb->Draw("HIST E SAME");
                // h_bb_bjet->Draw("HIST E SAME");

                // h_sumsig->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                // h_sumbkg->Draw("HIST E SAME"); // Combined sample for 1B bkg 
                // h_nob->Draw("HIST E SAME");
                // h_bkg_fit_nob->Draw("HIST E SAME");
                h_bkg_fit_1b->Draw("HIST E SAME");
                // h_sumbkg_0b_1b->Draw("HIST E SAME");
                h_bkg_fit->Draw("HIST E SAME");

                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
               



        // -- After fit: stacked, and normalized 
               
            THStack hstack_afterfit("hstack_afterfit","Mass stacked histogram");
                hstack_afterfit.SetTitle(";m_{2B} [GeV];");
                hstack_afterfit.Add(h_sig_fit);
                hstack_afterfit.Add(h_bkg_fit_1b);
                hstack_afterfit.SetMaximum(1.01* h_data_mb->GetMaximum());

            auto canva_stack_afterfit = new TCanvas("All_templates_Data_stacked_afterfit",Form(""), 800, 800 );
                canva_stack_afterfit->cd();
                hstack_afterfit.Draw("hist E"); 
                h_data_mb->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_afterfit->Modified();
                canva_stack_afterfit->Update();

            // normalized version: data and templates afterfit
            THStack hstack_norm_afterfit("hstack_norm_afterfit","Mass stacked histogram");
                hstack_norm_afterfit.SetTitle(";m_{2B} [GeV];");
    
                TH1D* h_sig_fit_normstack = (TH1D*)  h_sig_fit->Clone("h_sig_fit_normstack");  h_sig_fit_normstack ->Scale(1./integral_inputdata);
                TH1D* h_bkg_fit_1b_normstack = (TH1D*)  h_bkg_fit_1b->Clone("h_bkg_fit_1b_normstack");  h_bkg_fit_1b_normstack ->Scale(1./integral_inputdata);
            
                hstack_norm_afterfit.Add(h_sig_fit_normstack);
                hstack_norm_afterfit.Add(h_bkg_fit_1b_normstack);
                hstack_norm_afterfit.SetMaximum(1.5 * h_bkg_fit_1b_normstack->GetMaximum());
                
                auto canva_stack_norm_afterfit = new TCanvas("All_templates_Data_stacked_norm_afterfit",Form(""), 800, 800 );
                canva_stack_norm_afterfit->cd();
                hstack_norm_afterfit.Draw("hist E"); 
                hdata_selfnorm->Draw("PE same"); 

                gPad->Modified();   
                gPad->Update();
                canva_stack_norm_afterfit->Modified();
                canva_stack_norm_afterfit->Update();

            // -- for integrated bin 
                TString srangeDeltaR = Form(" %g < #DeltaR [bin %d] < %g", dr_binsVector[ibin_dr -1], ibin_dr ,dr_binsVector[ibin_dr]);
                if(!ibin_dr) srangeDeltaR = Form(" %g < #DeltaR [Full range] < %g", dr_binsVector[1], dr_binsVector[bins_dr]);
             
            //  -- Build legend and write plots     
                fout->cd();
                canva_afterfit->BuildLegend(0.7, 0.9, 0.7, 0.9, srangeDeltaR);
                // auto textlatex = new TLatex();
                // textlatex->DrawLatexNDC(0.7, 0.6, "pt range"); 
                canva_afterfit->Write();
                canva_afterfit->Print(Form("%s/%s_%d_%d.png", sDirname.Data(), canva_afterfit->GetName(), ibin_dr, ibin_pt));


                canva_stack_afterfit->BuildLegend(0.7, 0.9, 0.7, 0.9, srangeDeltaR);
                // textlatex->DrawLatexNDC(0.7, 0.6, "Some text here"); 
                canva_stack_afterfit->Write();
                canva_stack_afterfit->Print(Form("%s/%s_%d_%d.png", sDirname.Data(), canva_stack_afterfit->GetName(), ibin_dr, ibin_pt));


                canva_stack_norm_afterfit->BuildLegend(0.7, 0.9, 0.7, 0.9, srangeDeltaR);
                // textlatex->DrawLatexNDC(0.7, 0.6, "Some text here"); 
                canva_stack_norm_afterfit->Write();
                canva_stack_norm_afterfit->Print(Form("%s/%s_%d_%d.png", sDirname.Data(), canva_stack_norm_afterfit->GetName(), ibin_dr, ibin_pt));

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
        
    fout->Close();

}
// void do_template_fit_bjetonly(TString &dataset, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name, bool& alsoLowEG){
//     /*
//     This fit use sig and bkg templates from the input files.
//     */
//     // Dataset: for now HGtrigered data. 
//     //-- Samples: Dijet + bjet 
//     // -- Template fit pdfs : includes (0B, 1B) as bkg (added effectively), 2B as sig  

//     bool isSampled = false; // for templates 
    
//     // -- For output histograms 
//     TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");

//     // -- Using new naming convention
//     TString namehData = "h_count_data";
//     TString namehMC = "";

//     // -- Input Histogram name: signal and bkg templates
//     TString nameh2B = "h_count_bb";
//     TString nameh1B = "h_count_b";

//     TString nameh0B = "";
//     // -- Additional histograms
//     TString namehmore2B = "";

 
//     //-- Dijet sample:
//     TFile *file_dijet = new TFile(folder + templates, "read");
//         if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());}
// /*   
//     TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b"); // 
//     TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
//         if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
//         if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
//     // -- 0B and more B contribution for drawing only
//     TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone(nameh0B);
//        if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}
// */

//     //-- Bjet sample: Templates: Open file and get histograms from bjet: templates 
//     TFile *file_bjet = new TFile(templates_bjet, "read");
//          if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}
//         TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
//         TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
//         if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
//         if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}

// /*
//     // -- Additional contrinution to draw, not to fit 
//         TH3D *h3D_more2b = (TH3D*)file_dijet->Get(namehmore2B)->Clone(namehmore2B);
//             if(!h3D_more2b){Error("Get:", "histogram does not exist '%s' ",h3D_more2b->GetName());}
//         TH3D *h3D_more2b_bjet = (TH3D*) file_bjet->Get(namehmore2B)->Clone(Form("%s_bjet", namehmore2B.Data() ) );
//             if(!h3D_more2b_bjet){Error("Get:", "histogram does not exist '%s' ", h3D_more2b_bjet->GetName());}

// */
//     //Open dataset:
//     TFile *file_data = new TFile(dataset, "read");
//         if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
//         cout << "file data name " << file_data->GetName() << endl;
//         TH3D *h3D_data;
//         h3D_data = (TH3D*)file_data->Get(namehData)->Clone("h3D_data");
//         if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}

//         // -- Add also the lower triggerd data
//         if(alsoLowEG){
//             TString LowEGdata_name = "Prepare_Data/Data_Hist/withunder_overFlow_includedibins/hist_3d_4d_aggr_BDT_n1_data_LowEGJet_80_140.root";
//             TFile *file_data_LowEG = new TFile(LowEGdata_name, "read");
//             if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
//             cout << "file data name " << file_data_LowEG->GetName() << endl;
//             TH3D *h3D_dataLowEG;
//             h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
//             if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ",h3D_dataLowEG->GetName()); return;}
//             h3D_data->Add(h3D_dataLowEG);
//         }

//         h3D_data->SetMarkerStyle(20);
//         h3D_data->SetMarkerColor(0);
//         h3D_data->SetLineColor(0);

//         // -- For MC trivial test 
//         TH3D *h3D_MC_pseudo;
//         TH3D *h3D_MC_pseudo_bjet;
//         /*
//         if(isSampled){
//             h3D_MC_pseudo = (TH3D*)file_dijet->Get(namehMC)->Clone(namehMC);
//             h3D_MC_pseudo_bjet = (TH3D*)file_bjet->Get(namehMC)->Clone("h3D_data_bjet");
//         }
//         else {
//             // Trivial closure test
            
//             h3D_MC_pseudo = (TH3D*) h3D_b->Clone("h3D_pseudodata");
//             h3D_MC_pseudo->Add(h3D_bb);
//             h3D_MC_pseudo->Add(h3D_nob);
            
//             // then add bjet sample also
//             h3D_MC_pseudo_bjet = (TH3D*) h3D_b_bjet->Clone("h3D_data_bjet");
//             h3D_MC_pseudo_bjet->Add(h3D_bb_bjet);
//         }
//         if(!h3D_MC_pseudo){Error("Get:", "histogram does not exist '%s' ",h3D_MC_pseudo->GetName());}
//         if(!h3D_MC_pseudo_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_MC_pseudo_bjet->GetName());}
//         */

//     //Get number of bin entries
//     Int_t bins_pt = h3D_data->GetNbinsZ();
//     Int_t bins_dr = h3D_data->GetNbinsY(); 
//     Int_t mb_bins = h3D_data->GetNbinsX();
//     cout << "-- Data hist binning" << endl;
//     cout << "pt bins = "<< bins_pt << endl;
//     cout << "dr bins = "<< bins_dr << endl;
//     cout << "mb bins = "<< mb_bins << endl;

//     // -- Draw THstach for mb projection 
//     // Get projections X 
//     // hists names as for slices in bins (but outside {} so it is ok)
//         TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX("h_data_mb", 1,bins_dr, 1, bins_pt ); // all bins 
//             h_data_mb->SetTitle("Data");
//             h_data_mb->SetLineColor(kBlack);
//             h_data_mb->SetLineWidth(3);
//             h_data_mb->SetMarkerColor(kBlack);
//             h_data_mb->SetMarkerStyle(4);
//         // TH1D *h_b = (TH1D *) h3D_b->ProjectionX("h_b");
//         TH1D *h_b_bjet = (TH1D *) h3D_b_bjet->ProjectionX("h_b_bjet");
//         // TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX("h_bb");
//         TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet->ProjectionX("h_bb_bjet");
//         // TH1D *h_nob = (TH1D *) h3D_nob->ProjectionX("h_nob");
//         // TH1D *h_more2b = (TH1D *) h3D_more2b->ProjectionX("h_more2b");
//         // TH1D *h_more2b_bjet = (TH1D *) h3D_more2b_bjet ->ProjectionX("h_more2b_bjet");
//             /*h_more2b->SetDirectory(0);
//             h_more2b_bjet->SetDirectory(0);
//             h_more2b->SetFillColor(kOrange -2 ); 
//             h_more2b->SetLineColor(kOrange-2);
//             h_more2b_bjet->SetFillColor(kOrange -2 );
//             h_more2b_bjet->SetLineColor(kOrange-2);
//             h_more2b->SetTitle(">2B (Dijet)");
//             h_more2b_bjet->SetTitle(">2B (bjet)");
//             */

//         // -- Set titles for the hisograms  to be used in the default legends 
//         // h_b->SetTitle("1B (dijet)");
//         h_b_bjet->SetTitle("1B (bjet)");
//         // h_bb->SetTitle("2B (dijet)");
//         h_bb_bjet->SetTitle("2B (bjet)");
//         // h_nob->SetTitle("0B");

//         /*
//         // The trivial sum of sig+bkg t compare to full MC 
//         TH1D *hsum = (TH1D*) h_b->Clone("hsum");
//         hsum->SetTitle("BB signal + 1B bkg + 0B bkg");
//         hsum->Add(h_b_bjet);
//         hsum->Add(h_bb);
//         hsum->Add(h_bb_bjet);
//         hsum->Add(h_nob);
//         hsum->SetLineWidth(3);
//         hsum->SetLineColor(kMagenta+2);
//         hsum->SetDirectory(0);
//         fout->cd();
//         hsum->Write();
//         */
       

//        // -- Deattach hists from root file 
//         h_data_mb->SetDirectory(0);
//         // h_b->SetDirectory(0);
//         // h_bb->SetDirectory(0);
//         h_b_bjet->SetDirectory(0);
//         h_bb_bjet->SetDirectory(0);
//         // h_more2b->SetDirectory(0);
//         // h_nob->SetDirectory(0);
//         // Set styles of hist to be drawn:
//         // h_b->SetFillColor(kRed-10); h_b->SetLineColor(kRed-10);        
//         // h_bb->SetFillColor(kBlue-10); h_bb->SetLineColor(kBlue-10);
//         h_b_bjet->SetFillColor(kGreen-2); h_b_bjet->SetLineColor(kGreen-2);        
//         h_bb_bjet->SetFillColor(kOrange-6); h_bb_bjet->SetLineColor(kOrange-6); // brown 
//         // h_more2b->SetFillColor(kOrange); h_more2b->SetLineColor(kOrange);
//         // h_nob->SetFillColor(kRed); h_nob->SetLineColor(kRed);
       

//     THStack hstack_All("hstack_All","Mass stacked histogram");
//         hstack_All.SetTitle(";m_{2B} [GeV];");
//         // hstack_All.Add(h_more2b);
//         // hstack_All.Add(h_more2b_bjet);
//         // hstack_All.Add(h_nob);
//         // hstack_All.Add(h_b);
//         // hstack_All.Add(h_bb);
//         hstack_All.Add(h_b_bjet);
//         hstack_All.Add(h_bb_bjet);
//         // hstack_All.SetMaximum(20e+3);

//     auto canva = new TCanvas("All_templates_Data",Form("Templates(mb) and Full MC"), 800, 800 );
//         canva->cd();
//         canva->SetLogy();
//         hstack_All.Draw("hist E"); // templates stack
//         // -- try draw scaled down data histogram?
//         TH1D* hd_norm = (TH1D*) h_data_mb->Clone("hd_norm");
//             hd_norm->Scale(1./10000);
//             hd_norm->Draw("HIST E same"); 
//             hd_norm->SetTitle("data/10000");
//         // h_data_mb->Draw("HIST E same"); // data points 
//         // hsum->Draw(" Hist same"); // templates sum          
//         gPad->Modified();   
//         gPad->Update();
//         canva->Modified();
//         canva->Update();
//         fout->cd();
//         canva->BuildLegend(0.8, 0.9, 0.8, 0.9);
//         canva->Write();
//         canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));

// // return;
// /*  
// // -- Combined Dijet + bjet
//     // Another plot for the distributions  before slicing, but combined samples 
//         auto canvatotal = new TCanvas("All_templates_combined",Form(""), 800, 800 );
//         canvatotal->cd();
//         // canvatotal->SetLogy();
//         // Get combined versions + their normalization
//         TH1D* normhc_nob = (TH1D*) h_nob->Clone("normhc_nob"); normhc_nob->Scale(1./normhc_nob->Integral());
//         TH1D* hc_1b = (TH1D*) h_b->Clone("hc_1b"); hc_1b->Add(h_b_bjet); hc_1b->SetTitle("1B");
//             TH1D* normhc_1b = (TH1D*) hc_1b->Clone("normhc_1b");   normhc_1b->Scale(1./normhc_1b->Integral());
//         TH1D* hc_2b = (TH1D*) h_bb ->Clone("hc_2b"); hc_2b->Add(h_bb_bjet); hc_2b->SetTitle("2B");
//             TH1D* normhc_2b = (TH1D*) hc_2b->Clone("normhc_2b");   normhc_2b->Scale(1./normhc_2b->Integral());

//         // - draw stack 
//          THStack hstack_All_combined("hstack_All_combined","Mass stacked histogram");
//                 hstack_All_combined.SetTitle(";m_{2B} [GeV];");
//                 hstack_All_combined.Add(normhc_nob);
//                 hstack_All_combined.Add(normhc_2b);
//                 hstack_All_combined.Add(normhc_1b);
//         hstack_All_combined.Draw("hist E"); // templates stack
//         // -- normalizaed data 
//         TH1D* normhc_data = (TH1D*) h_data_mb->Clone("hd_norm");
//             normhc_data->Scale(1./normhc_data->Integral());
//             normhc_data->Draw("HIST E same"); 
//             normhc_data->SetTitle("data");
//         canvatotal->SetTitle("Normalized distributions");         
//         gPad->Modified();   
//         gPad->Update();
//         canvatotal->Modified();
//         canvatotal->Update();
//         fout->cd();
//         canvatotal->BuildLegend(0.7, 0.9, 0.7, 0.9);
//         canvatotal->Write();
//         canvatotal->Print(Form("%s/%s.png", sDirname.Data(), canvatotal->GetName()));
// */

// // return;

//     // -- Output hist declarations: empty histograms 
//         // -- For the fit result: per (deltaR, jtpt) bin
//         //Create the histograms storing the fit parameters in TH2D (the dr are stored in the X, the jtpt in the Y)
//         TH2D *h_sig_fraction = (TH2D*) h3D_data->Project3D("zy")->Clone("h_sig_fraction");// just to get the axes from hist
//         h_sig_fraction->SetTitle("x=#DeltaR, y=jtpt, fraction"); 
//         h_sig_fraction->Reset();
//         TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
//         TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
//         TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error");
//         //Store the true parameters and errors
//         TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
//         TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");
//         TH2D *h_bkg_frac_true = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true");
//         TH2D *h_bkg_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true_error");

//         // -- And for the integarted bin in deltaR, Fitparaemter Vs. ptbin 
//           TH1D *h_sig_fraction_DeltaRIntBin = (TH1D*) h3D_data->ProjectionZ("h_sig_fraction_DeltaRIntBin"); h_sig_fraction_DeltaRIntBin->Reset();
//           TH1D *h_bkg_fraction_DeltaRIntBin = (TH1D *) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_fraction_DeltaRIntBin");
//           TH1D *h_sig_true_fraction_DeltaRIntBin = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_sig_true_fraction_DeltaRIntBin");
//           TH1D *h_bkg_true_fraction_DeltaRIntBin = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_true_fraction_DeltaRIntBin");
//           TH1D *h_sig_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_sig_fraction_DeltaRIntBin_error");
//           TH1D *h_bkg_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_fraction_DeltaRIntBin_error");
//           TH1D *h_sig_true_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_sig_true_fraction_DeltaRIntBin_error");
//           TH1D *h_bkg_true_fraction_DeltaRIntBin_error = (TH1D*) h_sig_fraction_DeltaRIntBin->Clone("h_bkg_true_fraction_DeltaRIntBin_error");
          
//     // --- Vector to test the convergence
//     std::vector <std::pair<int, int>> non_converge_bins;


//     // Fitting - loop over dr and jtpt entries
//     for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
//     // for(Int_t ibin_pt = 1; ibin_pt <= 1; ibin_pt++){
//         // for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){

//         for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){
//         // for(Int_t ibin_dr = 0; ibin_dr <= 0; ibin_dr++){
//         // Bin = 0: integaretd deltaR range, other bins are the DeltaR bins 
            
//             // The slice range for projection is different for integarted Vs. one bin slice
//             Int_t SliceFirstbin_dr = ibin_dr;
//             Int_t SliceLastbin_dr =  ibin_dr;
           

//             // For integarted bin only: bin = 0
//             if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = bins_dr;}
//              // -- test rebinned dr in 2 
//             /*else{
               
//                 SliceLastbin_dr =  ibin_dr + 1;
//                 ibin_dr++; // to make step of 2  
//             }*/

//             // Make slices for data
//             TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
//                 h_data_mb->SetMarkerColor(kBlack);
//                 h_data_mb->SetMarkerStyle(20);
//             /*
//             // Make slices for dijet
//             TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
//             TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
//             TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
//            */
//             // set to zero for test 
//               TH1D *h_nob = (TH1D*) h_data_mb->Clone("h_nob_zeroinputfortest"); h_nob->Reset();

//            //  Make slices for Bjet 
//             TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
//             TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            

      
//             // -- Deattach hists from root file 
//             h_data_mb->SetDirectory(0);
//             // h_b->SetDirectory(0);
//             // h_bb->SetDirectory(0);
//             h_b_bjet->SetDirectory(0);
//             h_bb_bjet->SetDirectory(0);
//             // h_nob->SetDirectory(0);
//             // h_more2b->SetDirectory(0);

//             TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);

//             // Calculate true fractions to be used as initial values for the fit (the true fractions are the dijet qcd ones)            
//             // double int0 = h_bb->Integral(1, mb_bins, "width");// sig
//             // double int1 = h_b->Integral(1, mb_bins, "width");// bkg1
//             // double int2 = h_more2b->Integral(1, mb_bins, "width"); // bkg2
//             // double int3 = h_nob->Integral(1, mb_bins, "width"); // bkg 3

//             // From bjet sample 
//             double int0_bjet = h_bb_bjet->Integral(1, mb_bins, "width");// sig
//             double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");// bkg1

//             // For MC test case only
//             // double integral_inputdata1 = h_data_mb->Integral(1, mb_bins,"width");
//             // double integral_inputdata2 = h_data_bjet_mb->Integral(1, mb_bins,"width");
//                 // cout << "Before Samples additions:\n Dijet data integral = "<< integral_inputdata1 << endl;
//                 // cout << "Bjet data integral = "<< integral_inputdata2 << endl;

//             // Add Dijet and bjet sample for the pseudodata 
//                 // h_data_mb->Add(h_data_bjet_mb);

//             // For data: After the sum of dijet + bjet 
//             double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
                
//                 cout << "Total input data integral = "<< integral_inputdata << endl;
//                 // std::cout << "int0 h_bb=" << int0 << std::endl;
//                 // std::cout << "int1 h_b =" << int1 << std::endl;
//                 // cout << "2B/1B in Dijet sample = " << int0/int1 << endl;
//                 // cout << "Dijet: 2B/(1B + 2B) = " << int0/(int0 + int1) << endl;

//                 std::cout << "int0 h_bb_bjet=" << int0_bjet << std::endl;
//                 std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
//                 cout <<  "2B/1B in  Bjet sample = " << int0_bjet/int1_bjet << endl;
//                 // cout << "2B/(2B + 1B) in qcd sample = " << int0/(int0 + int1) << endl;
//                 // std::cout << "int2 h_more2b =" << int2 << std::endl;
//                 // std::cout << "int3 h_nob =" << int3 << std::endl;


//             // --- For now: the TRUE signal (2B) fraction is used from the Dijet sample.
//                 // the fraction of 2B/1B in Bjet and Dijet sample is exactly the same
//             // double sig_fraction_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int2 + int3 )); 
//             // double bkg_fraction_b_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int2 + int3 ));
            

//             // No contiubtion from hmore2B 
//             /*
//             double sig_fraction_true = (int0  + int1 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int3 )); 
//             double bkg_fraction_b_true = (int0  + int1 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int3 ));
//                 std::cout << "Dijet : 2B sig_fraction_true =" << sig_fraction_true << std::endl;
//                 std::cout << "Dijet: 1B bkg_fraction_b_true =" << bkg_fraction_b_true << std::endl;
//             // -- Compute the errors of true fractions for later fill
//             double True_bkg_b_err = 0.0;
//             double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");
//             double True_sig_err = 0.0;
//             double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");
//             */
            
//             // -- test bjet only as true signal 
//             double sig_fraction_true = (int0_bjet  + int1_bjet ) == 0 ? 0 : (int0_bjet / (int0_bjet + int1_bjet)); 
//             double bkg_fraction_b_true = (int0_bjet  + int1_bjet ) == 0 ? 0 : (int1_bjet / (int0_bjet + int1_bjet));
//                 std::cout << "bjet : 2B sig_fraction_true =" << sig_fraction_true << std::endl;
//                 std::cout << "bjet: 1B bkg_fraction_b_true =" << bkg_fraction_b_true << std::endl;  
//                 double True_bkg_b_err = 0.0;
//                 double True_bkg_b_integral = h_b_bjet->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");
//                 double True_sig_err = 0.0;
//                 double True_sig_integral = h_bb_bjet->IntegralAndError(1, mb_bins, True_sig_err, "width");

//                 // -- and their histograms for next normalization 
//                 TH1D* h_sig = (TH1D*) h_bb_bjet ->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); //  
//                 TH1D* h_bkg = (TH1D*) h_b_bjet ->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); //
//                     Int_t h_sig_bins = h_sig->GetNbinsX();
//                     Int_t h_bkg_bins = h_bkg->GetNbinsX();




//             /*
//             // -- Get Normalized templates as input  PDF to the fit(shape matters only): where PDF integral is 1 (Probability density function)
//             //signal: 2B from Dijet + Bjet sample 
//             TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
//             h_sumsig->Add(h_bb_bjet);
//             Int_t h_sig_bins = h_sumsig->GetNbinsX();
//             TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized
//             //Add background: here only from 1B 
//             TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
//             h_sumbkg->Add(h_b_bjet);
//             Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
//             TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized
//             // Total sum bkg = 0B + 1B  
//             TH1D* h_sumbkg_0b_1b = (TH1D*)  h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
//             h_sumbkg_0b_1b->Add(h_nob);
//             // write signal hists after combine Dijet and bjet sample 
//             fout->cd();
//             h_sumsig->Write();
//             h_sumbkg->Write();
//             h_sumbkg_0b_1b->Write();
//             // -- Write data hist 
//             h_data_mb->Write();
           
//             // -- Draw Stack for combined: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before ormalization 
//               THStack hstack_templatesforfit(Form("hstack_templatesforfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
//                     hstack_templatesforfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
//                     hstack_templatesforfit.Add(h_nob);
//                     hstack_templatesforfit.Add(h_sumbkg);
//                     hstack_templatesforfit.Add(h_sumsig);
//                     h_sumsig->SetTitle("2B: Dijet + bjet");
//                     h_sumbkg->SetTitle("1B: Dijet + bjet");
//                     h_nob->SetTitle("0B");
//                     // -- Set fill colors 
//                     h_nob->SetFillColor(kRed);
//                     h_sumsig->SetFillColor(kBlue -10);
//                     h_sumbkg->SetFillColor(kRed-10);
//                     h_sumbkg->SetFillStyle(3244);

//               auto canva = new TCanvas(Form("All_templates_%d_%d", ibin_dr, ibin_pt),Form("Templates(mb) and Full MC"), 800, 800 );
//                     canva->cd();
//                     canva->SetLogy();
//                     // hstack_templatesforfit.GetYaxis()->SetRange(0, 2e+08);
//                     hstack_templatesforfit.Draw("hist E"); 
//                     // h_data_mb->SetTitle("Data");
//                     // h_data_mb->Draw("PE same");
//                     gPad->Modified();   
//                     gPad->Update();
//                     canva->Modified();
//                     canva->Update();
//                     fout->cd();
//                     canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
//                     canva->Write();
//                     canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));
// */

//                 // ----- To avoid empty bins if exist!, set them to eps value
//                 const double eps = 1e-6; 
//                 for (int i = 1; i <= h_sig_bins; i++){
//                     if (h_sig->GetBinContent(i) <= 0){
//                         h_sig->SetBinContent(i, eps);
//                         cout << 
//                         Form("INFO: (ptbin %d, deltaRbin %d) has empty  signal Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i) 
//                         << endl; 
//                     }
//                 }

//                 for (int i = 1; i <= h_bkg_bins; i++){
//                     if (h_bkg->GetBinContent(i) <= 0){
//                         h_bkg->SetBinContent(i, eps);
//                         cout <<
//                         Form("INFO: (ptbin %d, deltaRbin %d) has empty bkg  Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i)  
//                         << endl; 
//                     }
//                 }
//                 /// Normalize safely for the Roofit
//                 if(h_sig->Integral(1, h_sig_bins, "width") == 0.0 || h_bkg->Integral(1, h_bkg_bins, "width") == 0.0)
//                 {
//                     cerr<< "Signal or bkg templates has zero integral"<< endl;
//                     continue;
//                 }

//                 h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));
//                 h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));
//                 cout << "After Normalization \n True 2B integral = " << h_sig->Integral(1, h_sig_bins, "width") << endl;
//                 cout << "True 1B integral with width option =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;
//                 cout << "True 1B integral without width option =  "<<  h_bkg->Integral(1, h_bkg_bins) << endl;

//             /// Hay que hacer esto ???
//             // ----------- Franc. added the hrest to the hsig template ----------- 
//             //Add signal and rest with the correct scaling (from Lida)
//             // add together light+c to sig with nominal MC ratio
//             // a sig + b bb + c light = n, a+b+c=1, a+c=1-b
//             // c'=c/(a+c)=(1-a-b)/(1-b)
//             // a'=a/(a+c)=a/(1-b)=1-c'
//             // -- Afnan: what are these param?
//             // -- Not used in my case: 
//              // double cl_frac = (1 - sig_fraction_true - bkg_fraction_b_true) / (1 - bkg_fraction_b_true); // nominal
//             // cl_frac = 2*cl_frac; // systematic
//             ////// h_sig->Add(h_sig, h_rest, sig_frac, 1-sig_frac);

            
//             /*
//             // If no other than 2B and 1B  are used: 
//             /////// For now: since i have pure signal, use charm-light contribution = 0
//             double cl_frac = 0.; // systematic
//             double sig_frac = 1 - cl_frac;
//             */

//             // -- Effective bkg PDF 
//             // Assume : a + b + c = 1; where a, b, c, are fractions in in nominal MC like 
//             // a: 2B fraction, b: 1B fraction, c: 0B fraction, in 
//             // and For Effective bkg: b`+ c` should = 1, so relate b`to b .. and simialr for c`
//             // This result in: b`= b/(b+c), and c`= c/(b+c); one can use b+c = 1-a;  
//             // I will write it in terms of the 2B sgnal, 1B bkg true fractions, before use effecive ones 
//             double eff_bkg0B = (1- sig_fraction_true - bkg_fraction_b_true)/(1- sig_fraction_true);// c`
//             double eff_bkg1B = 1. - eff_bkg0B;// b`
            
//             // -- test bjet, no 0B
//                eff_bkg0B = 0;
//                eff_bkg1B = 1; 

//             // Build effective signal hist: with new relaitve normalization, the integral should = 1
//             h_bkg->Add(h_bkg, h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist
            
//             double err_int;
//             double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins, err_int ,"width");
//             cout << "Effective Bkg normalized hist integral +/- uncertainity = "<< int_val << "+/-" << err_int << endl;
//             cout << "effective 1B fraction and 0B = "<< eff_bkg1B << ", " << eff_bkg0B << endl;
//             // -- Note: the effective bkg integral is Effective Bkg normalized hist integral = 0.997777
//             // if (h_bkg->Integral(1,  h_bkg_bins, "width") != 1.) {cout << "Effective BKG PDF is not normalized!"<< endl; return;}

//             ///// Fitting
//             // Create the observable
//             Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
//             Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
//             RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
//             mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
     
//             // Inputs 
//             // Create the RooDataHist object for the observed data + templates
//             RooDataHist *dh_data_mb = new RooDataHist(Form("dh_data_mb_%d_%d", ibin_dr, ibin_pt), "dh_data_mb", mb, RooFit::Import(*h_data_mb));
//             RooDataHist *dh_sig_mb = new RooDataHist(Form("dh_sig_mb_%d_%d", ibin_dr, ibin_pt), "dh_sig_mb", mb, RooFit::Import(*h_sig));
//             RooDataHist *dh_bkg_mb = new RooDataHist(Form("h_bkg_mb_%d_%d", ibin_dr, ibin_pt), "dh_bkg_mb", mb, RooFit::Import(*h_bkg));

//             // Create the RooHistPdf objects for the template PDFs
//             RooHistPdf sig_template(Form("sig_template_%d_%d", ibin_dr, ibin_pt), "sig_template", mb, *dh_sig_mb);
//             RooHistPdf bkg_template(Form("bkg_template_%d_%d", ibin_dr, ibin_pt), "bkg_template", mb, *dh_bkg_mb);

//             // Create list of templates
//             RooArgList template_list(sig_template, bkg_template, "template_list");

//             // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
//             // debug 
//             // sig_fraction_true = 0.5;
//             RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true , 0., 1);// , // 1-bkg_fraction_b_true

//             /*
//             // test fix fit par:
//             sig_fraction_val.setVal(sig_fraction_true);     // set it to zero
//             sig_fraction_val.setConstant(true); // fix it so it won't float in the fit
//             */

//             // Create the composite PDF using a linear combination of the template PDFs
//             // RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, true); 
//             //fit the scaling when adding the single b with the more b distribution

//             RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, false); // fasle: dont use extended fit: assume fraction not absolute yields  
//             RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name            
//             // -- test fit in a smaller  ragne? then fit 
//             /*
//             mb.setRange("fitRange", 4,  9); // limit the fit range 
//             RooFitResult* result = model0.fitTo(*dh_data_mb,
//                                                 RooFit::Range("fitRange"),
//                                                 RooFit::NormRange("fitRange"),
//                                                  RooFit::SumW2Error(true),
//                                                   RooFit::Save(),
//                                                    RooFit::CloneData(true),
//                                                     RooFit::PrintLevel(2),
//                                                      RooFit::Strategy(1),
//                                                       RooFit::Minos(false));

//             */

//                                                                                                                                                                                               //instead of sign histogram here we would put the data histogram
//             Int_t status = result->status();
//             result->Print();

//             std::cout << "covariance matrix:" << std::endl;
//             (result->covarianceMatrix().Print());

//             //Check if it converged for a dr and jtpt bin
//             if (status != 0) {
//                 std::cout << "\n\n\n\n!!!Fitting for ipt = " << ibin_pt 
//                          << ", ix = " << ibin_dr 
//                          << " did not converge\n\n\n\n" << std::endl;
//                 non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
//                 continue;
//             }
    
//             // Get the fitted parameter values
//             double a = sig_fraction_val.getValV();
//             double da = sig_fraction_val.getError();

//             //Print some check
//             std::cout << "RooFit result: \n Signal PDF fraction: a = " << a << " \n its uncertainity: da = " << da << std::endl;
//             /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
            
//             Double_t p0, p1, p2, errP0, errP1, errP2;
//             p0 = a; // signal 2B
//             // Total bkg: (1-a)
//             p1 = (1-a)*eff_bkg1B; // 1B bkg 
//             p2 = (1-a)*eff_bkg0B; // 0B bkg 
//             errP0 = da;
//             errP1 = da *eff_bkg1B;
//             errP2 = da *eff_bkg0B;

//             std::cout << "a: 2B =" << sig_fraction_true << ", a': after fit=" << p0 << std::endl;            
//             std::cout << "a'/a  for 2B = " << p0/sig_fraction_true << std::endl;


//             std::cout << "b: 1B =" << bkg_fraction_b_true << ", b': 1B =" << p1 << std::endl;
//             std::cout << "b'/b for 1B = " << p1/bkg_fraction_b_true << std::endl;


//             std::cout << "c: 0B =" << (1 - sig_fraction_true - bkg_fraction_b_true) << ", c': 0B after fit =" << p2 << std::endl;
//             std::cout << "c'/c = " << p2/(1-sig_fraction_true - bkg_fraction_b_true) << std::endl;



        
//             // save the fit
//              // rescaling for the signal only (not rest background)
//             h_sig_fraction->SetBinContent(ibin_dr, ibin_pt, p0);
//             h_sig_fraction->SetBinError(ibin_dr, ibin_pt, errP0);
//             h_sig_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP0);

//             h_bkg_fraction->SetBinContent(ibin_dr, ibin_pt, p1);
//             h_bkg_fraction->SetBinError(ibin_dr, ibin_pt, errP1);
//             h_bkg_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP1);

//             //save the true fraction
//             h_sig_frac_true->SetBinContent(ibin_dr, ibin_pt, sig_fraction_true);
//             h_bkg_frac_true->SetBinContent(ibin_dr, ibin_pt, bkg_fraction_b_true);

//             // -- Compute siganl and bkg fractions uncertainity 
//             // -- sigma True for S or B  = (B x sigmaS + SxsigmaB)/(S+B)²
//             // double err_true_frac = (int1 * True_sig_err + int0 * True_bkg_b_err)/(int0+ int1)/(int0+ int1);// int0 for signal, int1 for bkg 
//                 // when bjet oly 
//                      double err_true_frac = (int1_bjet * True_sig_err + int0_bjet * True_bkg_b_err)/(int0_bjet+ int1_bjet)/(int0_bjet+ int1_bjet);


//             // -- Save true fraction uncertainity
//             h_sig_frac_true_error->SetBinContent(ibin_dr, ibin_pt, err_true_frac);
//             h_bkg_frac_true_error->SetBinContent(ibin_dr, ibin_pt, err_true_frac);
            
//             // -- And for the integrated bin 
//            if(!ibin_dr){
//                 h_sig_fraction_DeltaRIntBin->SetBinContent(ibin_pt, p0);
//                 h_sig_fraction_DeltaRIntBin->SetBinError(ibin_pt, errP0);
//                 h_sig_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, errP0);

//                 h_bkg_fraction_DeltaRIntBin->SetBinContent(ibin_pt, p1);
//                 h_bkg_fraction_DeltaRIntBin->SetBinError(ibin_pt, errP1);
//                 h_bkg_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, errP1);

//                 //save the true fraction
//                 h_sig_true_fraction_DeltaRIntBin->SetBinContent(ibin_pt, sig_fraction_true);
//                 h_bkg_true_fraction_DeltaRIntBin->SetBinContent(ibin_pt, bkg_fraction_b_true);

//                 // -- Save true fraction uncertainity
//                 h_sig_true_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, err_true_frac);
//                 h_bkg_true_fraction_DeltaRIntBin_error->SetBinContent(ibin_pt, err_true_frac);
//            }
           


//             // -- After fits: save mass distribution re and post-fit 
//             TH1D *h_sig_fit = (TH1D*) h_sig->Clone(Form("h_sig_fit_%d_%d", ibin_dr, ibin_pt));
//                 h_sig_fit->Scale(p0 * integral_inputdata);
//                 h_sig_fit->SetTitle("After fit: Sig");
//                 h_sig_fit->SetMarkerColor(kRed+2);
//                 h_sig_fit->SetLineColor(kRed+2);
//                 cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

//             TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg 
//                 h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
//                 h_bkg_fit->SetTitle("After fit: bkg (1B+ 0B)");
//                 h_bkg_fit->SetMarkerColor(kCyan+2);
//                 h_bkg_fit->SetLineColor(kCyan+2);

//             TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only 
//                 h_bkg_fit_1b->Scale(integral_inputdata * p1);
//                 h_bkg_fit_1b->SetTitle("After fit: 1B");
//                 h_bkg_fit_1b->SetMarkerColor(kOrange+2);
//                 h_bkg_fit_1b->SetLineColor(kOrange+2);
//                 h_bkg_fit_1b->SetLineWidth(2);

//             TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only 
//                 h_bkg_fit_nob->Scale(integral_inputdata * p2);
//                 h_bkg_fit_nob->SetTitle("After fit: 0B");
//                 h_bkg_fit_nob->SetMarkerColor(kYellow+2);
//                 h_bkg_fit_nob->SetLineColor(kYellow+2);
//                 h_bkg_fit_nob->SetLineWidth(2);

//                 // And the post-fit template (fitted sig + bkg)
//             TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
//                   h_total_fit->Add(h_bkg_fit);
//                   h_total_fit->SetTitle("After fit: sig + bkg");
//                     h_total_fit->SetMarkerColor(kGreen+2);
//                     h_total_fit->SetLineColor(kGreen+2);
//                     h_total_fit->SetMarkerStyle(22);
//                     h_total_fit->SetLineWidth(2);


//             // Save Signal and bkg mass distributions after the fit 
//             fout->cd();
//             h_sig_fit->Write();
//             h_bkg_fit->Write();
//             h_total_fit->Write();

//             /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
//             TString sname_canvas_afterfit = sname_canvas + "_afterfit";
//             auto canva_afterfit = new TCanvas(sname_canvas_afterfit,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
//                 canva_afterfit->cd();
//                 h_data_mb->SetMaximum(1.2 * h_data_mb->GetMaximum());
//                 h_data_mb->SetTitle("Data");
//                 // h_sumbkg->SetTitle("1B (bkg) template (combined)");
//                 // h_sumsig->SetTitle("2B (sig) template (combined)");
//                 // h_sumbkg_0b_1b->SetTitle("Bkg template (0B + 1B) (combined)");
//                 // h_nob->SetTitle("0B template (Dijet only)");

//                 // h_sumbkg->SetLineColor(kOrange);
//                 // h_nob->SetLineColor(kMagenta);
//                 h_data_mb->SetLineWidth(2);


//                 // -- Make better visual
//                 // h_bb->SetLineWidth(3);
//                 h_bb_bjet->SetLineWidth(3);
//                 // h_sumsig->SetLineWidth(3);
//                 h_sig_fit->SetLineWidth(3);

//                 // h_b->SetLineStyle(9); // close dashes
//                 h_b_bjet->SetLineStyle(9);
//                 // h_sumbkg_0b_1b->SetLineStyle(7);
//                 // h_sumbkg_0b_1b->SetLineWidth(3);
//                 h_bkg_fit->SetLineWidth(2);
//                 h_bkg_fit->SetLineStyle(10);  
                
//                 // Add after fit templates
//                 h_data_mb->Draw("P E");
//                 h_total_fit->Draw("P E SAME");

//                 // h_bb->Draw("HIST E SAME");
//                 // h_bb_bjet->Draw("HIST E SAME");

//                 // h_sumsig->Draw("HIST E SAME");
//                 h_sig_fit->Draw("HIST E SAME");
//                 // h_sumbkg->Draw("HIST E SAME"); // Combined sample for 1B bkg 
//                 // h_nob->Draw("HIST E SAME");
//                 // h_bkg_fit_nob->Draw("HIST E SAME");
//                 h_bkg_fit_1b->Draw("HIST E SAME");
//                 // h_sumbkg_0b_1b->Draw("HIST E SAME");
//                 h_bkg_fit->Draw("HIST E SAME");

//                 gPad->Modified();   
//                 gPad->Update();
//                 canva_afterfit->Modified();
//                 canva_afterfit->Update();
//                 fout->cd();

//                 TString srangeDeltaR = Form(" %g < #DeltaR [bin %d] < %g", dr_binsVector[ibin_dr -1], ibin_dr ,dr_binsVector[ibin_dr]);
//                 if(!ibin_dr) srangeDeltaR = Form(" %g < #DeltaR [Full range] < %g", dr_binsVector[1], dr_binsVector[bins_dr]);
//                 canva_afterfit->BuildLegend(0.7, 0.9, 0.7, 0.9, srangeDeltaR);
//                 // Add another line for all bins 
//                 canva_afterfit->Write();


//                 canva_afterfit->Print(Form("%s/%s.png", sDirname.Data(), sname_canvas_afterfit.Data()));

//                 cout << "---------------------\n\n\n" << endl; 
//             } // loop over deltaR bins 
//     }

//     // Save histograms
//     // TH3D 
//     // for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
//     for (auto h : {h3D_data, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
        

        
//     // TH2D     
//     for (auto h : {h_sig_fraction, h_sig_fraction_error,
//                    h_bkg_fraction, h_bkg_fraction_error,
//                    h_sig_frac_true, h_sig_frac_true_error,
//                    h_bkg_frac_true, h_bkg_frac_true_error 
//                    }) {
//                     h->Write();
//     }

//     // TH1D              
//     for (auto h :{ 
//                 h_sig_fraction_DeltaRIntBin,
//                 h_sig_fraction_DeltaRIntBin_error,
//                 h_sig_true_fraction_DeltaRIntBin,
//                 h_sig_true_fraction_DeltaRIntBin_error,
//                 h_bkg_fraction_DeltaRIntBin,
//                 h_bkg_fraction_DeltaRIntBin_error,
//                 h_bkg_true_fraction_DeltaRIntBin,
//                 h_bkg_true_fraction_DeltaRIntBin_error
//                 }) {
//                 h->Write();
//     }


//     // //See if some bins did not converge
//     for (auto p : non_converge_bins) {
//         std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
//     }
        
//     fout->Close();

// }

void do_template_fit_Dijet(TString &dataset, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name){
    /*This fit use sign and bkg templates from the input files.*/
    //-- Samples: Dijet + bjet 
    // -- For now: Pseudodata is used as direct sum of sig + bkg 
    // -- Template fit pdfs : includes (0B, 1B) as bkg (added effectively), 2B as sig  

    bool isSampled = false; 
        // if the templates and pseudodata are seperated or not.
        // if False: pseudodata = direct sum of sig and bkg, if True: it has different template based on random sampling.
    

    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");


    // ---------------------------------------------not needed anymore ---------------------------------------------------
    // -- Temporare definiton: to be changed later to be common to the Dijet template 
    // TFile *file_bjet = new TFile("template_for_fit_histos_3D_makevtx_btag.root", "read");
        // if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());return;}
    // ---------------------------------------------------------------------------------------------------------

/*
    // -- using older Afnan templates
    // -- Input Histogram name: MC closure test  

    TString namehMC = "h3D_data";
    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";// h3D_bb
    TString nameh1B = "h3D_b"; // h3D_b
    // -- Add template of 0B case: ht_3d_trks_0Bbkgfor2Bsignal or h3D_nob
     TString nameh0B = "h3D_nob";// h3D_nob
     TString namehmore2B = "h3D_more2b"; // h3D_more2b
*/

// -- Using Zoe file names: without sampling! 
    TString namehMC = "h3D_data";
    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_2b";
    TString nameh1B = "h3D_1b";
    TString nameh0B = "h3D_0b";

    
    //-- Dijet sample: Templates: Open file and get histograms from Dijet: templates 
    TFile *file_dijet = new TFile(folder + templates, "read");
        if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());return;}
    TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b");
    TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
        if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
        if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
    // -- 0B and more B contribution for drawing only
    
    // TH3D *h3D_more2b = (TH3D*)file_dijet->Get(namehmore2B)->Clone(namehmore2B);
       // if(!h3D_more2b){Error("Get:", "histogram does not exist '%s' ",h3D_more2b->GetName()); return;}

    TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone(nameh0B);
       if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}


    //-- Bjet sample: Templates: Open file and get histograms from bjet: templates 
    TFile *file_bjet = new TFile(folder + templates_bjet, "read");
    TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
    TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
        if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
        if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}


    //Open dataset:
    TFile *file_data = new TFile(folder + dataset, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
        cout << "file dat name " << file_data->GetName() << endl;

        TH3D *h3D_data;
        TH3D *h3D_data_bjet;
        if(isSampled){
            h3D_data = (TH3D*)file_data->Get(namehMC)->Clone(namehMC);
            h3D_data_bjet = (TH3D*)file_bjet->Get(namehMC)->Clone("h3D_data_bjet");
        }
        else {
            // Trivial closure test

            // do sum: 1B and 2B!
            // from Dijet sample 
            h3D_data = (TH3D*) h3D_b->Clone(namehMC);
            h3D_data->Add(h3D_bb);
            h3D_data->Add(h3D_nob);

            // then add bjet sample also
            h3D_data_bjet = (TH3D*) h3D_b_bjet->Clone("h3D_data_bjet");
            h3D_data_bjet->Add(h3D_bb_bjet);
        }
        if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}
        if(!h3D_data_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_data_bjet->GetName()); return;}


    //Get number of bin entries
    Int_t bins_pt = h3D_data->GetNbinsZ();
    Int_t bins_dr = h3D_b->GetNbinsY(); // Re-defined bins(exists in the header file)
    Int_t mb_bins = h3D_data->GetNbinsX();
    

    // -- Draw THstach for mb projection 
    // Get projections X 
    // hists names as for slices in bins (but outside {} so it is ok)
        TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX("h_data_mb");
        TH1D *h_data_bjet_mb = (TH1D *) h3D_data_bjet->ProjectionX("h_data_bjet_mb");

        TH1D *h_b = (TH1D *) h3D_b->ProjectionX("h_b");
        TH1D *h_b_bjet = (TH1D *) h3D_b_bjet->ProjectionX("h_b_bjet");
        
        TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX("h_bb");
        TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet->ProjectionX("h_bb_bjet");
        
        // TH1D *h_more2b = (TH1D *) h3D_more2b->ProjectionX("h_more2b");
        TH1D *h_nob = (TH1D *) h3D_nob->ProjectionX("h_nob");
        
        // -- Set titles for the hisograms  to be used in the default legends 
        h_data_mb->SetTitle("0B+1B+2B (dijet)");
        h_data_bjet_mb->SetTitle("1B+2B (bjet)");
        h_b->SetTitle("1B (dijet)");
        h_b_bjet->SetTitle("1B (bjet)");
        h_bb->SetTitle("2B (dijet)");
        h_bb_bjet->SetTitle("2B (bjet)");
        h_nob->SetTitle("0B");


        // Get histogram of the residual difference ebtween FullMC - (BB signal + B bkg):
        TH1D *hRes = (TH1D*) h_data_mb->Clone("hRes");
            hRes->SetDirectory(0);
            hRes->SetTitle("Residual difference between MC and (sig + bkg)");
            hRes->Add(h_data_bjet_mb, +1);
            hRes->Add(h_b, -1);
            hRes->Add(h_b_bjet, -1);
            hRes->Add(h_bb, -1);
            hRes->Add(h_bb_bjet, -1);
            hRes->Add(h_nob, -1);
            // hRes->Add(h_more2b, -1);
            hRes->SetLineColor(kBlack);
            hRes->SetLineWidth(2);
            fout->cd();
            hRes->Write();

        // The trivial sum of sig+bkg t compare to full MC 
        TH1D *hsum = (TH1D*) h_b->Clone("hsum");
        hsum->SetTitle("BB signal + 1B bkg + 0B bkg");
        hsum->Add(h_b_bjet);
        hsum->Add(h_bb);
        hsum->Add(h_bb_bjet);
        hsum->Add(h_nob);
        // hsum->Add(h_more2b);
        hsum->SetDirectory(0);
        hsum->Write();

       // -- Deattach hists from root file 
        h_data_mb->SetDirectory(0);
        h_b->SetDirectory(0);
        h_bb->SetDirectory(0);
        h_data_bjet_mb->SetDirectory(0);
        h_b_bjet->SetDirectory(0);
        h_bb_bjet->SetDirectory(0);
        // h_more2b->SetDirectory(0);
        h_nob->SetDirectory(0);
       
        // Set styles of hist to be drawn:
        h_data_mb->SetLineColor(kBlack);
        h_data_mb->SetMarkerColor(kBlack);
        h_data_mb->SetMarkerStyle(4);

        h_b->SetFillColor(kRed-10); h_b->SetLineColor(kRed-10);        
        h_bb->SetFillColor(kBlue-10); h_bb->SetLineColor(kBlue-10);

        h_data_bjet_mb->SetLineColor(kGray+1);
        h_data_bjet_mb->SetMarkerColor(kGray+1);
        h_data_bjet_mb->SetMarkerStyle(8);

        h_b_bjet->SetFillColor(kGreen-2); h_b_bjet->SetLineColor(kGreen-2);        
        h_bb_bjet->SetFillColor(kOrange-6); h_bb_bjet->SetLineColor(kOrange-6); // brown 

        // h_more2b->SetFillColor(kOrange); h_more2b->SetLineColor(kOrange);

        h_nob->SetFillColor(kRed); h_nob->SetLineColor(kRed);

        hsum->SetLineWidth(3);
        hsum->SetLineColor(kMagenta+2);

    THStack hstack_All("hstack_All","Mass stacked histogram");
    hstack_All.SetTitle(";m_{2B} [GeV];");
    // hstack_All.Add(h_more2b);
    hstack_All.Add(h_nob);
    hstack_All.Add(h_b);
    hstack_All.Add(h_bb);
    hstack_All.Add(h_b_bjet);
    hstack_All.Add(h_bb_bjet);
    
    auto canva = new TCanvas("All_templates",Form("Templates(mb) and Full MC"), 800, 800 );
        canva->cd();
        canva->SetLogy();
        hstack_All.Draw("hist E"); // hist
        TH1D *h_data_mb_all = (TH1D*) h_data_mb->Clone("h_data_mb_all");
            h_data_mb_all->SetTitle("Dijet + bjet sample");
            h_data_mb_all->SetMarkerStyle(8);
            h_data_mb_all->Add(h_data_bjet_mb, 1);
            h_data_mb_all->Draw("PE same");  // full MC as points
            h_data_mb->Draw("HIST E same");
            h_data_bjet_mb->Draw("HIST E same");

            // For names to appear in the legend :
               h_data_mb->SetTitle("0B+1B+2B (dijet)");
                h_data_bjet_mb->SetTitle("1B+2B (bjet)");
                h_b->SetTitle("1B (dijet)");
                h_b_bjet->SetTitle("1B (bjet)");
                h_bb->SetTitle("2B (dijet)");
                h_bb_bjet->SetTitle("2B (bjet)");
                h_nob->SetTitle("0B");

        // hsum->Draw("histsame");
        gPad->Modified();   
        gPad->Update();
        canva->Modified();
        canva->Update();
        fout->cd();
        canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
        canva->Write();
        canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));

    // -- Output hist declarations: empty hists
    // -- For the fit result: per (deltaR, jtpt) bin
        //Create the histograms storing the fit parameters in TH2D (the dr are stored in the X, the jtpt in the Y)
        TH2D *h_sig_fraction = (TH2D*) h3D_data->Project3D("zy")->Clone("h_sig_fraction");// to get the axes from hist
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
          


    //Vector to test the convergence
    std::vector<std::pair<int, int>> non_converge_bins;


    // Fitting - loop over dr and jtpt entries
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){
        // for(Int_t ibin_dr = 0; ibin_dr <= 0; ibin_dr++){
        // Bin = 0: integaretd deltaR range, other bins are the DeltaR bins 
            
            // The slice range for projection is different for integarted Vs. one bin slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;

            // For integarted bin only: bin = 0
            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = bins_dr;}

            // Make slices for data
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_data_bjet_mb = (TH1D *) h3D_data_bjet->ProjectionX(Form("h_data_mb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);

            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            // TH1D *h_more2b = (TH1D *) h3D_more2b ->ProjectionX(Form("h3D_more2b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            //  Make slices for Bjet 
            TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            

            // -- Deattach hists from root file 
            h_data_mb->SetDirectory(0);
            h_b->SetDirectory(0);
            h_bb->SetDirectory(0);
            h_b_bjet->SetDirectory(0);
            h_bb_bjet->SetDirectory(0);
            h_nob->SetDirectory(0);
            // h_more2b->SetDirectory(0);

  TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);

            // Calculate true fractions to be used as initial values for the fit (the true fractions are the dijet qcd ones)            
            double int0 = h_bb->Integral(1, mb_bins, "width");// sig
            double int1 = h_b->Integral(1, mb_bins, "width");// bkg1
            // double int2 = h_more2b->Integral(1, mb_bins, "width"); // bkg2
            double int3 = h_nob->Integral(1, mb_bins, "width"); // bkg 3

            // From bjet sample 
            double int0_bjet = h_bb_bjet->Integral(1, mb_bins, "width");// sig
            double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");// bkg1

            // For MC test case only
            double integral_inputdata1 = h_data_mb->Integral(1, mb_bins,"width");
            double integral_inputdata2 = h_data_bjet_mb->Integral(1, mb_bins,"width");
                cout << "Before Samples additions:\n Dijet data integral = "<< integral_inputdata1 << endl;
                cout << "Bjet data integral = "<< integral_inputdata2 << endl;

            // Add Dijet and bjet sample for the pseudodata 
            h_data_mb->Add(h_data_bjet_mb);

            // For data: After the sum of dijet + bjet 
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
                
                cout << "Total input data integral = "<< integral_inputdata << endl;
                std::cout << "int0 h_bb=" << int0 << std::endl;
                std::cout << "int1 h_b =" << int1 << std::endl;
                cout << "2B/1B in Dijet sample = " << int0/int1 << endl;
                cout << "Dijet: 2B/(1B + 2B) = " << int0/(int0 + int1) << endl;

                std::cout << "int0 h_bb_bjet=" << int0_bjet << std::endl;
                std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                cout <<  "2B/1B in  Bjet sample = " << int0_bjet/int1_bjet << endl;
                cout << "2B/(2B + 1B) in Bjet sample = " << int0/(int0 + int1) << endl;
                // std::cout << "int2 h_more2b =" << int2 << std::endl;
                std::cout << "int3 h_nob =" << int3 << std::endl;


            // --- For now: the TRUE signal (2B) fraction is used from the Dijet sample.
                // the fraction of 2B/1B in Bjet and Dijet sample is exactly the same
            // double sig_fraction_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int2 + int3 )); 
            // double bkg_fraction_b_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int2 + int3 ));
            
            // No contiubtion from hmore2B 
            double sig_fraction_true = (int0  + int1 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int3 )); 
            double bkg_fraction_b_true = (int0  + int1 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int3 ));


                std::cout << "Dijet : 2B sig_fraction_true =" << sig_fraction_true << std::endl;
                std::cout << "Dijet: 1B bkg_fraction_b_true =" << bkg_fraction_b_true << std::endl;
            // -- Compute the errors of true fractions for later fill
            double True_bkg_b_err = 0.0;
            double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");
            double True_sig_err = 0.0;
            double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");


            // -- Get Normalized templates as input  PDF to the fit(shape matters only): where PDF integral is 1 (Probability density function)
            //signal: 2B from Dijet + Bjet sample 
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
            h_sumsig->Add(h_bb_bjet);
            Int_t h_sig_bins = h_sumsig->GetNbinsX();
            TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized


            //Add background: here only from 1B 
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
            h_sumbkg->Add(h_b_bjet);
            Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
            
            TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized
            
            // Total sum bkg = 0B + 1B  
            TH1D* h_sumbkg_0b_1b = (TH1D*)  h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
            h_sumbkg_0b_1b->Add(h_nob);


            // write signal hists after combine Dijet and bjet sample 
            fout->cd();
            h_sumsig->Write();
            h_sumbkg->Write();
            h_sumbkg_0b_1b->Write();
           

            // -- Draw Stack for combined: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before ormalization 
              THStack hstack_templatesfotfit(Form("hstack_templatesfotfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesfotfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesfotfit.Add(h_nob);
                    hstack_templatesfotfit.Add(h_sumbkg);
                    hstack_templatesfotfit.Add(h_sumsig);
                    h_sumsig->SetTitle("2B: Dijet + bjet");
                    h_sumbkg->SetTitle("1B: Dijet + bjet");
                    h_nob->SetTitle("0B");

                    // -- Set fill colors 
                    h_nob->SetFillColor(kRed);
                    h_sumsig->SetFillColor(kBlue -10);
                    h_sumbkg->SetFillColor(kRed-10);
                     h_sumbkg->SetFillStyle(3244);
                    h_data_mb->SetMarkerColor(kBlack); h_data_mb->SetMarkerStyle(20);

              auto canva = new TCanvas(Form("All_templates_%d_%d", ibin_dr, ibin_pt),Form("Templates(mb) and Full MC"), 800, 800 );
                    canva->cd();
                    canva->SetLogy();
                    hstack_templatesfotfit.Draw("hist E"); 
                    h_data_mb->SetTitle("All");
                    h_data_mb->Draw("PE same");
                    gPad->Modified();   
                    gPad->Update();
                    canva->Modified();
                    canva->Update();
                    fout->cd();
                    canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
                    canva->Write();
                    canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));




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
                cout << "True 1B integral without width option =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;


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
            
            // -- test 
               // eff_bkg0B = 0;
               // eff_bkg1B = 1; 

            // Build effective signal hist: with new relaitve normalization, the integral should = 1
            h_bkg->Add(h_bkg, h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist
            
            double err_int;
            double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins,err_int ,"width");
            cout << "Effective Bkg normalized hist integral +/- uncertainity = "<< int_val << "+/-" << err_int << endl;

            cout << "effective 1B fraction and 0B = "<< eff_bkg1B << ", " << eff_bkg0B << endl;
            // -- Note: the effective bkg integral is Effective Bkg normalized hist integral = 0.997777
            // if (h_bkg->Integral(1,  h_bkg_bins, "width") != 1.) {cout << "Effective BKG PDF is not normalized!"<< endl; return;}

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
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true, 0.,1.);// 1-bkg_fraction_b_true

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, true); //fit the scaling when adding the single b with the more b distribution
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name
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
            double err_true_frac = (int1 * True_sig_err + int0 * True_bkg_b_err)/(int0+ int1)/(int0+ int1);// int0 for signal, int1 for bkg 

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
                h_sig_fit->SetTitle("After fit: Sig");
                h_sig_fit->SetMarkerColor(kRed+2);
                h_sig_fit->SetLineColor(kRed+2);
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg 
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("After fit: bkg (1B+ 0B)");
                h_bkg_fit->SetMarkerColor(kCyan+2);
                h_bkg_fit->SetLineColor(kCyan+2);

            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only 
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("After fit: 1B");
                h_bkg_fit_1b->SetMarkerColor(kOrange+2);
                h_bkg_fit_1b->SetLineColor(kOrange+2);
                h_bkg_fit_1b->SetLineWidth(2);

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only 
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("After fit: 0B");
                h_bkg_fit_nob->SetMarkerColor(kYellow+2);
                h_bkg_fit_nob->SetLineColor(kYellow+2);
                h_bkg_fit_nob->SetLineWidth(2);

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("After fit: sig + bkg");
                    h_total_fit->SetMarkerColor(kGreen+2);
                    h_total_fit->SetLineColor(kGreen+2);
                    h_total_fit->SetMarkerStyle(20);
                    h_total_fit->SetLineWidth(2);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();
            h_bkg_fit->Write();
            h_total_fit->Write();

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
            auto canva_afterfit = new TCanvas(sname_canvas_afterfit,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetMaximum(1.2 * h_data_mb->GetMaximum());

                h_data_mb->SetTitle("Pseudodata (combined)");
                // h_b->SetTitle("B template (dijet)");
                // h_bb->SetTitle("2B template (dijet)");
                // h_b_bjet->SetTitle("B template (bjet)");
                // h_bb_bjet->SetTitle("2B template (bjet)");
                h_sumbkg->SetTitle("1B (bkg) template (combined)");
                h_sumsig->SetTitle("Signal template: 2B (combined)");
                h_sumbkg_0b_1b->SetTitle("Bkg template (0B + 1B) (combined)");
                h_nob->SetTitle("0B template (Dijet only)");

                h_sumbkg->SetLineColor(kOrange);
                h_nob->SetLineColor(kMagenta);
                h_data_mb->SetLineWidth(2);


                // -- Make better visual
                // h_bb->SetLineWidth(3);
                // h_bb_bjet->SetLineWidth(3);
                h_sumsig->SetLineWidth(3);
                h_sig_fit->SetLineWidth(3);

                // h_b->SetLineStyle(9); // close dashes
                // h_b_bjet->SetLineStyle(9);
                h_sumbkg_0b_1b->SetLineStyle(7);
                h_sumbkg_0b_1b->SetLineWidth(3);
                h_bkg_fit->SetLineWidth(2);
                h_bkg_fit->SetLineStyle(10);  
                
                // Add after fit templates
                h_data_mb->Draw("PE");
                h_total_fit->Draw("P E SAME");

                // h_bb->Draw("HIST E SAME");
                // h_bb_bjet->Draw("HIST E SAME");
                h_sumsig->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                

                h_sumbkg->Draw("HIST E SAME"); // Combined sample for 1B bkg 
                h_nob->Draw("HIST E SAME");
                h_bkg_fit_nob->Draw("HIST E SAME");
                h_bkg_fit_1b->Draw("HIST E SAME");
                h_sumbkg_0b_1b->Draw("HIST E SAME");
                h_bkg_fit->Draw("HIST E SAME");

                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
                fout->cd();

                TString srangeDeltaR = Form(" %g < #DeltaR [bin %d] < %g", dr_binsVector[ibin_dr -1], ibin_dr ,dr_binsVector[ibin_dr]);
                if(!ibin_dr) srangeDeltaR = Form(" %g < #DeltaR [Full range] < %g", dr_binsVector[1], dr_binsVector[bins_dr]);
                canva_afterfit->BuildLegend(0.7, 0.9, 0.7, 0.9, srangeDeltaR);
                // Add another line for all bins 
                canva_afterfit->Write();


                canva_afterfit->Print(Form("%s/%s.png", sDirname.Data(), sname_canvas_afterfit.Data()));

                cout << "---------------------\n\n\n" << endl; 
            } // loop over deltaR bins 
    }

    // Save histograms
    // TH3D 
    // for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_data_bjet, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
    for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
        

        
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
        
    fout->Close();

}

//Performs 3D template fit of mB for every dr and jtpt bin
void do_template_fit(TString &dataset, TString &templates, TString pT_selection, TString folder, TString &fout_name){
    /*This fit use the data and sign and bkg templates from the input files. No direct sum is done!*/

    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    // -- Input Histogram name:  Trivial MC closure test 
    TString namehMC = "h3D_data";
    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";
    TString nameh1B = "h3D_b";

    //Open dataset
    TFile *file_data = new TFile(folder + dataset, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
        cout << "file dat name " << file_data->GetName() << endl;
    TH3D *h3D_data = (TH3D*)file_data->Get(namehMC)->Clone("h3D_data");
        if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}

    //Open file and get histograms from bjet: templates 
    TFile *file_bjet = new TFile(folder + templates, "read");
        if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());return;}
    TH3D *h3D_b = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b");
    TH3D *h3D_bb = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb");
        if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
        if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}

      TH3D *h3D_more2b = (TH3D*)file_bjet->Get("h3D_more2b")->Clone("h3D_more2b");
       if(!h3D_more2b){Error("Get:", "histogram does not exist '%s' ",h3D_more2b->GetName()); return;}

    //Get number of bin entries
    Int_t bins_pt = h3D_data->GetNbinsZ();
    Int_t bins_dr = h3D_b->GetNbinsY(); // Re-defined bins(exists in the header file)
    Int_t mb_bins = h3D_data->GetNbinsX();
    

    // -- Draw THstach for mb projection 
        // Get projections X 
    // hists names as for slices in bins (but outside {} so it is ok)
        TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX("h_data_mb");
        TH1D *h_b_bjet = (TH1D *) h3D_b->ProjectionX("h_b_bjet");
        TH1D *h_bb_bjet = (TH1D *) h3D_bb->ProjectionX("h_bb_bjet");
        TH1D *h_more2b = (TH1D *) h3D_more2b->ProjectionX("h_more2b");

        // Get histogram of the residual difference ebtween FullMC - (BB signal + B bkg)
        TH1D *hRes = (TH1D*) h_data_mb->Clone("hRes");
            hRes->SetDirectory(0);
            hRes->SetTitle("Residual difference between MC and (sig + bkg)");
            hRes->Add(h_b_bjet, -1);
            hRes->Add(h_bb_bjet, -1);
            hRes->SetLineColor(kBlack);
            hRes->SetLineWidth(2);

            fout->cd();
            hRes->Write();

        // The trivial sum of sig+bkg t compare to full MC 
        TH1D *hsum = (TH1D*) h_b_bjet->Clone("hsum"); hsum->SetTitle("BB signal + 1B bkg");
        hsum->Add(h_bb_bjet);
        hsum->SetDirectory(0);

       // -- Deattach hists from root file 
        h_data_mb->SetDirectory(0);
        h_b_bjet->SetDirectory(0);
        h_bb_bjet->SetDirectory(0);
        h_more2b->SetDirectory(0);

        // Set styles of hist to be drawn:
         h_data_mb->SetLineColor(kBlack);
                // h_data_mb->SetLineWidth(3);
                h_data_mb->SetMarkerColor(kBlack);
                // h_data_mb->SetMarkerStyle(20);
                h_data_mb->SetMarkerStyle(4);

        h_b_bjet->SetFillColor(kRed-10); h_b_bjet->SetLineColor(kRed-10);
        
        h_bb_bjet->SetFillColor(kBlue-10); h_bb_bjet->SetLineColor(kBlue-10);

        h_more2b->SetFillColor(kRed); h_more2b->SetLineColor(kRed);

        hsum->SetLineWidth(3);
        hsum->SetLineColor(kMagenta+2);

    THStack hstack_All("hstack_All","Mass stacked histogram");
    hstack_All.Add(h_b_bjet);
    hstack_All.Add(h_bb_bjet);
    hstack_All.Add(h_more2b);
    
    auto canva = new TCanvas("All_templates",Form("Templates(mb) and Full MC (except zero B case)"), 800, 800 );
        canva->cd();
        hstack_All.Draw("histE"); // hist
        h_data_mb->Draw("PE same");  // full MC as points
        // hsum->Draw("histsame");

        gPad->Modified();   
        gPad->Update();
        canva->Modified();
        canva->Update();
        fout->cd();
        canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
        canva->Write();





    // -- Output hist
    // -- For the fit result: per (deltaR, jtpt) bin
        //Create the histograms storing the fit parameters in TH2D (the dr are stored in the X, the jtpt in the Y)
        TH2D *h_sig_fraction = (TH2D*) h3D_data->Project3D("zy")->Clone("h_sig_fraction");// to get the axes from hist
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
          


    //Vector to test the convergence
    std::vector<std::pair<int, int>> non_converge_bins;


    // Fitting - loop over dr and jtpt entries
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        // for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){
        for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){
        // Bin = 0: integaretd deltaR range, other bins are the DeltaR bins 
            
            // The slice range for projection is different for integarted Vs. one bin slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;

            // For integarted bin only: bin = 0
            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = bins_dr;}

            // Make slices for data
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);

            // Make slices for bjet
            TH1D *h_b_bjet = (TH1D *) h3D_b->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_bb_bjet = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            

            // -- Deattach hists from root file 
            h_data_mb->SetDirectory(0);
            h_b_bjet->SetDirectory(0);
            h_bb_bjet->SetDirectory(0);
            // Save mb projection in each deltaR, pt  for data, sig, bkg
            TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);
            auto canva = new TCanvas(sname_canvas,Form("Templates(mb) %s", sname_canvas.Data()), 800, 800 );
                canva->cd();
                h_data_mb->SetStats(0);
                h_b_bjet->SetStats(0);
                h_bb_bjet->SetStats(0);

                h_data_mb->SetLineColor(kBlack);
                h_data_mb->SetLineWidth(2);
                h_data_mb->SetMarkerColor(kBlack);
                h_data_mb->SetMarkerStyle(20);
                h_data_mb->Draw("PE ");

                h_b_bjet->SetMarkerColor(kMagenta+2);
                h_b_bjet->SetLineColor(kMagenta+2);
                // h_b_bjet->SetLineStyle(10);
                h_b_bjet->SetLineWidth(2);
                h_b_bjet->SetMarkerStyle(4);
                h_b_bjet->Draw("HIST E SAME");

                h_bb_bjet->SetMarkerStyle(25);
                h_bb_bjet->SetMarkerColor(kRed+2);
                h_bb_bjet->SetLineColor(kRed+2);
                h_bb_bjet->Draw("HIST E SAME");
                gPad->Modified();   
                gPad->Update();
                canva->Modified();
                canva->Update();
                fout->cd();
                canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
                canva->Write();


            // Calculate true fractions to be used as initial values for the fit (the true fractions are the dijet qcd ones)
            double int1 = h_b_bjet->Integral(1, mb_bins, "width");// bkg
            double int0 = h_bb_bjet->Integral(1, mb_bins, "width");// sig
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
        
            std::cout << "int0=" << int0 << std::endl;
            std::cout << "int1=" << int1 << std::endl;

            double sig_fraction_true = (int0 == 0 && int1 == 0) ? 0 : (int0 / (int0 + int1 )); 
            double bkg_fraction_b_true = (int0 == 0 && int1 == 0) ? 0 : (int1 / (int0 + int1));

            std::cout << "sig_fraction_true=" << sig_fraction_true << std::endl;
            std::cout << "bkg_fraction_b_true=" << bkg_fraction_b_true << std::endl;

            // Compute the errors of true fractions for later fill
            double True_bkg_err = 0.0;
            double True_bkg_integral = h_b_bjet->IntegralAndError(1, mb_bins, True_bkg_err, "width");
            double True_sig_err = 0.0;
            double True_sig_integral = h_bb_bjet->IntegralAndError(1, mb_bins, True_sig_err, "width");


            // -- Get Normalized templates as input  PDF to the fit(shape matters only): where PDF integral is 1 (Probability density function)
            //Add signal from bb jets 
            TH1D *h_sig = (TH1D*)h_bb_bjet->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt));
            Int_t h_sig_bins = h_sig->GetNbinsX();
            //Add background from b jets
            TH1D *h_bkg = (TH1D*)h_b_bjet->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt));
            Int_t h_bkg_bins = h_bkg->GetNbinsX();
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

            /// Hay que hacer esto ???
            //Add signal and rest with the correct scaling (from Lida)
            // add together light+c to sig with nominal MC ratio
            // a sig + b bb + c light = n, a+b+c=1, a+c=1-b
            // c'=c/(a+c)=(1-a-b)/(1-b)
            // a'=a/(a+c)=a/(1-b)=1-c'

            // -- Afnan: what are these param?
            // -- Not used in my case: 
             // double cl_frac = (1 - sig_fraction_true - bkg_fraction_b_true) / (1 - bkg_fraction_b_true); // nominal
           
            // cl_frac = 2*cl_frac; // systematic
            /////// For now: since i have pure signal, use charm-light contribution = 0
            double cl_frac = 0.; // systematic
            double sig_frac = 1 - cl_frac;
            ////// h_sig->Add(h_sig, h_rest, sig_frac, 1-sig_frac);
            

            //Fitting

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
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val",1-bkg_fraction_b_true, 0.,1.);//

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, true); //fit the scaling when adding the single b with the more b distribution
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name
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
            //std::cout << "a = " << a << " da = " << da << std::endl;
            /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
            Double_t p0, p1, p2, errP0, errP1, errP2;
            p0 = a*sig_frac;
            p1 = 1-a;
            p2 = a*(1-sig_frac);

            errP0 = da*sig_frac;
            errP1 = da;
            errP2 = da*(1-sig_frac);
            // std::cout << "errP0 = " << errP0 << std::endl;
            //std::cout << "a'/a = " << p0/sig_fraction_true << std::endl;
            //std::cout << "c'/c = " << p2/(1-sig_fraction_true-bkg_bb_fraction_true) << std::endl;
            //std::cout << "c=" << (1-sig_fraction_true-bkg_bb_fraction_true) << ", c'=" << p2 << std::endl;
            //std::cout << "a=" << sig_fraction_true << ", a'=" << p0 << std::endl;
            //std::cout << "b=" << bkg_bb_fraction_true << ", b'=" << p1 << std::endl;

        
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
            double err_true_frac = (int1 * True_sig_err + int0 * True_bkg_err)/(int0+ int1)/(int0+ int1);// int0 for signal, int1 for bkg 

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
                h_sig_fit->SetTitle("After fit: Sig");
                h_sig_fit->SetMarkerColor(kBlue+2);
                h_sig_fit->SetLineColor(kBlue+2);

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); 
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("After fit: bkg");
                h_bkg_fit->SetMarkerColor(kCyan+2);
                h_bkg_fit->SetLineColor(kCyan+2);

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("After fit: sig + bkg");
                    h_total_fit->SetMarkerColor(kGreen+2);
                    h_total_fit->SetLineColor(kGreen+2);
                    h_total_fit->SetMarkerStyle(20);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();
            h_bkg_fit->Write();
            h_total_fit->Write();

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
            auto canva_afterfit = new TCanvas(sname_canvas_afterfit,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();

                h_data_mb->SetTitle("Pseudodata");
                h_b_bjet->SetTitle("bkg template");
                h_bb_bjet->SetTitle("sig template");

                // -- Make better visual
                h_bb_bjet->SetLineWidth(3);
                h_sig_fit->SetLineWidth(3);
                h_b_bjet->SetLineStyle(9); // close dashes 
                h_bkg_fit->SetLineStyle(7); // dots 
                
                // Add after fit templates
                h_data_mb->Draw("PE");
                h_total_fit->Draw("P E SAME");

                h_bb_bjet->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                
                h_b_bjet->Draw("HIST E SAME");
                h_bkg_fit->Draw("HIST E SAME");

                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
                fout->cd();

                TString srangeDeltaR = Form(" %g < #DeltaR [bin %d] < %g", dr_binsVector[ibin_dr -1], ibin_dr ,dr_binsVector[ibin_dr]);
                if(!ibin_dr) srangeDeltaR = Form(" %g < #DeltaR [Full range] < %g", dr_binsVector[1], dr_binsVector[bins_dr]);
                canva_afterfit->BuildLegend(0.7, 0.9, 0.7, 0.9, srangeDeltaR);
                // Add another line for all bins 
                canva_afterfit->Write();

                TString sDirname = "TemplateFitOutputPlots";
                gSystem->mkdir(sDirname, kTRUE); // make directory if does not exist
                canva_afterfit->Print(Form("%s/%s.png", sDirname.Data(), sname_canvas_afterfit.Data()));

                cout << "---------------------\n\n\n" << endl; 
            } // loop over deltaR bins 
    }

    // Save histograms
    // TH3D 
    for (auto h : {h3D_data, h3D_bb, h3D_b}) {h->Write();}
        
    // TH2D     
    for (auto h : {h_sig_fraction, h_sig_fraction_error,
                   h_bkg_fraction, h_bkg_fraction_error,
                   h_sig_frac_true, h_sig_frac_true_error,
                   h_bkg_frac_true, h_bkg_frac_true_error 
                   }) {
                    h->Write();
    }

    // TH1D              
    for (auto h :{ h_sig_fraction_DeltaRIntBin,
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
        
    fout->Close();

}

void do_template_fit_trivialtest(TString &dataset, TString &templates, TString pT_selection, TString folder, TString &fout_name,  bool trivialMC = true){
/*This function assume full MC = SUM (templates from sign + Bkg). Regardless the MC template itself in the given input. */

    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");


    // -- Input Histogram name:  Trivial MC closure test 
    // TString namehMC = trivialMC ? "h3D_trivialMC": "h3D_data";
    TString namehMC = "h3D_data";
    
    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";
    TString nameh1B = "h3D_b";


    //Open dataset
    TFile *file_data = new TFile(folder + dataset, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
        cout << "file dat name " << file_data->GetName() << endl;

    // TH3D *h3D_data = (TH3D*)file_data->Get(namehMC)->Clone("h3D_data");
    //     if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}

    //Open file and get histograms from bjet: templates 
    TFile *file_bjet = new TFile(folder + templates, "read");
        if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());return;}
    TH3D *h3D_b = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b");
        h3D_b->SetDirectory(0);
    TH3D *h3D_bb = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb");
            h3D_bb->SetDirectory(0);
        if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
        if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}

    // Build the trivial  MC = signal + bkg 
    TH3D *h3D_data = (TH3D*) h3D_bb->Clone(namehMC); h3D_data->Add(h3D_b); h3D_data->SetTitle("Trivial sum of Sig & bkg");
        h3D_data->SetDirectory(0);

    // -- Draw THstach for mb projection 
        // Get projections X 
        TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX("h_data_mb");
        TH1D *h_b_bjet = (TH1D *) h3D_b->ProjectionX("h_b_bjet");
        TH1D *h_bb_bjet = (TH1D *) h3D_bb->ProjectionX("h_bb_bjet");
        // Set styles of hist to be drawn:
         h_data_mb->SetLineColor(kBlack);
                h_data_mb->SetLineWidth(2);
                h_data_mb->SetMarkerColor(kBlack);
                h_data_mb->SetMarkerStyle(20);
        h_b_bjet->SetFillColor(kRed-9);
        h_bb_bjet->SetFillColor(kBlue-2); 

    THStack hstack_All("hstack_All","Mass stacked histogram");
    hstack_All.Add(h_b_bjet);
    hstack_All.Add(h_bb_bjet);
    auto canva = new TCanvas("All_templates",Form("Templates(mb) and MC"), 800, 800 );
        canva->cd();
        hstack_All.Draw("hist");
        h_data_mb->Draw("PE same");  // full MC as points
        gPad->Modified();   
        gPad->Update();
        canva->Modified();
        canva->Update();
        fout->cd();
        canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
        canva->Write();




    //Get number of bin entries
    Int_t bins_pt = h3D_data->GetNbinsZ();
    Int_t bins_dr = h3D_b->GetNbinsY();
    Int_t mb_bins = h3D_data->GetNbinsX();
    
    // -- Output hist
    // -- For the fit result: per (deltaR, jtpt) bin
        //Create the histograms storing the fit parameters in TH2D (the dr are stored in the X, the jtpt in the Y)
        TH2D *h_sig_fraction = (TH2D*) h3D_data->Project3D("zy")->Clone("h_sig_fraction");// to get the axes from hist
        h_sig_fraction->SetTitle("x=#DeltaR, y=jtpt, fraction / error");
        h_sig_fraction->Reset();

        TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
        
        TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
        TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error");
       
        //Store the true parameters and errors
        TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
        TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");

        TH2D *h_bkg_true = (TH2D *) h_sig_fraction->Clone("h_bkg_true");
        TH2D *h_bkg_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_true_error");


    //Vector to test the convergence
    std::vector<std::pair<int, int>> non_converge_bins;


    // Fitting - loop over dr and jtpt entries
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){

            // Make slices for data
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            h_data_mb->SetDirectory(0);
            // Make slices for bjet
            TH1D *h_b_bjet = (TH1D *) h3D_b->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            h_b_bjet->SetDirectory(0);
            TH1D *h_bb_bjet = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            h_bb_bjet->SetDirectory(0);
            
            // Save mb projection in each deltaR, pt  for data, sig, bkg
            TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);
            auto canva = new TCanvas(sname_canvas,Form("Templates(mb) %s", sname_canvas.Data()), 800, 800 );
                canva->cd();
                h_data_mb->SetStats(0);
                h_b_bjet->SetStats(0);
                h_bb_bjet->SetStats(0);

                h_data_mb->SetLineColor(kBlack);
                h_data_mb->SetLineWidth(2);
                h_data_mb->SetMarkerColor(kBlack);
                h_data_mb->SetMarkerStyle(20);
                h_data_mb->Draw("PE ");
                h_b_bjet->SetMarkerColor(kBlue+2);
                h_b_bjet->SetLineColor(kBlue+2);
                // h_b_bjet->SetLineStyle(10);
                h_b_bjet->SetLineWidth(2);
                h_b_bjet->SetMarkerStyle(4);
                h_b_bjet->Draw("HIST E SAME");
                h_bb_bjet->SetMarkerStyle(25);
                h_bb_bjet->SetMarkerColor(kRed+2);
                h_bb_bjet->SetLineColor(kRed+2);
                h_bb_bjet->Draw("HIST E SAME");
                gPad->Modified();   
                gPad->Update();
                canva->Modified();
                canva->Update();
                fout->cd();
                canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
                canva->Write();


            // Calculate true fractions to be used as initial values for the fit (the true fractions are the dijet qcd ones)
            double int1 = h_b_bjet->Integral();// bkg
            double int0 = h_bb_bjet->Integral();// sig
        
            std::cout << "int0=" << int0 << std::endl;
            std::cout << "int1=" << int1 << std::endl;

            // -- if both sig and bkg true frac = 0, set intial fit value to zero 
            double sig_fraction_true = (int0 == 0 && int1 == 0) ? 0 : (int0 / (int0 + int1 )); 
            double bkg_fraction_b_true = (int0 == 0 && int1 == 0) ? 0 : (int1 / (int0 + int1));




            std::cout << "sig_fraction_true=" << sig_fraction_true << std::endl;
            std::cout << "bkg_fraction_b_true=" << bkg_fraction_b_true << std::endl;

            // -- Get Normalized templates as input  PDF to the fit(shape matters only): where PDF integral is 1 (Probability density function)
            //Add signal from bb jets 
            TH1D *h_sig = (TH1D*)h_bb_bjet->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt));
            Int_t h_sig_bins = h_sig->GetNbinsX();
            h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));

            //Add background from b jets
            TH1D *h_bkg = (TH1D*)h_b_bjet->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt));
            Int_t h_bkg_bins = h_bkg->GetNbinsX();
            h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));

        
            /// Hay que hacer esto ???
            //Add signal and rest with the correct scaling (from Lida)
            // add together light+c to sig with nominal MC ratio
            // a sig + b bb + c light = n, a+b+c=1, a+c=1-b
            // c'=c/(a+c)=(1-a-b)/(1-b)
            // a'=a/(a+c)=a/(1-b)=1-c'

            // -- Afnan: what are these param?
            // -- Not used in my case: 
             // double cl_frac = (1 - sig_fraction_true - bkg_fraction_b_true) / (1 - bkg_fraction_b_true); // nominal
           
            // cl_frac = 2*cl_frac; // systematic
            /////// For now: since i have pure signal, use charm-light contribution = 0
            double cl_frac = 0.; // systematic
            double sig_frac = 1 - cl_frac;
            ////// h_sig->Add(h_sig, h_rest, sig_frac, 1-sig_frac);
            

            //Fitting

            // Create the observable
            Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
            Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
            RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
            mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
     

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
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val",1-bkg_fraction_b_true,0.,1.);//

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, true); //fit the scaling when adding the single b with the more b distribution
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name
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
            //std::cout << "a = " << a << " da = " << da << std::endl;
            /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
            Double_t p0, p1, p2, errP0, errP1, errP2;
            p0 = a*sig_frac;
            p1 = 1-a;
            p2 = a*(1-sig_frac);

            errP0 = da*sig_frac;
            errP1 = da;
            errP2 = da*(1-sig_frac);
            // std::cout << "errP0 = " << errP0 << std::endl;
            //std::cout << "a'/a = " << p0/sig_fraction_true << std::endl;
            //std::cout << "c'/c = " << p2/(1-sig_fraction_true-bkg_bb_fraction_true) << std::endl;
            //std::cout << "c=" << (1-sig_fraction_true-bkg_bb_fraction_true) << ", c'=" << p2 << std::endl;
            //std::cout << "a=" << sig_fraction_true << ", a'=" << p0 << std::endl;
            //std::cout << "b=" << bkg_bb_fraction_true << ", b'=" << p1 << std::endl;

        
            //save the fit, rescaling for the signal only (not rest background)
            h_sig_fraction->SetBinContent(ibin_dr, ibin_pt, p0);
            h_sig_fraction->SetBinError(ibin_dr, ibin_pt, errP0);
            h_sig_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP0);
            h_bkg_fraction->SetBinContent(ibin_dr, ibin_pt, p1);
            h_bkg_fraction->SetBinError(ibin_dr, ibin_pt, errP1);
            h_bkg_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP1);
        

            //save the true fraction
            h_sig_frac_true->SetBinContent(ibin_dr, ibin_pt, sig_fraction_true);
            h_bkg_true->SetBinContent(ibin_dr, ibin_pt, bkg_fraction_b_true);

            }
    }

    // Save histograms


    for (auto h : {h3D_data, h3D_bb, h3D_b}) {h->Write();}
        

    for (auto h : {h_sig_fraction, h_sig_fraction_error,
                   h_bkg_fraction, h_bkg_fraction_error,
                   h_sig_frac_true, h_sig_frac_true_error,
                   h_bkg_true
                   }) 
    {
                    h->Write();
    }

    //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
        
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



/*//Draw a preliminary version of the EEC(dr) after template fit for checking
void draw_eec(TString fout_name, TString &dataset, TString &folder, TString &pT_selection, TString &pT_selection_label, bool &norm, bool &all, bool &eff_corr, Int_t &pt_bin, bool &data){
    //Define the canvas
    TCanvas *c = new TCanvas("c", " ",170,800,800,504);
    c->SetLogx();
    //c->SetLogy();
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetTitle("EEC corrected for single-b signal fraction, all mB bins");


    // Get signal and backgrounds fractions histograms
    TFile *file = new TFile(fout_name, "read");
    TH2D *h_sigfrac = (TH2D*)file->Get("h_sig_fraction");
    TH2D *h_bkg_bb_fraction = (TH2D*)file->Get("h_bkg_bb_fraction");
    
    // Get data
    TH3D *h3D = (TH3D*)file->Get("h_data");

    
    //Get full MC samples
    TH3D *h3D_mc = (TH3D*)file->Get("hmc_dijet");
    TH2D *hmc_2D = (TH2D*)h3D_mc->Project3D("zy")->Clone("hmc_2D");
    TH1D *hmc = (TH1D*)hmc_2D->ProjectionX("hmc", pt_bin, pt_bin);

    //Get 1 b samples
    TH3D *h3D_1b = (TH3D*)file->Get("hb_dijet");
    TH2D *h1b_2D = (TH2D*)h3D_1b->Project3D("zy")->Clone("h1b_2D");
    TH1D *h1b = (TH1D*)h1b_2D->ProjectionX("h1b", pt_bin, pt_bin);

    //Get more b samples
    TH3D *h3D_moreb = (TH3D*)file->Get("hmoreb_dijet");
    TH2D *hmoreb_2D = (TH2D*)h3D_moreb->Project3D("zy")->Clone("hmoreb_2D");
    TH1D *hmoreb = (TH1D*)hmoreb_2D->ProjectionX("hmoreb", pt_bin, pt_bin);

    //Get other background samples
    TH3D *h3D_other = (TH3D*)file->Get("hother_dijet");
    TH2D *hother_2D = (TH2D*)h3D_other->Project3D("zy")->Clone("hother_2D");
    TH1D *hother = (TH1D*)hother_2D->ProjectionX("hother", pt_bin, pt_bin);

    //Get the efficiency correction factor data
    TFile *file3D_1b_notag = new TFile(folder + "hist_3d_gen_aggr_n1_MC_dijet_notag_" + pT_selection + ".root", "read");
    //TFile *file3D_1b_notag = new TFile(folder + "hist_3d_gen_aggr_n1_b_dijet_notag_" + pT_selection + ".root", "read");
    TH3D *h3D_1b_notag = (TH3D*)file3D_1b_notag->Get("h3D_b1");
    TH2D *h1b_notag = (TH2D*)h3D_1b_notag->Project3D("zy")->Clone("h1b_notag");
    
    TFile *file3D_1b_tag = new TFile(folder + "hist_3d_gen_aggr_n1_MC_dijet_" + pT_selection + ".root", "read");
    //TFile *file3D_1b_tag = new TFile(folder + "hist_3d_gen_aggr_n1_b_dijet_" + pT_selection + ".root", "read");
    TH3D *h3D_1b_tag = (TH3D*)file3D_1b_tag->Get("h3D_b1");
    TH2D *h1b_tag = (TH2D*)h3D_1b_tag->Project3D("zy")->Clone("h1b_tag");

    //Prepare scaled EEC(dr) histograms for data
    TH2D *h2D_scaled = (TH2D*)h3D->Project3D("zy")->Clone("h2D_scaled");
    TH2D *h2D_bkg_bb_scaled = (TH2D*)h2D_scaled->Clone("h2D_bkg_bb_scaled");
    
    //save efficiency correction
    TH2D *h_eff = (TH2D*)h_sigfrac->Clone("h_eff");
    h_eff->Reset();

    //Get number of bins
    int pt_bins = h3D->GetNbinsZ();
    int dr_bins = h3D->GetNbinsY();
    int mB_bins = h3D->GetNbinsX();

    std::cout << pt_bins   << dr_bins   << mB_bins << std::endl; ////AQAA
    //Get the tagging efficiency
    for(int bin_pt = 1; bin_pt <= pt_bins; bin_pt++){
        for(int bin_dr = 1; bin_dr <= dr_bins; bin_dr++){

            Float_t tag_eff = h1b_tag->GetBinContent(bin_dr, bin_pt)/h1b_notag->GetBinContent(bin_dr, bin_pt);
            h_eff->SetBinContent(bin_dr, bin_pt, tag_eff);
        
        }
    }

    //Multiply data by signal and background fraction and efficiency correct
    h2D_scaled->Multiply(h_sigfrac);
    h2D_bkg_bb_scaled->Multiply(h_bkg_bb_fraction);
   
    //Correct already for b-tagging efficiency
    if(eff_corr){
        h2D_scaled->Divide(h_eff);
        h2D_bkg_bb_scaled->Divide(h_eff);
    }

    //Create slices in jtpt and scale them
    TH1D *heec = (TH1D*)h2D_scaled->ProjectionX("heec", pt_bin, pt_bin);

    TH1D *heec_bkg_bb = (TH1D*)h2D_bkg_bb_scaled->ProjectionX("heec_bkg_bb", pt_bin, pt_bin);

   
    //Normalise
    if(norm){
        heec->Scale(1/heec->Integral(), "width");
        heec_bkg_bb->Scale(1/heec_bkg_bb->Integral(), "width");
        h1b->Scale(1/h1b->Integral(), "width");
        hmoreb->Scale(1/hmoreb->Integral(), "width");}

    //Plot
    
    TLegend *leg = new TLegend(0.6,0.7,0.9,0.9, "Legend");

    //Fitted histograms
    heec->SetStats(0);
    if(data) dataset = "data";
    if(!eff_corr) heec->SetTitle("EEC distribution integrated, " + dataset + " (not efficiency corrected) "+ " (MC from dijet), pt bin = " + pt_bin);
    else heec->SetTitle("EEC distribution integrated, " + dataset + ", after efficiency correction, "+  + " (MC from dijet), pt_bin = " + pt_bin);
    heec->GetXaxis()->SetTitle("\\Delta r");
    heec->GetXaxis()->CenterTitle(true);
    if(norm) heec->GetYaxis()->SetTitle("eec (norm)");
    else heec->GetYaxis()->SetTitle("eec");
    heec->GetYaxis()->CenterTitle(true);
    //if (!norm) heec->GetYaxis()->SetRangeUser(10000, 7000000);
    if(norm) heec->GetYaxis()->SetRangeUser(0, 20);
    //heec->GetXaxis()->SetRangeUser(0.001, 1);
    heec->SetMarkerStyle(20);
    heec->SetMarkerColor(kRed);
    heec->SetLineColor(kRed);
    leg->AddEntry(heec, "signal fitted (data)");
    heec->Draw("P0EHIST");
    
    //Plot all the data and MC distributions
    if(all){
        heec_bkg_bb->SetMarkerStyle(20);
        heec_bkg_bb->SetMarkerColor(kBlue);
        heec_bkg_bb->SetLineColor(kBlue);
        leg->AddEntry(heec_bkg_bb, "more B background fitted (data)");
        heec_bkg_bb->Draw("P0EHIST SAME");
       

        //Comparisons
        h1b->Draw("P0EHIST SAME");
        h1b->SetMarkerStyle(24);
        h1b->SetMarkerColor(kRed+1);
        h1b->SetLineColor(kRed+1);
        hmoreb->Draw("P0EHIST SAME");
        hmoreb->SetMarkerStyle(24);
        hmoreb->SetMarkerColor(kBlue+1);
        hmoreb->SetLineColor(kBlue+1);
        hother->Draw("P0EHIST SAME");
        hother->SetMarkerStyle(24);
        hother->SetMarkerColor(kOrange);
        hother->SetLineColor(kOrange);
        hmc->Draw("P0EHIST SAME");
        hmc->SetMarkerStyle(3);
        hmc->SetMarkerColor(kMagenta);
        hmc->SetLineColor(kMagenta);
    

        leg->AddEntry(hmc, "MC not fitted");
        leg->AddEntry(h1b, "1 B");
        leg->AddEntry(hmoreb, "more B");
        leg->AddEntry(hother, "other");

        leg->Draw("SAME");
    }

    TString plot_name = "eec_template_fit_check_" + dataset + "_" + pT_selection;
    if(eff_corr) plot_name += "_effcorr_";
    if(norm) plot_name += "norm.pdf";
    c->Print(folder + plot_name);

    TFile *file_efficiency = new TFile(folder + "file_efficiency_" + dataset + "_" + pT_selection + ".root", "recreate");
    h_eff->Write();
    file_efficiency->Close();
}*/

void do_template_fit_data_Afnan(TString &dataset, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name){
    /*This fit use the data and sign and bkg templates from the input files.*/
    // Templates from Combined sample: Bjet+dijet: simple addition of hist; the relative ratio 1b/(1b + 2b) or 2b/(sum) is same in both samples.
    // Data his are Real data.


    bool isSampled = false; 

    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");


    // ---------------------------------------------------------------------------------------------------------
    // -- Temporare definiton: to be changed later to be common to the Dijet template 
    // TFile *file_bjet = new TFile("template_for_fit_histos_3D_makevtx_btag.root", "read");
        // if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());return;}
    // ---------------------------------------------------------------------------------------------------------

/*
    // -- Input Histogram name: MC closure test  
    TString namehMC = "h3D_data";
    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";// h3D_bb
    TString nameh1B = "h3D_b"; // h3D_b
    // -- Add template of 0B case: ht_3d_trks_0Bbkgfor2Bsignal or h3D_nob
     TString nameh0B = "h3D_nob";// h3D_nob
     TString namehmore2B = "h3D_more2b"; // h3D_more2b
*/

// -- Using Zoe file names: without sampling! 
    TString namehMC = "h3D_data";
    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_2b";
    TString nameh1B = "h3D_1b";
    TString nameh0B = "h3D_0b";

    
    //-- Dijet sample: Templates: Open file and get histograms from Dijet: templates 
    TFile *file_dijet = new TFile(folder + templates, "read");
        if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());return;}
    TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b");
    TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
        if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
        if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
    // -- 0B and more B contribution for drawing only
    
    // TH3D *h3D_more2b = (TH3D*)file_dijet->Get(namehmore2B)->Clone(namehmore2B);
       // if(!h3D_more2b){Error("Get:", "histogram does not exist '%s' ",h3D_more2b->GetName()); return;}

    TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone(nameh0B);
       if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}


    //-- Bjet sample: Templates: Open file and get histograms from bjet: templates 
    TFile *file_bjet = new TFile(folder + templates_bjet, "read");
    TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
    TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
        if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
        if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}


    //Open dataset:
    TFile *file_data = new TFile(folder + dataset, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
        cout << "file dat name " << file_data->GetName() << endl;

        TH3D *h3D_data;
        TH3D *h3D_data_bjet;
        if(isSampled){
            h3D_data = (TH3D*)file_data->Get(namehMC)->Clone(namehMC);
            h3D_data_bjet = (TH3D*)file_bjet->Get(namehMC)->Clone("h3D_data_bjet");
        }
        else {
            // Trivial closure test

            // do sum: 1B and 2B!
            // from Dijet sample 
            h3D_data = (TH3D*) h3D_b->Clone(namehMC);
            h3D_data->Add(h3D_bb);
            h3D_data->Add(h3D_nob);

            // then add bjet sample also
            h3D_data_bjet = (TH3D*) h3D_b_bjet->Clone("h3D_data_bjet");
            h3D_data_bjet->Add(h3D_bb_bjet);
        }
        if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}
        if(!h3D_data_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_data_bjet->GetName()); return;}


    //Get number of bin entries
    Int_t bins_pt = h3D_data->GetNbinsZ();
    Int_t bins_dr = h3D_b->GetNbinsY(); // Re-defined bins(exists in the header file)
    Int_t mb_bins = h3D_data->GetNbinsX();
    

    // -- Draw THstach for mb projection 
    // Get projections X 
    // hists names as for slices in bins (but outside {} so it is ok)
        TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX("h_data_mb");
        TH1D *h_data_bjet_mb = (TH1D *) h3D_data_bjet->ProjectionX("h_data_bjet_mb");

        TH1D *h_b = (TH1D *) h3D_b->ProjectionX("h_b");
        TH1D *h_b_bjet = (TH1D *) h3D_b_bjet->ProjectionX("h_b_bjet");
        
        TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX("h_bb");
        TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet->ProjectionX("h_bb_bjet");
        
        // TH1D *h_more2b = (TH1D *) h3D_more2b->ProjectionX("h_more2b");
        TH1D *h_nob = (TH1D *) h3D_nob->ProjectionX("h_nob");
        
        // -- Set titles for the hisograms  to be used in the default legends 
        h_data_mb->SetTitle("0B+1B+2B (dijet)");
        h_data_bjet_mb->SetTitle("1B+2B (bjet)");
        h_b->SetTitle("1B (dijet)");
        h_b_bjet->SetTitle("1B (bjet)");
        h_bb->SetTitle("2B (dijet)");
        h_bb_bjet->SetTitle("2B (bjet)");
        h_nob->SetTitle("0B");


        // Get histogram of the residual difference ebtween FullMC - (BB signal + B bkg):
        TH1D *hRes = (TH1D*) h_data_mb->Clone("hRes");
            hRes->SetDirectory(0);
            hRes->SetTitle("Residual difference between MC and (sig + bkg)");
            hRes->Add(h_data_bjet_mb, +1);
            hRes->Add(h_b, -1);
            hRes->Add(h_b_bjet, -1);
            hRes->Add(h_bb, -1);
            hRes->Add(h_bb_bjet, -1);
            hRes->Add(h_nob, -1);
            // hRes->Add(h_more2b, -1);
            hRes->SetLineColor(kBlack);
            hRes->SetLineWidth(2);
            fout->cd();
            hRes->Write();

        // The trivial sum of sig+bkg t compare to full MC 
        TH1D *hsum = (TH1D*) h_b->Clone("hsum");
        hsum->SetTitle("BB signal + 1B bkg + 0B bkg");
        hsum->Add(h_b_bjet);
        hsum->Add(h_bb);
        hsum->Add(h_bb_bjet);
        hsum->Add(h_nob);
        // hsum->Add(h_more2b);
        hsum->SetDirectory(0);
        hsum->Write();

       // -- Deattach hists from root file 
        h_data_mb->SetDirectory(0);
        h_b->SetDirectory(0);
        h_bb->SetDirectory(0);
        h_data_bjet_mb->SetDirectory(0);
        h_b_bjet->SetDirectory(0);
        h_bb_bjet->SetDirectory(0);
        // h_more2b->SetDirectory(0);
        h_nob->SetDirectory(0);
       
        // Set styles of hist to be drawn:
        h_data_mb->SetLineColor(kBlack);
        h_data_mb->SetMarkerColor(kBlack);
        h_data_mb->SetMarkerStyle(4);

        h_b->SetFillColor(kRed-10); h_b->SetLineColor(kRed-10);        
        h_bb->SetFillColor(kBlue-10); h_bb->SetLineColor(kBlue-10);

        h_data_bjet_mb->SetLineColor(kGray+1);
        h_data_bjet_mb->SetMarkerColor(kGray+1);
        h_data_bjet_mb->SetMarkerStyle(8);

        h_b_bjet->SetFillColor(kGreen-2); h_b_bjet->SetLineColor(kGreen-2);        
        h_bb_bjet->SetFillColor(kOrange-6); h_bb_bjet->SetLineColor(kOrange-6); // brown 

        // h_more2b->SetFillColor(kOrange); h_more2b->SetLineColor(kOrange);

        h_nob->SetFillColor(kRed); h_nob->SetLineColor(kRed);

        hsum->SetLineWidth(3);
        hsum->SetLineColor(kMagenta+2);

    THStack hstack_All("hstack_All","Mass stacked histogram");
    hstack_All.SetTitle(";m_{2B} [GeV];");
    // hstack_All.Add(h_more2b);
    hstack_All.Add(h_nob);
    hstack_All.Add(h_b);
    hstack_All.Add(h_bb);
    hstack_All.Add(h_b_bjet);
    hstack_All.Add(h_bb_bjet);
    
    auto canva = new TCanvas("All_templates",Form("Templates(mb) and Full MC"), 800, 800 );
        canva->cd();
        canva->SetLogy();
        hstack_All.Draw("hist E"); // hist
        TH1D *h_data_mb_all = (TH1D*) h_data_mb->Clone("h_data_mb_all");
            h_data_mb_all->SetTitle("Dijet + bjet sample");
            h_data_mb_all->SetMarkerStyle(8);
            h_data_mb_all->Add(h_data_bjet_mb, 1);
            h_data_mb_all->Draw("PE same");  // full MC as points
            h_data_mb->Draw("HIST E same");
            h_data_bjet_mb->Draw("HIST E same");

            // For names to appear in the legend :
               h_data_mb->SetTitle("0B+1B+2B (dijet)");
                h_data_bjet_mb->SetTitle("1B+2B (bjet)");
                h_b->SetTitle("1B (dijet)");
                h_b_bjet->SetTitle("1B (bjet)");
                h_bb->SetTitle("2B (dijet)");
                h_bb_bjet->SetTitle("2B (bjet)");
                h_nob->SetTitle("0B");

        // hsum->Draw("histsame");
        gPad->Modified();   
        gPad->Update();
        canva->Modified();
        canva->Update();
        fout->cd();
        canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
        canva->Write();
        canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));

    // -- Output hist declarations: empty hists
    // -- For the fit result: per (deltaR, jtpt) bin
        //Create the histograms storing the fit parameters in TH2D (the dr are stored in the X, the jtpt in the Y)
        TH2D *h_sig_fraction = (TH2D*) h3D_data->Project3D("zy")->Clone("h_sig_fraction");// to get the axes from hist
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
          


    //Vector to test the convergence
    std::vector<std::pair<int, int>> non_converge_bins;


    // Fitting - loop over dr and jtpt entries
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){
        // for(Int_t ibin_dr = 0; ibin_dr <= 0; ibin_dr++){
        // Bin = 0: integaretd deltaR range, other bins are the DeltaR bins 
            
            // The slice range for projection is different for integarted Vs. one bin slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;

            // For integarted bin only: bin = 0
            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = bins_dr;}

            // Make slices for data
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_data_bjet_mb = (TH1D *) h3D_data_bjet->ProjectionX(Form("h_data_mb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);

            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            // TH1D *h_more2b = (TH1D *) h3D_more2b ->ProjectionX(Form("h3D_more2b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            //  Make slices for Bjet 
            TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            

            // -- Deattach hists from root file 
            h_data_mb->SetDirectory(0);
            h_b->SetDirectory(0);
            h_bb->SetDirectory(0);
            h_b_bjet->SetDirectory(0);
            h_bb_bjet->SetDirectory(0);
            h_nob->SetDirectory(0);
            // h_more2b->SetDirectory(0);

  TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);

            // Calculate true fractions to be used as initial values for the fit (the true fractions are the dijet qcd ones)            
            double int0 = h_bb->Integral(1, mb_bins, "width");// sig
            double int1 = h_b->Integral(1, mb_bins, "width");// bkg1
            // double int2 = h_more2b->Integral(1, mb_bins, "width"); // bkg2
            double int3 = h_nob->Integral(1, mb_bins, "width"); // bkg 3

            // From bjet sample 
            double int0_bjet = h_bb_bjet->Integral(1, mb_bins, "width");// sig
            double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");// bkg1

            // For MC test case only
            double integral_inputdata1 = h_data_mb->Integral(1, mb_bins,"width");
            double integral_inputdata2 = h_data_bjet_mb->Integral(1, mb_bins,"width");
                cout << "Before Samples additions:\n Dijet data integral = "<< integral_inputdata1 << endl;
                cout << "Bjet data integral = "<< integral_inputdata2 << endl;

            // Add Dijet and bjet sample for the pseudodata 
            h_data_mb->Add(h_data_bjet_mb);

            // For data: After the sum of dijet + bjet 
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
                
                cout << "Total input data integral = "<< integral_inputdata << endl;
                std::cout << "int0 h_bb=" << int0 << std::endl;
                std::cout << "int1 h_b =" << int1 << std::endl;
                cout << "2B/1B in Dijet sample = " << int0/int1 << endl;
                cout << "Dijet: 2B/(1B + 2B) = " << int0/(int0 + int1) << endl;

                std::cout << "int0 h_bb_bjet=" << int0_bjet << std::endl;
                std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                cout <<  "2B/1B in  Bjet sample = " << int0_bjet/int1_bjet << endl;
                cout << "2B/(2B + 1B) in Bjet sample = " << int0/(int0 + int1) << endl;
                // std::cout << "int2 h_more2b =" << int2 << std::endl;
                std::cout << "int3 h_nob =" << int3 << std::endl;


            // --- For now: the TRUE signal (2B) fraction is used from the Dijet sample.
                // the fraction of 2B/1B in Bjet and Dijet sample is exactly the same
            // double sig_fraction_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int2 + int3 )); 
            // double bkg_fraction_b_true = (int0  + int1 + int2 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int2 + int3 ));
            
            // No contiubtion from hmore2B 
            double sig_fraction_true = (int0  + int1 + int3 ) == 0 ? 0 : (int0 / (int0 + int1 + int3 )); 
            double bkg_fraction_b_true = (int0  + int1 + int3 ) == 0 ? 0 : (int1 / (int0 + int1 + int3 ));


                std::cout << "Dijet : 2B sig_fraction_true =" << sig_fraction_true << std::endl;
                std::cout << "Dijet: 1B bkg_fraction_b_true =" << bkg_fraction_b_true << std::endl;
            // -- Compute the errors of true fractions for later fill
            double True_bkg_b_err = 0.0;
            double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");
            double True_sig_err = 0.0;
            double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");


            // -- Get Normalized templates as input  PDF to the fit(shape matters only): where PDF integral is 1 (Probability density function)
            //signal: 2B from Dijet + Bjet sample 
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
            h_sumsig->Add(h_bb_bjet);
            Int_t h_sig_bins = h_sumsig->GetNbinsX();
            TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized


            //Add background: here only from 1B 
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
            h_sumbkg->Add(h_b_bjet);
            Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
            
            TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized
            
            // Total sum bkg = 0B + 1B  
            TH1D* h_sumbkg_0b_1b = (TH1D*)  h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
            h_sumbkg_0b_1b->Add(h_nob);


            // write signal hists after combine Dijet and bjet sample 
            fout->cd();
            h_sumsig->Write();
            h_sumbkg->Write();
            h_sumbkg_0b_1b->Write();
           

            // -- Draw Stack for combined: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before ormalization 
              THStack hstack_templatesfotfit(Form("hstack_templatesfotfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesfotfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesfotfit.Add(h_nob);
                    hstack_templatesfotfit.Add(h_sumbkg);
                    hstack_templatesfotfit.Add(h_sumsig);
                    h_sumsig->SetTitle("2B: Dijet + bjet");
                    h_sumbkg->SetTitle("1B: Dijet + bjet");
                    h_nob->SetTitle("0B");

                    // -- Set fill colors 
                    h_nob->SetFillColor(kRed);
                    h_sumsig->SetFillColor(kBlue -10);
                    h_sumbkg->SetFillColor(kRed-10);
                     h_sumbkg->SetFillStyle(3244);
                    h_data_mb->SetMarkerColor(kBlack); h_data_mb->SetMarkerStyle(20);

              auto canva = new TCanvas(Form("All_templates_%d_%d", ibin_dr, ibin_pt),Form("Templates(mb) and Full MC"), 800, 800 );
                    canva->cd();
                    canva->SetLogy();
                    hstack_templatesfotfit.Draw("hist E"); 
                    h_data_mb->SetTitle("All");
                    h_data_mb->Draw("PE same");
                    gPad->Modified();   
                    gPad->Update();
                    canva->Modified();
                    canva->Update();
                    fout->cd();
                    canva->BuildLegend(0.7, 0.9, 0.7, 0.9);
                    canva->Write();
                    canva->Print(Form("%s/%s.png", sDirname.Data(), canva->GetName()));




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
                cout << "True 1B integral without width option =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;


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
            
            // -- test 
               // eff_bkg0B = 0;
               // eff_bkg1B = 1; 

            // Build effective signal hist: with new relaitve normalization, the integral should = 1
            h_bkg->Add(h_bkg, h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist
            
            double err_int;
            double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins,err_int ,"width");
            cout << "Effective Bkg normalized hist integral +/- uncertainity = "<< int_val << "+/-" << err_int << endl;

            cout << "effective 1B fraction and 0B = "<< eff_bkg1B << ", " << eff_bkg0B << endl;
            // -- Note: the effective bkg integral is Effective Bkg normalized hist integral = 0.997777
            // if (h_bkg->Integral(1,  h_bkg_bins, "width") != 1.) {cout << "Effective BKG PDF is not normalized!"<< endl; return;}

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
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true, 0.,1.);// 1-bkg_fraction_b_true

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, true); //fit the scaling when adding the single b with the more b distribution
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name
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
            double err_true_frac = (int1 * True_sig_err + int0 * True_bkg_b_err)/(int0+ int1)/(int0+ int1);// int0 for signal, int1 for bkg 

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
                h_sig_fit->SetTitle("After fit: Sig");
                h_sig_fit->SetMarkerColor(kRed+2);
                h_sig_fit->SetLineColor(kRed+2);
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg 
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("After fit: bkg (1B+ 0B)");
                h_bkg_fit->SetMarkerColor(kCyan+2);
                h_bkg_fit->SetLineColor(kCyan+2);

            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only 
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("After fit: 1B");
                h_bkg_fit_1b->SetMarkerColor(kOrange+2);
                h_bkg_fit_1b->SetLineColor(kOrange+2);
                h_bkg_fit_1b->SetLineWidth(2);

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only 
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("After fit: 0B");
                h_bkg_fit_nob->SetMarkerColor(kYellow+2);
                h_bkg_fit_nob->SetLineColor(kYellow+2);
                h_bkg_fit_nob->SetLineWidth(2);

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("After fit: sig + bkg");
                    h_total_fit->SetMarkerColor(kGreen+2);
                    h_total_fit->SetLineColor(kGreen+2);
                    h_total_fit->SetMarkerStyle(20);
                    h_total_fit->SetLineWidth(2);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();
            h_bkg_fit->Write();
            h_total_fit->Write();

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
            auto canva_afterfit = new TCanvas(sname_canvas_afterfit,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetMaximum(1.2 * h_data_mb->GetMaximum());

                h_data_mb->SetTitle("Pseudodata (combined)");
                // h_b->SetTitle("B template (dijet)");
                // h_bb->SetTitle("2B template (dijet)");
                // h_b_bjet->SetTitle("B template (bjet)");
                // h_bb_bjet->SetTitle("2B template (bjet)");
                h_sumbkg->SetTitle("1B (bkg) template (combined)");
                h_sumsig->SetTitle("Signal template: 2B (combined)");
                h_sumbkg_0b_1b->SetTitle("Bkg template (0B + 1B) (combined)");
                h_nob->SetTitle("0B template (Dijet only)");

                h_sumbkg->SetLineColor(kOrange);
                h_nob->SetLineColor(kMagenta);
                h_data_mb->SetLineWidth(2);


                // -- Make better visual
                // h_bb->SetLineWidth(3);
                // h_bb_bjet->SetLineWidth(3);
                h_sumsig->SetLineWidth(3);
                h_sig_fit->SetLineWidth(3);

                // h_b->SetLineStyle(9); // close dashes
                // h_b_bjet->SetLineStyle(9);
                h_sumbkg_0b_1b->SetLineStyle(7);
                h_sumbkg_0b_1b->SetLineWidth(3);
                h_bkg_fit->SetLineWidth(2);
                h_bkg_fit->SetLineStyle(10);  
                
                // Add after fit templates
                h_data_mb->Draw("PE");
                h_total_fit->Draw("P E SAME");

                // h_bb->Draw("HIST E SAME");
                // h_bb_bjet->Draw("HIST E SAME");
                h_sumsig->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                

                h_sumbkg->Draw("HIST E SAME"); // Combined sample for 1B bkg 
                h_nob->Draw("HIST E SAME");
                h_bkg_fit_nob->Draw("HIST E SAME");
                h_bkg_fit_1b->Draw("HIST E SAME");
                h_sumbkg_0b_1b->Draw("HIST E SAME");
                h_bkg_fit->Draw("HIST E SAME");

                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
                fout->cd();

                TString srangeDeltaR = Form(" %g < #DeltaR [bin %d] < %g", dr_binsVector[ibin_dr -1], ibin_dr ,dr_binsVector[ibin_dr]);
                if(!ibin_dr) srangeDeltaR = Form(" %g < #DeltaR [Full range] < %g", dr_binsVector[1], dr_binsVector[bins_dr]);
                canva_afterfit->BuildLegend(0.7, 0.9, 0.7, 0.9, srangeDeltaR);
                // Add another line for all bins 
                canva_afterfit->Write();


                canva_afterfit->Print(Form("%s/%s.png", sDirname.Data(), sname_canvas_afterfit.Data()));

                cout << "---------------------\n\n\n" << endl; 
            } // loop over deltaR bins 
    }

    // Save histograms
    // TH3D 
    // for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_data_bjet, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
    for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
        

        
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
        
    fout->Close();

}
void template_fit_test(){

    //Get data and mc labels
    TString pT_selection = "80_140";// full range 
    TString folder = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/TestTemplateFits_newcode/";
   
    //Normalise the eec plot
    bool norm = true;

    //Plot also the mc eec
    bool all = true;

    //Correct for tagging efficiency
    bool eff_corr = false;

    // Add LowEG data 
    bool alsoLowEG = true;

    // -- Output folder to save the result of the tests 
    gSystem->mkdir(sDirname, kTRUE); // make directory if does not exist for the output files to be all in one place 

    TString sfoutputPlots_dijet = "Summary_histo_templatefit_test_26mars.root";
    TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "RECREATE"); // For canvas output
    if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    // --  HG data and bjet with HLT80: using result of Zoe code and condor result 
    alsoLowEG = false;
    TString dataset_HG = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/condor/condor_output/template_for_fit_histos_3D_HighEG_btag_HLT80.root"; 
        TString dataname = "HighEG_HLT80";
        TString templates_dijet = ""; // not used now 
        TString templates_bjet = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/condor/condor_output/template_for_fit_histos_3D_bjet_btag_HLT80.root";
        TString fout_name = "TemplateFits_histos_3d_HLT80_HighEGdata_bjet_" + pT_selection +  ".root";
        do_template_fit_bjetonly(dataset_HG, templates_dijet, templates_bjet,  pT_selection, folder, fout_name, alsoLowEG);


   
    foutputPlots_dijet->Print();
    foutputPlots_dijet->Close();
    delete foutputPlots_dijet;
} 