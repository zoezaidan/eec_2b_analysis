

#include "tTree.h"
#include "binning_histos_all.h"
#include "Help_Functions.h"



void do_template_fit_combined(const TString &HighEGdata_name, const TString &LowEGdata_name, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name, bool& alsoLowEG){
    /*
        // ----- WORK in PROGRESS ---- 14/04/2026
        // In the new version of files: 2B contribution has both 2B + morethan2B in the same histogram. No anymore seprated histogram for it.
        // 0B onl from qcd sample. bjet sample 0B contribution is not physical (it is already filtered to be bjets so these light quark templates are not real but mistaken).
    */

    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");

    // -- Using new naming convention
    TString namehData = "h3D_data"; // h3D_data with eec on the weight 
    TString namehMC = "";

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";
    TString nameh1B = "h3D_b";

    TString nameh0B = "h3D_0b";
    // -- Additional histograms
    TString namehmore2B = "";

 
    //-- Dijet sample:
    TFile *file_dijet = new TFile(templates, "read");
            if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());}
        TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b");
        TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
            if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
            if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
        TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone(nameh0B);
            if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}

    //-- Bjet sample:
    TFile *file_bjet = new TFile(templates_bjet, "read");
         if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}
        TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
        TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
        if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
        if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}

    //Open dataset:
    TFile *file_data = new TFile(HighEGdata_name, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
            cout << "file data name " << file_data->GetName() << endl;
        TH3D *h3D_data;
            h3D_data = (TH3D*)file_data->Get(namehData)->Clone("h3D_data");
            if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}
            if(alsoLowEG){
                TFile *file_data_LowEG = new TFile(LowEGdata_name, "read");
                    if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
                    cout << "file data name " << file_data_LowEG->GetName() << endl;
                TH3D *h3D_dataLowEG;
                    h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
                    if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ",h3D_dataLowEG->GetName()); return;}
                h3D_data->Add(h3D_dataLowEG);
            }
        // -- data: set style 
        h3D_data->SetTitle("Data");
        h3D_data->SetMarkerStyle(20);
        h3D_data->SetMarkerColor(kBlack);
        h3D_data->SetLineColor(kBlack);
        h3D_data->SetLineWidth(3);

        // bjet 
        h3D_b_bjet->SetTitle("1B (bjet)");h3D_bb_bjet->SetTitle("2B (bjet)");
        h3D_b_bjet->SetFillColor(kMagenta-5); h3D_b_bjet->SetLineColor(kMagenta-5);        
        h3D_bb_bjet->SetFillColor(kMagenta-3); h3D_bb_bjet->SetLineColor(kMagenta-3);


        // qcd
        h3D_b->SetTitle("1B (qcd)");h3D_bb->SetTitle("2B (qcd)"); h3D_nob ->SetTitle("0B (qcd)");
        h3D_b->SetFillColor(kBlue-4); h3D_b->SetLineColor(kBlue-4);        
        h3D_bb->SetFillColor(kBlue-6); h3D_bb->SetLineColor(kBlue-6); 
        // h3D_nob->SetFillColor(kYellow-6); h3D_nob->SetLineColor(kYellow-6); // not visible 
        h3D_nob->SetFillColor(kRed); h3D_nob->SetLineColor(kRed); 

    // -- general style: remove stats box 
    gStyle->SetOptStat(0);


    //define used bins 
    Int_t bins_pt = h3D_data->GetNbinsZ();
    Int_t bins_dr = h3D_data->GetNbinsY(); 
    Int_t mb_bins = h3D_data->GetNbinsX();
        cout << "-- Data hist binning" << endl;
        cout << "pt bins = "<< bins_pt << endl;
        cout << "dr bins = "<< bins_dr << endl;
        cout << "mb bins = "<< mb_bins << endl;

    /// -- For later drawing of S/B fractions, store true and fit result in the following histograms
        // Note that: since inetgaretd bins are pt 0 and dr 0, the hist of fractions should have #bins + 1 size 
        //  // x = dr, y = jetpt
        // Axis here is for the bin number instead of the values 
          TH2D *h_sig_fraction = new TH2D("h_sig_fraction", ";dr; jet pt",h3D_data->GetNbinsY() + 1,  1,  h3D_data->GetNbinsY() + 2,
                                                                          h3D_data->GetNbinsZ() + 1 , 1, h3D_data->GetNbinsZ()+ 2 ); // yes GetXmin() fpr Y an Z-axis
                                                                                                                                                 
            h_sig_fraction->Reset();
            TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
            TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
            TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error"); 
            TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
            TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");
            TH2D *h_bkg_frac_true = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true");
            TH2D *h_bkg_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true_error");


    // --- Vector to test the convergence
    std::vector <std::pair<int, int>> non_converge_bins;


    // Fitting - loop over dr and jtpt entries
    // Bin0: is integarted over the range 
    for(Int_t ibin_pt = 0; ibin_pt <= bins_pt; ibin_pt++){
    // for(Int_t ibin_pt = 0; ibin_pt <= 0; ibin_pt++){
        for(Int_t ibin_dr = 0; ibin_dr <= bins_dr; ibin_dr++){
        // for(Int_t ibin_dr = 0; ibin_dr <= 0; ibin_dr++){
            
            // define slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;
            Int_t SliceFirstbin_pt = ibin_pt;
            Int_t SliceLastbin_pt =  ibin_pt;
            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = bins_dr;}
            if (!ibin_pt){ SliceFirstbin_pt = 1;  SliceLastbin_pt = bins_pt;}

            // Make projections 
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
                h_data_mb->GetXaxis()->SetTitle("m_{2B} [GeV]");
                h_data_mb->SetTitle(h3D_data->GetTitle()); // upadte projection title

            
            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
                h_bb->SetTitle( h3D_bb->GetTitle() );
                h_b->SetTitle( h3D_b->GetTitle() );
                h_nob->SetTitle( h3D_nob->GetTitle() );

            //  Make slices for bjet 
            TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
            TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr,  SliceFirstbin_pt, SliceLastbin_pt);
                h_b_bjet->SetTitle(h3D_b_bjet->GetTitle());
                h_bb_bjet->SetTitle(h3D_bb_bjet->GetTitle());
      
            // -- Deattach hists from root file 
            h_data_mb->SetDirectory(0);
            h_b->SetDirectory(0);
            h_bb->SetDirectory(0);
            h_b_bjet->SetDirectory(0);
            h_bb_bjet->SetDirectory(0);
            h_nob->SetDirectory(0);

            // string (dr, pt) 
            TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);

            // --  Calculate true fractions to be used as initial values for the fit (the true fractions are the qcd ones)            
            double int2 = h_bb->Integral(1, mb_bins, "width");
            double int1 = h_b ->Integral(1, mb_bins, "width");
            double int0 = h_nob->Integral(1, mb_bins, "width"); 
                // compute the true fractions and their errors (for later comparison)           
                double sig_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : (int2 / (int0 + int1 + int2)); 
                double bkg_fraction_b_true = (int0  + int1 + int2 ) == 0 ? 0 : (int1 / (int0 + int1 + int2));
                double bkg_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : ( (int0 + int1)/(int0 + int1 + int2) );

                double True_bkg_b_err = 0.0;
                double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");

                double True_sig_err = 0.0;
                double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");

                // total bkg fraction: 0B + 1B (qcd only)
                TH1D* h_b_nob_dijet = (TH1D*) h_b ->Clone("h_b_nob_dijet");  h_b_nob_dijet->Add(h_nob);
                double True_bkg_err = 0.0;
                double True_bkg_integral = h_b_nob_dijet->IntegralAndError(1, mb_bins, True_bkg_err, "width");


            // -- Compute other useful integrals 
            // From bjet sample 
            double int2_bjet = h_bb_bjet->Integral(1, mb_bins, "width"); 
            double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");
            // data:
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
                
                // -- check integrals in bjet and qcd sample
                // relative fraction of 1B: 2B in bjet and qcd samples is SAME. So you can combine the two samples with simple + (without reweighting).
                // cout << "Total input data integral = "<< integral_inputdata << endl;
                // std::cout << "int2 h_bb=" << int2 << std::endl;
                // std::cout << "int1 h_b =" << int1 << std::endl;
                // cout << "Dijet: int0 of 0B" << int0 << endl;

                // cout << "Dijet: 2B/(1B + 2B) = " << int2/(int2 + int1) << endl;
                // std::cout << "int2 h_bb_bjet=" << int2_bjet << std::endl;
                // std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                // cout << "2B/(2B + 1B) in qcd sample = " << int2_bjet/(int2_bjet + int1_bjet) << endl;
                // cout << " Dijet: 0B/(1B+2B) = " << int0/(int2 + int1) << endl;

            // --- Prepare PDFs for template fit

        
            // 1- Combine qcd + bjet samples 
            // Signal: 2B   
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
                h_sumsig->Add(h_bb_bjet);
                h_sumsig->SetTitle("2B: qcd+bjet");
                Int_t h_sig_bins = h_sumsig->GetNbinsX();
                TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

            // Bkg: first 1B, and later added effectively 0B
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
                h_sumbkg->Add(h_b_bjet);
                h_sumbkg->SetTitle("1B: qcd+bjet");
                Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
                TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

                // for drawings
                    // Total sum bkg = 0B + 1B  
                    TH1D* h_sumbkg_0b_1b = (TH1D*) h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
                    h_sumbkg_0b_1b->Add(h_nob);
                    h_sumbkg_0b_1b->SetTitle("0B +1B (qcd+bjet)");


                // -- set the combined samples styles before fit 
                h_sumbkg->SetFillColor(kGreen-2); h_sumbkg->SetLineColor(kGreen-2);        
                h_sumsig->SetFillColor(kOrange-6); h_sumsig->SetLineColor(kOrange-6);
                    // and for the bkg 1B + 0B 
                    h_sumbkg_0b_1b->SetFillColor(kGreen); h_sumbkg_0b_1b->SetLineColor(kGreen);  

                // write to rootfile
                fout->cd();
                h_data_mb->Write(); // data 
                h_sumsig->Write(); // 2B (qcd + bjet)
                h_sumbkg->Write(); // 1B (qcd + bjet)
                h_sumbkg_0b_1b->Write(); // 1B (qcd + bjet) + 0B 


    // -------- Draw prefits 
        // contirbutions seperated  
                // absolute yields: seperated contributions: qcd, bjet 
            THStack hstack_all_beforefit (Form("hstack_all_beforefit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms before fit");
                hstack_all_beforefit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                hstack_all_beforefit.Add(h_bb_bjet);
                hstack_all_beforefit.Add(h_bb);
                hstack_all_beforefit.Add(h_b_bjet); 
                hstack_all_beforefit.Add(h_b);
                hstack_all_beforefit.Add(h_nob);
                hstack_all_beforefit.SetMaximum(10e+03);
                auto canva_beforefit = new TCanvas(Form("All_contributions_beforefit_%d_%d", ibin_dr, ibin_pt),"", 800, 800 );
                        canva_beforefit->cd();
                        // canva_beforefit->SetLogy();
                        hstack_all_beforefit.Draw("hist E");
                        TH1D* hd_norm_slice = (TH1D*) h_data_mb->Clone("hd_norm_slice");
                            hd_norm_slice->Scale(1./10000);
                            hd_norm_slice->Draw("HIST E same"); 
                            hd_norm_slice->SetTitle("data/1e+04");
                            gPad->Modified();   
                            gPad->Update();
                            canva_beforefit->Modified();
                            canva_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_beforefit, ibin_dr, ibin_pt);
                            canva_beforefit->Write();
                            canva_beforefit->Print(Form("%s/%s_allcontributions_beforefit.png", sDirname.Data(), sname_canvas.Data()));

        // -- Draw Stack for combined qcd + bjet: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before normalization
            // absolute yields of Signal and bkg  
            THStack hstack_templatesforfit(Form("hstack_templatesforfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesforfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesforfit.Add(h_sumsig);
                    hstack_templatesforfit.Add(h_sumbkg_0b_1b);                    
                    auto canva_sum_beforefit = new TCanvas(Form("templates_beforefit_%d_%d", ibin_dr, ibin_pt), "", 800, 800 );
                        canva_sum_beforefit->cd();
                        hstack_templatesforfit.Draw("hist E");
                        hd_norm_slice->Draw("HIST E same"); 
                            gPad->Modified();   
                            gPad->Update();
                            canva_sum_beforefit->Modified();
                            canva_sum_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_sum_beforefit, ibin_dr, ibin_pt);
                            canva_sum_beforefit->Write();
                            canva_sum_beforefit->Print(Form("%s/%s_templates_beforefit.png", sDirname.Data(), sname_canvas.Data()));

        // -- Safety for empty bins 
            if ((int0 + int1 + int2) < 1  ) { cout << " ----------- empty bin -----------  "; continue;}   
                          

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
                // normalize signal and bkg 
                h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));
                h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));
                    // cout << "After Normalization (qcd +bjet) True 2B integral = " << h_sig->Integral(1, h_sig_bins, "width") << endl;
                    // cout << "After Normalization (qcd+bjet)  True 1B integral =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;

            // -- Effective bkg PDF 
                // Assume : a + b + c = 1; where a, b, c, are fractions in in nominal MC like 
                // a: 2B fraction, b: 1B fraction, c: 0B fraction, in 
                // and For Effective bkg: b`+ c` should = 1, so relate b`to b .. and simialr for c`
                // This result in: b`= b/(b+c), and c`= c/(b+c); one can use b+c = 1-a;  
                // I will write it in terms of the 2B sig, 1B bkg true fractions, before use effecive ones 
                    double eff_bkg0B = (1 - sig_fraction_true - bkg_fraction_b_true)/(1- sig_fraction_true);// c`
                    double eff_bkg1B = 1. - eff_bkg0B;// b`
                    // -- when test bjet, no 0B
                       // eff_bkg0B = 0;
                       // eff_bkg1B = 1; 

            // Build effective bkg hist: with new relaitve normalization, the integral should = 1
            // Normlaize h_nob to be added effectively to the normalized bkg distribution
                TH1D* norm_h_nob = (TH1D*) h_nob->Clone("norm_h_nob"); norm_h_nob->Scale(1./norm_h_nob->Integral(1,  h_bkg_bins, "width"));
                h_bkg->Add(h_bkg, norm_h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist            
                // update combined normalized bkg integral and error  
                double err_int;
                double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins, err_int ,"width");
                    // cout << "Effective Bkg normalized hist integral +/- uncertainity = ? (should be 1):  "<< int_val << "+/-" << err_int << endl;
                    // cout << "effective 1B and 0B fractions = (sum should be 1): "<< eff_bkg1B << ", " << eff_bkg0B << endl;
                    // avoid float point failur 
                    if ( std::abs(int_val - 1)  > 1e-06) {cout << "Effective BKG PDF is not normalized to 1!"<< endl; return;}
                    if ( std::abs( h_sig->Integral(1,  h_sig_bins, "width") - 1)  > 1e-06) {cout << "Signal PDF is not normalized to 1!"<< endl; return;}

            // -- Draw normalized PDFs of signal and effective background as inputs: prefit
            // set PDF titles 
            h_bkg->SetTitle("1B + 0B");
            h_sig->SetTitle("2B");
            h_sig->SetFillStyle(3244);
            h_bkg->SetFillColor(kGreen+2);
            h_sig->SetFillColor(kOrange+4);
            h_bkg->SetLineColor(kGreen+2);
            h_sig->SetLineColor(kOrange+4);


            // normlaized data for comparison 
            TH1D* hnorm_data_self = (TH1D*)h_data_mb ->Clone("hnorm_data_self");
                hnorm_data_self->Scale(1./hnorm_data_self->Integral(1, mb_bins, "width"));
                hnorm_data_self->SetTitle("data self normalized");

                    auto canva_pdf_norm_beforefit = new TCanvas("Pdfs_norm_beforefit","", 800, 800 );
                        canva_pdf_norm_beforefit->cd();
                        h_bkg->Draw("hist E");
                        h_sig->Draw("hist E same");
                        hnorm_data_self->Draw("PE same");
                        canva_pdf_norm_beforefit->SetTitle("PDFs before fit");
                            gPad->Modified();   
                            gPad->Update();
                            canva_pdf_norm_beforefit->Modified();
                            canva_pdf_norm_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_pdf_norm_beforefit, ibin_dr, ibin_pt);
                            canva_pdf_norm_beforefit->Write();
                            canva_pdf_norm_beforefit->Print(Form("%s/%s_PDF_norm_beforefit.png", sDirname.Data(), sname_canvas.Data()));
            

            // what about a stack of PDFs before fit 
                // these are not normalized to 1, but normlaized such that the total PDFs are 1, using their qcd fractions  
                THStack hstack_pfds_scaledtoqcd (Form("hstack_pfds_scaledtoqcd_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true fractions, before fit");
                    hstack_pfds_scaledtoqcd.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true fractions in qcd 
                        TH1D* h_sig_scaled = (TH1D*) h_sig ->Clone("h_sig_scaled"); h_sig_scaled->Scale(sig_fraction_true);
                        TH1D* h_bkg_scaled = (TH1D*) h_bkg ->Clone("h_bkg_scaled"); h_bkg_scaled->Scale(1. - sig_fraction_true);
                            hstack_pfds_scaledtoqcd.Add(h_sig_scaled); 
                            hstack_pfds_scaledtoqcd.Add(h_bkg_scaled);
                            auto canva_pdfs_scaledtoqcd = new TCanvas("canva_pdfs_scaledtoqcd","", 800, 800 );
                                canva_pdfs_scaledtoqcd->cd();
                                hstack_pfds_scaledtoqcd.Draw("Hist E");
                                hnorm_data_self->Draw("HIST E same"); 
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd->Modified();
                                canva_pdfs_scaledtoqcd->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd, ibin_dr, ibin_pt);
                                canva_pdfs_scaledtoqcd->Write();
                                canva_pdfs_scaledtoqcd->Print(Form("%s/%s_pdfs_scaledtoqcdfractions_beforefit.png", sDirname.Data(), sname_canvas.Data()));
            

                // -- Prefit PDF: Signal, BKG
                THStack hstack_pfds_scaledtoqcd_int (Form("hstack_pfds_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        TH1D* h_sig_scaledint = (TH1D*) h_sig ->Clone("h_sig_scaledint"); h_sig_scaledint->Scale(int2);
                        TH1D* h_bkg_scaledint = (TH1D*) h_bkg ->Clone("h_bkg_scaledint"); h_bkg_scaledint->Scale(int1 + int0);
                            hstack_pfds_scaledtoqcd_int.Add(h_sig_scaledint); 
                            hstack_pfds_scaledtoqcd_int.Add(h_bkg_scaledint);
                            auto canva_pdfs_scaledtoqcd_int = new TCanvas("canva_pdfs_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_scaledtoqcd_int->cd();
                                hstack_pfds_scaledtoqcd_int.Draw("Hist E");
                                TH1D* hd_scaledtoqcdint = (TH1D*) hnorm_data_self->Clone("hd_scaledtoqcdint"); hd_scaledtoqcdint->Scale(int0+int1+int2);
                                        hd_scaledtoqcdint->SetTitle("data scaled to qcd integral");
                                        hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd_int->Modified();
                                canva_pdfs_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd_int, ibin_dr, ibin_pt);
                                canva_pdfs_scaledtoqcd_int->Write();
                                canva_pdfs_scaledtoqcd_int->Print(Form("%s/%s_pdfs_scaledtoqcd_int_beforefit.png", sDirname.Data(), sname_canvas.Data()));
            
                                    // test 0B: (1B + 2B) after scaling to dijet fractions: is it like the true fractions?
                                    cout << "----------------" << endl;
                                    cout << "After combine Dijet+ bjet, and scaling to qcd integrals " << endl;
                                    cout << "Integarl of scaled PDF(2B) = " <<   h_sig_scaledint->Integral() << endl;
                                    cout << "Integral of scaled PDF(1B+0B) = " << h_bkg_scaledint->Integral() << endl;

                // -- Prefit PDF: signal, BKG (seperated 0B and 1B)
                 THStack hstack_pfds_seperated_scaledtoqcd_int (Form("hstack_pfds_seperated_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_seperated_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        // 1B (bjet and dijet)
                        TH1D* hnorm_sumbkg_1B_scaledint = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledint");
                                 hnorm_sumbkg_1B_scaledint->Scale(1./hnorm_sumbkg_1B_scaledint->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledint->Scale(int1);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_sig_scaledint); 
                        hstack_pfds_seperated_scaledtoqcd_int.Add(hnorm_sumbkg_1B_scaledint);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_nob);

                            auto canva_pdfs_seperated_scaledtoqcd_int = new TCanvas("canva_pdfs_seperated_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_seperated_scaledtoqcd_int->cd();
                                hstack_pfds_seperated_scaledtoqcd_int.Draw("Hist E");
                                hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_seperated_scaledtoqcd_int->Modified();
                                canva_pdfs_seperated_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_seperated_scaledtoqcd_int, ibin_dr, ibin_pt);
                                canva_pdfs_seperated_scaledtoqcd_int->Write();
                                canva_pdfs_seperated_scaledtoqcd_int->Print(Form("%s/%s_pdfs_seperated_scaledtoqcd_int_beforefit.png", sDirname.Data(), sname_canvas.Data()));
                                 

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
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true , 0., 1);// , // 1-bkg_fraction_b_true

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, false); // fasle: dont use extended fit: assume fraction not absolute yields  
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

            // std::cout << "a: 2B =" << sig_fraction_true << ", a': after fit=" << p0 << std::endl;            
            // std::cout << "a'/a  for 2B = " << p0/sig_fraction_true << std::endl;


            // std::cout << "b: 1B =" << bkg_fraction_b_true << ", b': 1B =" << p1 << std::endl;
            // std::cout << "b'/b for 1B = " << p1/bkg_fraction_b_true << std::endl;


            // std::cout << "c: 0B =" << (1 - sig_fraction_true - bkg_fraction_b_true) << ", c': 0B after fit =" << p2 << std::endl;
            // std::cout << "c'/c = " << p2/(1-sig_fraction_true - bkg_fraction_b_true) << std::endl;


            // -- updated for the new hsit binning for the true and fitresult S/B fractions (for later drawings)
            // Fit result 
            h_sig_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, p0);
            h_sig_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0);
            h_sig_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP0);

            h_bkg_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, 1- p0);
            h_bkg_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0); // errorp1 = error (1-p0)
            h_bkg_fraction_error->SetBinContent(ibin_dr +1, ibin_pt +1, errP0);


            // -- Compute siganl and bkg fractions uncertainity 
            //save the true fraction
            // -- sigma True for S or B  = (B x sigmaS + SxsigmaB)/(S+B)²
            double err_true_frac = ( (int0 + int1) * True_sig_err + int2 * True_bkg_err )/TMath::Power(int0 + int1+ int2, 2);// 0B + 1B and 2B 
            h_sig_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, sig_fraction_true);
            h_sig_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            
            h_bkg_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, bkg_fraction_true);
            h_bkg_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            h_bkg_frac_true_error->SetBinContent(ibin_dr +1, ibin_pt +1, err_true_frac);
            

            // -- After fits: save mass distribution re and post-fit 
            TH1D *h_sig_fit = (TH1D*) h_sig->Clone(Form("h_sig_fit_%d_%d", ibin_dr, ibin_pt));
                h_sig_fit->Scale(p0 * integral_inputdata);
                h_sig_fit->SetTitle("2B");
                h_sig_fit->GetYaxis()->SetTitle("Counts/[GeV^{2}]");

                h_sig_fit->SetMarkerColor(kRed+2);
                h_sig_fit->SetLineColor(kRed+2);
                h_sig_fit->SetFillColor(kRed+2);
                h_sig_fit->SetFillStyle(3001); // solid
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg 
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("1B+ 0B");
                h_bkg_fit->SetMarkerColor(kCyan+2);
                h_bkg_fit->SetLineColor(kCyan+2);
                h_bkg_fit->SetFillColor(kCyan+2);


            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only 
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("1B");
                h_bkg_fit_1b->SetMarkerColor(kOrange+2);
                h_bkg_fit_1b->SetLineColor(kOrange+2);
                h_bkg_fit_1b->SetLineWidth(2);
                h_bkg_fit_1b->SetFillColor(kOrange+2);

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only 
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("0B");
                h_bkg_fit_nob->SetMarkerColor(kYellow+2);
                h_bkg_fit_nob->SetLineColor(kYellow+2);
                h_bkg_fit_nob->SetFillColor(kYellow+2);
                h_bkg_fit_nob->SetLineWidth(2);

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("Total fit");
                    h_total_fit->SetFillColor (0); // remove color
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
            auto canva_afterfit = new TCanvas(Form("ALLHist_%s", sname_canvas_afterfit.Data()) ,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetMaximum(1.2 * h_data_mb->GetMaximum());
                h_data_mb->SetTitle("Data");
                h_data_mb->SetLineWidth(2);
                h_bb_bjet->SetLineWidth(2);
                h_sig_fit->SetLineWidth(2);
                h_b_bjet->SetLineStyle(9);
                h_bkg_fit->SetLineWidth(2);
                h_data_mb->Draw("P E");
                h_total_fit->Draw("P E SAME");
                h_bkg_fit_1b->Draw("HIST E SAME");
                h_bkg_fit_nob->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
               

        // -- After fit: stacked: seperated contibutions 
            THStack hstack_afterfit("hstack_afterfit","Mass stacked histogram");
                hstack_afterfit.SetTitle(";m_{2B} [GeV];");
                hstack_afterfit.Add(h_sig_fit);
                hstack_afterfit.Add(h_bkg_fit_1b);
                hstack_afterfit.Add(h_bkg_fit_nob);
                hstack_afterfit.SetMaximum(1.2 * hstack_afterfit.GetMaximum());

            auto canva_stack_afterfit = new TCanvas(Form("All_templates_Data_stacked_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_afterfit->cd();
                hstack_afterfit.Draw("hist E"); 
                h_data_mb->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_afterfit->Modified();
                canva_stack_afterfit->Update();

            // normalized version: data and templates afterfit: Data Vs. pdf signal Vs. bkg 
            THStack hstack_norm_afterfit("hstack_norm_afterfit","Mass stacked histogram");
                hstack_norm_afterfit.SetTitle(";m_{2B} [GeV];");    
                TH1D* h_sig_fit_normstack = (TH1D*)  h_sig_fit->Clone("h_sig_fit_normstack");  h_sig_fit_normstack ->Scale(1./integral_inputdata);
                TH1D* h_bkg_fit_normstack = (TH1D*)  h_bkg_fit->Clone("h_bkg_fit_normstack");  h_bkg_fit_normstack ->Scale(1./integral_inputdata);
                hstack_norm_afterfit.Add(h_sig_fit_normstack);
                hstack_norm_afterfit.Add(h_bkg_fit_normstack);
                hstack_norm_afterfit.SetMaximum(1.5 * h_bkg_fit_normstack->GetMaximum());
                auto canva_stack_norm_afterfit = new TCanvas(Form("PDFs_Data_stacked_norm_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_norm_afterfit->cd();
                hstack_norm_afterfit.Draw("hist E"); 
                hnorm_data_self->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_norm_afterfit->Modified();
                canva_stack_norm_afterfit->Update();

            //  -- Build legend and write plots     
                fout->cd();
                DrawCommonTextTopRight(canva_afterfit, ibin_dr, ibin_pt);
                canva_afterfit->Write();
                canva_afterfit->Print(Form("%s/%s.png", sDirname.Data(), canva_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_afterfit,ibin_dr, ibin_pt);
                canva_stack_afterfit->Write();
                canva_stack_afterfit->Print(Form("%s/%s.png", sDirname.Data(), canva_stack_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_norm_afterfit, ibin_dr, ibin_pt);
                canva_stack_norm_afterfit->Write();
                canva_stack_norm_afterfit->Print(Form("%s/%s.png", sDirname.Data(), canva_stack_norm_afterfit->GetName()));

            //-- ratio plot: data /total fit
                TCanvas* c = new TCanvas(Form("RatioPlot_%s", sname_canvas_afterfit.Data()), "", 900, 1100); // 800, 1100
                    TPad* pad1 = new TPad("pad1","",0,0.2,1,1);
                    TPad* pad2 = new TPad("pad2","",0,0,1,0.24);
                    // pad1->SetBottomMargin(0.13);
                    pad1->SetBottomMargin(0.06);
                    pad2->SetTopMargin(0.02);     // bottom pad (very small)
                    pad2->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad1->Draw();
                    pad2->Draw(); 
                    pad1->cd();
                    // Draw frame to control y axis name 
                    TH1F *frame = pad1->DrawFrame(hstack_afterfit.GetXaxis()->GetXmin(), 0, hstack_afterfit.GetXaxis()->GetXmax(), hstack_afterfit.GetMaximum()*1.);
                        frame->GetYaxis()->SetTitle("Counts/[GeV]");
                        hstack_afterfit.Draw("hist E same"); 
                        h_data_mb->Draw("PE same"); 
                        // h_total_fit->Draw("P E SAME");
                        DrawCommonTextTopRight(pad1, ibin_dr, ibin_pt);
                        pad1->Modified(); // force refresh 
                        pad1->Update();
                        // gPad->Modified();
                        // gPad->Update();
                    pad2->cd(); 
                    AddRatioPlot(h_data_mb, h_total_fit);
                    pad2->SetTickx(1);// → draws ticks on both bottom and top


                        fout->cd();
                        c ->Write();
                        c ->Print(Form("%s/%s.png", sDirname.Data(), c->GetName()));


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

    // //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
        
    fout->Close();

}


void template_fit(){

    //Get data and mc labels
    TString pT_selection = "80_140";// full range 
    TString folder = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/TemplateFit_test_combined/";
   
    
    // Add LowEG data 
    bool alsoLowEG = true;

    // -- Output folder to save the result of the tests 
    gSystem->mkdir(sDirname, kTRUE); // make directory if does not exist for the output files to be all in one place 

    TString sfoutputPlots_dijet = "Summary_histo_templatefit_test.root";
    TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "RECREATE"); // For canvas output
    if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    // --  using result of Zoe code and condor result 
    alsoLowEG = true;
    TString dataset_HG = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/data_MC_samples_Zoe_14April/template_for_fit_histos_3D_HighEG_btag.root"; 
    TString dataset_LG = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/data_MC_samples_Zoe_14April/template_for_fit_histos_3D_LowEG_btag.root"; 
        TString dataname = "All";
        TString templates_dijet = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/data_MC_samples_Zoe_14April/template_for_fit_histos_3D_qcd_btag.root";
        TString templates_bjet = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/data_MC_samples_Zoe_14April/template_for_fit_histos_3D_bjet_btag.root";
        TString fout_name = "TemplateFits_histos_3d_" + pT_selection +  ".root";
        do_template_fit_combined (dataset_HG,dataset_LG,templates_dijet, templates_bjet,  pT_selection, folder, fout_name, alsoLowEG);

    // -- next draw S, B farctions true Vs. after fit (draw_template_fit_result)
    // --- 



    foutputPlots_dijet->Print();
    foutputPlots_dijet->Close();
    delete foutputPlots_dijet;
}


