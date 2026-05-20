//Applies 1D and 2D unfolding to data (both have 2D respomse matrices, which unfolding is used depends on the bool matched)
#include "binning_histos_small.h"


void apply_unfolding(TString &dataset, TString &label, TString &folder, bool matched, bool btag, Int_t n, TString pT_selection)
{   
    //Select unfolding options
    bool unfoldBayes =  false;
    bool multiply_sigfrac = true;
    bool correct_tageff = false;
    bool reduced = false;
    TString dataset_response = "bjet";
    TString filename_template_fit = folder + "FRAN_histos_3d_from_templ_" + dataset + "_" + pT_selection + ".root";
    //histos_3d_from_templ_data_80_140.root

    TString fin_name;
    if(matched) fin_name = "histos_response_2D_";
    else fin_name = "histos_response_1D_";
    if(!btag) label += "_notag"; 
    fin_name += TString(Form("n%i_", n)) + dataset_response + "_" + label + ".root";
    TString filename_response = folder + fin_name;
    std::cout << "Using response file: " << filename_response << std::endl;
    
    //Select central pT bin
    int ibin_pt = 2;

    //Print options
    std::cout << "Options:"
              << "\n\tunfoldBayes:" << unfoldBayes
              << "\n\tibin_pt:" << ibin_pt
              << "\n\tmultiply_sigfrac:" << multiply_sigfrac
              << "\n\t\tcorrect_tageff:" << correct_tageff
              << std::endl;


    if(unfoldBayes) label += "_bayesian";
    if(correct_tageff) label += "_tageff_corrected";
    TString fout_name;
    if(matched) fout_name = folder + "histos_" + label + "_after_unfolding_2D.root";
    else fout_name = folder + "histos_" + label + "_after_unfolding_1D.root";

    // ---------------- Plotting setup ------------
    gSystem->Load("libRooUnfold.so");
    gStyle->SetErrorX(0.5);

    Float_t text_size = 20.;
    /*gStyle->SetTextSize(text_size);
    gStyle->SetLegendTextSize(text_size);
    gStyle->SetLabelSize(text_size, "XYZ");
    gStyle->SetTitleSize(text_size, "XYZ");*/
    // --------------------------------------------


    // ----------- Grab data ----------- 

    TString fname_data = filename_template_fit;
    std::cout << "Getting data from " << fname_data << std::endl;
    TFile *fin_data = new TFile(fname_data);
    TH3D *h_data_reco_3D = (TH3D*)fin_data->Get("h_data")->Clone("h_data_reco_3D");
    TH2D *h_data_reco = (TH2D*)h_data_reco_3D->Project3D("zy");

    //Dimension of the response matrix
    int dim = bins_pt*bins_dr;

    int ibin_dr_min = 1;
    int ibin_dr_max = bins_dr;
    

    // Multiply histograms by signal fraction
    TH2D *h_data_after_fit = (TH2D *) h_data_reco->Clone("h_data_after_fit");
    if (multiply_sigfrac) {
        std::cout << "\t---->Multiplying by signal fraction" << std::endl;
        // Grab signal fraction from template fit
        TString fname_fit = filename_template_fit;
        std::cout << "Getting signal fraction from " << fname_fit << std::endl;
        TFile *fin_fit = new TFile(fname_fit);
        TH2D *h_sig_fraction = (TH2D *) fin_fit->Get("h_sig_fraction");    
        h_data_after_fit->Multiply(h_sig_fraction);
    }
    
    // Note: Result = unfold(raw * purity) * 1 / (efficiency)
    //       fakes are negligible

    // ---- Grab response matrix + corrections
    TString fname_unfolding = filename_response;
    std::cout << "Getting response + corrections from : " << fname_unfolding << std::endl;
    TFile *fin_unfolding = new TFile(fname_unfolding);

    TH2D *h_full_purity;
    TH2D *h_full_efficiency;
    TH2D *h_mc_reco;
    RooUnfoldResponse *response;
    TH2D *h_mc_true_no_eff;

    // get purity and efficiency
    h_full_purity = (TH2D *) fin_unfolding->Get("h_full_purity_eecpt"); // reconstruction purity correction
    h_full_efficiency = (TH2D *) fin_unfolding->Get("h_full_efficiency_eecpt"); 
    
    // get MC for comparison and response matrix
    h_mc_reco = (TH2D *) fin_unfolding->Get("h_full_purity_numerator_eecpt"); // reco MC to compare w/ data after purity correction
    response = (RooUnfoldResponse *) fin_unfolding->Get("response_full_eecpt"); // response 
    h_mc_true_no_eff = (TH2D *) fin_unfolding->Get("h_full_efficiency_numerator_eecpt");
   

    // ---- Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    std::cout << "\t---->Condition number nominal = " << cond_number
              << std::endl;

    // ---- Grab the truth level MC ---- 
    TString fname_response_truth = "/data_CMS/cms/zaidan/eec_trees/fran_bins/hist_3d_gen_aggr_n1_MC_bjet_80_140.root";
    std::cout << "Getting truth from : " << fname_response_truth << std::endl;
    TFile *fin_response_truth = new TFile(fname_response_truth);
    TH3D *h_mc_true_3D = (TH3D*)fin_response_truth->Get("h3D_b1"); //ATTENTION 
    TH2D *h_mc_true = (TH2D *)h_mc_true_3D->Project3D("zy")->Clone("h_mc_true");
    //TH2D *h_mc_true = (TH2D *) fin_response_truth->Get("h_full_efficiency_denominator_eecpt"); // true MC to compare w/ data after BOTH efficiency corrections
    

    //------- Apply purity correction
    std::cout << "\t---->Multiplying data by purity" << std::endl;
    TH2D *h_data_purity_corrected = (TH2D *) h_data_after_fit->Clone("h_data_purity_corrected");
    h_data_purity_corrected->Multiply(h_full_purity);

    // ---- Unfold
    // ACA!!
    std::cout << "\t---->Unfolding" << std::endl;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    TH2D *h_data_unfolded;
    TMatrixD covariance_matrix_before_unfolding(dim,dim);
    TMatrixD covariance_matrix_after_unfolding(dim,dim);
    if (unfoldBayes) {
        Int_t niter =  100;
        RooUnfoldBayes unfold(response, h_data_purity_corrected, niter);
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment);
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    } else {
        RooUnfoldInvert unfold(response, h_data_purity_corrected);
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment);
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    }    
    h_data_unfolded->SetName("h_data_unfolded");

    // ---- Fold back
    std::cout << "\t---->Refolding" << std::endl;
    TH2D *h_data_refolded = (TH2D *) response->ApplyToTruth(h_data_unfolded, "h_data_refolded");

    // ---- Apply efficiency correction 
    std::cout << "\t---->Dividing by recostruction efficiency" << std::endl;
    TH2D *h_data_efficiency_corrected = (TH2D *) h_data_unfolded->Clone("h_data_efficiency_corrected");
    h_data_efficiency_corrected->Divide(h_full_efficiency);

    // ---- Final corrections
    TH2D *h_data_fully_corrected = (TH2D *) h_data_efficiency_corrected->Clone("h_data_fully_corrected");
    TH2D *h_eff;
    if (correct_tageff) {
        // ---- Grab b tagging efficiency correction 
      TString fname_b_tag_eff = folder + "file_efficiency_" + dataset + "_" + pT_selection + ".root";
        std::cout << "Getting b tagging efficiency from: " << fname_b_tag_eff << std::endl;
        TFile *fin_b_tag_eff = new TFile(fname_b_tag_eff);
        h_eff = (TH2D *) fin_b_tag_eff->Get("h_eff");
        std::cout << "\t---->Dividing by b tagging efficiency" << std::endl;
        h_data_fully_corrected->Divide(h_eff);
        //h_mc_true->Divide(h_eff);
    } 

    // ---- Graphical bottomline test

    std::cout << "Performing graphical bottomline test" << std::endl;
    TH1D *h_mc_reco_1d = (TH1D *) h_mc_reco->ProjectionX("h_mc_reco_1d", ibin_pt, ibin_pt);
    TH1D *h_mc_true_1d = (TH1D *) h_mc_true->ProjectionX("h_mc_true_1d", ibin_pt, ibin_pt);
    TH1D *h_data_purity_corrected_1d = (TH1D *) h_data_purity_corrected->ProjectionX("h_data_purity_corrected_1d", ibin_pt, ibin_pt);
    TH1D *h_data_fully_corrected_1d = (TH1D *) h_data_fully_corrected->ProjectionX("h_data_fully_corrected_1d", ibin_pt, ibin_pt);
    TH1D *h_data_refolded_1d = (TH1D *) h_data_refolded->ProjectionX("h_data_refolded_1d", ibin_pt, ibin_pt);
    TH1D* h_data_after_fit_1d = (TH1D*) h_data_after_fit->ProjectionX("h_data_after_fit_1d", ibin_pt, ibin_pt);
    TH1D* h_data_unfolded_1d = (TH1D*) h_data_unfolded->ProjectionX("h_data_unfolded_1d", ibin_pt, ibin_pt);
    
    if (true) {
        
        double ymax = 0.;
        for (auto h : {
                    h_mc_reco_1d, 
                    h_mc_true_1d,
                    h_data_purity_corrected_1d,
                    h_data_fully_corrected_1d,
                    h_data_refolded_1d,
                    h_data_after_fit_1d,
                    h_data_unfolded_1d
                    }) {
                        h->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
                        h->Scale(1/h->Integral(), "width");
                        if (h->GetMaximum()>ymax) ymax =  h->GetMaximum();
                    }
        
        Float_t pt_min_plot = 100;
        Float_t pt_max_plot = 120;

        TLegend *leg = new TLegend(0.35, 0.55, 0.65, 0.85);
        if (false) leg = new TLegend(0.2, 0.4, 0.6, 0.8);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetMargin(0.15);
        leg->SetHeader(Form("%.0f < p_{T}^{jet} < %.0f GeV", pt_min_plot, pt_max_plot));

        if(correct_tageff) h_data_purity_corrected_1d->SetTitle("Data, " + label + " " + dataset_response + " response matrix, b-tag efficiency corrected");
        else h_data_purity_corrected_1d->SetTitle("");//Data, " + label + " " + dataset_response + " response matrix");
        h_data_purity_corrected_1d->SetStats(0);
        h_data_purity_corrected_1d->SetMarkerColor(kBlack);
        h_data_purity_corrected_1d->SetMarkerStyle(kFullCircle);
        h_data_purity_corrected_1d->SetMarkerSize(1);
        h_data_purity_corrected_1d->GetYaxis()->SetRangeUser(0., ymax*1.1);
        h_data_purity_corrected_1d->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
        leg->AddEntry(h_data_purity_corrected_1d, "Detector level data (tagged single B jets)", "pe1");
        //leg->AddEntry(h_data_purity_corrected_1d, "before unfolding", "pe1");
        

        /*h_data_after_fit_1d->SetMarkerColor(kOrange);
        h_data_after_fit_1d->SetLineColor(kOrange);
        h_data_after_fit_1d->SetMarkerStyle(kFullTriangleUp);
        h_data_after_fit_1d->SetMarkerSize(1);
        leg->AddEntry(h_data_after_fit_1d, "Before purity correction", "pe1");

        h_data_unfolded_1d->SetMarkerColor(kGreen);
        h_data_unfolded_1d->SetLineColor(kGreen);
        h_data_unfolded_1d->SetMarkerStyle(kFullTriangleUp);
        h_data_unfolded_1d->SetMarkerSize(1);
        leg->AddEntry(h_data_unfolded_1d, "After unfolding (no efficiency correction)", "pe1");*/

        h_mc_reco_1d->SetStats(0);
        h_mc_reco_1d->SetMarkerColor(kRed);
        h_mc_reco_1d->SetLineColor(kRed);
        h_mc_reco_1d->SetMarkerStyle(kFullTriangleUp);
        h_mc_reco_1d->SetMarkerSize(1);
        h_mc_reco_1d->GetYaxis()->SetRangeUser(0., ymax*1.1);
        h_mc_reco_1d->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
        leg->AddEntry(h_mc_reco_1d, "Detector level MC (tagged single B jets)", "pe1");

        h_data_fully_corrected_1d->SetMarkerColor(kBlue);
        h_data_fully_corrected_1d->SetLineColor(kBlue);
        h_data_fully_corrected_1d->SetMarkerStyle(kOpenCross);
        h_data_fully_corrected_1d->SetMarkerSize(1);
        if(correct_tageff) leg->AddEntry(h_data_fully_corrected_1d, "Unfolded data (single B jets)", "pe1");
        else leg->AddEntry(h_data_fully_corrected_1d, "Unfolded data (tagged single B jets)", "pe1");

        h_mc_true_1d->SetMarkerColor(kRed);
        h_mc_true_1d->SetLineColor(kRed);
        h_mc_true_1d->SetMarkerStyle(kOpenTriangleUp);
        h_mc_true_1d->SetMarkerSize(1);
        leg->AddEntry(h_mc_true_1d, "Particle level MC (tagged single B jets)", "pe1");

        h_data_refolded_1d->SetMarkerColor(kBlue);
        h_data_refolded_1d->SetLineColor(kBlue);
        h_data_refolded_1d->SetMarkerStyle(kFullCross);
        h_data_refolded_1d->SetMarkerSize(1);
        //if(!reduced) leg->AddEntry(h_data_refolded_1d, "Refolded data", "pe1");

        TCanvas *c_unfold = new TCanvas("c_unfold", "", 800, 600);
        TPad *pad1 = new TPad("pad1", "", 0., 0., 1., 0.3);
        TPad *pad2 = new TPad("pad2", "", 0., 0.3, 1., 1.);
        pad1->SetTopMargin(0.01);
        pad1->SetBottomMargin(0.3);
        pad2->SetBottomMargin(0.01);
        pad1->SetLogx();
        pad2->SetLogx();

        pad2->cd();

        h_data_purity_corrected_1d->Draw("pe1");
        h_data_purity_corrected_1d->SetMaximum(9);
        //h_data_unfolded_1d->Draw("pe1 same");
        h_mc_reco_1d->Draw("pe1 same");
        h_data_fully_corrected_1d->Draw("pe1 same");
        h_mc_true_1d->Draw("pe1 same");
        //if(!reduced) h_data_refolded_1d->Draw("pe1 same");
        leg->Draw();
        TLatex *test_info_text = new TLatex;
        test_info_text->SetNDC();
        test_info_text->SetTextSize(0.03);
        test_info_text->DrawLatex(0.69, 0.7, "Data HighEGJet + LowEGJet");
        test_info_text->DrawLatex(0.69, 0.75, "MC bjet response matrix");
        if(matched) test_info_text->DrawLatex(0.69, 0.65, "2D jet p_{T} and #Deltar unfolding");
        else test_info_text->DrawLatex(0.69, 0.65, "1D jet p_{T} unfolding");
        test_info_text->Draw("same");
        //drawHeader();    

        TLine *line = new TLine(dr_min, 1., dr_max, 1.);
        line->SetLineWidth(2.); 
        line->SetLineStyle(kDashed);
        line->SetLineColor(kGray);

        TLegend *leg_ratio = new TLegend(0.7, 0.8, 0.9, 0.99);//0.5, 0.3, 0.85, 0.5);
        leg_ratio->SetBorderSize(1);

        TH1D *h_data_mc_reco_ratio = (TH1D *) h_data_purity_corrected_1d->Clone("h_data_mc_reco_ratio");
        h_data_mc_reco_ratio->SetTitle("");
        h_data_mc_reco_ratio->SetStats(0);
        h_data_mc_reco_ratio->Divide(h_mc_reco_1d);
        h_data_mc_reco_ratio->SetMarkerStyle(kFullCircle);
        h_data_mc_reco_ratio->SetMarkerColor(kBlack);
        h_data_mc_reco_ratio->SetLineColor(kBlack);
        h_data_mc_reco_ratio->SetMarkerSize(1);
        h_data_mc_reco_ratio->GetYaxis()->SetRangeUser(0.,2.);
        h_data_mc_reco_ratio->GetYaxis()->SetTitle("ratio");
        h_data_mc_reco_ratio->GetYaxis()->SetLabelSize(0.08);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleSize(0.08);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleOffset(1.5);
        h_data_mc_reco_ratio->GetXaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
        h_data_mc_reco_ratio->GetXaxis()->SetLabelSize(0.1);
        //h_data_mc_reco_ratio->GetXaxis()->SetTitleSize(text_size);
        //h_data_mc_reco_ratio->GetXaxis()->SetTitleOffset(3.5);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleSize(0.08);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleOffset(0.4);
        h_data_mc_reco_ratio->GetYaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetYaxis()->SetNdivisions(8);
        leg_ratio->AddEntry(h_data_mc_reco_ratio, "reco data / reco mc", "pe1");

        TH1D *h_data_mc_true_ratio = (TH1D *) h_data_fully_corrected_1d->Clone("h_data_mc_true_ratio");
        h_data_mc_true_ratio->Divide(h_mc_true_1d);
        h_data_mc_true_ratio->SetMarkerStyle(kOpenCross);
        h_data_mc_true_ratio->SetMarkerColor(kBlue);
        h_data_mc_true_ratio->SetLineColor(kBlue);
        h_data_mc_true_ratio->SetMarkerSize(1);
        leg_ratio->AddEntry(h_data_mc_true_ratio, "unfolded data / true mc", "pe1");

        pad1->cd();
        h_data_mc_reco_ratio->Draw("pe1 same");
        h_data_mc_true_ratio->Draw("pe1 same");
        leg_ratio->Draw();
        line->Draw();

        c_unfold->cd();
        pad1->Draw();
        pad2->Draw();
        c_unfold->Draw();
        if(matched) dataset += "_matched";
        if(reduced) dataset += "_reduced";
        if(unfoldBayes) dataset += "_bayesian";
        if(matched) c_unfold->Print(folder + "unfolding_plot"+dataset+"_bottomline_test_eec_2D.pdf");
        else c_unfold->Print(folder + "unfolding_plot"+dataset+"_bottomline_test_eec_1D.pdf");
        

    }

    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    h_data_fully_corrected->SetName("h_data_unfolded");
    h_data_fully_corrected->Write();

    h_data_after_fit->SetName("h_data_singleb");
    h_data_after_fit->Write();

    h_mc_reco->SetName("h_mc_reco_singleb");
    h_mc_reco->Write();

    //Write 1D histograms for plotting
    h_mc_reco_1d->Write();
    h_mc_true_1d->Write();
    h_data_purity_corrected_1d->Write();
    h_data_fully_corrected_1d->Write();
    h_data_refolded_1d->Write();
    
   


    fout->Close();
    delete fout;

    // gApplication -> Terminate(0);
}

void apply_unfolding_2d(){
    std::vector<TString> datasets{"data"};   //{"data"};//, "dijet"};
    
    TString folder = "/data_CMS/cms/zaidan/eec_trees/fran_bins/";

    TString pT_selection = "80_140";

    bool btag = true;

    bool matched = false;

    Int_t n = 1;

    //Create labels
    std::vector<TString> labels_vec{"b1"};//, "moreb", "other","mc"};

    for(Int_t i = 0; i < datasets.size(); i++){
        for(Int_t j = 0; j < labels_vec.size(); j++){
	  apply_unfolding(datasets.at(i), labels_vec.at(j), folder, matched, btag, n, pT_selection);
        }

    }
}
