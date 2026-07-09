//Applies 2D and 2D unfolding to data (both have 2D respomse matrices, which unfolding is used depends on the bool matched)
#include "binning_histos_small.h"

namespace ROCColor {
  Color_t blue()   { return TColor::GetColor("#4E79A7"); }
  Color_t red()    { return TColor::GetColor("#E15759"); }
  Color_t green()  { return TColor::GetColor("#59A14F"); }
  Color_t purple() { return TColor::GetColor("#B07AA1"); }
  Color_t orange() { return TColor::GetColor("#F28E2B"); }
}
void apply_unfolding(TString &label, TString &folder, bool btag, Int_t n, TString pT_selection)
{   
    //Select unfolding options
    bool unfoldBayes =  true;
    bool multiply_sigfrac = true;
    bool split_test = true;
    
    const Color_t blue = ROCColor::blue();
    const Color_t red = ROCColor::red();
    const Color_t green = ROCColor::green();
    const Color_t purple = ROCColor::purple();
    const Color_t orange = ROCColor::orange();

    TString filename_template_fit = "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/TemplateFit_Run3/TemplateFits_Run3_minHLT60_LinearBin/nominal_Run3_TemplateFits_histos_3d_80_200.root"; 
    std::cout << "Using template file: " << filename_template_fit << std::endl;
    
    TString filename_response = "/data_CMS/cms/zaidan/analysis_lise/Run3/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_split_test_f.root";  //RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root"; //folder + "RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
    std::cout << "Using response file: " << filename_response << std::endl;
    
    TString filename_data =  "/data_CMS/cms/zaidan/analysis_lise/Run3/Run3_btagWP868_template_for_fit_histos_3D_qcd_split_test.root"; //Run3_btagWP868_template_for_fit_histos_3D_data_f.root";
    std::cout << "Getting data from " << filename_data << std::endl;
    
    //Select central pT bin
    int ibin_pt = 2;
    //Print options
    std::cout << "Options:"
              << "\n\tunfoldBayes:" << unfoldBayes
              << "\n\tibin_pt:" << ibin_pt
              << "\n\tmultiply_sigfrac:" << multiply_sigfrac
              << std::endl;


   
    if(unfoldBayes) label += "_bayesian";
    TString fout_name;
    fout_name = folder + "histos_" + label + "_after_unfolding_2D_split_test.root";

    // ---------------- Plotting setup ------------
    gSystem->Load("libRooUnfold.so");
    gStyle->SetErrorX(0.5);
    Float_t text_size = 20.;

    // ----------- Grab data ----------- 
    TFile *fin_data = new TFile(filename_data);
    TString histname = "h3D_data";
    if (split_test) histname = "h3D_pseudodata";
    TH3D *h_data_reco_3D = (TH3D*)fin_data->Get(histname)->Clone("h_data_reco_3D"); //att
    TH2D *h_data_reco = (TH2D*)h_data_reco_3D->Project3D("zy");
    TH2D *h_data_after_fit = (TH2D*) h_data_reco->Clone("h_data_after_fit");
        
    

    // Multiply histograms by signal fraction
    
    if (multiply_sigfrac) {
        std::cout << "\t---->Multiplying by signal fraction" << std::endl;
        // Grab signal fraction from template fit
        TString fname_fit = filename_template_fit;
        std::cout << "Getting signal fraction from " << fname_fit << std::endl;
        TFile *fin_fit = new TFile(fname_fit);
        TH2D *h_sig_fraction = (TH2D *) fin_fit->Get("h_sig_fraction_fit");    
        h_data_after_fit->Multiply(h_sig_fraction);
    
    }
    else {
        std::cout << "\t---->Not multiplying by signal fraction" << std::endl;
    }

    

     //Dimension of the response matrix
    int dim = bins_pt*bins_dr;
    int ibin_dr_min = 1;
    int ibin_dr_max = bins_dr;
    
    std::cout << "dim = " << dim << std::endl;
    std::cout << "bins_pt = " << bins_pt << std::endl;
    std::cout << "bins_dr = " << bins_dr << std::endl;

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

     if (split_test){
        std::cout << "\t----> Doing split test" << std::endl;
        // get purity and efficiency
        h_full_purity = (TH2D *) fin_unfolding->Get("h_full_pseudo_purity_tf"); // reconstruction purity correction
        h_full_efficiency = (TH2D *) fin_unfolding->Get("h_full_pseudo_efficiency_tf"); 
        
        // get MC for comparison and response matrix
        h_mc_reco = (TH2D *) fin_unfolding->Get("h_full_pseudo_numerator_tf"); // reco MC to compare w/ data after purity correction
        response = (RooUnfoldResponse *) fin_unfolding->Get("response_tf_pseudo_full"); // response 
        h_mc_true_no_eff = (TH2D *) fin_unfolding->Get("h_full_pseudo_efficiency_numerator_tf");
     }

    else {
        // get purity and efficiency
        h_full_purity = (TH2D *) fin_unfolding->Get("h_full_purity_tf"); // reconstruction purity correction
        h_full_efficiency = (TH2D *) fin_unfolding->Get("h_full_efficiency_tf"); 
        
        // get MC for comparison and response matrix
        h_mc_reco = (TH2D *) fin_unfolding->Get("h_full_purity_numerator_tf"); // reco MC to compare w/ data after purity correction
        response = (RooUnfoldResponse *) fin_unfolding->Get("response_tf_full"); // response 
        h_mc_true_no_eff = (TH2D *) fin_unfolding->Get("h_full_efficiency_numerator_tf");
    }

    // ---- Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    std::cout << "\t---->Condition number nominal = " << cond_number
              << std::endl;

    // ---- Grab the truth level MC ---- 
    TString fname_response_truth =  "/data_CMS/cms/zaidan/analysis_lise/Run3/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
    std::cout << "Getting truth from : " << fname_response_truth << std::endl;
    TFile *fin_response_truth = new TFile(fname_response_truth);
    TH2D *h_mc_true = (TH2D*)fin_response_truth->Get("hgenjet_2b_passbtag"); //ATTENTION 

    //------- Apply purity correction
    std::cout << "\t---->Multiplying data by purity" << std::endl;
    TH2D *h_data_purity_corrected = (TH2D *) h_data_after_fit->Clone("h_data_purity_corrected");
    h_data_purity_corrected->Multiply(h_full_purity);

    //double underflowY = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1, 0, 0);
    //double overflowY  = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1, 
                          //h_data_purity_corrected->GetNbinsY()+1, h_data_purity_corrected->GetNbinsY()+1);
    //std::cout << "underflowY=" << underflowY << " overflowY=" << overflowY << std::endl;
    
    // ---- Unfold
    std::cout << "\t---->Unfolding" << std::endl;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    TH2D *h_data_unfolded;
    TMatrixD covariance_matrix_before_unfolding(dim,dim);
    TMatrixD covariance_matrix_after_unfolding(dim,dim);
    if (unfoldBayes) {
        Int_t niter = 4;
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


    // ---- Graphical bottomline test

    std::cout << "Performing graphical bottomline test" << std::endl;
    TH2D *h_mc_reco_2D = (TH2D *) h_mc_reco->ProjectionX("h_mc_reco_2D", ibin_pt, ibin_pt);
    TH2D *h_mc_true_2D = (TH2D *) h_mc_true->ProjectionX("h_mc_true_2D", ibin_pt, ibin_pt);
    TH2D *h_data_purity_corrected_2D = (TH2D *) h_data_purity_corrected->ProjectionX("h_data_purity_corrected_2D", ibin_pt, ibin_pt);
    TH2D *h_data_fully_corrected_2D = (TH2D *) h_data_fully_corrected->ProjectionX("h_data_fully_corrected_2D", ibin_pt, ibin_pt);
    TH2D *h_data_refolded_2D = (TH2D *) h_data_refolded->ProjectionX("h_data_refolded_2D", ibin_pt, ibin_pt);
    TH2D* h_data_after_fit_2D = (TH2D*) h_data_after_fit->ProjectionX("h_data_after_fit_2D", ibin_pt, ibin_pt);
    TH2D* h_data_unfolded_2D = (TH2D*) h_data_unfolded->ProjectionX("h_data_unfolded_2D", ibin_pt, ibin_pt);
    
    if (true) {
        
        double ymax = 0.;
        for (auto h : {
                    h_mc_reco_2D, 
                    h_mc_true_2D,
                    h_data_purity_corrected_2D,
                    h_data_fully_corrected_2D,
                    h_data_refolded_2D,
                    h_data_after_fit_2D,
                    h_data_unfolded_2D
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

        h_data_purity_corrected_2D->SetTitle("");//Data, " + label + " " + " response matrix");
        h_data_purity_corrected_2D->SetStats(0);
        h_data_purity_corrected_2D->SetMarkerColor(purple);
        h_data_purity_corrected_2D->SetMarkerStyle(kFullCircle);
        h_data_purity_corrected_2D->SetMarkerSize(1);
        h_data_purity_corrected_2D->GetYaxis()->SetRangeUser(0., ymax*1.1);
        h_data_purity_corrected_2D->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
        leg->AddEntry(h_data_purity_corrected_2D, "Detector level data (tagged BB jets)", "pe1");
    

        //h_data_after_fit_2D->SetMarkerColor(orange);
        //h_data_after_fit_2D->SetLineColor(orange);
        //h_data_after_fit_2D->SetMarkerStyle(kFullTriangleUp);
        //h_data_after_fit_2D->SetMarkerSize(1);
        //leg->AddEntry(h_data_after_fit_2D, "Before purity correction", "pe1");

        //h_data_unfolded_2D->SetMarkerColor(green);
        //h_data_unfolded_2D->SetLineColor(green);
        //h_data_unfolded_2D->SetMarkerStyle(kFullTriangleUp);
        //h_data_unfolded_2D->SetMarkerSize(1);
        //leg->AddEntry(h_data_unfolded_2D, "After unfolding (no efficiency correction)", "pe1");

        h_mc_reco_2D->SetStats(0);
        h_mc_reco_2D->SetMarkerColor(red);
        h_mc_reco_2D->SetLineColor(red);
        h_mc_reco_2D->SetMarkerStyle(kFullTriangleUp);
        h_mc_reco_2D->SetMarkerSize(1);
        h_mc_reco_2D->GetYaxis()->SetRangeUser(0., ymax*1.1);
        h_mc_reco_2D->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
        leg->AddEntry(h_mc_reco_2D, "Detector level MC (tagged BB jets)", "pe1");

        h_data_fully_corrected_2D->SetMarkerColor(blue);
        h_data_fully_corrected_2D->SetLineColor(blue);
        h_data_fully_corrected_2D->SetMarkerStyle(kOpenCross);
        h_data_fully_corrected_2D->SetMarkerSize(1);
        leg->AddEntry(h_data_fully_corrected_2D, "Unfolded data (tagged BB jets)", "pe1");

        h_mc_true_2D->SetMarkerColor(orange);
        h_mc_true_2D->SetLineColor(orange);
        h_mc_true_2D->SetMarkerStyle(kOpenTriangleUp);
        h_mc_true_2D->SetMarkerSize(1);
        leg->AddEntry(h_mc_true_2D, "Particle level MC (tagged BB jets)", "pe1");

        h_data_refolded_2D->SetMarkerColor(blue);
        h_data_refolded_2D->SetLineColor(blue);
        h_data_refolded_2D->SetMarkerStyle(kFullCross);
        h_data_refolded_2D->SetMarkerSize(1);
        
        
        TCanvas *c_unfold = new TCanvas("c_unfold", "", 800, 800);
        TPad *pad1 = new TPad("pad1", "", 0., 0., 1., 0.3);
        TPad *pad2 = new TPad("pad2", "", 0., 0.3, 1., 1.);
        pad1->SetTopMargin(0.01);
        pad1->SetBottomMargin(0.3);
        pad2->SetBottomMargin(0.01);
        //pad1->SetLogx();
        //pad2->SetLogx();

        pad2->cd();

        h_data_purity_corrected_2D->Draw("pe1");
        h_data_purity_corrected_2D->SetMaximum(9);
        //h_data_unfolded_2D->Draw("pe1 same");
        h_mc_reco_2D->Draw("pe1 same");
        h_data_fully_corrected_2D->Draw("pe1 same");
        //h_data_after_fit_2D->Draw("pe1 same");
        h_mc_true_2D->Draw("pe1 same");

        leg->Draw();
        TLatex *test_info_text = new TLatex;
        test_info_text->SetNDC();
        test_info_text->SetTextSize(0.03);
        test_info_text->DrawLatex(0.69, 0.7, "Data Run 3");
        test_info_text->DrawLatex(0.69, 0.75, "MC bjet response matrix");
        test_info_text->DrawLatex(0.69, 0.65, "2D jet p_{T} unfolding");
        test_info_text->Draw("same");
        //drawHeader();    

        TLine *line = new TLine(dr_min, 1., dr_max, 1.);
        line->SetLineWidth(2.); 
        line->SetLineStyle(kDashed);
        line->SetLineColor(kGray);

        TLegend *leg_ratio = new TLegend(0.7, 0.8, 0.9, 0.99);//0.5, 0.3, 0.85, 0.5);
        leg_ratio->SetBorderSize(1);

        TH2D *h_data_mc_reco_ratio = (TH2D *) h_data_purity_corrected_2D->Clone("h_data_mc_reco_ratio");
        h_data_mc_reco_ratio->SetTitle("");
        h_data_mc_reco_ratio->SetStats(0);
        h_data_mc_reco_ratio->Divide(h_mc_reco_2D);
        h_data_mc_reco_ratio->SetMarkerStyle(kFullCircle);
        h_data_mc_reco_ratio->SetMarkerColor(purple);
        h_data_mc_reco_ratio->SetLineColor(purple);
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

        TH2D *h_data_mc_true_ratio = (TH2D *) h_data_fully_corrected_2D->Clone("h_data_mc_true_ratio");
        h_data_mc_true_ratio->Divide(h_mc_true_2D);
        h_data_mc_true_ratio->SetMarkerStyle(kOpenCross);
        h_data_mc_true_ratio->SetMarkerColor(blue);
        h_data_mc_true_ratio->SetLineColor(blue);
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
        label += "_reduced??";
        if(unfoldBayes) label += "_bayesian";
        c_unfold->Print(folder + "unfolding_plot"+label+"_bottomline_test_eec_2D.pdf");
        

    }

    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    h_data_fully_corrected->SetName("h_data_unfolded");
    h_data_fully_corrected->Write();
    h_data_after_fit->SetName("h_data_singleb");
    h_data_after_fit->Write();
    h_mc_reco->SetName("h_mc_reco_singleb");
    h_mc_reco->Write();

    //Write 2D histograms for plotting
    h_mc_reco_2D->Write();
    //h_mc_true_2D->Write();
    h_data_purity_corrected_2D->Write();
    h_data_fully_corrected_2D->Write();
    h_data_refolded_2D->Write();
    
    fout->Close();
    delete fout;

    // gApplication -> Terminate(0);
}

void apply_unfolding_2d(){    
    TString dataset = "data";  
    TString folder = "/data_CMS/cms/zaidan/analysis_lise/Run3/";
    TString pT_selection = "80_200";
    bool btag = true;
    Int_t n = 1;
	apply_unfolding(dataset, folder, btag, n, pT_selection);
     
}
