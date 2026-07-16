#include "binning_histos_small.h"
#include <algorithm>
#include <vector>

namespace ROCColor {
  Color_t blue()   { return TColor::GetColor("#4C72B0"); }
  Color_t red()    { return TColor::GetColor("#C44E52"); }
  Color_t green()  { return TColor::GetColor("#4F8F52"); }
  Color_t purple() { return TColor::GetColor("#8C6BB1"); }
  Color_t orange() { return TColor::GetColor("#C97430"); }
  // Neutral highlight for the ratio-pad gap bands: a hue used by none of the five curves above,
  // so a filled region never reads as one of the plotted quantities. Drawn at two alphas.
  Color_t teal()   { return TColor::GetColor("#1F8A8A"); }
}

template <typename T>
T *getOrWarn(TFile *f, const TString &name)
{
    T *obj = dynamic_cast<T *>(f->Get(name.Data()));
    if (!obj)
        std::cerr << "ERROR: '" << name << "' (" << T::Class_Name() << ") not found in "
                  << f->GetName() << std::endl;
    return obj;
}

TFile *openOrWarn(const TString &name)
{
    TFile *f = TFile::Open(name);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << name << std::endl;
        return nullptr;
    }
    return f;
}

bool normalizeToUnitArea(TH1D *h)
{
    if (!h) return false;
    const Double_t integral = h->Integral();
    if (integral <= 0.) {
        std::cerr << "WARNING: " << h->GetName() << " has integral " << integral
                  << ", skipping normalisation" << std::endl;
        return false;
    }
    h->Scale(1. / integral, "width");
    return true;
}


void apply_unfolding(TString &label, TString &folder, bool btag, Int_t n, TString pT_selection,
                     int test_mode, bool unfoldBayes)
{
    //Select unfolding options. test_mode and unfoldBayes are passed in as arguments:
    //   test_mode 0 = FULL-MC closure : input h3D_bb (full sample), full-sample response/purity/eff,
    //                         truth = h_full_efficiency_denominator_tf. Same events in and out
    //                         -> a perfect/technical closure (NOT statistically independent).
    //   test_mode 1 = SPLIT test      : input h3D_pseudodata_bb (even half), pseudo (odd-half) corrections,
    //                         truth = h_pseudodata_truth_tf (even half) -> real, independent closure.
    //   test_mode 2 = DATA            : input h3D_data + template-fit signal fraction, full-sample corrections,
    //                         truth = hgenjet_2b_passbtag (reference).
    //   unfoldBayes: true = Bayesian unfolding, false = matrix inversion.
    const bool split_test       = (test_mode == 1);  // -> pseudo (odd-half) corrections + even-half truth
    const bool is_data          = (test_mode == 2);
    const bool multiply_sigfrac = is_data;           // only real data needs the bb template fit
    // Matching purity = (reco_pass && gen_pass) / (reco_pass), both already restricted to
    // true-2b jets. It is NOT a bb-purity, so h3D_pseudodata_bb still needs it: that input
    // contains jets whose gen jet fails gen_pass, and the response has no truth row for them.
    bool apply_purity = true;
    // scan_niter = true  -> run niter = 1..100, save the multi-page PDF + one PNG per iteration.
    // scan_niter = false -> single unfolding with niter = 4, save the usual single bottomline plot.
    bool scan_niter = true;


    const Color_t blue = ROCColor::blue();
    const Color_t red = ROCColor::red();
    const Color_t green = ROCColor::green();
    const Color_t purple = ROCColor::purple();
    const Color_t orange = ROCColor::orange();
    const Color_t teal = ROCColor::teal();

    TString filename_template_fit = "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/TemplateFit_Run3/TemplateFits_Run3_minHLT60_LinearBin/nominal_Run3_TemplateFits_histos_3d_80_200.root";
    std::cout << "Using template file: " << filename_template_fit << std::endl;

    TString filename_response = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";  //RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root"; //folder + "RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
    std::cout << "Using response file: " << filename_response << std::endl;

    // Modes 0/1 unfold MC (h3D_bb / h3D_pseudodata_bb), which live in the MCGEN file.
    // Mode 2 unfolds real data from the data file.
    TString filename_data = is_data
        ? "/data_CMS/cms/zaidan/analysis_lise/Run3/Run3_btagWP868_template_for_fit_histos_3D_data_f.root"
        : "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP868_template_for_fit_histos_3D_qcd_fMCGEN.root";
    std::cout << "Getting data from " << filename_data << std::endl; //
    //Select central pT bin
    int ibin_pt = 2;
    //Print options
    std::cout << "Options:"
              << "\n\tunfoldBayes:" << unfoldBayes
              << "\n\tibin_pt:" << ibin_pt
              << "\n\ttest_mode:" << test_mode
              << "\n\tmultiply_sigfrac:" << multiply_sigfrac
              << "\n\tsplit_test:" << split_test
              << "\n\tapply_purity:" << apply_purity
              << std::endl;



    label += unfoldBayes ? "_bayesian" : "_MI";
    if(test_mode==0)  label += "_full_closure";
    if(test_mode==1)  label += "_split_test";
    if(test_mode==2)  label += "_data";
    if(!apply_purity) label += "_nopurity";
    TString fout_name;
    fout_name = folder + "histos_" + label + "_after_unfolding_2D.root";

    // ---------------- Plotting setup ------------
    gSystem->Load("libRooUnfold.so");
    gStyle->SetErrorX(0.5);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    TGaxis::SetMaxDigits(3);   // "x10^-3" on the axis instead of 0.000123 on every label

    const Float_t  font_scale = 1200. / 800.;
    const Style_t  font_code  = 43;
    const Float_t  label_size = 15. * font_scale;
    const Float_t  title_size = 15. * font_scale;
    const Float_t  legend_size = 15. * font_scale;
    // One offset for every axis title, so they all sit the same distance from their numbers.
    const Float_t  title_offset = 1.0;

    // ---- Grab response matrix + corrections
    TString fname_unfolding = filename_response;
    std::cout << "Getting response + corrections from : " << fname_unfolding << std::endl;
    TFile *fin_unfolding = openOrWarn(fname_unfolding);
    if (!fin_unfolding) return;

    // ----------- Grab data -----------
    TFile *fin_data = openOrWarn(filename_data);
    if (!fin_data) return;
    TString histname = (test_mode == 0) ? "h3D_bb"
                     : (test_mode == 1) ? "h3D_pseudodata_bb"
                     :                     "h3D_data";
    std::cout << "Using input histogram: " << histname << std::endl;
    TH3D *h_data_reco_3D_in = getOrWarn<TH3D>(fin_data, histname);
    if (!h_data_reco_3D_in) {
        std::cerr << "       re-run create_files_for_template_fit.cpp to produce it." << std::endl;
        return;
    }
    TH3D *h_data_reco_3D = (TH3D*) h_data_reco_3D_in->Clone("h_data_reco_3D"); //att
    TH2D *h_data_reco = (TH2D*)h_data_reco_3D->Project3D("zy");
    TH2D *h_data_after_fit = (TH2D*) h_data_reco->Clone("h_data_after_fit");



    // Multiply histograms by signal fraction

    if (multiply_sigfrac) {
        std::cout << "\t---->Multiplying by signal fraction" << std::endl;
        // Grab signal fraction from template fit
        TString fname_fit = filename_template_fit;
        std::cout << "Getting signal fraction from " << fname_fit << std::endl;
        TFile *fin_fit = openOrWarn(fname_fit);
        if (!fin_fit) return;
        TH2D *h_sig_fraction = getOrWarn<TH2D>(fin_fit, "h_sig_fraction_fit");
        if (!h_sig_fraction) return;
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

    


    TH2D *h_full_purity = nullptr;
    TH2D *h_full_efficiency = nullptr;
    TH2D *h_mc_reco = nullptr;
    RooUnfoldResponse *response = nullptr;
    TH2D *h_mc_true_no_eff = nullptr;

    // The reco-level MC comparison must match what the data curve has had done to it:
    // purity-corrected data -> matched reco (numerator); uncorrected data -> all reco (denominator).
    if (split_test){
        std::cout << "\t----> Doing split test" << std::endl;
        h_full_purity = getOrWarn<TH2D>(fin_unfolding, "h_full_pseudo_purity_tf");
        h_full_efficiency = getOrWarn<TH2D>(fin_unfolding, "h_full_pseudo_efficiency_tf");
        h_mc_reco = getOrWarn<TH2D>(fin_unfolding, apply_purity ? "h_full_pseudo_purity_numerator_tf"
                                                                : "h_full_pseudo_purity_denominator_tf");
        response = getOrWarn<RooUnfoldResponse>(fin_unfolding, "response_tf_pseudo_full");
        h_mc_true_no_eff = getOrWarn<TH2D>(fin_unfolding, "h_full_pseudo_efficiency_numerator_tf");
    }
    else {
        h_full_purity = getOrWarn<TH2D>(fin_unfolding, "h_full_purity_tf");
        h_full_efficiency = getOrWarn<TH2D>(fin_unfolding, "h_full_efficiency_tf");
        h_mc_reco = getOrWarn<TH2D>(fin_unfolding, apply_purity ? "h_full_purity_numerator_tf"
                                                                : "h_full_purity_denominator_tf");
        response = getOrWarn<RooUnfoldResponse>(fin_unfolding, "response_tf_full");
        h_mc_true_no_eff = getOrWarn<TH2D>(fin_unfolding, "h_full_efficiency_numerator_tf");
    }
    if (!h_full_purity || !h_full_efficiency || !h_mc_reco || !response) return;



    // ---- Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    std::cout << "\t---->Condition number nominal = " << cond_number
              << std::endl;

    // ---- Grab the truth level MC ----
    TString fname_response_truth =  "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root"; //"/data_CMS/cms/zaidan/analysis_lise/Run3/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
    std::cout << "Getting truth from : " << fname_response_truth << std::endl;
    TFile *fin_response_truth = openOrWarn(fname_response_truth);
    if (!fin_response_truth) return;
    TH2D *h_mc_true = nullptr;
    if (test_mode == 1) {
        // Split test: the even half's gen distribution -- same gen_pass gate and same w_reco
        // weight the correction chain outputs. Independent from the odd-half corrections.
        h_mc_true = getOrWarn<TH2D>(fin_response_truth, "h_pseudodata_truth_tf");
    } else if (test_mode == 0) {
        // Full-MC closure: the full-sample all-gen distribution the chain recovers by construction.
        h_mc_true = getOrWarn<TH2D>(fin_response_truth, "h_full_efficiency_denominator_tf");
    } else {
        // Data: gen reference.
        h_mc_true = getOrWarn<TH2D>(fin_response_truth, "hgenjet_2b_passbtag");
    }
    if (!h_mc_true) return;


    //------- Apply purity correction
    TH2D *h_data_purity_corrected = (TH2D *) h_data_after_fit->Clone("h_data_purity_corrected");
    if (apply_purity) {
        std::cout << "\t---->Multiplying data by purity" << std::endl;
        h_data_purity_corrected->Multiply(h_full_purity);
    } else {
        std::cout << "\t---->NOT multiplying data by purity" << std::endl;
    }

    //double underflowY = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1, 0, 0);
    //double overflowY  = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1,
                          //h_data_purity_corrected->GetNbinsY()+1, h_data_purity_corrected->GetNbinsY()+1);
    //std::cout << "underflowY=" << underflowY << " overflowY=" << overflowY << std::endl;

    // ---- Scan the Bayesian iteration count, one PDF page per value (only when scan_niter).
    TString home_llr = "/home/llr/cms/zaidan/";                       // LLR home directory
    // label already encodes the option flags (_bayesian / _split_test / _nopurity), so toggling
    // those bools writes to distinct files instead of overwriting the previous run's outputs.
    TString iter_pdf = home_llr + "bayesian_unfolding_iterations_" + label + ".pdf";  // multi-page PDF of all iterations
    TString png_dir  = home_llr + "bayesian_unfolding_iterations_" + label + "_png/"; // one PNG per iteration lands here
    if (scan_niter) gSystem->mkdir(png_dir, kTRUE);                     // create the PNG folder if it does not exist
    const int niter_min = 0;
    const int niter_max = scan_niter ? 99 : 0;   // one pass only when not scanning
    const Int_t prev_ignore = gErrorIgnoreLevel;
    if (scan_niter) gErrorIgnoreLevel = kError;  // silence "Replacing existing histogram" spam from re-cloning each pass

    for (int iter = niter_min; iter <= niter_max; ++iter) {

    // ---- Unfold
    std::cout << "\t---->Unfolding (niter = " << (scan_niter ? iter + 1 : 4) << ")" << std::endl;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    TH2D *h_data_unfolded = nullptr;
    TMatrixD covariance_matrix_before_unfolding(dim,dim);
    TMatrixD covariance_matrix_after_unfolding(dim,dim);
    if (unfoldBayes) {
        Int_t niter = scan_niter ? iter + 1 : 4;   // scan: 1..100 ; otherwise fixed 4
        RooUnfoldBayes unfold(response, h_data_purity_corrected, niter);
        // Clone before `unfold` leaves scope: some RooUnfold versions hand back a cached
        // histogram that the unfolder owns and deletes with itself.
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment)->Clone("h_data_unfolded");
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    } else {
        RooUnfoldInvert unfold(response, h_data_purity_corrected);
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment)->Clone("h_data_unfolded");
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    }

    // ---- Fold back
    std::cout << "\t---->Refolding" << std::endl;
    TH2D *h_data_refolded = (TH2D *) response->ApplyToTruth(h_data_unfolded, "h_data_refolded");

    // ---- Apply efficiency correction
    std::cout << "\t---->Dividing by recostruction efficiency" << std::endl;
    TH2D *h_data_efficiency_corrected = (TH2D *) h_data_unfolded->Clone("h_data_efficiency_corrected");
    h_data_efficiency_corrected->Divide(h_full_efficiency);

    // ---- Final corrections
    TH2D *h_data_fully_corrected = (TH2D *) h_data_efficiency_corrected->Clone("h_data_fully_corrected");


    // ---- Graphical bottomline test

    std::cout << "Performing graphical bottomline test" << std::endl;
    TH1D *h_mc_reco_2D = h_mc_reco->ProjectionX("h_mc_reco_2D", ibin_pt, ibin_pt);
    TH1D *h_mc_true_2D = h_mc_true->ProjectionX("h_mc_true_2D", ibin_pt, ibin_pt);
    TH1D *h_data_purity_corrected_2D = h_data_purity_corrected->ProjectionX("h_data_purity_corrected_2D", ibin_pt, ibin_pt);
    TH1D *h_data_fully_corrected_2D = h_data_fully_corrected->ProjectionX("h_data_fully_corrected_2D", ibin_pt, ibin_pt);
    TH1D *h_data_refolded_2D = h_data_refolded->ProjectionX("h_data_refolded_2D", ibin_pt, ibin_pt);
    TH1D *h_data_after_fit_2D = h_data_after_fit->ProjectionX("h_data_after_fit_2D", ibin_pt, ibin_pt);
    TH1D *h_data_unfolded_2D = h_data_unfolded->ProjectionX("h_data_unfolded_2D", ibin_pt, ibin_pt);

    // A Divide() across incompatible binnings returns silently and leaves the ratio empty.
    if (h_mc_true_2D->GetNbinsX() != h_data_fully_corrected_2D->GetNbinsX()) {
        std::cerr << "ERROR: truth has " << h_mc_true_2D->GetNbinsX() << " dr bins but the unfolded result has "
                  << h_data_fully_corrected_2D->GetNbinsX() << std::endl;
        return;
    }

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
                        if (!h) continue;
                        h->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
                        normalizeToUnitArea(h);
                    }
        // The y range has to clear the error bars of the curves that are actually drawn, not
        // just their central values, and not the ones left out of the plot.
        for (auto h : {h_mc_reco_2D, h_mc_true_2D, h_data_purity_corrected_2D, h_data_fully_corrected_2D}) {
            for (int i = ibin_dr_min; i <= ibin_dr_max; ++i)
                ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
        }

        Float_t pt_min_plot = h_data_reco->GetYaxis()->GetBinLowEdge(ibin_pt);
        Float_t pt_max_plot = h_data_reco->GetYaxis()->GetBinUpEdge(ibin_pt);

        TLegend *leg = new TLegend(0.15, 0.55, 0.55, 0.80);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetMargin(0.15);
        leg->SetTextFont(font_code);
        leg->SetTextSize(legend_size);
        leg->SetHeader(Form("%.0f < p_{T}^{jet} < %.0f GeV", pt_min_plot, pt_max_plot));

        // Name every curve by what it is and, in the split test, which half it came from.
        TString lbl_reco     = split_test ? "Reco pseudodata"  : "Reco data";
        TString lbl_mc_reco  = split_test ? (apply_purity ? "Reco MC, matched"
                                                          : "Reco MC, all reco")
                                          : "Reco MC";
        TString lbl_unfolded = split_test ? "Unfolded pseudodata"  : "Unfolded data";
        // Same gen distribution, two disjoint samples. Only the half differs.
        TString lbl_mc_true  = split_test ? "Gen MC"  : "Gen MC";

        // Only the first histogram drawn owns the pad's axes; style them there and nowhere else.
        h_data_purity_corrected_2D->SetTitle("");//Data, " + label + " " + " response matrix");
        h_data_purity_corrected_2D->SetStats(0);
        h_data_purity_corrected_2D->SetMarkerColor(purple);
        h_data_purity_corrected_2D->SetLineColor(purple);
        h_data_purity_corrected_2D->SetMarkerStyle(kFullCircle);
        h_data_purity_corrected_2D->SetMarkerSize(1);
        // 1.75x headroom keeps the legend clear of the EEC peak.
        h_data_purity_corrected_2D->GetYaxis()->SetRangeUser(0., ymax*1.75);
        // TLatex (#Delta), not TMathText (\Delta\mbox{r}): TMathText axis titles do not get
        // placed reliably — the x-axis title was silently dropped altogether.
        h_data_purity_corrected_2D->GetYaxis()->SetTitle("EEC(#Delta r)");
        h_data_purity_corrected_2D->GetYaxis()->CenterTitle(true);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleFont(font_code);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleSize(title_size);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleOffset(1.55);
        h_data_purity_corrected_2D->GetYaxis()->SetLabelFont(font_code);
        h_data_purity_corrected_2D->GetYaxis()->SetLabelSize(label_size);
        h_data_purity_corrected_2D->GetXaxis()->SetLabelSize(0); // x labels belong to the ratio pad
        h_data_purity_corrected_2D->GetXaxis()->SetTitleSize(0);
        leg->AddEntry(h_data_purity_corrected_2D, lbl_reco, "pe1");


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
        leg->AddEntry(h_mc_reco_2D, lbl_mc_reco, "pe1");

        h_data_fully_corrected_2D->SetMarkerColor(blue);
        h_data_fully_corrected_2D->SetLineColor(blue);
        h_data_fully_corrected_2D->SetMarkerStyle(kOpenCross);
        h_data_fully_corrected_2D->SetMarkerSize(1);
        leg->AddEntry(h_data_fully_corrected_2D, lbl_unfolded, "pe1");

        h_mc_true_2D->SetMarkerColor(orange);
        h_mc_true_2D->SetLineColor(orange);
        h_mc_true_2D->SetMarkerStyle(kOpenTriangleUp);
        h_mc_true_2D->SetMarkerSize(1);
        leg->AddEntry(h_mc_true_2D, lbl_mc_true, "lp");   // line + open-triangle marker

        h_data_refolded_2D->SetMarkerColor(blue);
        h_data_refolded_2D->SetLineColor(blue);
        h_data_refolded_2D->SetMarkerStyle(kFullCross);
        h_data_refolded_2D->SetMarkerSize(1);


      
        
        TCanvas *c_unfold = new TCanvas("c_unfold", "", 800, 800);

        TPad *pad_main  = new TPad("pad_main",  "", 0., 0.3, 1., 1.);   // placeholder coords;
        TPad *pad_ratio = new TPad("pad_ratio", "", 0., 0.,  1., 0.3);  // applySquareLayout owns them
        
    
        
        
        for (TPad *p : {pad_main, pad_ratio}) {
            p->SetTicks(1, 0);   // mirrored ticks on all four sides
            p->SetFillColor(0);
        
        }
        //pad_ratio->SetLogx();
        //pad_main->SetLogx();

        pad_main->SetMargin(0.1, 0.1, 0.0, 0.1);
        pad_ratio->SetMargin(0.1, 0.1, 0.23, 0.0);
        
    
    
        // Attach the pads to the canvas before drawing into them, so gPad and the pad's
        // primitive list agree at every Draw().
        c_unfold->cd();
        pad_main->Draw();
        pad_ratio->Draw();

        pad_main->cd();
        h_data_purity_corrected_2D->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleOffset(title_offset);
        h_data_purity_corrected_2D->Draw("PE X0");
        //h_data_unfolded_2D->Draw("pe1 same");
        h_mc_reco_2D->Draw("PE X0 same");
        h_data_fully_corrected_2D->Draw("PE X0 same");
        //h_data_after_fit_2D->Draw("pe1 same");
        h_mc_true_2D->Draw("hist E same");

        leg->Draw();
        TLatex *test_info_text = new TLatex;
        
        test_info_text->SetNDC();
        test_info_text->SetTextFont(font_code);
        test_info_text->SetTextSize(label_size);
        test_info_text->DrawLatex(0.57, 0.75,
            test_mode == 0 ? "MC full-sample closure" :
            test_mode == 1 ? "MC split closure test"  :
                             "Data Run 3");
        test_info_text->DrawLatex(0.57, 0.70, "2D unfolding");
        if(unfoldBayes) {
            test_info_text->DrawLatex(0.57, 0.65, "Bayesian unfolding");
            test_info_text->DrawLatex(0.57, 0.60, Form("N iter = %d", scan_niter ? iter + 1 : 4));
        }
        else {
            test_info_text->DrawLatex(0.57, 0.65, "Matrix inversion unfolding");
        }
        //drawHeader();

        pad_main->RedrawAxis();

        
        TH1D *h_data_mc_reco_ratio = (TH1D *) h_data_purity_corrected_2D->Clone("h_data_mc_reco_ratio");
        
        h_data_mc_reco_ratio->SetTitle("");
        h_data_mc_reco_ratio->SetStats(0);
        h_data_mc_reco_ratio->Divide(h_mc_reco_2D);
        h_data_mc_reco_ratio->SetMarkerStyle(kFullCircle);
        h_data_mc_reco_ratio->SetMarkerColor(purple);
        h_data_mc_reco_ratio->SetLineColor(purple);
        h_data_mc_reco_ratio->SetMarkerSize(1.1);

        TH1D *h_data_mc_true_ratio = (TH1D *) h_data_fully_corrected_2D->Clone("h_data_mc_true_ratio");
        h_data_mc_true_ratio->Divide(h_mc_true_2D);
        h_data_mc_true_ratio->SetMarkerStyle(kOpenCross);
        h_data_mc_true_ratio->SetMarkerColor(blue);
        h_data_mc_true_ratio->SetLineColor(blue);
        h_data_mc_true_ratio->SetMarkerSize(1.3);

        TH1D *h_mc_gen_reco_ratio = (TH1D *) h_mc_true_2D->Clone("h_mc_gen_reco_ratio");


        h_mc_gen_reco_ratio->Divide(h_mc_reco_2D);
        h_mc_gen_reco_ratio->SetLineColor(red);
        h_mc_gen_reco_ratio->SetLineWidth(1);
        h_mc_gen_reco_ratio->SetMarkerSize(0);
                h_mc_gen_reco_ratio->GetYaxis()->SetTitle("ratio");
                h_mc_gen_reco_ratio->GetXaxis()->SetTitle("#Delta r");
                        h_mc_gen_reco_ratio->GetYaxis()->CenterTitle(true);
                        h_mc_gen_reco_ratio->GetXaxis()->CenterTitle(true);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleSize(0.2); 
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleOffset(0.28);

        h_mc_gen_reco_ratio->GetXaxis()->SetLabelSize(0.06);
        h_mc_gen_reco_ratio->GetYaxis()->SetLabelSize(0.1);
        h_mc_gen_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleOffset(0.24);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleSize(0.2);
        h_mc_gen_reco_ratio->GetYaxis()->SetLimits(0.5, 2);
        h_mc_gen_reco_ratio->SetLineStyle(2);
        h_mc_gen_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);

        // Match the ratio-pad axes to the main pad: pixel (precision-43) fonts with the same
        // title size, label size and title offset, so every axis name is the same size and the
        // same distance from its numbers. Overrides the fractional sizes set just above.
        for (TAxis *ax : { h_mc_gen_reco_ratio->GetXaxis(), h_mc_gen_reco_ratio->GetYaxis() }) {
            ax->SetTitleFont(font_code);
            ax->SetTitleSize(title_size);
            ax->SetTitleOffset(title_offset);
            ax->SetLabelFont(font_code);
            ax->SetLabelSize(label_size);
        }


        // Green, not blue: blue already belongs to the unfolded points, and two blue curves in
        // one pad read as the same quantity.
        TH1D *h_unfolded_reco_ratio = (TH1D *) h_data_fully_corrected_2D->Clone("h_unfolded_reco_ratio");
        h_unfolded_reco_ratio->Divide(h_data_purity_corrected_2D);
        h_unfolded_reco_ratio->SetLineColor(green);
        h_unfolded_reco_ratio->SetLineWidth(1);
        h_unfolded_reco_ratio->SetMarkerSize(0);
        h_unfolded_reco_ratio->SetLineStyle(2);
        

        h_data_mc_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);
        h_data_mc_reco_ratio->GetYaxis()->SetTitle("ratio");
        h_data_mc_reco_ratio->GetYaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleFont(font_code);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleSize(title_size);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleOffset(0.05);
        h_data_mc_reco_ratio->GetYaxis()->SetLabelFont(font_code);
        h_data_mc_reco_ratio->GetYaxis()->SetLabelSize(label_size);
        h_data_mc_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_data_mc_reco_ratio->GetXaxis()->SetTitle("#Delta r");
        h_data_mc_reco_ratio->GetXaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleFont(font_code);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleSize(title_size);
        
        // With precision-3 fonts this offset scales off the (short) pad height: 3.2 pushed the
        // title clean off the canvas. 1.35 lands it just under the tick labels.
        h_data_mc_reco_ratio->GetXaxis()->SetTitleOffset(0.05);
        h_data_mc_reco_ratio->GetXaxis()->SetLabelFont(12);
        h_data_mc_reco_ratio->GetXaxis()->SetLabelSize(16);
        h_data_mc_reco_ratio->GetXaxis();
       

        // The unity line and the gap bands span the axis actually plotted, not the full dr
        // range: the two disagree as soon as ibin_dr_min/ibin_dr_max restrict the view.
        const Double_t x_lo = h_mc_reco_2D->GetXaxis()->GetBinLowEdge(ibin_dr_min);
        const Double_t x_hi = h_mc_reco_2D->GetXaxis()->GetBinUpEdge(ibin_dr_max);

        pad_ratio->cd();
        //h_data_mc_reco_ratio->Draw("axis");   // axes first, so the bands land behind everything

      
        TLine *line = new TLine(x_lo, 1., x_hi, 1.);
        line->SetLineWidth(1.);
        line->SetLineStyle(kSolid);   // solid, so the unity line reads apart from the dashed red/green curves
        line->SetLineColor(kGray+2);

        // "Before" band: for each dr bin, fill the full bin width between the gen MC/reco MC line
        // (red) and the unfolded/reco line (green) -- the size of the correction. Stepped closed
        // polygon: forward along red at the bin edges, back along green, so each bin is filled edge
        // to edge. Drawn behind, translucent light grey.
        const int nb_band = ibin_dr_max - ibin_dr_min + 1;
        TGraph *band_before_after = new TGraph(4 * nb_band);
        int ip = 0;
        for (int i = 0; i < nb_band; ++i) {              // top edge: along red, left -> right
            const int bin = ibin_dr_min + i;
            const double xlo = h_mc_gen_reco_ratio->GetXaxis()->GetBinLowEdge(bin);
            const double xhi = h_mc_gen_reco_ratio->GetXaxis()->GetBinUpEdge(bin);
            const double yr  = h_mc_gen_reco_ratio->GetBinContent(bin);
            band_before_after->SetPoint(ip++, xlo, yr);
            band_before_after->SetPoint(ip++, xhi, yr);
        }
        for (int i = nb_band - 1; i >= 0; --i) {         // bottom edge: along green, right -> left
            const int bin = ibin_dr_min + i;
            const double xlo = h_unfolded_reco_ratio->GetXaxis()->GetBinLowEdge(bin);
            const double xhi = h_unfolded_reco_ratio->GetXaxis()->GetBinUpEdge(bin);
            const double yg  = h_unfolded_reco_ratio->GetBinContent(bin);
            band_before_after->SetPoint(ip++, xhi, yg);
            band_before_after->SetPoint(ip++, xlo, yg);
        }
        // Solid light grey, NOT SetFillColorAlpha: the PDF backend does not render alpha
        // transparency (it shows on screen in XQuartz but disappears in the PDF), whereas a
        // solid fill renders identically everywhere. Curves are redrawn on top, so opaque is fine.
        band_before_after->SetFillColor(TColor::GetColor(224, 224, 224));
        band_before_after->SetFillStyle(1001);
        band_before_after->SetLineWidth(0);

        h_mc_gen_reco_ratio->Draw("hist");     // axes + red line
        band_before_after->Draw("f same");     // band behind the curves
        h_mc_gen_reco_ratio->Draw("hist same");// redraw red crisp over the band
        line->Draw("same");

        h_unfolded_reco_ratio->Draw("hist same");
        h_data_mc_reco_ratio->Draw("PE X0 same");
        h_data_mc_true_ratio->Draw("PE X0 same");

        // Curves grouped by family (points = closure, lines = size of correction), then the two
        // band swatches, so the legend reads top-down as the story: correction, closure, and the
        // two gaps that connect them.
        TLegend *leg_ratio = new TLegend(0.16, 0.66, 0.96, 0.985);
        leg_ratio->SetFillStyle(0);
        leg_ratio->SetBorderSize(0);
        leg_ratio->SetTextFont(font_code);
        leg_ratio->SetTextSize(legend_size - 2.);
        leg_ratio->SetNColumns(2);
        leg_ratio->AddEntry(h_data_mc_reco_ratio, "reco data/reco MC", "pe1");
        leg_ratio->AddEntry(h_mc_gen_reco_ratio, "gen MC/reco MC", "l");
        leg_ratio->AddEntry(h_data_mc_true_ratio, "unfolded data/gen MC", "pe1");
        leg_ratio->AddEntry(h_unfolded_reco_ratio, "unfolded/reco", "l");
        leg_ratio->Draw();
        pad_ratio->RedrawAxis();
        c_unfold->cd();
        c_unfold->Update();
        
        h_data_mc_reco_ratio->GetXaxis();
        c_unfold->Modified();
        c_unfold->Update();
       
        if (scan_niter) {
            // One page per iteration: "(" opens the PDF on the first page, ")" closes it on the last.
            TString iter_page = iter_pdf;
            if (iter == niter_min)      iter_page += "(";
            else if (iter == niter_max) iter_page += ")";
            c_unfold->Print(iter_page, "pdf");

            // Also drop a standalone PNG per iteration into png_dir (zero-padded so they sort 001..100).
            c_unfold->Print(Form("%siteration_%03d.png", png_dir.Data(), iter + 1), "png");
        } else {
            // Single run: the usual one-off bottomline plot next to the other outputs.
            TString plot_stem = folder + "unfolding_plot" + label + "_bottomline_test_eec_2D";
            c_unfold->Print(plot_stem + ".pdf");
            c_unfold->Print(plot_stem + ".png");
        }

        // Only delete the canvas when another pass will recreate it (avoids the name clash during
        // the scan). On the final/only pass keep it, so the window stays open to view.
        if (iter < niter_max) delete c_unfold;


    }

    // Write the corrected histograms once, from the final iteration's result.
    if (iter == niter_max) {
    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    h_data_fully_corrected->SetName("h_data_unfolded_corrected");
    h_data_fully_corrected->Write();
    h_data_after_fit->SetName("h_data");
    h_data_after_fit->Write();
    h_mc_reco->SetName("h_mc_reco");
    h_mc_reco->Write();

    //Write 2D histograms for plotting
    h_mc_reco_2D->Write();
    //h_mc_true_2D->Write();
    h_data_purity_corrected_2D->Write();
    h_data_fully_corrected_2D->Write();
    h_data_refolded_2D->Write();

    fout->Close();
    delete fout;
    }

    }   // end loop over iteration counts
    gErrorIgnoreLevel = prev_ignore;

    // gApplication -> Terminate(0);
}

// test_mode: 0 = full-MC closure, 1 = split test, 2 = data
// unfoldBayes: true = Bayesian, false = matrix inversion
void apply_unfolding_2d(int test_mode = 1, bool unfoldBayes = true){
    TString dataset = "data";
    TString folder = "/data_CMS/cms/zaidan/analysis_lise/Run3/";
    TString pT_selection = "80_200";
    bool btag = true;
    Int_t n = 1;
	apply_unfolding(dataset, folder, btag, n, pT_selection, test_mode, unfoldBayes);

}
