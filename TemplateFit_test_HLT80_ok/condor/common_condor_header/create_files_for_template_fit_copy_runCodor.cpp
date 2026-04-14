#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2.h>
#include "TAxis.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <iostream>
#include <Math/Vector4D.h>
#include <map>
#include <unordered_map>
#include <string>
#include <random>
#include <vector>
#include "TString.h"
#include "Math/VectorUtil.h"
#include "tTree.h"
#include "binning_histos_all.h" 
#include "Aggregation.h"



// Build templates for the template fit.
// MC:   fills h3D_b (jtNbHad==1) and h3D_bb (jtNbHad==2) using reco-level cuts and reco SV reconstruction.
// Data: fills h3D_data with the same reco logic — no truth classification.
void make_templates(TString filename, TString output_folder, TString output_hist, TString domain,
                    Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC, Int_t dataType,
                    Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1, bool applyHLT80 = false, bool applyHLT80100 = false, bool applyHLT40 = false, bool applyHLT4060 = false,) {

  /// test flag 
  // bool myflag = false;

  // -- prescale factor (only for 40 GeV trigger)
    double prescale_pf40 = 33.917210;


  tTree t;
  t.Init(filename, isMC);
  t.SetBranchStatus("*", 0);

  double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0;

  std::vector<TString> active_branches = {
    // reco branches — identical for data and MC
    "jtpt", "jteta", "nref", "jtNtrk", "jtNsvtx", "discr_particleNet_BvsAll",
    "ntrk", "trkJetId", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
    "trkSvtxId", "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",
    "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
    "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1"};
  if (isMC) {
    // MC-only: event weight, pthat for skipMC, and truth b-hadron count for template classification
    std::vector<TString> mc_branches = {"weight", "pthat", "jtNbHad"};
    active_branches.insert(active_branches.end(), mc_branches.begin(), mc_branches.end());
  }
  t.SetBranchStatus(active_branches, 1);

  // MC: separate 0b, b and bb templates
  TH3D *h3D_0b = new TH3D("h3D_0b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b  = new TH3D("h3D_b",  "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_bb = new TH3D("h3D_bb", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  // Data: single distribution to be fit
  TH3D *h3D_data = new TH3D("h3D_data", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);

  h3D_0b->Sumw2();   h3D_0b->SetCanExtend(TH1::kNoAxis);
  h3D_b->Sumw2();    h3D_b->SetCanExtend(TH1::kNoAxis);
  h3D_bb->Sumw2();   h3D_bb->SetCanExtend(TH1::kNoAxis);
  h3D_data->Sumw2(); h3D_data->SetCanExtend(TH1::kNoAxis);

  // Jet counts (no EEC weight): 3D (mB, dr, jtpt) — same axes as the EEC histograms
  TH3D *h_count_0b   = new TH3D("h_count_0b",   "jet counts 0b;m_{B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_b    = new TH3D("h_count_b",    "jet counts 1b;m_{B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_bb   = new TH3D("h_count_bb",   "jet counts 2b;m_{B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_data = new TH3D("h_count_data", "jet counts data;m_{B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  
  h_count_0b->Sumw2();   h_count_0b->SetCanExtend(TH1::kNoAxis);
  h_count_b->Sumw2();    h_count_b->SetCanExtend(TH1::kNoAxis);
  h_count_bb->Sumw2();   h_count_bb->SetCanExtend(TH1::kNoAxis);
  h_count_data->Sumw2(); h_count_data->SetCanExtend(TH1::kNoAxis);

  Long64_t n_events = t.GetEntries();
  if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
  if (ev_first < 0) ev_first = 0;
  std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

    cout << "n used entries = " << ev_last << endl;


  for (Long64_t ient = ev_first; ient < ev_last; ient++) { 

    // cout << "this is event " << ient << endl;
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
    t.GetEntry(ient);

    double weight_tree = isMC ? t.weight : 1.0;
    double prescale = 1; 


    // trigger selection
    // data triggers 
    if (!isMC && dataType == 0) { // High EG data 
      if (applyHLT80 && !(t.HLT_HIAK4PFJet80_v1 == 1)) {continue;}
      // default configuration 
      if (applyHLT80100 && !(t.HLT_HIAK4PFJet80_v1 == 1 || t.HLT_HIAK4PFJet100_v1 == 1)) continue; 
    }

    if (!isMC && dataType == -1) { // low EG data 
      //  
      if (!((t.HLT_HIAK4PFJet60_v1 == 1 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) ||
            (t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0))) continue;
    
    }

    // MC triggers 
    if (isMC) {
      if ( applyHLT80 && !(t.HLT_HIAK4PFJet80_v1 == 1)) continue;
      // default configuration
      if ( applyHLT40 && !(t.HLT_HIAK4PFJet40_v1 == 1)) continue;
    }




    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // if (myflag) { cout << "test event is FOUND! EXIT now :D "<< endl;return;}
      // cout  << "this is ijet: "<< ijet << endl;
      // std::cout << "jet pt: " << t.jtpt[ijet] << std::endl;


      // reco-level cuts — identical for data and MC
      if (std::abs(t.jteta[ijet]) > 1.9) { continue;}
      if (isMC && skipMC(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < pT_low || t.jtpt[ijet] > pT_high) { continue;}
      if (btag && t.discr_particleNet_BvsAll[ijet] <= 0.99) {continue;}

      // reco SV reconstruction — same for data and MC
      vector<ROOT::Math::PtEtaPhiMVector> reco_sv = makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr);
      // cout << "Before checking #aggreagted , reco vector size  " << reco_sv.size() << endl; // << " and " << reco_sv [1] << endl;
           

      if (reco_sv.size() != 2) {continue;}

      double dr   = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(), reco_sv[1].Eta(), reco_sv[1].Phi());
      double pt1  = reco_sv[0].Pt();
      double pt2  = reco_sv[1].Pt();
      double eec  = std::pow(pt1 * pt2, n);
      double jtpt = t.jtpt[ijet];
      double mB   = reco_sv[0].M() + reco_sv[1].M();

      if (mB > mb_max_fill) mB = mb_max_fill;  // fold overflow into last bin
      /// The overflow fr dr is missing, but hsould be zero anyway :) 
/*
      // std::cout << "-------------------------------" << endl;
      // myflag = true;
      std::cout << "First reco event with 2Bs, change flag to true, print info and break, flag is " << myflag << endl; 
    
      std::cout << "event number: " << ient << std::endl;
      std::cout << "jet: " << ijet << std::endl;
      std::cout << "tree weight: " << weight_tree << std::endl;
      cout << "reco Bs are " << reco_sv[0] << " and " << reco_sv [1] << endl;
      std::cout << "eec: " << eec << std::endl;
      std::cout << "pt1 , pt2: " << pt1  << ", "<<  pt2 << std::endl;
      std::cout << "mass: " << mB << std::endl;
      std::cout << "dr: " << dr << std::endl;
      std::cout << "jet pt: " << jtpt << std::endl;
      std::cout << "-------------------------------" << endl;
*/
      
      if (isMC) {
        // use truth to classify: fill separate 0b, b and bb templates
        if      (t.jtNbHad[ijet] == 0) { h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree); h_count_0b->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 1) { h3D_b ->Fill(mB, dr, jtpt, eec * weight_tree); h_count_b ->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 2) { h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree); h_count_bb->Fill(mB, dr, jtpt, weight_tree); }
      } else {
        h3D_data->Fill(mB, dr, jtpt, eec * weight_tree);
        h_count_data->Fill(mB, dr, jtpt, weight_tree);
      }
    }
  }
  std::cout << std::endl;

  TString label = btag ? "_btag" : "_nobtag";
  TString label_trigger = applyHLT80 ? "_HLT80": "";
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
  TFile outFile((output_folder + output_hist + label_trigger + label + job_suffix + domain).Data(), "RECREATE");
  if (isMC) {
    h3D_0b->Write();
    h3D_b->Write();
    h3D_bb->Write();
    h_count_0b->Write();
    h_count_b->Write();
    h_count_bb->Write();

    // for fast test 
    // with eec 
    TH1D* h3dproj_mass_0b = (TH1D*) h3D_0b->ProjectionX("h3dproj_mass_0b", 1, h3D_0b->GetNbinsY(), 1, h3D_0b->GetNbinsZ());
          h3dproj_mass_0b->SetTitle("eec*weight (when needed)"); h3dproj_mass_0b->Write(); 
    TH1D* h3dproj_mass_b = (TH1D*) h3D_b->ProjectionX("h3dproj_mass_b", 1, h3D_b->GetNbinsY(), 1, h3D_b->GetNbinsZ());
          h3dproj_mass_b->SetTitle("eec*weight (when needed)"); h3dproj_mass_b->Write(); 
    TH1D* h3dproj_mass_bb = (TH1D*) h3D_bb->ProjectionX("h3dproj_mass_bb", 1, h3D_bb->GetNbinsY(), 1, h3D_bb->GetNbinsZ());
          h3dproj_mass_bb->SetTitle("eec*weight (when needed)"); h3dproj_mass_bb->Write(); 

    // without eec 
    TH1D* h3dproj_mass_count_0b = (TH1D*) h_count_0b->ProjectionX("h3dproj_mass_count_0b", 1, h_count_0b->GetNbinsY(), 1, h_count_0b->GetNbinsZ());
          h3dproj_mass_count_0b->SetTitle("weight (when needed)"); h3dproj_mass_count_0b->Write(); 
    
    TH1D* h3dproj_mass_count_b = (TH1D*) h_count_b->ProjectionX("h3dproj_mass_count_b", 1, h_count_b->GetNbinsY(), 1, h_count_b->GetNbinsZ());
          h3dproj_mass_count_b->SetTitle("weight (when needed)"); h3dproj_mass_count_b->Write(); 

    TH1D* h3dproj_mass_count_bb = (TH1D*) h_count_bb->ProjectionX("h3dproj_mass_count_bb", 1, h_count_bb->GetNbinsY(), 1, h_count_bb->GetNbinsZ());
          h3dproj_mass_count_bb->SetTitle("weight (when needed)"); h3dproj_mass_count_bb->Write(); 
  } else {
    h3D_data->Write();
    h_count_data->Write();
  
    // for fast test 
    // with eec 
    TH1D* h3dproj_mass = (TH1D*) h3D_data->ProjectionX("h3dproj_mass", 1, h3D_data->GetNbinsY(), 1, h3D_data->GetNbinsZ());
    h3dproj_mass->SetTitle("eec*weight (when needed)");
    h3dproj_mass->Write(); 

    // without eec 
    TH1D* h3dproj_mass_counts = (TH1D*) h_count_data->ProjectionX("h3dproj_mass_counts", 1, h_count_data->GetNbinsY(), 1, h_count_data->GetNbinsZ());
    h3dproj_mass_counts->SetTitle("weight (when needed)");
    h3dproj_mass_counts->Write();
  }

 

  outFile.Close();
}


void create_files_for_template_fit(Int_t dataType = 0, Float_t pT_low = 80, Float_t pT_high = 140, Int_t n = 1, bool btag = true, bool isMC = true, bool applyHLT80 = false, bool applyHLT40 = false){
 
TString filename;
TString output_hist;
TString output_folder = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/Test_vs_Zoe_code/Zoe_version/eec_2b_analysis-main/";
TString domain = ".root";

// -- test result of high EG triggered data  
// dataType = 0; // HG 
// isMC = false; // data 

// -- test: bjet MC 
dataType = 1;
isMC = true;
applyHLT80 = true;

//sanity check
if (isMC && dataType < 1) {
  std::cerr << "Invalid data type for MC sample" << std::endl;
  return;}

if (!isMC && dataType > 1) {
  std::cerr << "Invalid data type for data sample" << std::endl;
  return;}

 


if(dataType == -1){//________________________________data______________________________
  filename = "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root";
  output_hist = "template_for_fit_histos_3D_LowEG";
  isMC = false;
  cout<<"you chose data Low" <<endl;
  }

else if(dataType == 0) {
  filename = "/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
  output_hist = "template_for_fit_histos_3D_HighEG";
  isMC = false;
  cout<<"you chose data High" <<endl;       
  }      
                                                                                                                                                                                                                                                                        
else if(dataType == 1){//________________________________bjet______________________________
  filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
  output_hist = "template_for_fit_histos_3D_bjet";
  std::cout << "Creating files for template fit for bjet sample" << std::endl;
  cout<<"you chose bjet MC" <<endl;
  }

else if(dataType == 2){//________________________________dijet______________________________
  filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root"; 
  output_hist = "template_for_fit_histos_3D_qcd";
  std::cout << "Creating files for template fit for qcd sample" << std::endl;
  cout<<"you chose qcd MC" <<endl;
  }

else{
  cout<<" undefined data type"<<endl;
  return; 
  }


  // 
  make_templates(filename, output_folder, output_hist, domain, pT_low, pT_high, n, btag, isMC, dataType, applyHLT80);


  std::cout << "finished :)" << std::endl;
}

