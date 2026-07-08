
/*
  Studies b-decay track BDT vs jet pT as well as (full, gen) b-hadron pT
  b-decay tracks are associated only to _their_ parent b-hadron
  This code is a useful example for correlating a track to a given b-hadron
  -Matt
 */

/*
- Using Matt. code looping over full B info for the b-decay track BDT vs jet pT and full, gen bhadron pT.
- Afnan: Study SVx finding eff.  
Study in the nominal Jet pt interval(detector level MC selection): 100- 120 GeV ONLY 
Note that no min pt ct on track is applied.

Get_1Svx_eff() : run to get Svx eff. for single and double B. where for each sv, only one track is required to have same sattus as True B.
                this function gives Eff. as function of Full B pt (GeV).
                Eff. is fitted with pol5 function, from 0-140 GeV, while Full B pt drawn fo 0-200 GeV.

test_closure(): test applying Eff. for events pass the Svx cut, corrected for this eff. event by event (value/Eff.). 
                Then compared to the Total histogram without cut for each case.
                A faster closure is done using histpgrams (bin by bin correction).


*/ 

// Legend 
#include"binning_histos_small.h"
#include "TemplateFit_Run3/Help_Functions.h"
#include <vector>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

// #include "CMSStyle.C"
// #include "MyRootColors.h" // didnt work

// -- Global tree variables: 
  const int MAXJETS = 20;
  const int MAXB = 20;
  const int MAXTRACKS = 100;
  Float_t weight;
  Int_t nref;
  Float_t jtpt[MAXJETS];
  Float_t refpt[MAXJETS];
  Float_t jteta[MAXJETS];
  Int_t jtNbHad[MAXJETS];

  Int_t jtNsvtx[MAXJETS]; 

  Int_t ntrk;
  Int_t trkJetId[MAXTRACKS];
  Int_t trkMatchSta[MAXTRACKS];
  Float_t trkBdtScore[MAXTRACKS];
  Int_t trkSvtxId[MAXTRACKS];
  
  Int_t nfullB;
  Int_t fullBJetId[MAXB];
  Int_t fullBSta[MAXB];
  Float_t fullBPt[MAXB];
  Float_t fullBPhi[MAXB];
  Float_t fullBEta[MAXB];
  Float_t fullBM[MAXB];

  Int_t nrefTrk;
  Int_t refTrkJetId[MAXTRACKS];
  Float_t refTrkPt[MAXTRACKS];
  Float_t  refTrkPhi[MAXTRACKS];
  Float_t  refTrkEta[MAXTRACKS];
  Int_t refTrkPdgId[MAXTRACKS];
  Int_t refTrkSta[MAXTRACKS];




// Functions declarations:

TH1D* DivideByFunction(TH1D* h, TF1* f, const char* name);
void SetBranches(TTree* t);
void PartialBsAggregation_globaltree(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, Int_t ijet);


void Get_1Svx_eff(){
  TFile *fin = new TFile("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root");
  TTree *t = (TTree*)fin->Get("ak4PFJetAnalyzer/t");
  TTree *hie = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");  
  // const int MAXJETS = 20;
  // const int MAXB = 20;
  // const int MAXTRACKS = 100;
 
  // Float_t weight;

  // Int_t nref;
  // Float_t jtpt[MAXJETS];
  // Float_t refpt[MAXJETS];
  // Float_t jteta[MAXJETS];
  // Int_t jtNbHad[MAXJETS];

  // Int_t ntrk;
  // Int_t trkJetId[MAXTRACKS];
  // Int_t trkMatchSta[MAXTRACKS];
  // Float_t trkBdtScore[MAXTRACKS];
  // Int_t trkSvtxId[MAXTRACKS];
  
  // Int_t nfullB;
  // Int_t fullBJetId[MAXB];
  // Int_t fullBSta[MAXB];
  // Float_t fullBPt[MAXB];
  
  hie->SetBranchAddress("weight",  &weight);
  SetBranches(t);

  // t->SetBranchAddress("nref",  &nref);
  // t->SetBranchAddress("jtpt",  jtpt);
  // t->SetBranchAddress("refpt",  refpt);
  // t->SetBranchAddress("jteta",  jteta);
  // t->SetBranchAddress("jtNbHad", jtNbHad);

  // t->SetBranchAddress("ntrk",  &ntrk);
  // t->SetBranchAddress("trkBdtScore",  trkBdtScore);
  // t->SetBranchAddress("trkJetId", trkJetId);
  // t->SetBranchAddress("trkMatchSta", trkMatchSta);
  // t->SetBranchAddress("trkSvtxId", trkSvtxId);
  // t->SetBranchAddress("nfullB",  &nfullB);
  // t->SetBranchAddress("fullBJetId",  fullBJetId);
  // t->SetBranchAddress("fullBSta",  fullBSta);
  // t->SetBranchAddress("fullBPt",  fullBPt);
  
  // --- MC weight used 
  TH1D* hBpt_w1svxCut_1b = new TH1D("hBpt_w1svxCut_1b", "At least 1 particle from B belongs to 1Svx;p_{T, fullB} [GeV]; counts", 100, 0, 200);
  TH1D* hBpt_1b = new TH1D("hBpt_1b", "All 1b;p_{T, fullB} [GeV]; counts", 100, 0, 200);
  
  TH1D* hBpt_w2svxCut_2b = new TH1D("hBpt_w2svxCut_2b", "At least 2Svx, each has at least 1 particle from a B ;p_{T, fullB} [GeV]; counts", 100, 0, 200);
  TH1D* hBpt_2b = new TH1D("hBpt_2b", "All 2b;p_{T, fullB} [GeV]; counts", 100, 0, 200);
  
  Long64_t nentries = t->GetEntries();
  // nentries = 1e+05;
  cout << "#entries = " << nentries << endl;

  for (Long64_t i = 0; i < nentries; i++) {
    
    if(i%100000==0) cout<<" i "<<i<<endl;
    
    t->GetEntry(i);
    hie->GetEntry(i);

    for (int j = 0; j < nref; j++) { // Jet loop 

      float jetPt = jtpt[j];
      // Jet Kineamtic cuts 
      if(jetPt < 100 || jetPt > 120) continue; // use nominal jet pt 

      if(jteta[j]<-2 || jteta[j]>2)  continue;
      if (refpt[j] <= 0) continue;
      
      // -- Double B jets and single B case 
      // loop over each true B
      for( int ib = 0; ib < nfullB; ib++) // per event
      {
        // For jet, have 2B --> each ID will be j
        if(fullBJetId[ib] != j) continue; // B in jet
        // Before SVX cut
        if (jtNbHad[j] == 1) {hBpt_1b->Fill(fullBPt[ib], weight); } // All single B
        if (jtNbHad[j] >= 2) {hBpt_2b->Fill(fullBPt[ib], weight); }// All Double B
        
        Int_t currentBSta = fullBSta[ib];

          // Find if at least one particle of same status as B belong to Svx in this jet
          for (int k = 0; k < ntrk; k++) // over tracks
          {
              if(trkJetId[k]!=j) continue;
              if(trkMatchSta[k]<1) continue;
              // Shall I cut trk pt > 1 
              // to pass: pass at least one svx requriment with one particle at least coming from the mother B
              if(trkMatchSta[k] == currentBSta)
              {
                if(trkSvtxId[k] >= 0 )
                {
                  // jtNsvtx [nref]
                  // svtxJetId [nsvtx] = jet index 
                  // where the svx id is trkSvtxId[k] which is from 0 to the totoal number of SVx: based on #svx per jet + counter on jets

                  if (jtNbHad[j] == 1) {hBpt_w1svxCut_1b->Fill(fullBPt[ib], weight); } 
                  if (jtNbHad[j] >= 2) {hBpt_w2svxCut_2b->Fill(fullBPt[ib], weight); } 
                  break; // go to next B.
                } 
              }
          }
      }// end B loop
  
    }// end jet loop
  } // end event loop


  // -- Compute SVX eff. for single B and Double-B jets   
    auto divide = [](TH1D* num, TH1D* den, const char* name) -> TH1D* {
        TH1D *h = (TH1D*) num->Clone(name);
        h->Divide(num, den, 1., 1., "b"); // Bionmial error propogation
        return h;
    };

    TH1D*  hSvxEff_1B  = divide(hBpt_w1svxCut_1b, hBpt_1b, "1B"); hSvxEff_1B->GetYaxis()->SetTitle("Efficiency"); hSvxEff_1B ->SetTitle("");
    TH1D*  hSvxEff_2B  = divide(hBpt_w2svxCut_2b, hBpt_2b, "2B"); hSvxEff_2B->GetYaxis()->SetTitle("Efficiency"); hSvxEff_2B ->SetTitle("");

    // -- To Fit Eff. 
    // Fit functions are saved in root file for ater use as func = file->get("functionname") and  func->Eval(100.0);
    TF1* f1 = new TF1("f1", "pol5", 0, 200);
      f1->SetLineColor(kBlue-2);
      f1->SetLineWidth(2);
      hSvxEff_1B->Fit(f1, "R");
      f1->SetRange(0, 200);// To draw full range
      cout << "---------------" << endl;
      cout << "Fit paraemters for single B case" << endl; 
      for (int i = 0; i < 5; i++) cout << f1->GetParameter(i) << std::endl;

    TF1* f2 = new TF1("f2", "pol5", 0, 140);
      f2->SetLineColor(kRed-2);
      f2->SetLineWidth(2);
      hSvxEff_2B->Fit(f2, "R");
      f2->SetRange(0, 200);// To draw full range 
      cout << "---------------" << endl;
      cout << "Fit paraemters for double B case" << endl; 
      for (int i = 0; i < 5; i++) cout << f2->GetParameter(i) << std::endl;



    auto canvas = new TCanvas("canvas", "",1200, 1200);
      gStyle->SetOptStat(0);
      // hSvxEff_1B->SetLineColor(MyColors::kblue);
      // hSvxEff_2B->SetLineColor(MyColors::kred);
      hSvxEff_1B->SetLineColor(kBlue-2);
      hSvxEff_2B->SetLineColor(kRed-2);

      hSvxEff_1B->Draw();
      hSvxEff_2B->Draw("SAME");
      f1->Draw("same");
      f2->Draw("same");

      TLegend* leg = CreateLegend(0.35, 0.35, 0.8, 0.6,
                                      {hSvxEff_1B, hSvxEff_2B},
                                      {"LE", "LE"},
                                      {"1B", "2B"}
                                      );
                                      leg->SetHeader(Form(" %d < p_{T} < %d GeV", 100, 120), "L"); //centered 
                                      leg->Draw("same");
      canvas->SaveAs("Eff_SV/SVEff_1B_2B_atleast1trkin1SvxMatchfullB.png");
            canvas->SaveAs("Eff_SV/SVEff_1B_2B_atleast1trkin1SvxMatchfullB.root");


  TFile *fout = new TFile("Eff_SV/Eff_SV.root","recreate");
    // MC weights
    hBpt_1b->Write();
    hBpt_w1svxCut_1b->Write();
    hSvxEff_1B->Write();
    f1->Write("SvxEff_1B_pol5fit");

    hBpt_2b->Write();
    hBpt_w2svxCut_2b->Write();
    hSvxEff_2B->Write();
    f2->Write("SvxEff_2B_pol5fit");

    canvas->Write();

  fout->Close();
}

// What change in Eff. If I apply b tagging 
// Fit the historgams --> OK
// USe the Eff. to closure test on cutted histogra --> Does it give you bakc the other ones?
// Check if you use Matt. purity on SV labeling here and in the closure test? 



// Now test the fit fnction of Ef. when applied to hist w cut/Eff. --> does it give the total histogram ? 
void test_closure()
{
  gStyle->SetOptStat(0);
  ///////////  ---------------- Fast test: Test Eff. correction on the histogram level 
  TFile* finput = new TFile("Eff_SV/Eff_SV.root");
  // Until is produced: 
  TF1* f1 = (TF1*) finput->Get("SvxEff_1B_pol5fit");
  TF1* f2 = (TF1*) finput->Get("SvxEff_2B_pol5fit");
    // Use hard coded values: until the code end the run
        // TF1* f1 = new TF1("f1", "pol5", 0, 200);
        //   f1->SetParameters( -0.0558516, 0.0334072 , -0.000226583 , -3.16479e-06, 4.44619e-08,-1.48192e-10  );
          // with new fit to 200 : // f1->SetParameters(-0.0670951, 0.0387897, -0.000530821,  2.9222e-06, -6.24571e-09);
        // TF1* f2 = new TF1("f2", "pol5", 0, 200);
        //   f2->SetParameters( -0.0454185, 0.0450178, -0.000834719 , 7.58105e-06, -3.55534e-08 , 6.60356e-11  );

  // Test: correction oh hist(aftercut)/Eff. ? hist(beforecut)
    TH1D* hTotal_1B = (TH1D*) finput->Get("hBpt_1b");
    TH1D* hCut_1B   = (TH1D*) finput->Get("hBpt_w1svxCut_1b");
    TH1D* hTotal_2B = (TH1D*) finput->Get("hBpt_2b");
    TH1D* hCut_2B   = (TH1D*) finput->Get("hBpt_w2svxCut_2b");
    // hcut/Eff ?= hTotal? 
    TH1D* hClosure_1B = DivideByFunction(hCut_1B, f1, "hClosure_1B"); // hcut /Eff.
        hClosure_1B->SetTitle("hClosure (1B): 1b w SVx cut / Eff.");
    TH1D* hClosure_2B = DivideByFunction(hCut_2B, f2, "hClosure_2B");
       hClosure_2B->SetTitle("hClosure (2B): 2b w SVx cut / Eff.");

    // Compare the hist/Eff ?= hTtoal ? 
    auto canvas =  new TCanvas("canvas","Bin by bin level closure test", 1000, 1000);
      hTotal_1B->SetLineColor(kBlue-2); hTotal_1B ->SetMarkerStyle(20); hTotal_1B->SetMarkerColor(kBlue-2); hTotal_1B ->SetLineWidth(2);
      hTotal_2B->SetLineColor(kRed +1); hTotal_2B->SetMarkerStyle(20);  hTotal_2B ->SetMarkerColor(kRed +1); hTotal_2B ->SetLineWidth(2);
      hClosure_1B ->SetLineColor(kBlue-2); hClosure_1B->SetMarkerStyle(24); hClosure_1B->SetMarkerColor(kBlue-2);
      hClosure_2B ->SetLineColor(kRed+ 1); hClosure_2B->SetMarkerStyle(24); hClosure_2B->SetMarkerColor(kRed+ 1);
      canvas->cd();
      hTotal_1B->Draw("E");
      hTotal_2B->Draw("E SAME");
      hClosure_1B->Draw("E SAME");
      hClosure_2B->Draw("E SAME");
      canvas ->BuildLegend(0.35, 0.35, 0.65, 0.60);
      canvas->SetTitle("Bin level correction");
      canvas->Update();
      canvas->SaveAs("Eff_SV/FasctClosure_SvxEff.png");



   /////////  ---------------- Event level Closure test:
  TFile *fin = new TFile("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root");
  TTree *t = (TTree*)fin->Get("ak4PFJetAnalyzer/t");
  TTree *hie = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  // const int MAXJETS = 20;
  // const int MAXB = 20;
  // const int MAXTRACKS = 100;
  // Float_t weight;
  // Int_t nref;
  // Float_t jtpt[MAXJETS];
  // Float_t refpt[MAXJETS];
  // Float_t jteta[MAXJETS];
  // Int_t jtNbHad[MAXJETS];
  // Int_t ntrk;
  // Int_t trkJetId[MAXTRACKS];
  // Int_t trkMatchSta[MAXTRACKS];
  // Float_t trkBdtScore[MAXTRACKS];
  // Int_t trkSvtxId[MAXTRACKS];
  // Int_t nfullB;
  // Int_t fullBJetId[MAXB];
  // Int_t fullBSta[MAXB];
  // Float_t fullBPt[MAXB];
  hie->SetBranchAddress("weight",  &weight);
      SetBranches(t);
  // t->SetBranchAddress("nref",  &nref);
  // t->SetBranchAddress("jtpt",  jtpt);
  // t->SetBranchAddress("refpt",  refpt);
  // t->SetBranchAddress("jteta",  jteta);
  // t->SetBranchAddress("jtNbHad", jtNbHad);
  // t->SetBranchAddress("ntrk",  &ntrk);
  // t->SetBranchAddress("trkBdtScore",  trkBdtScore);
  // t->SetBranchAddress("trkJetId", trkJetId);
  // t->SetBranchAddress("trkMatchSta", trkMatchSta);
  // t->SetBranchAddress("trkSvtxId", trkSvtxId);
  // t->SetBranchAddress("nfullB",  &nfullB);
  // t->SetBranchAddress("fullBJetId",  fullBJetId);
  // t->SetBranchAddress("fullBSta",  fullBSta);
  // t->SetBranchAddress("fullBPt",  fullBPt);
  
  // --- corrected histograms with cut: does it recover the full histogram ? 
    // Correction with self eff. of same type 
  TH1D* hBpt_wcut_1b_self = new TH1D("hBpt_wcut_1b_self", "True 1B: self closure, event level ;p_{T, fullB} [GeV]; counts", 100, 0, 200);
  TH1D* hBpt_wcut_2b_self = new TH1D("hBpt_wcut_2b_self", "True 2B: self closure, event level ;p_{T, fullB} [GeV]; counts", 100, 0, 200);
    // correction with cross eff. of the other type 
  TH1D* hBpt_wcut_1b_cross = new TH1D("hBpt_wcut_1b_cross", "True 1B: cross closure, event level;p_{T, fullB} [GeV]; counts", 100, 0, 200);
  TH1D* hBpt_wcut_2b_cross = new TH1D("hBpt_wcut_2b_cross", "True 2B: cross closure, event level;p_{T, fullB} [GeV]; counts", 100, 0, 200);

  Long64_t nentries = t->GetEntries();
    // nentries = 1e+05;
  cout << "#entries = " << nentries << endl;
  for (Long64_t i = 0; i < nentries; i++) {
    if(i%100000==0) cout<<" i "<<i<<endl;
    t->GetEntry(i);
    hie->GetEntry(i);

    for (int j = 0; j < nref; j++) { // Jet loop 
      float jetPt = jtpt[j];
      // Jet Kineamtic cuts 
      if(jetPt < 100 || jetPt > 120) continue; // use nominal jet pt 
      if(jteta[j]<-2 || jteta[j]>2)  continue;
      if (refpt[j] <= 0) continue;
      // -- Double B jets and single B case 
      // loop over each true B
      for( int ib = 0; ib < nfullB; ib++) // per event
      {
        // For jet, have 2B --> each ID will be j
        if(fullBJetId[ib] != j) continue; // B in jet
        // Before SVX cut
        // if (jtNbHad[j] == 1) {hBpt_1b->Fill(fullBPt[ib], weight); } // All single B
        // if (jtNbHad[j] >= 2) {hBpt_2b->Fill(fullBPt[ib], weight); }// All Double B
        
        Int_t currentBSta = fullBSta[ib];

          // Find if at least one particle of same status as B belong to Svx in this jet
          for (int k = 0; k < ntrk; k++) // over tracks
          {
              if(trkJetId[k]!=j) continue;
              if(trkMatchSta[k]<1) continue;
              // Shall I cut trk pt > 1 
              // to pass: pass at least one svx requriment with one particle at least coming from the mother B
              if(trkMatchSta[k] == currentBSta)
              {
                if(trkSvtxId[k] >= 0 )
                {
                  // jtNsvtx [nref]
                  // svtxJetId [nsvtx] = jet index 
                  // where the svx id is trkSvtxId[k] which is from 0 to the totoal number of SVx: based on #svx per jet + counter on jets

                  ///--- APPLY EFF. correction here for closure test 
                  // if (jtNbHad[j] == 1) {hBpt_w1svxCut_1b->Fill(fullBPt[ib], weight); } 
                  // if (jtNbHad[j] >= 2) {hBpt_w2svxCut_2b->Fill(fullBPt[ib], weight); } 

                  double eff_true1b = f1->Eval(fullBPt[ib]);
                  double eff_true2b = f2->Eval(fullBPt[ib]);

                  if (jtNbHad[j] == 1)
                  {
                    if (eff_true1b > 0)   hBpt_wcut_1b_self->Fill(fullBPt[ib], weight/eff_true1b);
                    if (eff_true2b > 0)  hBpt_wcut_1b_cross ->Fill(fullBPt[ib], weight/eff_true2b);
                  } 
                  if (jtNbHad[j] >= 2)
                  {
                    if (eff_true2b > 0)  hBpt_wcut_2b_self->Fill(fullBPt[ib], weight/eff_true2b);
                    if (eff_true1b > 0)  hBpt_wcut_2b_cross ->Fill(fullBPt[ib], weight/eff_true1b);
                  } 

                  break; // go to next B.
                } 
              }
          }
      }// end B loop
  
    }// end jet loop
  } // end event loop


  ///- Draw Comparison: httoal Vs. the corrected histograms after cut
   auto canvas_eventl = new TCanvas("canvas_eventl","Event level closure test", 1000, 1000);
      hBpt_wcut_1b_self ->SetLineColor(kBlue-2); hBpt_wcut_1b_self->SetMarkerStyle(24); hBpt_wcut_1b_self->SetMarkerColor(kBlue-2); hBpt_wcut_1b_self->SetLineWidth(2);
      hBpt_wcut_2b_self ->SetLineColor(kRed+ 1); hBpt_wcut_2b_self->SetMarkerStyle(24); hBpt_wcut_2b_self->SetMarkerColor(kRed+ 1); hBpt_wcut_2b_self ->SetLineWidth(2);
      hBpt_wcut_1b_cross->SetMarkerStyle(25); hBpt_wcut_1b_cross->SetMarkerColor(kBlue-2); hBpt_wcut_1b_cross->SetLineColor(kBlue-2); hBpt_wcut_1b_cross ->SetLineWidth(2);
      hBpt_wcut_2b_cross->SetMarkerStyle(25); hBpt_wcut_2b_cross->SetMarkerColor(kRed+ 1); hBpt_wcut_2b_cross->SetLineColor(kRed+ 1); hBpt_wcut_2b_cross ->SetLineWidth(2);

      canvas_eventl->cd();
      hTotal_1B->Draw("E");
      hTotal_2B->Draw("E SAME");
      hBpt_wcut_1b_self->Draw("E SAME");
      hBpt_wcut_2b_self->Draw("E SAME");
      hBpt_wcut_1b_cross->Draw("E SAME");
      hBpt_wcut_2b_cross->Draw("E SAME");
      canvas_eventl ->BuildLegend(0.35, 0.35, 0.65, 0.60);
      canvas_eventl->SetTitle("Event Level correction");
      canvas_eventl->Update();
      canvas_eventl->SaveAs("Eff_SV/EventLevelClosure_SvxEff.png");

    //////// Write Output
    TFile *fout = new TFile("Eff_SV/Closure_Eff_SV.root","recreate"); 
        // Used inputs 
        hTotal_1B->Write();
        hTotal_2B->Write();
        f1->Write();
        f2->Write();

        // Fast closure
        hClosure_1B->Write();
        hClosure_2B->Write();
        canvas->Write();

        // Event Level Closure and Cross test
        hBpt_wcut_1b_self->Write();
        hBpt_wcut_2b_self->Write();
        hBpt_wcut_1b_cross->Write();
        hBpt_wcut_2b_cross->Write();
        canvas_eventl->Write();



}


//-- Now test Eff. correction of EEC at true level
void test_factorization()
{
  /////  Test Svx Eff. factorization at Gen level using EEC(2B): our analysis observable

  // -- Read Eff. (1B) 
  TFile* finput = new TFile("Eff_SV/Eff_SV.root");
  // Until is produced: 
  TF1* f1 = (TF1*) finput->Get("SvxEff_1B_pol5fit");
  // TF1* f2 = (TF1*) finput->Get("SvxEff_2B_pol5fit"); // will not be used 

    //test 
    // TF1* f1 = new TF1("f1", "pol5", 0, 200);
    // f1->SetParameters(-0.0670951, 0.0387897, -0.000530821,  2.9222e-06, -6.24571e-09); // new fit



  //// For EEC 
  TFile *fin = new TFile("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root");
  TTree *t = (TTree*)fin->Get("ak4PFJetAnalyzer/t");
  TTree *hie = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  hie->SetBranchAddress("weight",  &weight);
  SetBranches(t);
  
  // EEC histograms 
  // at gen level, without Svx cut, and using True B directly (no aggregation)
  TH1D *heec_true_2b = new TH1D("heec_2b", "Best pair;#DeltaR; EEC", bins_dr, dr_binsVector);

  // True level + NjetSvx >=2: this is not ok since this cut is reco level (#svx is reconstructed number)
  TH1D *heec_true_2b_jtsvxcut = new TH1D("heec_true_2b_jtsvxcut", ";#DeltaR; EEC", bins_dr, dr_binsVector);
  TH1D *heec_true_2b_jtsvxcut_corr = new TH1D("heec_true_2b_jtsvxcut_corr", ";#DeltaR; EEC", bins_dr, dr_binsVector);// doent go back t the original True level

  // Aggregated 2Bs at gen level, fullB used to get the Eff.(of SVx as function of its fullBpt)
    // Here Dr from aggragted Bs 
  TH1D *heec_aggreagted_gen_2b = new TH1D("heec_aggreagted_gen_2b", "Aggreagted 2Bs gen level;#DeltaR; EEC", bins_dr, dr_binsVector);
  TH1D *heec_aggreagted_gen_2b_corr = new TH1D("heec_aggreagted_gen_2b_corr", "Aggreagted 2Bs gen level (best pair) corrected w svx Eff 1B from fullB info ;#DeltaR; EEC", bins_dr, dr_binsVector);
    // What if you also use the corresponding True B ? 
      TH1D *heec_Afteraggreagted_trueDr_aggEEC_2b = new TH1D("heec_Afteraggreagted_trueDr_aggEEC_2b", "pass aggregtion condition (2Svx equivalent);#DeltaR; EEC", bins_dr, dr_binsVector); // true DR, agrregated eec 
      TH1D *heec_Afteraggreagted_trueDr_aggEEC_2b_corr = new TH1D("heec_Afteraggreagted_trueDr_aggEEC_2b_corr", "pass aggregtion condition (2Svx equivalent);#DeltaR; EEC", bins_dr, dr_binsVector); // correcred 
      // add another one for true Dr, true eec ? 
      TH1D *heec_Afteraggreagted_trueDrEEC_2b = new TH1D("heec_Afteraggreagted_trueDrEEC_2b", "pass aggregtion condition (2Svx equivalent);#DeltaR; EEC", bins_dr, dr_binsVector); // true DR, true eec 
      TH1D *heec_Afteraggreagted_trueDrEEC_2b_corr = new TH1D("heec_Afteraggreagted_trueDrEEC_2b_corr", "pass aggregtion condition (2Svx equivalent);#DeltaR; EEC", bins_dr, dr_binsVector); // correcred 
     


  Long64_t nentries = t->GetEntries();
  // nentries = 50e+06;
  // nentries = 10e+06;
  cout << "#entries = " << nentries << endl;

  for (Long64_t i = 0; i < nentries; i++) {
    
    if(i%100000==0) cout<<" i "<<i<<endl;
    
    t->GetEntry(i);
    hie->GetEntry(i);

    // cout << "nref = " << nref << endl;

    for (int j = 0; j < nref; j++) { // Jet loop 

      float jetPt = jtpt[j];
      // Jet Kineamtic cuts 
      if(jetPt < 100 || jetPt > 120) continue; // use nominal jet pt 
      if(jteta[j]<-2 || jteta[j]>2)  continue;
      if (refpt[j] <= 0) continue;
      
      // Use only 2B true level
      if ( jtNbHad[j]  < 2 ) continue; // can have >=2 true B

      // Get the kinematics of True Bs in vectors 
      std::vector <ROOT::Math::PtEtaPhiMVector> TrueB_vec;
      std::vector<Int_t> TrueB_status;
      TrueB_status.clear();
      TrueB_vec.clear();   // at the beginning of each event
      for( int ib = 0; ib < nfullB; ib++)
      {
        // Keep only Bs in this jet
        if (fullBJetId[ib] != j) continue;

         ROOT::Math::PtEtaPhiMVector v1;
          v1.SetPt(fullBPt[ib]);
          v1.SetEta(fullBEta[ib]);
          v1.SetPhi(fullBPhi[ib]);

          TrueB_vec.push_back(v1);
          TrueB_status.push_back(fullBSta[ib]);
      }
      // By default: I hsould have 2Bs 
      if (TrueB_vec.size() < 2) continue;  // why true plot is empty ?

      // cout << "Also full B info are consistent" << endl;

      // Compute eec 
      // pick best 2B pair
      // choose best 2Bs if more than 2 
      int best_true_i = 0, best_true_j = 1;
      double best_true_pt_prod = -1;
      for (size_t ti = 0; ti < TrueB_vec.size(); ti++)
          for (size_t tj = ti+1; tj < TrueB_vec.size(); tj++) {
              double pp = TrueB_vec[ti].Pt() * TrueB_vec[tj].Pt();
              if (pp > best_true_pt_prod) { best_true_pt_prod = pp; best_true_i = ti; best_true_j = tj; }
      }

      // -- EEC at gen level, without aggregation and no svx cut
        double dr_true   = ROOT::Math::VectorUtil::DeltaR(TrueB_vec[best_true_i], TrueB_vec[best_true_j]);
        double eec_true  = std::pow(TrueB_vec[best_true_i].Pt() * TrueB_vec[best_true_j].Pt(), 1);
        heec_true_2b ->Fill( dr_true, weight * eec_true);

        // hTrue --> Apply Eff. -> Get EEC as t supposed to be after SVx cut!  
        double eff1 = f1->Eval(TrueB_vec[best_true_i].Pt());
        double eff2 = f1->Eval(TrueB_vec[best_true_j].Pt());

      // -- Add explict cut on Njet SVx on True level
        if (jtNsvtx[j] >=2 )
        {
            heec_true_2b_jtsvxcut ->Fill( dr_true, weight * eec_true);
            heec_true_2b_jtsvxcut_corr ->Fill( dr_true, weight * eec_true / (eff1 * eff2)); // to be compared with True EEC
        }

      // -- Aggregated True Bs: call it gen level aggregation --> already need 2svx each of differet status. (No 2-0 topology at this level).
        std::vector<ROOT::Math::PtEtaPhiMVector> gen_bh;
          std::vector<Int_t> gen_bh_sta;
          PartialBsAggregation_globaltree(gen_bh, gen_bh_sta, j);
          if (gen_bh.size() < 2) continue;
          // choose best 2Bs if more than 2 
            int best_i = 0, best_j = 1;
            double best_pt_prod = -1;
            for (size_t gi = 0; gi < gen_bh.size(); gi++)
                for (size_t gj = gi+1; gj < gen_bh.size(); gj++) {
                    double pp = gen_bh[gi].Pt() * gen_bh[gj].Pt();
                    if (pp > best_pt_prod) { best_pt_prod = pp; best_i = gi; best_j = gj; }
            }
            double eec_gen_agg = std::pow(gen_bh[best_i].Pt() * gen_bh[best_j].Pt(), 1);
            double dr_gen_agg  = ROOT::Math::VectorUtil::DeltaR(gen_bh[best_i],  gen_bh[best_j] );
            // The Full B Pt and Eff. corresponding to choosen 2B pair: via status 
            double status_fullB1 = gen_bh_sta[best_i];
            double status_fullB2 = gen_bh_sta[best_j];
              // Get the fullB index that have the target aggreagted status
              int index1 = -1, index2 = -1;
              for (size_t k = 0; k < TrueB_vec.size(); k++)
              {
                if (TrueB_status[k] == status_fullB1 )  index1 = k;
                else if (TrueB_status[k] == status_fullB2 )  index2 = k;
                if (index1 != -1 && index2 != -1) break; // both found 
              } 
              // if (index1 > -1 && index2 > -1) { cout << "full B found" << endl;}
              double genB_eff1 = f1->Eval(TrueB_vec[index1].Pt());
              double genB_eff2 = f1->Eval(TrueB_vec[index2].Pt());
             
              // // what if using true kineamtics correponding of this aggreagted pair ?
               double eec_trueof_agg = std::pow(TrueB_vec[index1].Pt() * TrueB_vec[index2].Pt(), 1);
               double dr_trueof_agg  =  ROOT::Math::VectorUtil::DeltaR(TrueB_vec[index1], TrueB_vec[index2]);

              // Fill EEC(2B): aggreaged gen level, eff from full B
              heec_aggreagted_gen_2b ->Fill( dr_gen_agg, weight * eec_gen_agg);
              heec_aggreagted_gen_2b_corr ->Fill( dr_gen_agg, weight * eec_gen_agg / (genB_eff1 * genB_eff2) );
              // using True kinematics after aggreagtion (here aggreagtion as a SVx condition)
              heec_Afteraggreagted_trueDr_aggEEC_2b  ->Fill( dr_trueof_agg, weight * eec_gen_agg);
              heec_Afteraggreagted_trueDr_aggEEC_2b_corr  ->Fill( dr_trueof_agg, weight * eec_gen_agg  / (genB_eff1 * genB_eff2) );
              heec_Afteraggreagted_trueDrEEC_2b  ->Fill( dr_trueof_agg, weight * eec_trueof_agg);
              heec_Afteraggreagted_trueDrEEC_2b_corr  ->Fill( dr_trueof_agg, weight * eec_trueof_agg / (genB_eff1 * genB_eff2));


    }// end jet loop
  } // end event loop


  // Write output
  TFile *fout = new TFile("Eff_SV/Test_SvxEffFactorization_EEC_trueAggregTruecorr.root","recreate"); 
      heec_true_2b->Write();
      
      heec_true_2b_jtsvxcut->Write();
      heec_true_2b_jtsvxcut_corr ->Write();

      heec_aggreagted_gen_2b ->Write();
      heec_aggreagted_gen_2b_corr ->Write();

      heec_Afteraggreagted_trueDr_aggEEC_2b->Write();
      heec_Afteraggreagted_trueDr_aggEEC_2b_corr->Write();

      heec_Afteraggreagted_trueDrEEC_2b->Write();
      heec_Afteraggreagted_trueDrEEC_2b_corr->Write();
    
}



void ClosureTest_SV_eff()
{
  // main function
  // step1: compute Eff.
  // Get_1Svx_eff();
  // step2: test closure self and cross closures
  // test_closure();

  // step3: test fatorization
    test_factorization();
}

////////////////////////////////////////////////////////////////////
// Help function for hist divide by function bin by bin 
TH1D* DivideByFunction(TH1D* h, TF1* f, const char* name)
{
    TH1D* hnew = (TH1D*) h->Clone(name);

    for (int i = 1; i <= hnew->GetNbinsX(); i++) {

        double x   = hnew->GetBinCenter(i);
        double eff = f->Eval(x);

        if (eff <= 0) continue;

        hnew->SetBinContent(i,
                            hnew->GetBinContent(i)/eff);

        hnew->SetBinError(i,
                          hnew->GetBinError(i)/eff);
    }

    return hnew;
}
  
void SetBranches(TTree* t){

  t->SetBranchAddress("nref", &nref);
  t->SetBranchAddress("jtpt",  jtpt);
  t->SetBranchAddress("refpt",  refpt);
  t->SetBranchAddress("jteta",  jteta);
  t->SetBranchAddress("jtNbHad", jtNbHad);

  t->SetBranchAddress("ntrk",  &ntrk);
  t->SetBranchAddress("trkBdtScore",  trkBdtScore);
  t->SetBranchAddress("trkJetId", trkJetId);
  t->SetBranchAddress("trkMatchSta", trkMatchSta);
  t->SetBranchAddress("trkSvtxId", trkSvtxId);

  
  t->SetBranchAddress("nfullB",  &nfullB);
  t->SetBranchAddress("fullBJetId",  fullBJetId);
  t->SetBranchAddress("fullBSta",  fullBSta);
  t->SetBranchAddress("fullBPt",  fullBPt);
  // add fullB kinematics
  t->SetBranchAddress("fullBEta",  fullBEta);
  t->SetBranchAddress("fullBPhi",  fullBPhi);
  t->SetBranchAddress("fullBM",  fullBM);

  // Ref Tracks
  t->SetBranchAddress("nrefTrk", &nrefTrk);
  t->SetBranchAddress("refTrkJetId", refTrkJetId);
  t->SetBranchAddress("refTrkPdgId", refTrkPdgId);
  t->SetBranchAddress("refTrkSta", refTrkSta);

  t->SetBranchAddress("refTrkPt", refTrkPt);
  t->SetBranchAddress("refTrkPhi", refTrkPhi);
  t->SetBranchAddress("refTrkEta", refTrkEta);
  }

  void PartialBsAggregation_globaltree(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, Int_t ijet){
  ///// Difference wrt: PartialBsAggregation() --> this is using Global tree varaibles instead of tTree objects 

  // Only using gen level info
  // Add tacks of status >= 100 to make Bs.

  hadrons_4vec.clear();
  hadrons_stat.clear();   

  for (Int_t itrk = 0; itrk < nrefTrk; itrk++) {
    // Track must belong to this jet
    if (refTrkJetId[itrk] != ijet) continue;
    // pT cut
    if (refTrkPt[itrk] < 1) continue;                                                                                                  // Assign mass                                                                                                                  
    double mass = 0.0;                                                                                                               
    int pid = std::abs(refTrkPdgId[itrk]);
    if      (pid == 211)  mass = 0.139570;                                                                                           
    else if (pid == 13)   mass = 0.105658;                                                                                           
    else if (pid == 11)   mass = 0.000510;                                                                                           
    else if (pid == 2212) mass = 0.938272;                                                                                           
    else if (pid == 321)  mass = 0.493677;                                                                                           
    else if (pid == 3112) mass = 1.19744;                                                                                            
    else if (pid == 3222) mass = 1.18937;                                                                                           
    else if (pid == 3312) mass = 1.32171;                                                                                            
    else if (pid == 3334) mass = 1.67245;                                                                                            
    else std::cout << "PDG:" << pid << std::endl;
    // 4 vector                                                                                                                      
    ROOT::Math::PtEtaPhiMVector v(refTrkPt [itrk], refTrkEta[itrk], refTrkPhi[itrk], mass);
    
    Int_t status = refTrkSta[itrk];                                                                                                
    if (status < 100) continue;                                                                                                     
    else {
      auto it = std::find(hadrons_stat.begin(), hadrons_stat.end(), status);
      if (it == hadrons_stat.end()) {                                                                                                
        hadrons_stat.push_back(status);
        hadrons_4vec.push_back(v);                                                                                                   
      }                                                                                                                              
      else {                                                                                                                         
        size_t index = std::distance(hadrons_stat.begin(), it);                                                                      
        hadrons_4vec[index] += v;}
    }                                                                                                                                
  }                                                                                                                                  
}

/*
//-- Eff Fit result: f1 and f2 
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
Array of paraemters: ( -0.0558516, 0.0334072 , -0.000226583 , -3.16479e-06, 4.44619e-08,-1.48192e-10  )
****************************************
Minimizer is Linear / Migrad
Chi2                      =      204.074
NDf                       =           64
p0                        =   -0.0558516   +/-   0.00509306  
p1                        =    0.0334072   +/-   0.00150163  
p2                        = -0.000226583   +/-   7.27896e-05 
p3                        = -3.16479e-06   +/-   1.34395e-06 
p4                        =  4.44619e-08   +/-   1.06949e-08 
p5                        = -1.48192e-10   +/-   3.08108e-11 
****************************************

Paraemters array: ( -0.0454185, 0.0450178, -0.000834719 , 7.58105e-06, -3.55534e-08 , 6.60356e-11  )
Minimizer is Linear / Migrad
Chi2                      =      126.135
NDf                       =           64
p0                        =   -0.0454185   +/-   0.00368619  
p1                        =    0.0450178   +/-   0.00145107  
p2                        = -0.000834719   +/-   8.3187e-05  
p3                        =  7.58105e-06   +/-   1.76137e-06 
p4                        = -3.55534e-08   +/-   1.57243e-08 
p5                        =  6.60356e-11   +/-   4.97734e-11 
*/

/*
// Updated fit for 1b
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
****************************************
Minimizer is Linear / Migrad
Chi2                      =      272.706
NDf                       =           94
p0                        =   -0.0670951   +/-   0.00447904  
p1                        =    0.0387897   +/-   0.000971801 
p2                        = -0.000530821   +/-   3.62135e-05 
p3                        =   2.9222e-06   +/-   5.1396e-07  
p4                        = -6.24571e-09   +/-   3.1283e-09  
p5                        =  1.76937e-12   +/-   6.82053e-12 
---------------
Fit paraemters for single B case
-0.0670951
0.0387897
-0.000530821
2.9222e-06
-6.24571e-09
****************************************
Minimizer is Linear / Migrad
Chi2                      =      126.135
NDf                       =           64
p0                        =   -0.0454185   +/-   0.00368619  
p1                        =    0.0450178   +/-   0.00145107  
p2                        = -0.000834719   +/-   8.3187e-05  
p3                        =  7.58105e-06   +/-   1.76137e-06 
p4                        = -3.55534e-08   +/-   1.57243e-08 
p5                        =  6.60356e-11   +/-   4.97734e-11 
---------------
Fit paraemters for double B case
-0.0454185
0.0450178
-0.000834719
7.58105e-06
-3.55534e-08

*/