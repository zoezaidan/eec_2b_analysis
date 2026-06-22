#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include <iostream>

void drawTrigEff(){


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = TFile::Open("merged_HiForestMiniAOD_v2.root");

   TTree *HltTree = (TTree*) f->Get("hltanalysis/HltTree");
   TTree *skimTree = (TTree*) f->Get("skimanalysis/HltTree");
   TTree *t = (TTree*) f->Get("ak4PFJetAnalyzer/t");
   TTree *hie = (TTree*) f->Get("hiEvtAnalyzer/HiTree");

//Declaration of leaves types
   Int_t jet40, jet60, jet80, jet100, jet120;

   // Set branch addresses.
   HltTree->SetBranchAddress("HLT_AK4PFJet40_v8",&jet40);
   HltTree->SetBranchAddress("HLT_AK4PFJet60_v8",&jet60);
   HltTree->SetBranchAddress("HLT_AK4PFJet80_v8",&jet80);
   HltTree->SetBranchAddress("HLT_AK4PFJet100_v8",&jet100);
   HltTree->SetBranchAddress("HLT_AK4PFJet120_v8",&jet120);

   Int_t nref;
   Float_t jtpt[20], jteta[20], btagb[20], btagbb[20], btaglepb[20];  

   t->SetBranchAddress("nref",&nref);
   t->SetBranchAddress("jtpt",jtpt);
   t->SetBranchAddress("jteta",jteta);
   t->SetBranchAddress("discr_unifiedParticleTransformer_probb",btagb);
   t->SetBranchAddress("discr_unifiedParticleTransformer_probbb",btagbb);
   t->SetBranchAddress("discr_unifiedParticleTransformer_problepb",btaglepb);
   
   Int_t           pprimaryVertexFilter;
   skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);

   Float_t vz;
   hie->SetBranchAddress("vz",&vz);

   
   Long64_t nentries = HltTree->GetEntries();


   TH1F *h40 = new TH1F("h40","h40",500,0,500);
   TH1F *h60 = new TH1F("h60","h60",500,0,500);
   TH1F *h80 = new TH1F("h80","h80",500,0,500);
   TH1F *h100 = new TH1F("h100","h100",500,0,500);
   TH1F *h120 = new TH1F("h120","h120",500,0,500);

   TH1F *hb40 = new TH1F("hb40","hb40",500,0,500);
   TH1F *hb60 = new TH1F("hb60","hb60",500,0,500);
   TH1F *hb80 = new TH1F("hb80","hb80",500,0,500);
   TH1F *hb100 = new TH1F("hb100","hb100",500,0,500);
   TH1F *hb120 = new TH1F("hb120","hb120",500,0,500);

   
   // # of jets firing 120
   Long64_t n120=0;
   
   // # of jets firing other triggers AND 120
   Long64_t nc40=0;
   Long64_t nc60=0;
   Long64_t nc80=0;
   Long64_t nc100=0;
   
   
   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
     if(i%1000000==0) cout<<" i "<<i<<" nentries "<<nentries<<endl;     
     
     nbytes += HltTree->GetEntry(i);
     if(jet120){
       n120++;       
       if(jet100)nc100++;
       if(jet80)nc80++;
       if(jet60)nc60++;
       if(jet40)nc40++;
     }       
     
     t->GetEntry(i);

     skimTree->GetEntry(i);
     if(!pprimaryVertexFilter) continue;

     if(vz < -24. || vz > 24.) continue;
     
     for(int j=0; j< nref;j++){
       if(jteta[j] < -2 ||jteta[j] > 2) continue;

       if(jet120)h120->Fill(jtpt[j]);
       else if(jet100)h100->Fill(jtpt[j]);
       else if(jet80)h80->Fill(jtpt[j]);
       else if(jet60)h60->Fill(jtpt[j]);
       else if(jet40)h40->Fill(jtpt[j]);

       float btagScore = btagb[j]+btagbb[j]+btaglepb[j];
       if(btagScore<0.872) continue;

       if(jet120)hb120->Fill(jtpt[j]);
       else if(jet100)hb100->Fill(jtpt[j]);
       else if(jet80)hb80->Fill(jtpt[j]);
       else if(jet60)hb60->Fill(jtpt[j]);
       else if(jet40)hb40->Fill(jtpt[j]);

       
     }
     
     
   }

   cout<<" n120 "<<n120<<endl;
   cout<<" nc100 "<<nc100<<" nc100 / n120 "<<n120/(float)nc100<<endl;
   cout<<" nc80 "<<nc80<<" nc80 / n120 "<<n120/(float)nc80<<endl;
   cout<<" nc60 "<<nc60<<" nc60 / n60 "<<n120/(float)nc60<<endl;
   cout<<" nc40 "<<nc40<<" nc40 / n120 "<<n120/(float)nc40<<endl;


   TFile *fout = new TFile("trigEffHistos.root","recreate");
   h40->Write();
   h60->Write();
   h80->Write();
   h100->Write();
   h120->Write();
   hb40->Write();
   hb60->Write();
   hb80->Write();
   hb100->Write();
   hb120->Write();
   fout->Close();
}
