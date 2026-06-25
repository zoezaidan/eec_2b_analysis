
/*
  Studies b-decay track BDT vs jet pT as well as (full, gen) b-hadron pT
  b-decay tracks are associated only to _their_ parent b-hadron
  This code is a useful example for correlating a track to a given b-hadron
  -Matt

 */

void drawBDTvsPt(){

  
  TFile *fin = new TFile("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root");
  TTree *t = (TTree*)fin->Get("ak4PFJetAnalyzer/t");
  TTree *hie = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  

  const int MAXJETS = 20;
  const int MAXB = 20;
  const int MAXTRACKS = 100;
 
  Float_t weight;

  Int_t nref;
  Float_t jtpt[MAXJETS];
  Float_t refpt[MAXJETS];
  Float_t jteta[MAXJETS];
  Int_t jtNbHad[MAXJETS];

  Int_t ntrk;
  Int_t trkJetId[MAXTRACKS];
  Int_t trkMatchSta[MAXTRACKS];
  Float_t trkBdtScore[MAXTRACKS];
  
  Int_t nfullB;
  Int_t fullBJetId[MAXB];
  Int_t fullBSta[MAXB];
  Float_t fullBPt[MAXB];
  
  hie->SetBranchAddress("weight",  &weight);
  t->SetBranchAddress("nref",  &nref);
  t->SetBranchAddress("jtpt",  jtpt);
  t->SetBranchAddress("refpt",  refpt);
  t->SetBranchAddress("jteta",  jteta);
  t->SetBranchAddress("jtNbHad", jtNbHad);

  t->SetBranchAddress("ntrk",  &ntrk);
  t->SetBranchAddress("trkBdtScore",  trkBdtScore);
  t->SetBranchAddress("trkJetId", trkJetId);
  t->SetBranchAddress("trkMatchSta", trkMatchSta);
  
  t->SetBranchAddress("nfullB",  &nfullB);
  t->SetBranchAddress("fullBJetId",  fullBJetId);
  t->SetBranchAddress("fullBSta",  fullBSta);
  t->SetBranchAddress("fullBPt",  fullBPt);
  

  TH3F *h2B_sig_bdt_jtpt_bpt = new TH3F("h2B_sig_bdt_jtpt_bpt",";BDT;p_{T,jet};p_{T,B}",10, 0., 1., 6,80,200, 100,0.,200.);
  TH3F *h1B_sig_bdt_jtpt_bpt = new TH3F("h1B_sig_bdt_jtpt_bpt",";BDT;p_{T,jet};p_{T,B}",10, 0., 1., 6,80,200, 100,0.,200.);
  TH2F *h2B_bkd_bdt_jtpt = new TH2F("h2B_bkd_bdt_jtpt",";BDT;p_{T,jet}",10, 0., 1.,  6,80,200);
  TH2F *h1B_bkd_bdt_jtpt = new TH2F("h1B_bkd_bdt_jtpt",";BDT;p_{T,jet}",10, 0., 1.,  6,80,200);
  TH2F *h0B_bkd_bdt_jtpt = new TH2F("h0B_bkd_bdt_jtpt",";BDT;p_{T,jet}",10, 0., 1.,  6,80,200);
  
  
  Long64_t nentries = t->GetEntries();

  for (Long64_t i = 0; i < nentries; i++) {
    
    if(i%100000==0) cout<<" i "<<i<<endl;
    
    t->GetEntry(i);
    hie->GetEntry(i);

    for (int j = 0; j < nref; j++) {

      //if(jtNbHad[j]==0) continue;

      float jetPt = jtpt[j];
      
      if(jetPt < 80 || jetPt > 200) continue;
      if(jteta[j]<-2 || jteta[j]>2)  continue;
      if (refpt[j] <= 0) continue;
      
      for (int k = 0; k < ntrk; k++) {
	if(trkJetId[k]!=j) continue;
	if(trkMatchSta[k]<1) continue;
	if(jtNbHad[j] == 0 && trkMatchSta[k]>1) continue;

	float discr = trkBdtScore[k];	

	if(trkMatchSta[k]==1){
	  if(jtNbHad[j]==0)h0B_bkd_bdt_jtpt->Fill(discr,jetPt,weight); 
	  else if(jtNbHad[j]==1)h1B_bkd_bdt_jtpt->Fill(discr,jetPt,weight);
	  else h2B_bkd_bdt_jtpt->Fill(discr,jetPt,weight);
	}
	else { // means sta > 1	  
	  for( int ib = 0; ib < nfullB; ib++){
	    if(fullBJetId[ib] != j) continue;
	    if(fullBSta[ib] != trkMatchSta[k]) continue;
	    
	    float bpt = fullBPt[ib];
	    
	    if(jtNbHad[j]==1)h1B_sig_bdt_jtpt_bpt->Fill(discr,jetPt,bpt,weight);
	    else h2B_sig_bdt_jtpt_bpt->Fill(discr,jetPt,bpt,weight);
	  }
	}
      }
            
    }
  }	 
  
  TFile *fout = new TFile("BDTvsPt.root","recreate");
  h2B_sig_bdt_jtpt_bpt->Write();
  h1B_sig_bdt_jtpt_bpt->Write();
  h2B_bkd_bdt_jtpt->Write();
  h1B_bkd_bdt_jtpt->Write();
  h0B_bkd_bdt_jtpt->Write();
  fout->Close();
}
