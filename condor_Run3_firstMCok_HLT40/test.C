void testRead(){                                                                                                                                                                                                                                       
  TFile *_file0 = TFile::Open("root://cms-xrd-global.cern.ch//store/user/mnguyen/bJetAggRun3/QCD_pThat-15to1200_TuneCP5_5p36TeV_pythia8/bJetAgg_2024PPRef_QCD/0/merged_HiForestMiniAOD.root");                                                         
  TTree *t = (TTree*)_file0->Get("ak4PFJetAnalyzer/t");                                                                                                                                                                                                
  cout<<"hello "<<t->GetEntries()<<endl;                                                                                                                                                                                                               
}
