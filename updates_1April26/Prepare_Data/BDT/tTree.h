#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <TH3D.h>
#include <Math/Vector4D.h>
#include <Math/Vector4Dfwd.h>
#include <Math/VectorUtil.h>

// Header file for the classes stored in the TTree if any.

class tTree {
public :
   TString         fname;
   TTree           *tree;

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           lumi;
   Int_t           nref;
   Int_t           nvtx;
   Float_t         rawpt[500];   //[nref]
   Float_t         jtpt[500];   //[nref]
   Float_t         jtpt_gen[500];   //[nref]
   Float_t         jteta[500];   //[nref]
   Float_t         jty[500];   //[nref]
   Float_t         jtphi[500];   //[nref]
   Float_t         jtpu[500];   //[nref]
   Float_t         jtm[500];   //[nref]
   Float_t         jtarea[500];   //[nref]
   Float_t         jtPfCHF[500];   //[nref]
   Float_t         jtPfNHF[500];   //[nref]
   Float_t         jtPfCEF[500];   //[nref]
   Float_t         jtPfNEF[500];   //[nref]
   Float_t         jtPfMUF[500];   //[nref]
   Int_t           jtPfCHM[500];   //[nref]
   Int_t           jtPfNHM[500];   //[nref]
   Int_t           jtPfCEM[500];   //[nref]
   Int_t           jtPfNEM[500];   //[nref]
   Int_t           jtPfMUM[500];   //[nref]
   Float_t         jttau1[500];   //[nref]
   Float_t         jttau2[500];   //[nref]
   Float_t         jttau3[500];   //[nref]
   Int_t           jtNtrk[500];   //[nref]
   Float_t         discr_particleNet_BvsAll[500];
   Int_t           ntrk;
   Int_t           trkJetId[500];   //[ntrk]
   Int_t           trkSvtxId[500];   //[ntrk]
   Int_t           trkPdgId[500];    //[ntrk]
   Float_t         trkPt[500];   //[ntrk]
   Float_t         trkEta[500];   //[ntrk]
   Float_t         trkPhi[500];   //[ntrk]
   Float_t         trkY[500];   //[ntrk]
   Float_t         trkIp3d[500];   //[ntrk]
   Float_t         trkIp3dSig[500];   //[ntrk]
   Float_t         trkDistToAxisSig[500];   //[ntrk]
   Float_t         trkDistToAxis[500];   //[ntrk]
   Int_t           trkMatchSta[500];   //[ntrk]
   Float_t         trkBdtScore[500]; 
   Int_t           jtNsvtx[500];   //[nref]
   Float_t         trkMass[500];
   Int_t           nsvtx;
   Int_t           svtxJetId[50];   //[nsvtx]
   Int_t           svtxNtrk[50];   //[nsvtx]
   Float_t         svtxdl[50];   //[nsvtx]
   Float_t         svtxdls[50];   //[nsvtx]
   Float_t         svtxdl2d[50];   //[nsvtx]
   Float_t         svtxdls2d[50];   //[nsvtx]
   Float_t         svtxm[50];   //[nsvtx]
   Float_t         svtxpt[50];   //[nsvtx]
   Float_t         svtxmcorr[50];   //[nsvtx]
   Float_t         svtxnormchi2[50];   //[nsvtx]

   //new
   Int_t           nrefTrk;
   Int_t           refTrkJetId[500];   //[ntrk]
   Int_t           refTrkPdgId[500];    //[ntrk]
   Float_t         refTrkPt[500];   //[ntrk]
   Float_t         refTrkEta[500];   //[ntrk]
   Float_t         refTrkPhi[500];   //[ntrk]
   Float_t         refTrkY[500];   //[ntrk]
   Int_t           refTrkSta[500];   //[ntrk]
   Float_t         refTrkMass[500];

   //HLT selection
   Int_t           HLT_HIAK4PFJet100_v1;
   Int_t           HLT_HIAK4PFJet80_v1;
   Int_t           HLT_HIAK4PFJet60_v1;
   Int_t           HLT_HIAK4PFJet40_v1;
   Int_t           HLT_HIAK4PFJet30_v1;

   //Prescale
   //Double_t        prescale_pf40;


   Int_t           ntrkInSvtxNotInJet;
   Int_t           trkInSvtxNotInJetSvId[500];
   Int_t           trkInSvtxNotInJetOtherJetId[500];
   Int_t           trkInSvtxNotInJetMatchSta[500];
   Float_t         trkInSvtxNotInJetPt[500];
   Float_t         trkInSvtxNotInJetEta[500];
   Float_t         trkInSvtxNotInJetPhi[500];

   // aod compatibility
   Float_t         jtDiscDeepFlavourB[500];   //[nref]
   Float_t         jtDiscDeepFlavourBB[500];   //[nref]
   Float_t         jtDiscDeepFlavourLEPB[500];   //[nref]

   Float_t         discr_deepCSV[500];   //[nref]
   Float_t         discr_pfJP[500];   //[nref]
   Float_t         discr_deepFlavour_b[500];   //[nref]
   Float_t         discr_deepFlavour_bb[500];   //[nref]
   Float_t         discr_deepFlavour_lepb[500];   //[nref]
   Float_t         pthat;
   Float_t         refpt[500];   //[nref]
   Float_t         refeta[500];   //[nref]
   Float_t         refy[500];   //[nref]
   Float_t         refphi[500];   //[nref]
   Float_t         refm[500];   //[nref]
   Float_t         refarea[500];   //[nref]
   Float_t         refdphijt[500];   //[nref]
   Float_t         refdrjt[500];   //[nref]
   Float_t         refparton_pt[500];   //[nref]
   Int_t           refparton_flavor[500];   //[nref]
   Int_t           refparton_flavorForB[500];   //[nref]
   Float_t         genChargedSum[500];   //[nref]
   Float_t         genHardSum[500];   //[nref]
   Float_t         signalChargedSum[500];   //[nref]
   Float_t         signalHardSum[500];   //[nref]
   Int_t           ngen;
   Int_t           genmatchindex[100];   //[ngen]
   Float_t         genpt[100];   //[ngen]
   Float_t         geneta[100];   //[ngen]
   Float_t         geny[100];   //[ngen]
   Float_t         genphi[100];   //[ngen]
   Float_t         genm[100];   //[ngen]
   Float_t         gendphijt[100];   //[ngen]
   Float_t         gendrjt[100];   //[ngen]

   // True flavour
   Int_t           jtHadFlav[500];   //[nref]
   Int_t           jtParFlav[500];   //[nref]
   Int_t           jtNbHad[500];   //[nref]
   Int_t           jtNcHad[500];   //[nref]
   Int_t           jtNbPar[500];   //[nref]
   Int_t           jtNcPar[500];   //[nref]

   // aod compatibility
//    Float_t jtHadFlav[500]; //[nref]
//    Float_t jtParFlav[500]; //[nref]

   // Subjets
   Float_t         sjt1Pt[500];
   Float_t         sjt1Eta[500];
   Float_t         sjt1Phi[500];
   Float_t         sjt1Y[500];

   Float_t         sjt2Pt[500];
   Float_t         sjt2Eta[500];
   Float_t         sjt2Phi[500];
   Float_t         sjt2Y[500];

   Float_t         rsjt1Pt[500];
   Float_t         rsjt1Eta[500];
   Float_t         rsjt1Phi[500];
   Float_t         rsjt1Y[500];

   Float_t         rsjt2Pt[500];
   Float_t         rsjt2Eta[500];
   Float_t         rsjt2Phi[500];
   Float_t         rsjt2Y[500];

   Float_t         jtmB[500]; //[nref]
   Float_t         jtBpt[500]; //[nref]
   Float_t         jtptCh[500]; //[nref]
   Float_t         refmB[500]; //[nref]
   Float_t         refBpt[500]; //[nref]
   Float_t         refptCh[500]; //[nref]
   Int_t           refNtrk[500]; //[nref]

   Float_t         weight;


   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jtpt_gen;   //!
   TBranch        *b_discr_particleNet_BvsAll;
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_jtarea;   //!
   TBranch        *b_jtPfCHF;   //!
   TBranch        *b_jtPfNHF;   //!
   TBranch        *b_jtPfCEF;   //!
   TBranch        *b_jtPfNEF;   //!
   TBranch        *b_jtPfMUF;   //!
   TBranch        *b_jtPfCHM;   //!
   TBranch        *b_jtPfNHM;   //!
   TBranch        *b_jtPfCEM;   //!
   TBranch        *b_jtPfNEM;   //!
   TBranch        *b_jtPfMUM;   //!
   TBranch        *b_jttau1;   //!
   TBranch        *b_jttau2;   //!
   TBranch        *b_jttau3;   //!
   TBranch        *b_jtNtrk;   //!
   TBranch        *b_ntrk;   //!
   TBranch        *b_trkJetId;   //!
   TBranch        *b_trkSvtxId;   //!
   TBranch        *b_trkPdgId;
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkY;   //!
   TBranch        *b_trkIp3d;   //!
   TBranch        *b_trkIp3dSig;   //!
   TBranch        *b_trkDistToAxisSig;   //!
   TBranch        *b_trkDistToAxis;   //!
   TBranch        *b_trkMatchSta;   //!
   TBranch        *b_trkBdtScore;
   TBranch        *b_jtNsvtx;   //!
   TBranch        *b_nsvtx;   //!
   TBranch        *b_svtxJetId;   //!
   TBranch        *b_svtxNtrk;   //!
   TBranch        *b_svtxdl;   //!
   TBranch        *b_svtxdls;   //!
   TBranch        *b_svtxdl2d;   //!
   TBranch        *b_svtxdls2d;   //!
   TBranch        *b_svtxm;   //!
   TBranch        *b_svtxpt;   //!
   TBranch        *b_svtxmcorr;   //!
   TBranch        *b_svtxnormchi2;   //!

  
   TBranch        *b_refTrkJetId;
   TBranch        *b_refTrkPdgId;
   TBranch        *b_refTrkPt;
   TBranch        *b_refTrkEta;
   TBranch        *b_refTrkPhi;
   TBranch        *b_refTrkY;
   TBranch        *b_refTrkSta;
   TBranch        *b_refTrkMass;
   TBranch        *b_trkMass;
   TBranch        *b_nrefTrk;

   //HLT branches
   TBranch        *b_HLT_HIAK4PFJet100_v1;
   TBranch        *b_HLT_HIAK4PFJet80_v1;
   TBranch        *b_HLT_HIAK4PFJet60_v1;
   TBranch        *b_HLT_HIAK4PFJet40_v1;
   TBranch        *b_HLT_HIAK4PFJet30_v1;

   //Prescale
   //TBranch        *b_prescale_pf40;


   TBranch        *b_ntrkInSvtxNotInJet;   //!
   TBranch        *b_trkInSvtxNotInJetSvId;   //!
   TBranch        *b_trkInSvtxNotInJetOtherJetId;   //!
   TBranch        *b_trkInSvtxNotInJetMatchSta;   //!
   TBranch        *b_trkInSvtxNotInJetPt;   //!
   TBranch        *b_trkInSvtxNotInJetEta;   //!
   TBranch        *b_trkInSvtxNotInJetPhi;   //!

   // aod compatibily
   TBranch        *b_jtDiscDeepFlavourB;   //!
   TBranch        *b_jtDiscDeepFlavourBB;   //!
   TBranch        *b_jtDiscDeepFlavourLEPB;   //!

   TBranch        *b_discr_deepCSV;   //!
   TBranch        *b_discr_pfJP;   //!
   TBranch        *b_discr_deepFlavour_b;   //!
   TBranch        *b_discr_deepFlavour_bb;   //!
   TBranch        *b_discr_deepFlavour_lepb;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refm;   //!
   TBranch        *b_refarea;   //!
   TBranch        *b_refdphijt;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refparton_pt;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_refparton_flavorForB;   //!
   TBranch        *b_genChargedSum;   //!
   TBranch        *b_genHardSum;   //!
   TBranch        *b_signalChargedSum;   //!
   TBranch        *b_signalHardSum;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genmatchindex;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_geny;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_genm;   //!
   TBranch        *b_gendphijt;   //!
   TBranch        *b_gendrjt;   //!

   TBranch        *b_jtHadFlav;   //!
   TBranch        *b_jtParFlav;   //!
   TBranch        *b_jtNbHad;   //!
   TBranch        *b_jtNcHad;   //!
   TBranch        *b_jtNbPar;   //!
   TBranch        *b_jtNcPar;   //!

   TBranch        *b_sjt1Pt;
   TBranch        *b_sjt1Eta;
   TBranch        *b_sjt1Phi;
   TBranch        *b_sjt1Y;

   TBranch        *b_sjt2Pt;
   TBranch        *b_sjt2Eta;
   TBranch        *b_sjt2Phi;
   TBranch        *b_sjt2Y;

   TBranch        *b_rsjt1Pt;
   TBranch        *b_rsjt1Eta;
   TBranch        *b_rsjt1Phi;
   TBranch        *b_rsjt1Y;

   TBranch        *b_rsjt2Pt;
   TBranch        *b_rsjt2Eta;
   TBranch        *b_rsjt2Phi;
   TBranch        *b_rsjt2Y;


   TBranch        *b_jtmB;   //!
   TBranch        *b_jtBpt;   //!
   TBranch        *b_jtptCh;   //!
   TBranch        *b_refmB;   //!
   TBranch        *b_refBpt;   //!
   TBranch        *b_refptCh;   //!
   TBranch        *b_refNtrk;   //!

   TBranch        *b_weight;   //!

   tTree(TString rootf);
   ~tTree();
   Int_t GetEntry(Long64_t entry);
   Long64_t GetEntries();
   void Init(TTree *tree);
   void SetBranchStatus(TString branchName, Int_t status);
   void SetBranchStatus(std::vector<TString> branchNames, Int_t status);
   void plot_rgzgkt(TString foutname, Float_t bTagWP);
   Float_t calc_dr(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);
   Float_t calc_rg(Float_t y1, Float_t phi1, Float_t y2, Float_t phi2);
};

tTree::tTree(TString rootf)
{
   TFile *fin = TFile::Open(rootf);
   tree = (TTree*) fin->Get("ak4PFJetAnalyzer/t");
   tree->AddFriend("hiEvtAnalyzer/HiTree");
   tree->AddFriend("hltanalysis/HltTree");

   //Data
   /*tree = (TTree*) fin->Get("akCs4PFJetAnalyzer/t");
   tree->AddFriend("hltanalysis/HltTree");
   tree->AddFriend("hiEvtAnalyzer/HiTree");*/
   
   Init(tree);
}

tTree::~tTree()
{
   if (!tree) return;
   delete tree;
}

Int_t tTree::GetEntry(Long64_t entry)
{
   if (!tree) return 0;
   return tree->GetEntry(entry);
}

Long64_t tTree::GetEntries()
{
   if (!tree) return 0;
   return tree->GetEntries();
}

void tTree::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;

   tree->SetBranchAddress("run", &run, &b_run);
   tree->SetBranchAddress("evt", &evt, &b_evt);
   tree->SetBranchAddress("lumi", &lumi, &b_lumi);
   tree->SetBranchAddress("nref", &nref, &b_nref);
   tree->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   tree->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   tree->SetBranchAddress("jtpt", jtpt, &b_jtpt);
   tree->SetBranchAddress("jtpt_gen", jtpt_gen, &b_jtpt_gen);
   tree->SetBranchAddress("discr_particleNet_BvsAll", discr_particleNet_BvsAll, &b_discr_particleNet_BvsAll);
   tree->SetBranchAddress("jteta", jteta, &b_jteta);
   tree->SetBranchAddress("jty", jty, &b_jty);
   tree->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   tree->SetBranchAddress("jtpu", jtpu, &b_jtpu);
   tree->SetBranchAddress("jtm", jtm, &b_jtm);
   tree->SetBranchAddress("jtarea", jtarea, &b_jtarea);
   tree->SetBranchAddress("jtPfCHF", jtPfCHF, &b_jtPfCHF);
   tree->SetBranchAddress("jtPfNHF", jtPfNHF, &b_jtPfNHF);
   tree->SetBranchAddress("jtPfCEF", jtPfCEF, &b_jtPfCEF);
   tree->SetBranchAddress("jtPfNEF", jtPfNEF, &b_jtPfNEF);
   tree->SetBranchAddress("jtPfMUF", jtPfMUF, &b_jtPfMUF);
   tree->SetBranchAddress("jtPfCHM", jtPfCHM, &b_jtPfCHM);
   tree->SetBranchAddress("jtPfNHM", jtPfNHM, &b_jtPfNHM);
   tree->SetBranchAddress("jtPfCEM", jtPfCEM, &b_jtPfCEM);
   tree->SetBranchAddress("jtPfNEM", jtPfNEM, &b_jtPfNEM);
   tree->SetBranchAddress("jtPfMUM", jtPfMUM, &b_jtPfMUM);
   tree->SetBranchAddress("jttau1", jttau1, &b_jttau1);
   tree->SetBranchAddress("jttau2", jttau2, &b_jttau2);
   tree->SetBranchAddress("jttau3", jttau3, &b_jttau3);
   tree->SetBranchAddress("jtNtrk", jtNtrk, &b_jtNtrk);
   tree->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
   tree->SetBranchAddress("trkJetId", trkJetId, &b_trkJetId);
   tree->SetBranchAddress("trkSvtxId", trkSvtxId, &b_trkSvtxId);
   tree->SetBranchAddress("trkPdgId", trkPdgId, &b_trkPdgId);
   tree->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   tree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   tree->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   tree->SetBranchAddress("trkY", trkY, &b_trkY);
   tree->SetBranchAddress("trkIp3d", trkIp3d, &b_trkIp3d);
   tree->SetBranchAddress("trkIp3dSig", trkIp3dSig, &b_trkIp3dSig);
   tree->SetBranchAddress("trkDistToAxisSig", trkDistToAxisSig, &b_trkDistToAxisSig);
   tree->SetBranchAddress("trkDistToAxis", trkDistToAxis, &b_trkDistToAxis);
   tree->SetBranchAddress("trkMatchSta", trkMatchSta, &b_trkMatchSta);
   tree->SetBranchAddress("trkBdtScore", trkBdtScore, &b_trkBdtScore); //
   tree->SetBranchAddress("jtNsvtx", jtNsvtx, &b_jtNsvtx);
   tree->SetBranchAddress("nsvtx", &nsvtx, &b_nsvtx);
   tree->SetBranchAddress("svtxJetId", svtxJetId, &b_svtxJetId);
   tree->SetBranchAddress("svtxNtrk", svtxNtrk, &b_svtxNtrk);
   tree->SetBranchAddress("svtxdl", svtxdl, &b_svtxdl);
   tree->SetBranchAddress("svtxdls", svtxdls, &b_svtxdls);
   tree->SetBranchAddress("svtxdl2d", svtxdl2d, &b_svtxdl2d);
   tree->SetBranchAddress("svtxdls2d", svtxdls2d, &b_svtxdls2d);
   tree->SetBranchAddress("svtxm", svtxm, &b_svtxm);
   tree->SetBranchAddress("svtxpt", svtxpt, &b_svtxpt);
   tree->SetBranchAddress("svtxmcorr", svtxmcorr, &b_svtxmcorr);
   tree->SetBranchAddress("svtxnormchi2", svtxnormchi2, &b_svtxnormchi2);

   tree->SetBranchAddress("ntrkInSvtxNotInJet", &ntrkInSvtxNotInJet, &b_ntrkInSvtxNotInJet);  
   tree->SetBranchAddress("trkInSvtxNotInJetSvId", trkInSvtxNotInJetSvId, &b_trkInSvtxNotInJetSvId);  
   tree->SetBranchAddress("trkInSvtxNotInJetOtherJetId", trkInSvtxNotInJetOtherJetId, &b_trkInSvtxNotInJetOtherJetId);  
   tree->SetBranchAddress("trkInSvtxNotInJetMatchSta", trkInSvtxNotInJetMatchSta, &b_trkInSvtxNotInJetMatchSta);  
   tree->SetBranchAddress("trkInSvtxNotInJetPt", trkInSvtxNotInJetPt, &b_trkInSvtxNotInJetPt);  
   tree->SetBranchAddress("trkInSvtxNotInJetEta", trkInSvtxNotInJetEta, &b_trkInSvtxNotInJetEta);  
   tree->SetBranchAddress("trkInSvtxNotInJetPhi", trkInSvtxNotInJetPhi, &b_trkInSvtxNotInJetPhi);  

   // aod compatibility
   tree->SetBranchAddress("jtDiscDeepFlavourB", jtDiscDeepFlavourB, &b_jtDiscDeepFlavourB);
   tree->SetBranchAddress("jtDiscDeepFlavourBB", jtDiscDeepFlavourBB, &b_jtDiscDeepFlavourBB);
   tree->SetBranchAddress("jtDiscDeepFlavourLEPB", jtDiscDeepFlavourLEPB, &b_jtDiscDeepFlavourLEPB);

   tree->SetBranchAddress("discr_deepCSV", discr_deepCSV, &b_discr_deepCSV);
   tree->SetBranchAddress("discr_pfJP", discr_pfJP, &b_discr_pfJP);
   tree->SetBranchAddress("discr_deepFlavour_b", discr_deepFlavour_b, &b_discr_deepFlavour_b);
   tree->SetBranchAddress("discr_deepFlavour_bb", discr_deepFlavour_bb, &b_discr_deepFlavour_bb);
   tree->SetBranchAddress("discr_deepFlavour_lepb", discr_deepFlavour_lepb, &b_discr_deepFlavour_lepb);
   tree->SetBranchAddress("pthat", &pthat, &b_pthat);
   tree->SetBranchAddress("refpt", refpt, &b_refpt);
   tree->SetBranchAddress("refeta", refeta, &b_refeta);
   tree->SetBranchAddress("refy", refy, &b_refy);
   tree->SetBranchAddress("refphi", refphi, &b_refphi);
   tree->SetBranchAddress("refm", refm, &b_refm);
   tree->SetBranchAddress("refarea", refarea, &b_refarea);
   tree->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
   tree->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
   tree->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
   tree->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
   tree->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
   tree->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
   tree->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
   tree->SetBranchAddress("signalChargedSum", signalChargedSum, &b_signalChargedSum);
   tree->SetBranchAddress("signalHardSum", signalHardSum, &b_signalHardSum);
   tree->SetBranchAddress("ngen", &ngen, &b_ngen);
   tree->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
   tree->SetBranchAddress("genpt", genpt, &b_genpt);
   tree->SetBranchAddress("geneta", geneta, &b_geneta);
   tree->SetBranchAddress("geny", geny, &b_geny);
   tree->SetBranchAddress("genphi", genphi, &b_genphi);
   tree->SetBranchAddress("genm", genm, &b_genm);
   tree->SetBranchAddress("gendphijt", gendphijt, &b_gendphijt);
   tree->SetBranchAddress("gendrjt", gendrjt, &b_gendrjt);

   tree->SetBranchAddress("jtHadFlav", jtHadFlav, &b_jtHadFlav);
   tree->SetBranchAddress("jtParFlav", jtParFlav, &b_jtParFlav);
   tree->SetBranchAddress("jtNbHad", jtNbHad, &b_jtNbHad);
   tree->SetBranchAddress("jtNcHad", jtNcHad, &b_jtNcHad);
   tree->SetBranchAddress("jtNbPar", jtNbPar, &b_jtNbPar);
   tree->SetBranchAddress("jtNcPar", jtNcPar, &b_jtNcPar);

   tree->SetBranchAddress("sjt1Pt", sjt1Pt, &b_sjt1Pt);
   tree->SetBranchAddress("sjt1Eta", sjt1Eta, &b_sjt1Eta);
   tree->SetBranchAddress("sjt1Phi", sjt1Phi, &b_sjt1Phi);
   tree->SetBranchAddress("sjt1Y", sjt1Y, &b_sjt1Y);

   tree->SetBranchAddress("sjt2Pt", sjt2Pt, &b_sjt2Pt);
   tree->SetBranchAddress("sjt2Eta", sjt2Eta, &b_sjt2Eta);
   tree->SetBranchAddress("sjt2Phi", sjt2Phi, &b_sjt2Phi);
   tree->SetBranchAddress("sjt2Y", sjt2Y, &b_sjt2Y);

   tree->SetBranchAddress("rsjt1Pt", rsjt1Pt, &b_rsjt1Pt);
   tree->SetBranchAddress("rsjt1Eta", rsjt1Eta, &b_rsjt1Eta);
   tree->SetBranchAddress("rsjt1Phi", rsjt1Phi, &b_rsjt1Phi);
   tree->SetBranchAddress("rsjt1Y", rsjt1Y, &b_rsjt1Y);

   tree->SetBranchAddress("rsjt2Pt", rsjt2Pt, &b_rsjt2Pt);
   tree->SetBranchAddress("rsjt2Eta", rsjt2Eta, &b_rsjt2Eta);
   tree->SetBranchAddress("rsjt2Phi", rsjt2Phi, &b_rsjt2Phi);
   tree->SetBranchAddress("rsjt2Y", rsjt2Y, &b_rsjt2Y);

   tree->SetBranchAddress("jtmB", jtmB, &b_jtmB);
   tree->SetBranchAddress("jtBpt", jtBpt, &b_jtBpt);
   tree->SetBranchAddress("jtptCh", jtptCh, &b_jtptCh);
   tree->SetBranchAddress("refmB", refmB, &b_refmB);
   tree->SetBranchAddress("refBpt", refBpt, &b_refBpt);
   tree->SetBranchAddress("refptCh", refptCh, &b_refptCh);
   tree->SetBranchAddress("refNtrk", refNtrk, &b_refNtrk);

   tree->SetBranchAddress("weight", &weight, &b_weight);

   //new
   tree->SetBranchAddress("refTrkJetId", refTrkJetId, &b_refTrkJetId);
   tree->SetBranchAddress("refTrkPdgId", refTrkPdgId, &b_refTrkPdgId);
   tree->SetBranchAddress("refTrkPt", refTrkPt, &b_refTrkPt);
   tree->SetBranchAddress("refTrkEta", refTrkEta, &b_refTrkEta);
   tree->SetBranchAddress("refTrkPhi", refTrkPhi, &b_refTrkPhi);
   tree->SetBranchAddress("refTrkY", refTrkY, &b_refTrkY);
   tree->SetBranchAddress("refTrkSta", refTrkSta, &b_refTrkSta);
   tree->SetBranchAddress("refTrkMass", refTrkMass, &b_refTrkMass);
   tree->SetBranchAddress("trkMass", trkMass, &b_trkMass);
   tree->SetBranchAddress("nrefTrk", &nrefTrk, &b_nrefTrk);

   //HLT
   tree->SetBranchAddress("HLT_HIAK4PFJet100_v1", &HLT_HIAK4PFJet100_v1, &b_HLT_HIAK4PFJet100_v1);
   tree->SetBranchAddress("HLT_HIAK4PFJet80_v1", &HLT_HIAK4PFJet80_v1, &b_HLT_HIAK4PFJet80_v1);
   tree->SetBranchAddress("HLT_HIAK4PFJet60_v1", &HLT_HIAK4PFJet60_v1, &b_HLT_HIAK4PFJet60_v1);
   tree->SetBranchAddress("HLT_HIAK4PFJet40_v1", &HLT_HIAK4PFJet40_v1, &b_HLT_HIAK4PFJet40_v1);
   tree->SetBranchAddress("HLT_HIAK4PFJet30_v1", &HLT_HIAK4PFJet30_v1, &b_HLT_HIAK4PFJet30_v1);

   //Prescale
   //tree->SetBranchAddress("prescale_pf40", &prescale_pf40, &b_prescale_pf40);

}

void tTree::SetBranchStatus(TString branchName, Int_t status)
{
    tree->SetBranchStatus(branchName, status);
}

void tTree::SetBranchStatus(vector<TString> branchNames, Int_t status)
{
    for (TString branchName : branchNames) {
        tree->SetBranchStatus(branchName, status);
    }
}

void tTree::plot_rgzgkt(TString foutname, Float_t bTagWP = 0.9)
{
    // Activate relevant branches
    tree->SetBranchStatus("*", 0);
    for (auto activeBranchName : {"nref", "refeta", "refpt", "jtHadFlav",
                                  "rsjt1Pt", "rsjt1Eta", "rsjt1Phi", "rsjt1Y",
                                  "rsjt2Pt", "rsjt2Eta", "rsjt2Phi", "rsjt2Y",
                                  "jtpt", "jteta",
                                  "sjt1Pt", "sjt1Eta", "sjt1Phi", "sjt1Y",
                                  "sjt2Pt", "sjt2Eta", "sjt2Phi", "sjt2Y",
                                  "discr_deepFlavour_b", "discr_deepFlavour_bb", "discr_deepFlavour_lepb",
                                  "weight"
                                  }) {
        tree->SetBranchStatus(activeBranchName, 1);
    }

    // Create histograms

    // ln(0.4/rg)
    Int_t x1bins = 10;
    Float_t x1min = 0.;
    Float_t x1max = 2.5;

    // rg
    // Int_t x1bins = 10;
    // Float_t x1min = 0.;
    // Float_t x1max = 0.4;

    // ln(kt)
    Int_t y1bins = 40;
    Float_t y1min = -5.;
    Float_t y1max = 5.;

    // zg
    Int_t y2bins = 10;
    Float_t y2min = 0.1;
    Float_t y2max = 0.5;

    // refpt
    Int_t z1bins = 27*2;
    Float_t z1min = 30.;
    Float_t z1max = 300.;

    // reco level
    TH3D *hB_drkt = new TH3D("hB_drkt", "dr, kt, pt, bjets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hB_rgkt = new TH3D("hB_rgkt", "rg, kt, pt, bjets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hB_rgzg = new TH3D("hB_rgzg", "rg, zg, pt, bjets", x1bins, x1min, x1max, y2bins, y2min, y2max, z1bins, z1min, z1max);
    TH3D *hB_zgkt = new TH3D("hB_zgkt", "zg, kt, pt, bjets", y2bins, y2min, y2max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    TH3D *hBtag_drkt = new TH3D("hBtag_drkt", "dr, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hBtag_rgkt = new TH3D("hBtag_rgkt", "rg, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hBtag_rgzg = new TH3D("hBtag_rgzg", "rg, zg, pt, b tagged jets", x1bins, x1min, x1max, y2bins, y2min, y2max, z1bins, z1min, z1max);
    TH3D *hBtag_zgkt = new TH3D("hBtag_zgkt", "zg, kt, pt, b tagged jets", y2bins, y2min, y2max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    TH3D *hBtagNoBB_drkt = new TH3D("hBtagNoBB_drkt", "dr, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hBtagNoBB_rgkt = new TH3D("hBtagNoBB_rgkt", "rg, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    // gen level
    TH3D *hB_drkt_gen = new TH3D("hB_drkt_gen", "dr, kt, pt, bjets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hB_rgkt_gen = new TH3D("hB_rgkt_gen", "rg, kt, pt, bjets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hB_rgzg_gen = new TH3D("hB_rgzg_gen", "rg, zg, pt, bjets", x1bins, x1min, x1max, y2bins, y2min, y2max, z1bins, z1min, z1max);
    TH3D *hB_zgkt_gen = new TH3D("hB_zgkt_gen", "zg, kt, pt, bjets", y2bins, y2min, y2max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    TH3D *hBtag_drkt_gen = new TH3D("hBtag_drkt_gen", "dr, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hBtag_rgkt_gen = new TH3D("hBtag_rgkt_gen", "rg, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hBtag_rgzg_gen = new TH3D("hBtag_rgzg_gen", "rg, zg, pt, b tagged jets", x1bins, x1min, x1max, y2bins, y2min, y2max, z1bins, z1min, z1max);
    TH3D *hBtag_zgkt_gen = new TH3D("hBtag_zgkt_gen", "zg, kt, pt, b tagged jets", y2bins, y2min, y2max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    TH3D *hL_drkt_gen = new TH3D("hL_drkt_gen", "dr, kt, pt, guds jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hL_rgkt_gen = new TH3D("hL_rgkt_gen", "rg, kt, pt, guds jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hL_rgzg_gen = new TH3D("hL_rgzg_gen", "rg, zg, pt, guds jets", x1bins, x1min, x1max, y2bins, y2min, y2max, z1bins, z1min, z1max);
    TH3D *hL_zgkt_gen = new TH3D("hL_zgkt_gen", "zg, kt, pt, guds jets", y2bins, y2min, y2max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    TH3D *hC_drkt_gen = new TH3D("hC_drkt_gen", "dr, kt, pt, cjets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hC_rgkt_gen = new TH3D("hC_rgkt_gen", "rg, kt, pt, cjets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hC_rgzg_gen = new TH3D("hC_rgzg_gen", "rg, zg, pt, cjets", x1bins, x1min, x1max, y2bins, y2min, y2max, z1bins, z1min, z1max);
    TH3D *hC_zgkt_gen = new TH3D("hC_zgkt_gen", "zg, kt, pt, cjets", y2bins, y2min, y2max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    TH3D *hSingleBtag_rgkt_gen = new TH3D("hSingleBtag_rgkt_gen", "rg, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);
    TH3D *hSingleBtag_drkt_gen = new TH3D("hSingleBtag_drkt_gen", "dr, kt, pt, b tagged jets", x1bins, x1min, x1max, y1bins, y1min, y1max, z1bins, z1min, z1max);

    const Float_t jetR = 0.4;

    std::cout << "Creating histograms ..." << std::endl;
    Long64_t nentries = tree->GetEntries();
    for (Long64_t ient = 0; ient < nentries; ient++) {
        // Print progress
        if (ient % 1000000 == 0) {
            std::cout << "ient = " << ient << std::endl;
        }

        // Choose nb of events
        const Long64_t total_events = 1000000;
        if (ient > total_events) break;
		// if (ient > 100) break;
 
        tree->GetEntry(ient);

        for (Int_t ijet = 0; ijet < nref; ijet++) {
            // universal eta cut
            if (std::abs(refeta[ijet]) > 2) continue;

            Float_t dr = -1.;
            Float_t rg = -1.;
            Float_t kt = -1.;
            Float_t zg = -1.;

            Float_t logdr = -1.;
            Float_t logrg = -1.;
            Float_t logkt = -10.;

            Float_t dr_gen = -1.;
            Float_t rg_gen = -1.;
            Float_t kt_gen = -1.;
            Float_t zg_gen = -1.;

            Float_t logdr_gen = -1.;
            Float_t logrg_gen = -1.;
            Float_t logkt_gen = -10.;

            // Calculate rg, kt only for 2 prong jets
            if (sjt2Pt[ijet] > 0.) {
                dr = calc_dr(sjt1Eta[ijet], sjt1Phi[ijet], sjt2Eta[ijet], sjt2Phi[ijet]);
                rg = calc_rg(sjt1Y[ijet], sjt1Phi[ijet], sjt2Y[ijet], sjt2Phi[ijet]);
                kt = sjt2Pt[ijet] * rg;
                zg = sjt2Pt[ijet] / (sjt1Pt[ijet] + sjt2Pt[ijet]);
                
                // calculate logs
                logdr = std::log(jetR/dr);
                logrg = std::log(jetR/rg);
                logkt = std::log(kt);

                // for poster
                // logrg = rg;
            }

            if (rsjt2Pt[ijet] > 0.) {
                dr_gen = calc_dr(rsjt1Eta[ijet], rsjt1Phi[ijet], rsjt2Eta[ijet], rsjt2Phi[ijet]);
                rg_gen = calc_rg(rsjt1Y[ijet], rsjt1Phi[ijet], rsjt2Y[ijet], rsjt2Phi[ijet]);
                kt_gen = rsjt2Pt[ijet] * rg_gen;
                zg_gen = rsjt2Pt[ijet] / (rsjt1Pt[ijet] + rsjt2Pt[ijet]);
                
                // calculate logs
                logdr_gen = std::log(jetR/dr_gen);
                logrg_gen = std::log(jetR/rg_gen);
                logkt_gen = std::log(kt_gen);

                // for poster
                // logrg_gen = rg_gen;
            }

            // Fill true-flavour histograms
            bool isBjet = (jtHadFlav[ijet] == 5);
            if (isBjet) {
               hB_drkt->Fill(logdr, logkt, refpt[ijet], weight);
               hB_rgkt->Fill(logrg, logkt, refpt[ijet], weight);
               hB_rgzg->Fill(logrg, zg, refpt[ijet], weight);
               hB_zgkt->Fill(zg, logkt, refpt[ijet], weight);
            
               hB_drkt_gen->Fill(logdr_gen, logkt_gen, refpt[ijet], weight);
               hB_rgkt_gen->Fill(logrg_gen, logkt_gen, refpt[ijet], weight);
               hB_rgzg_gen->Fill(logrg_gen, zg_gen, refpt[ijet], weight);
               hB_zgkt_gen->Fill(zg_gen, logkt_gen, refpt[ijet], weight);
            

               // Fill the b-tag histogram
               bool passWP = ((discr_deepFlavour_b[ijet] + discr_deepFlavour_bb[ijet] + discr_deepFlavour_lepb[ijet]) > bTagWP);
            //    bool passWP = ((discr_deepFlavour_b[ijet] + discr_deepFlavour_lepb[ijet]) > 0.7);
               if (passWP) {
                    hBtag_drkt->Fill(logdr, logkt, refpt[ijet], weight);
                    hBtag_rgkt->Fill(logrg, logkt, refpt[ijet], weight);
                    hBtag_rgzg->Fill(logrg, zg, refpt[ijet], weight);
                    hBtag_zgkt->Fill(zg, logkt, refpt[ijet], weight);

                    hBtag_drkt_gen->Fill(logdr_gen, logkt_gen, refpt[ijet], weight);
                    hBtag_rgkt_gen->Fill(logrg_gen, logkt_gen, refpt[ijet], weight);
                    hBtag_rgzg_gen->Fill(logrg_gen, zg_gen, refpt[ijet], weight);
                    hBtag_zgkt_gen->Fill(zg_gen, logkt_gen, refpt[ijet], weight);  
               }

               if (passWP && (discr_deepFlavour_bb[ijet] < 0.2)) {
                    hBtagNoBB_drkt->Fill(logdr, logkt, refpt[ijet], weight);
                    hSingleBtag_drkt_gen->Fill(logdr_gen, logkt_gen, refpt[ijet], weight);

                    hBtagNoBB_rgkt->Fill(logrg, logkt, refpt[ijet], weight);
                    hSingleBtag_rgkt_gen->Fill(logrg_gen, logkt_gen, refpt[ijet], weight);
               }
            } // end if is b jet

            bool isLightJet = (jtHadFlav[ijet] == 0);
            if (isLightJet) {
                hL_drkt_gen->Fill(logdr_gen, logkt_gen, refpt[ijet], weight);
                hL_rgkt_gen->Fill(logrg_gen, logkt_gen, refpt[ijet], weight);
                hL_rgzg_gen->Fill(logrg_gen, zg_gen, refpt[ijet], weight);
                hL_zgkt_gen->Fill(zg_gen, logkt_gen, refpt[ijet], weight);
            }

            bool isCJet = (jtHadFlav[ijet] == 4);
            if (isCJet) {
                hC_drkt_gen->Fill(logdr_gen, logkt_gen, refpt[ijet], weight);
                hC_rgkt_gen->Fill(logrg_gen, logkt_gen, refpt[ijet], weight);
                hC_rgzg_gen->Fill(logrg_gen, zg_gen, refpt[ijet], weight);
                hC_zgkt_gen->Fill(zg_gen, logkt_gen, refpt[ijet], weight);
            }
        } // jet loop
    } // entry loop

    // Save histograms
    std::cout << "\n(Re)creating file " << foutname << std::endl;
    TFile *fout = new TFile(foutname, "recreate");

    for (auto h : {hB_drkt, hB_rgkt, hB_rgzg, hB_zgkt, 
                   hBtag_drkt, hBtag_rgkt, hBtag_rgzg, hBtag_zgkt,
                   hB_drkt_gen, hB_rgkt_gen, hB_rgzg_gen, hB_zgkt_gen, 
                   hBtag_drkt_gen, hBtag_rgkt_gen, hBtag_rgzg_gen, hBtag_zgkt_gen,
                   hL_drkt_gen, hL_rgkt_gen, hL_rgzg_gen, hL_zgkt_gen, 
                   hC_drkt_gen, hC_rgkt_gen, hC_rgzg_gen, hC_zgkt_gen, 
                   hBtagNoBB_drkt,
                   hSingleBtag_drkt_gen,
                   hBtagNoBB_rgkt,
                   hSingleBtag_rgkt_gen}) {
        h->Write();
    }

    fout->Close();
}

Float_t tTree::calc_dr(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
    ROOT::Math::PtEtaPhiMVector v1;
    v1.SetPt(1.);
    v1.SetEta(eta1);
    v1.SetPhi(phi1);

    ROOT::Math::PtEtaPhiMVector v2;
    v2.SetPt(1.);
    v2.SetEta(eta2);
    v2.SetPhi(phi2);

    Float_t dr = ROOT::Math::VectorUtil::DeltaR(v1, v2);
    return dr;
}

Float_t tTree::calc_rg(Float_t y1, Float_t phi1, Float_t y2, Float_t phi2) {
    ROOT::Math::PtEtaPhiMVector v1;
    v1.SetPhi(phi1);

    ROOT::Math::PtEtaPhiMVector v2;
    v2.SetPhi(phi2);

    Float_t dphi = ROOT::Math::VectorUtil::DeltaPhi(v1, v2);
    Float_t dy = y1 - y2;

    // Float_t dphi_test = std::acos(std::cos(phi1-phi2));

    // std::cout << "dphi = " << dphi << std::endl;
    // std::cout << "dphi_test = " << dphi_test << std::endl;
    Float_t rg = std::sqrt((dy*dy) + (dphi*dphi));
    return rg;
}