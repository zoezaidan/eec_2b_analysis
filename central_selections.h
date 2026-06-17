
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------START: PART CAN GO TO COMMON HEADER --------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// -- Jet kinematics 
struct KinematicConfig {
  float ptLow;
  float ptHigh;
  float etaMax;
};

// -- Dataset
struct DatasetConfig {
  int RunN;        // 2 or 3
  int dataType;    // data / MC type (0: Data, 1: bjet MC, 2: qcd MC), and for Run2: -1 is for lowEG, and 0 for HighEG
  bool isMC;

  TString filename;
  TString output_folder = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/test_unfoldingcodewithTemplates/";
  TString output_hist;
  TString domain = ".root";

  // prescale
  double data_prescale;


  int fileindex = 0; // to be checked
};

// -- Physics selections: trigger and btagging 
struct PhysicsConfig {
  bool useBtag;
  double btagWP; // use -1 if you want default Run2 and Run3 values, otherwise set a custom WP value!  
};

// -- Full Analysis configuration structure 
struct AnalysisConfig {

  DatasetConfig dataset;
  KinematicConfig kin;
  PhysicsConfig physics;

  int n;

};

// -- Helper function for reco/gen MC choice of branches
float reco_jet_pt(const tTree& t, int ijet) {
  return t.jtpt[ijet];
}

float reco_jet_eta(const tTree& t, int ijet) {
  return t.jteta[ijet];
}
float gen_jet_pt(const tTree& t, int ijet) {
  return t.refpt[ijet];
}

float gen_jet_eta(const tTree& t, int ijet) {
  return t.refeta[ijet];
}

// ----------- Selection function for each configuration choice 
// -- event selection based on trigger 
bool passEventSelection(const tTree& t,
                        const AnalysisConfig& cfg) { 

  int RunN = cfg.dataset.RunN;
  bool isMC = cfg.dataset.isMC;
  int dataType = cfg.dataset.dataType;

  if (RunN == 2) {

    if (!isMC && dataType == 0) // HighEG data 
      return (t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1);

    if (!isMC && dataType == -1) // LowEG data 
      return !(t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1) &&
             (t.HLT_HIAK4PFJet40_v1 || t.HLT_HIAK4PFJet60_v1);

    if (isMC)
      return t.HLT_HIAK4PFJet40_v1; // minHLT in Run2
  }

  if (RunN == 3) { // minHLT is used at 60 GeV. 40GeV is heavily prescaled --> not used.

    if (!isMC)
      return (t.HLT_AK4PFJet60_v8 ||
              t.HLT_AK4PFJet80_v8 ||
              t.HLT_AK4PFJet100_v8);

    if (isMC)
      return t.HLT_AK4PFJet60_v8;
  }

  return true;
}

// -- Jet kinematic selection: seperated for Gen/Reco + cleaned events of large weights based on jetpt
bool passGenJetKinematics(const tTree& t,
                          int ijet,
                          const AnalysisConfig& cfg) {
  
  // This function have check for refpt > 0
  // It return false if the event pthat cut is not satisfied (which is based on reco jtpt, yes, reco pt not gen pt!).
  
  // ---------------- MC EVENT CLEANING ----------------
  if (cfg.dataset.isMC) {
    if (t.pthat < 0.40 * t.jtpt[ijet]) // reco jet is used for this cut
      return false;
  }

  float pt  = gen_jet_pt(t, ijet);
  float eta = gen_jet_eta(t, ijet);

  // -- Safety:
  if (pt < 0) return false;

  if (fabs(eta) > cfg.kin.etaMax) return false;
  if (pt < cfg.kin.ptLow || pt > cfg.kin.ptHigh) return false;

  return true;
}

bool passRecoJetKinematics(const tTree& t,
                           int ijet,
                           const AnalysisConfig& cfg) {

  // ---------------- MC EVENT CLEANING ----------------
  if (cfg.dataset.isMC) {
    if (t.pthat < 0.40 * t.jtpt[ijet]) // reco jet is used for this cut
      return false;
  }

  float pt  = reco_jet_pt(t, ijet);
  float eta = reco_jet_eta(t, ijet);

  if (fabs(eta) > cfg.kin.etaMax) return false;
  if (pt < cfg.kin.ptLow || pt > cfg.kin.ptHigh) return false;

  return true;
}

// -- Jet Btagging selection 
bool passBtag(const tTree& t,
              int ijet,
              const AnalysisConfig& cfg) {

  if (!cfg.physics.useBtag) return true;// safety if btag = false, and this function is used.

  double btagVar = -999.;

  if (cfg.dataset.RunN == 2) {
    btagVar = t.discr_particleNet_BvsAll[ijet];
  }
  else if (cfg.dataset.RunN == 3) {
    btagVar =
      t.discr_unifiedParticleTransformer_probb[ijet] +
      t.discr_unifiedParticleTransformer_problepb[ijet] +
      t.discr_unifiedParticleTransformer_probbb[ijet];
  }
  else{
      std::cerr << "Unknown RunN for b tagging= "
                << cfg.dataset.RunN << std::endl;
      return false;
  }

  return (btagVar > cfg.physics.btagWP);
}

// -- Function: set the Dataset and output templates names 
DatasetConfig buildDataset(int RunN, int dataType, bool isMC, const PhysicsConfig& physics) {

  // Set the rest of DatasetConfig information (Input file name + output hist for templates)

  DatasetConfig d;
  d.RunN = RunN;
  d.dataType = dataType;
  d.isMC = isMC;

  // -- Set pre-calculated prescales 
  if (RunN == 2){
      d.data_prescale = 33.917210; // prescale 40-60 GeV
    }
    else if (RunN == 3){
      d.data_prescale = 6.2336493; // prescale 60-80 GeV
    }

  TString add_BtagWP = physics.useBtag ? Form("_btagWP%d", static_cast<int>(std::round(1000 * physics.btagWP))):""; // round to int without floating digits
  TString sample = "";

  if (RunN == 2) {

    if (dataType == -1) {
      d.filename =  "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root";
      sample = "LowEG";
    }

    if (dataType == 0) {
      d.filename =  "/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      sample = "HighEG";
    }

    if (dataType == 1) {
      d.filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      sample = "bjet";
    }
    if (dataType == 2) {
      d.filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      sample = "qcd";
    }
  }

  if (RunN == 3) {

    if (dataType == 0) {
      d.filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/HiForestMiniAOD_v2_TChains.root"; // Chaine (available working as tree)
      sample = "data";
    }

    if (dataType == 1) {
      d.filename = ""; // Not available yet
      sample = "bjet";
    }
    if (dataType == 2) {
      Int_t fileindex = 0; // temporare set
      d.filename = Form("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_recalJP_chunks_20260617/merged_block_000%d_Pythia8_recalJP_20260617.root", fileindex); // [0-9]
      sample = "qcd";
    }

    d.output_hist = Form("Run%d%s_template_for_fit_histos_3D_%s_f",RunN, add_BtagWP.Data(), sample.Data());
  }


  return d;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------END: PART CAN GO TO COMMON HEADER --------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// --- HERE: USER: Factory function of analysis configurations: set all configurations to use ---
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AnalysisConfig buildConfig(
    int RunN,
    int dataType,
    float ptLow,
    float ptHigh,
    float etaCut,
    int n,
    bool btag,
    bool isMC,
    double btagWP
){
  // -- USER: SET HERE YOUR DESIRED ANALYSIS CONFIG -------

  AnalysisConfig cfg;
  // Kinematics
  cfg.kin.ptLow = ptLow;
  cfg.kin.ptHigh = ptHigh;
  cfg.kin.etaMax = etaCut;

  // Physics 
  cfg.physics.useBtag = btag;
  if (btagWP > 0) cfg.physics.btagWP = btagWP; // custom value 
  else{ // default value 
    cfg.physics.btagWP = (RunN == 2) ? 0.898 : 0.872;
  }

  cfg.n = n;

  // Dataset setting
  cfg.dataset = buildDataset(RunN, dataType, isMC, cfg.physics); 

  return cfg;
}

/////////////////////////////////////
// ----------- Centralize activated branches based on need from configuration ----------------
std::vector<TString> getActiveBranches(const AnalysisConfig& cfg)
{
    std::vector<TString> branches = {

        "jtpt",
        "jteta",
        "nref",

        "jtNtrk",
        "jtNsvtx",

        "ntrk",
        "trkJetId",
        "trkBdtScore",
        "trkPdgId",
        "trkMatchPdgId",
        "trkMatchSta",
        "trkPt",
        "trkEta",
        "trkPhi",

        "trkSvtxId",

        "nsvtx",
        "svtxJetId",
        "svtxNtrk",
        "svtxm",
        "svtxmcorr",
        "svtxpt",

        "svtxdl",
        "svtxdls",
        "svtxdl2d",
        "svtxdls2d",
        "svtxnormchi2"
    };

    if (cfg.dataset.RunN == 2) {

        branches.insert(branches.end(), {

            "HLT_HIAK4PFJet40_v1",
            "HLT_HIAK4PFJet60_v1",
            "HLT_HIAK4PFJet80_v1",
            "HLT_HIAK4PFJet100_v1",

            "discr_particleNet_BvsAll"
        });
    }

    else if (cfg.dataset.RunN == 3) {

        branches.insert(branches.end(), {

            "HLT_AK4PFJet60_v8",
            "HLT_AK4PFJet80_v8",
            "HLT_AK4PFJet100_v8",

            "discr_unifiedParticleTransformer_probb",
            "discr_unifiedParticleTransformer_problepb",
            "discr_unifiedParticleTransformer_probbb"
        });
    }

    if (cfg.dataset.isMC) {

        branches.insert(branches.end(), {
            "weight",
            "pthat",
            "jtNbHad"
        });
    }

    return branches;
}




