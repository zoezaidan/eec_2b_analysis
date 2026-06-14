// Entry point called by each condor job.
// Computes its event range from job_idx and n_jobs, then calls make_templates.
//
// Usage (from run_job.sh):
//   root -b -q -l 'run_condor_job.C+(job_idx, n_jobs, dataType, pT_low, pT_high, n, btag_int, isMC_int, fileindex)'

#include "../create_files_for_template_fit_Run3.cpp"

void run_condor_job_Run3_wproxy(Int_t job_idx, Int_t n_jobs, Int_t dataType,
                    Float_t pT_low, Float_t pT_high, Int_t n,
                    Int_t btag_int, Int_t isMC_int, Int_t fileindex) {


  std::cout << "ENTERED run_condor_job_Run3()" << std::endl;
  std::cout << "JOBID = " << job_idx << std::endl;


  bool btag = (btag_int != 0);
  bool isMC = (isMC_int != 0);

  TString filename, output_hist;
  TString output_folder = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3/JobResult/";
  TString domain = ".root";

  if (dataType == -1) {
    filename     = "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root";
    output_hist  = "template_for_fit_histos_3D_LowEG";
    isMC = false;
  } else if (dataType == 0) {
    filename     = "/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    output_hist  = "template_for_fit_histos_3D_HighEG";
    isMC = false;
  } else if (dataType == 1) {
    filename     = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    output_hist  = "template_for_fit_histos_3D_bjet_test";
  } else if (dataType == 2) {

    // -- Run 3 MC 
    // TFile *f = TFile::Open("root://cms-xrd-global.cern.ch//store/user/mnguyen/.../merged_HiForestMiniAOD.root");
    // filename = Form("root://cms-xrd-global.cern.ch//store/user/mnguyen//bJetAggRun3/QCD_pThat-15to1200_TuneCP5_5p36TeV_pythia8/bJetAgg_2024PPRef_QCD/%d/merged_HiForestMiniAOD.root", fileindex); 
    // test read locally 
    // From Mat.. path
    // filename = Form("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/%d/merged_HiForestMiniAOD.root", fileindex);
    filename = Form("root://cms-xrd-global.cern.ch//store/user/mnguyen//bJetAggRun3/QCD_pThat-15to1200_TuneCP5_5p36TeV_pythia8/bJetAgg_2024PPRef_QCD/%d/merged_HiForestMiniAOD.root", fileindex); 
    output_hist = Form("Run3_WP90_template_for_fit_histos_3D_qcd_file_%d", fileindex);
  } else {
    std::cerr << "Unknown dataType " << dataType << std::endl;
    return;
  }

  // Get total events to compute this job's range
  tTree tmp;
  tmp.Init(filename, isMC);
  Long64_t n_total = tmp.GetEntries();

  // job_idx is global 
  // Compute local one
  int local_job_idx = job_idx % n_jobs; 
  Long64_t chunk   = (n_total + n_jobs - 1) / n_jobs;
  Long64_t ev_first = (Long64_t) local_job_idx * chunk;
  Long64_t ev_last  = std::min(ev_first + chunk, n_total);

  // TEST few entries 
  // -- test with one job
  n_total = 1000; // total events
  ev_first = 0; 
  ev_last = 1000;

  std::cout
  << "  global job = " << job_idx 
   << "  Job " << local_job_idx << "/" << n_jobs
            << "  events [" << ev_first << ", " << ev_last << ")"
            << "  total=" << n_total << std::endl;


  make_templates(filename, output_folder, output_hist, domain,
                 pT_low, pT_high, n, btag, isMC, dataType,
                 ev_first, ev_last, local_job_idx);
}
