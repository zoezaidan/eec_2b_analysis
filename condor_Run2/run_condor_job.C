// Entry point called by each condor job.
// Computes its event range from job_idx and n_jobs, then calls make_templates.
//
// Usage (from run_job.sh):
//   root -b -q -l 'run_condor_job.C+(job_idx, n_jobs, dataType, pT_low, pT_high, n, btag_int, isMC_int)'

#include "../create_files_for_template_fit.cpp"

void run_condor_job(Int_t job_idx, Int_t n_jobs, Int_t dataType,
                    Float_t pT_low, Float_t pT_high, Int_t n,
                    Int_t btag_int, Int_t isMC_int) {

  bool btag = (btag_int != 0);
  bool isMC = (isMC_int != 0);

  TString filename, output_hist;
  TString output_folder = "/data_CMS/cms/zaidan/analysis_lise/";
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
    filename     = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    output_hist  = "template_for_fit_histos_3D_qcd";
  } else {
    std::cerr << "Unknown dataType " << dataType << std::endl;
    return;
  }

  // Get total events to compute this job's range
  tTree tmp;
  tmp.Init(filename, isMC);
  Long64_t n_total = tmp.GetEntries();

  Long64_t chunk   = (n_total + n_jobs - 1) / n_jobs;
  Long64_t ev_first = (Long64_t)job_idx * chunk;
  Long64_t ev_last  = std::min(ev_first + chunk, n_total);

  std::cout << "Job " << job_idx << "/" << n_jobs
            << "  events [" << ev_first << ", " << ev_last << ")"
            << "  total=" << n_total << std::endl;

  make_templates(filename, output_folder, output_hist, domain,
                 pT_low, pT_high, n, btag, isMC, dataType,
                 ev_first, ev_last, job_idx);
}
