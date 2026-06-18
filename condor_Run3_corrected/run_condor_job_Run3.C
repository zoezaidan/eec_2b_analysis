// Entry point called by each condor job.
// Computes its event range from job_idx and n_jobs, then calls make_templates.
//
// Usage (from run_job.sh):
//   root -b -q -l 'run_condor_job.C+(job_idx, n_jobs, dataType, pT_low, pT_high, n, btag_int, isMC_int, fileindex)'


// ---- Notes -------
// --- Condor job for Run3 has input root files in directories from 0-9 for MC, and from 2-4 for data
// -- Output directory to be reviewd before run code. 
// -------------------


#include "../create_files_for_template_fit_AllRun2Run3.cpp"

void run_condor_job_Run3(Int_t job_idx, Int_t n_jobs,
                    Int_t RunN, 
                    Int_t dataType,
                    Float_t pT_low, Float_t pT_high, Int_t n,
                    Int_t btag_int, Int_t isMC_int
                    , Int_t fileindex // fileindex to choose input file only
                    ){


  std::cout << "ENTERED run_condor_job_Run3()" << std::endl;
  std::cout << "JOBID = " << job_idx << std::endl;


  bool btag = (btag_int != 0);
  bool isMC = (isMC_int != 0);

  TString filename, output_hist;
  TString output_folder = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3_corrected/JobResult/";
  TString domain = ".root";

  TString RunN_str = (RunN == 2) ? "Run2" : (RunN == 3) ? "Run3" : "UnknownRun";

  if(RunN == 2) {output_folder = "/data_CMS/cms/zaidan/analysis_lise/Run2/";
    //sanity check

    cout<<"---->>>> RUN 2" <<endl;

    if (isMC && dataType < 1) {
      std::cerr << "Invalid data type for MC sample" << std::endl;
      return;}

    if (!isMC && dataType > 1) {
      std::cerr << "Invalid data type for data sample" << std::endl;
      return;}

    if(dataType == -1){//________________________________data______________________________
      filename = "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root";
      output_hist = RunN_str + "secondbinsplitting_MAY_WP0898_template_for_fit_histos_3D_LowEG_f";
      isMC = false;
      cout<<"you chose data Low" <<endl;
      }

    else if(dataType == 0) {
      filename = "/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      output_hist = RunN_str + "secondbinsplitting_MAY_WP0898_template_for_fit_histos_3D_HighEG_f";
      isMC = false;
      cout<<"you chose data High" <<endl;       
      }      
                                                                                                                                                                                                                                                                            
    else if(dataType == 1){//________________________________bjet______________________________
      filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      output_hist = RunN_str + "secondbinsplitting_MAY_WP0898_template_for_fit_histos_3D_bjet_f";
      std::cout << "Creating files for template fit for bjet sample" << std::endl;
      cout<<"you chose bjet MC" <<endl;
      }

    else if(dataType == 2){//________________________________dijet______________________________
      filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root"; 
      output_hist = RunN_str + "secondbinsplitting_MAY_WP0898_template_for_fit_histos_3D_qcd_f";
      std::cout << "Creating files for template fit for qcd sample" << std::endl;
      cout<<"you chose qcd MC" <<endl;
      }

    else{
      cout<<"undefined data type"<<endl;
      return; 
      }
    }


  if(RunN == 3) {
    // output_folder = "/data_CMS/cms/zaidan/analysis_lise/Run3/";
    output_folder = "/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/condor_Run3_corrected/JobResult/";
    
    cout<<"---->>>> RUN 3" <<endl;
    //sanity check
    if (isMC && dataType < 1) {
      std::cerr << "Invalid data type for MC sample" << std::endl;
      return;}

    if (!isMC && dataType > 1) {
      std::cerr << "Invalid data type for data sample" << std::endl;
      return;}


    else if(dataType == 0) { //________________________________data______________________________
      // -- Corrected data: all triggers starting from 60 GeV
      // filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/HiForestMiniAOD_v2_TChains.root";

      filename = Form("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/%d/merged_HiForestMiniAOD_v2.root", fileindex);
      output_hist = Form("%s_secondbinsplitting_WP0872_template_for_fit_histos_3D_data_%d", RunN_str.Data(), fileindex);

      isMC = false;
      cout<<"you chose data" <<endl;   
      cout << "file name = " << filename <<endl;    
      }      
                                                                                                                                                                                                                                                                            
    else if(dataType == 1){//________________________________bjet______________________________
      filename = "";
      output_hist = RunN_str + "secondbinsplitting_MAY_WP0898_template_for_fit_histos_3D_bjet_f";
      std::cout << "Creating files for template fit for bjet sample" << std::endl;
      cout<<"you chose bjet MC" <<endl;
      }

    else if(dataType == 2){//________________________________MC: dijet______________________________
      // -- Corrected file: QCD
        // filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root";// TChain but needs Matt. path in reading events
        filename = Form("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/%d/merged_HiForestMiniAOD_v2.root", fileindex);
        output_hist = Form("%s_secondbinsplitting_June_WP0872_template_for_fit_histos_3D_qcd_%d", RunN_str.Data(), fileindex);
        isMC = true;

        std::cout << "Creating files for template fit for qcd sample" << std::endl;
        cout<<"you chose qcd MC" <<endl;
      }
    
    else{
          cout<<"undefined data type"<<endl;
          return; 
          }
        }

    cout << "input file: "<< filename << endl;


  // Get total events to compute this job's range
  tTree tmp;
  cout << "Before Init() in Run condor() " << endl;
  
    // Tchain root file reads relative paths from where Matt. creted the TChain! --> need to read it from there 
    // TString oldDir = gSystem->WorkingDirectory();
    // gSystem->ChangeDirectory("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes");

    tmp.Init(filename, isMC, RunN);
    cout << "After init of tree inside run condor()" << endl;
    Long64_t n_total = tmp.GetEntries();
    cout << "n total = " << n_total << endl;

  // job_idx is global 
  // Compute local one
  int local_job_idx = job_idx % n_jobs; 
  Long64_t chunk   = (n_total + n_jobs - 1) / n_jobs;
  Long64_t ev_first = (Long64_t) local_job_idx * chunk;
  Long64_t ev_last  = std::min(ev_first + chunk, n_total);

// TEST
//   cout << "using Few events for test only" <<endl;
// ev_first = 0;
// ev_last = 100;

  std::cout
  << "  global job = " << job_idx 
   << "  Job " << local_job_idx << "/" << n_jobs
            << "  events [" << ev_first << ", " << ev_last << ")"
            << "  total=" << n_total << std::endl;


  make_templates(RunN,
                filename, output_folder, output_hist, domain,
                pT_low, pT_high, n, btag, isMC,
                dataType,
                ev_first, ev_last, local_job_idx);
}
