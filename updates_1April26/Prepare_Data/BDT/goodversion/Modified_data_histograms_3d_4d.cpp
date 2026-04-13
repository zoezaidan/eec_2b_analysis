// ---- Local HEADERS ---------
#include "tTree_data.h" //If using data use "tTree_data.h"
#include "binning_histos_data.h"

// --- ROOT/c++  includes 
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <Math/Vector4D.h>
#include <map>
#include <string>
#include <random>
#include <vector>
#include "TString.h"
#include "Math/VectorUtil.h"


//Skip MC events that have a too large weight to stabilize the distribution
bool skipMC(double jtpt, double pthat) {//double refpt
    //if (!(refpt>0)) return true;    
    if (pthat<0.35*jtpt) return true;
    return false;
}

//Print a vector (for debugging)
void print_vector(std::vector<Int_t> &vec){
    for(Int_t i = 0; i < vec.size(); i++){
        std::cout << vec.at(i) << " ";
    }
    std::cout << std::endl;
}

//Check that the same track is not twice in the list of track vectors (for debugging)
bool check_twice(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors){

    Int_t ntrk = trackVectors.size();

    for(Int_t i = 0; i < ntrk; i++){

        ROOT::Math::PtEtaPhiMVector v;
        v = trackVectors.at(i);

        for(Int_t j = 0; j < ntrk; j++){
            if((i != j) && (v == trackVectors.at(j))) return true;
        }
    }

    return false;
}

auto deltaR2 = [](const ROOT::Math::PtEtaPhiMVector &a,
                  const ROOT::Math::PtEtaPhiMVector &b) {
    //given two vectors, computes deltaR squared 
    double dEta = a.Eta() - b.Eta();
    double dPhi = std::acos(std::cos(a.Phi() - b.Phi()));
    return dEta*dEta + dPhi*dPhi;
};

struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    std::vector<int> trkMatchSta;
  vector<int> trkSvtxId;   
};



void groupVertexes1_v2(std::vector<Vertex>& vertices)
{
    if (vertices.size() < 2) return;

    while (vertices.size() > 2) {

        double min_distance = std::numeric_limits<double>::infinity();
        size_t index1 = 0;
        size_t index2 = 1;

        // find closest pair in (eta, phi)
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = i + 1; j < vertices.size(); ++j) {

                double dEta = vertices[i].p4.Eta() - vertices[j].p4.Eta();
                double dPhi = std::acos(std::cos(
                    vertices[i].p4.Phi() - vertices[j].p4.Phi()));

                double dist = dEta*dEta + dPhi*dPhi;

                if (dist < min_distance) {
                    min_distance = dist;
                    index1 = i;
                    index2 = j;
                }
            }
        }

    // ---- merge vertex index2 into index1 ----
    vertices[index1].p4 += vertices[index2].p4;
    
    // merge track lists
    auto &t1 = vertices[index1].tracks;
    auto &t2 = vertices[index2].tracks;
    t1.insert(t1.end(), t2.begin(), t2.end());
    
    //merge vertex id
    vertices[index1].trkSvtxId.insert(
    vertices[index1].trkSvtxId.end(),
        vertices[index2].trkSvtxId.begin(),
        vertices[index2].trkSvtxId.end());
    
    // remove merged vertex
    vertices.erase(vertices.begin() + index2);

    }
}

vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs_BDTcut_v2(
  tTree& t,     //event
  Int_t& ijet,  //jet nb
  Long64_t& ient, 
  double BDTcutTh // 0.365 
){

/* 
// Modified version of makeSvtxs_witBDT to remove all mathcing info (was used for quality checks,
// but not availble for DATA)
// Removed all non necessary histograms and counters 
// Added argument for the BDT cut threshold 
*/
  std::map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;// map from vertex id to list of tracks :  secVtxs[svId] = list of track 4 vec in that vertex 
  vector<ROOT::Math::PtEtaPhiMVector> empty;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;


    // Loop over tracks
    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //standard cuts for tracks 
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;

    // -- Add trkbdt score
    if (t.trkBdtScore[itrk]<= BDTcutTh) continue;
    
    //build v1, the 4 vec of the track 
        ROOT::Math::PtEtaPhiMVector v1;
        v1.SetEta(t.trkEta[itrk]);
        v1.SetPt(t.trkPt[itrk]);
        v1.SetPhi(t.trkPhi[itrk]);
        int pid = std::abs(t.trkPdgId[itrk]);
        if      (pid == 211)  v1.SetM(0.139570);
        else if (pid == 13)   v1.SetM(0.105658);
        else if (pid == 11)   v1.SetM(0.000510);
        else if (pid == 2212) v1.SetM(0.938272);
        else if (pid == 321)  v1.SetM(0.493677);
        else if (pid == 3112) v1.SetM(1.19744);                                                                                             
        else if (pid == 3222) v1.SetM(1.18937);                                                                                             
        else if (pid == 3312) v1.SetM(1.32171);                                                                                             
        else if (pid == 3334) v1.SetM( 1.67245);
        else                  v1.SetM(0.139570); //  safety

       // no_sv_tracks  
        if (t.trkSvtxId[itrk] < 0) {
            no_sv_list.push_back(v1);
        } 
        //map vertex id to list of tracks
        else if (t.trkSvtxId[itrk] >= 0 ) {
          secVtxs[t.trkSvtxId[itrk]].push_back(v1);
        } 
    }// end tracks loop
    
    // Build Vertex objects
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      // auto& v = secVtxsMatchSta[originalSvId];

      Vertex vtx;
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      for (size_t i = 0; i < it.second.size(); ++i) 
      {      //loop over tracks in this vertex (it.second) 
            vtx.p4 += it.second[i];                       //sums the tracks' 4 vector into vtx.p4   (it.second[i] = 4 vec of track i 
            vtx.tracks.push_back(it.second[i]);           //copies the track itself in vtx.tracks 
            vtx.trkSvtxId.push_back(originalSvId);  //store provenance 
        }
      vertices.push_back(vtx);  //adds the Vertex vtx to the vector vertices 
    }// end for 
    // what we have created is tracks[i] ↔ trkMatchSta[i] ↔ trkSvtxId[i] correspondance :  originalSvId = vertex.trkSvtxId[i]

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1_v2(vertices);
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl; 
       return empty; // try to overcome break segementation 
    }

    //vertices names
    auto &vertex0 = vertices[0];
    auto &vertex1 = vertices[1];
    // -- Add in the particles that are not part of a sv
    for (const auto& v1 : no_sv_list) {
      //std::cout<<"add tracks "<<no_sv_list.size()<<std::endl;
      double d0 = deltaR2 (v1, vertex0.p4);
      double d1 = deltaR2 (v1, vertex1.p4);
       if (d0 < d1 ) {
        vertex0.p4 += v1;
      }
      else if ( d0 > d1 ) {
            vertex1.p4 += v1 ;
      }
    }
    no_sv_list.clear();

    std::set<int> originalVertices0(  vertices[0].trkSvtxId.begin(), vertices[0].trkSvtxId.end());   // keeps only one copy of each sv nb 
    int nOriginalVertices0 = originalVertices0.size();
    // vertex 1
    std::set<int> originalVertices1(  vertices[1].trkSvtxId.begin(), vertices[1].trkSvtxId.end() );
    int nOriginalVertices1 = originalVertices1.size();


    // Return summed 4-vectors of the 2 vertices
    vector<ROOT::Math::PtEtaPhiMVector> vecFinalSecVtxs;
    vecFinalSecVtxs.push_back(vertices[0].p4);
    vecFinalSecVtxs.push_back(vertices[1].p4);


return vecFinalSecVtxs;
}
   
//_____________________________________________________________________________________________
//_______________________________Get EEC histograms at reco____________________________________
//_____________________________________________________________________________________________

//Create a 4D 2-point EEC histogram from tree
void do_hist(TString &filename,  TString &dataset, Float_t &pT_low, Float_t &pT_high, bool &aggregated, double &BDTcutTh, bool &btag, Int_t &n, Int_t &lowEn, Int_t &highEn, bool &data, TString folder){
    
    TString fin_name = filename;
    tTree t(fin_name);

    /* 
    // Function modifed mainly for data files. 
    // Note that TH3D are also added (considering the Sumw2() option t propoerly assign the stats. uncert.)
    // Which is not the case of 4D hist since the entries for each eec are filled in a seperate EEC axis, 
    //  and not as histogram content assuming Sumw2().
    //////////////////////////////////////////////////////////
    //// The difference on the stats. uncert. between 4D and 3D histograms should be checked after output is read. /// 
    //////////////////////////////////////////////////////////
    */

    // Turn off all branches and turn on only the interesting ones
    // Attention! If a branch is off, it will return bs without crashing 
    t.SetBranchStatus("*", 0);
    std::vector<TString> active_branches = {
        // "weight",
        "jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB",
        "ntrk", "trkPt", "trkJetId",
        "trkEta", "trkPhi", "jtNbHad", "jtHadFlav",
        "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId",
        "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
    };
    t.SetBranchStatus(active_branches, 1);
    
    
    //Define the number of bins 4d histogram (the 4 axes are: mB, dr, eec weight, jet pT)
    Int_t nbins[4] = {bins_mb, bins_dr, bins_eec, jtpt_bins}; // #Bins per axis 
    
    //Define 4D histogram and set bin edges (variable bin size)
    THnD *h4D = new THnD("h4D", "h4D", 4, nbins, NULL, NULL);
    h4D->SetBinEdges(mb_dim, mb_binsVector);
    h4D->SetBinEdges(dr_dim, dr_binsVector);
    h4D->SetBinEdges(eec_dim, eec_binsVector);
    h4D->SetBinEdges(pt_dim, jtpt_binsVector);

    //create a 3D histogram for the mB, dr and pT distributions

    TH3D *h3D = new TH3D("h3D", "h3D", bins_mb, mb_binsVector,
                                       bins_dr, dr_binsVector,
                                       jtpt_bins, jtpt_binsVector);
    h3D->Sumw2();

    // -- test another h3D without EEC value, fill counts only 
    TH3D *h3D_noweights = new TH3D("h3D_noweights", "h3D_noweights", bins_mb, mb_binsVector,
                                       bins_dr, dr_binsVector,
                                       jtpt_bins, jtpt_binsVector);
    h3D_noweights->Sumw2();



    // -- hist for test: 1D mass distribution 
    TH1D *hm = new TH1D("hm", ";mass;counts",200, 0, 20);
        TH1D *hm_witheecweight = new TH1D("hm_witheecweight", ";mass;counts",200, 0, 20);
            TH1D *hm_out = new TH1D("hm_out", ";mass;counts",200, 0, 20);

    TH1D *hm1 = new TH1D("hm1", ";mass;counts",200, 0, 20);
    TH1D *hm2 = new TH1D("hm2", ";mass;counts",200, 0, 20);

    //Save the prescale factor (only for 40 GeV trigger)
    double prescale_pf40 = 33.917210;

    std::cout << "Dataset = " << dataset << std::endl;
    std::cout << "Events = " << t.GetEntries() << std::endl;


    //looping over events
    // std::cout << "Looping over events" << std::endl;
    Long64_t  n_entries = t.GetEntries(); // 365 706 800
    //////////
    //// For tests
    n_entries = 50e+05;

    ///////// 
    for (Long64_t ient = 0; ient < n_entries; ient++) { //
        // Print progress
        if (ient % 50000 == 0) {                                                                                                         
              float percent = 100.0 * ient / n_entries;
              std::cout << "\rProcessing: "  << percent << " %" << std::flush;
        }

        //get tree entry ient
        t.GetEntry(ient); 

        //select HLT events in data depending if we are using the HighEG dataset or the LowEG dataset
        if(data && highEn){
            if(!(t.HLT_HIAK4PFJet100_v1 || t.HLT_HIAK4PFJet80_v1)) {continue;}
        }
         // cout <<"skipped event by HG trigger" << endl;
        
        if(data && lowEn){
            if(!((t.HLT_HIAK4PFJet60_v1 == 1 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) ||
                 (t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0))) continue;
        }

        // Loop over jets 
        for (Int_t ijet = 0; ijet < t.nref; ijet++) {

            // Skip jets outside tracker 
            if (std::abs(t.jteta[ijet]) > 1.9) continue;
            //Select jet pt - low limit
            if (std::abs(t.jtpt[ijet]) < pT_low) continue;
            //Select jet pt - high limit
            if (std::abs(t.jtpt[ijet]) > pT_high) continue;
            
            // Select jet flavour and/or select on the number of b hadrons: none for data 
            
            //Select jets passing the b-jet tagging
            if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;


            // -- test reading branches properly ?
            // only one event 
            // cout << "here is some branches result for a choosen event"  << endl;
            // cout <<  "Jet Eta, pt  = " << t.jteta[ijet] << ", "  << t.jtpt[ijet] << endl;

            //Define 2B hadron mass
            double mB = 0, dr = 0, pt1 = 0, pt2 = 0, eec = 0, jtpt = 0;

            //Aggregate tracks into 2Bs (with BDT information)
            if (aggregated){
                //Aggregate tracks into 2Bs (with BDT information)            
                vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs_BDTcut_v2(t, ijet, ient, BDTcutTh);
                if (reco_hadrons_4vec.size() !=  2) continue;
                mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();
           
                cout << "mass =  " << mB << endl;

                // -- Compute eec between aggregated Bs and fill histograms
                  dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
                  pt1 = reco_hadrons_4vec[0].Pt();
                  pt2 = reco_hadrons_4vec[1].Pt();
                  eec = std::pow(pt1 * pt2, n);
                  jtpt = t.jtpt[ijet];

                // // -- test hist 
                hm ->Fill(mB);
                hm_witheecweight->Fill(mB, eec);
                hm1->Fill(reco_hadrons_4vec[0].M());
                hm2->Fill(reco_hadrons_4vec[1].M());


                // Filling TH3D here : sape is almost wrong? 
                  // Double_t filling[4] = {mB, dr, eec, jtpt};     
                  //   h4D->Fill(filling);

                  //   h3D->Fill(mB, dr, jtpt, eec); // 3D hist with bin content = eec  value 
                  //   h3D_noweights ->Fill(mB, dr, jtpt);


              }
              else{
                cout << "aggregated condition not true!!!! ESCAPE" << endl; return;
              }



              // While filling here is ok! 
              // test fill TH3D first! 
                    // Double_t filling[4] = {mB, dr, eec, jtpt};     
                    // h4D->Fill(filling);

                    // h3D->Fill(mB, dr, jtpt, eec); // 3D hist with bin content = eec  value 
                    // h3D_noweights ->Fill(mB, dr, jtpt);

                // cout << "After aggregation " << endl;

                // -- what happwens If I fil lthe mass here, is it different mass !! ?? 
                hm_out ->Fill(mB, eec);


/*
                   //Add prescale weight if necessary: for the jet40 trigger --> apply a weight 
                   if(t.HLT_HIAK4PFJet40_v1 == 1 &&
                         t.HLT_HIAK4PFJet60_v1 == 0 &&
                         t.HLT_HIAK4PFJet80_v1 == 0 &&
                         t.HLT_HIAK4PFJet100_v1 == 0){ eec*=prescale_pf40; cout << "prescale40 applied!!! " << endl;}


                    //Fix the under/overflow
                    if(dr < dr_min) dr = dr_min_fill;
                    if(dr >= dr_max) dr = dr_max_fill;
                    if(mB >= mb_max) mB = mb_max_fill;
                    if(eec >= eec_max) eec = eec_max_fill;
  */              

                    // cout << "now fill " << endl;

                    //Fill ThnD 
                    Double_t filling[4] = {mB, dr, eec, jtpt};     
                    h4D->Fill(filling);// no MC weight 
                    h3D->Fill(mB, dr, jtpt, eec); // 3D hist with bin content = eec  value 
                    h3D_noweights ->Fill(mB, dr, jtpt);


                    // cout << "After fill " << endl;
        } // jet loop 
    }// event loop

//Create the fout name depending on the selection
TString fout_name = "hist_3d_4d_";

TString label = "";
if (data) label+= "data";
if(aggregated) fout_name += "aggr_BDT_";
if(!btag) label += "_notag"; 
fout_name += TString(Form("n%i_",n)) + label + "_" + dataset + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high))) + ".root";

//Save histograms and close file
std::cout << "Creating file: " << fout_name << std::endl;
TFile *fhist = new TFile(folder+fout_name, "recreate");

h3D->Write();
h4D->Write();

// -- for test 
hm->Write();
hm1->Write();
 hm2->Write();
hm_witheecweight->Write();
hm_out->Write();

// -- test projection 
TH1D* h3dproj_mass = (TH1D*) h3D->ProjectionX("h3dproj_mass", 1, h3D->GetNbinsY(), 1, h3D->GetNbinsZ());
h3dproj_mass->Write();
h3D_noweights->Write();

fhist->Close();

}

void Modified_data_histograms_3d_4d(){

    //________________________________data______________________________
    //Create vectors for more datasets
    std::vector<TString> filenames{"/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root"
                                    ,"/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root"};
    std::vector<TString> datasets{"LowEGJet","HighEGJet"};
    
    //Folder name
    // TString folder = "/home/llr/cms/shatat/CMSAnalysis/EECs/Prepare_Data/Data_Hist/";

    TString folder = "/home/llr/cms/shatat/CMSAnalysis/EECs/Prepare_Data/Test_Data_Hist/";


    bool data = true;
    
    //Select pT range (later divided into 3 bins, central bin is the nominal one)
    Float_t pT_low = 80;
    Float_t pT_high = 140;

    //select eec weight exponent
    Int_t n=1;

    //apply b-tagging
    bool btag = true;

    //Type of datafile (true only for data)
    std::vector<Int_t> lowEns{1, 0}; // first dataset true
    std::vector<Int_t> highEns{0, 1}; // second dataset true 

    //Aggregate tracks coming from B hadrons
    bool aggregated = true;
    double BDTcutTh = 0.365;

    //Loop over datasets:
    // test HG set only 
    for(Int_t i = 1; i < filenames.size(); i++){
        do_hist(filenames.at(i), datasets.at(i),
          pT_low, pT_high,
          aggregated, 
          BDTcutTh,
          btag,
          n,
          lowEns.at(i), highEns.at(i),
          data, 
          folder);
    }
}