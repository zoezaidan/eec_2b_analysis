#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <iostream>
#include <Math/Vector4D.h>
#include <map>
#include <unordered_map>
#include <string>
#include <random>
#include <vector>
#include "TString.h"
#include "Math/VectorUtil.h"
#include "tTree.h"
#include "binning_histos_small.h"
struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    std::vector<int> trkMatchSta;
  std::vector<int> trkSvtxId;
};


//given two vectors, computes deltaR squared 
auto deltaR2 = [](const ROOT::Math::PtEtaPhiMVector &a,
                  const ROOT::Math::PtEtaPhiMVector &b) {
    double dEta = a.Eta() - b.Eta();
    double dPhi = std::remainder(a.Phi() - b.Phi(), 2*M_PI);
    return dEta*dEta + dPhi*dPhi;
};


//printing function for a Vertex object 
void printvtx (const Vertex& vertex,
              const TString& txt,
	       int vtxnb,  Int_t& ijet, Long64_t& ient) {
  for (size_t itrk = 0; itrk < vertex.trkSvtxId.size(); itrk++) {
    std::cout<<"ient : "<<ient <<"  ijet : "<<ijet << txt << "vertex number: "<<vtxnb<<"  track original vertex : "<< vertex.trkSvtxId[itrk] <<"    and track original status : "<<vertex.trkMatchSta[itrk] << " pT = " << vertex.tracks[itrk].Pt() <<  std::endl;
  }
}
//Skip MC events that have a too large weight to stabilize the distribution                                                          
bool skipMC(double pt, double pthat) {
  if (pthat<0.35*pt) return true;                                                                                                  
  return false;
}                                                                                                                                  
                                                                                                                                     
void PartialBsAggregation(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, tTree& t, Int_t ijet){
  hadrons_4vec.clear();
  hadrons_stat.clear();         
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
    // Track must belong to this jet
    if (t.refTrkJetId[itrk] != ijet) continue;
    // pT cut
    if (t.refTrkPt[itrk] < 1) continue;                                                                                                  // Assign mass                                                                                                                  
    double mass = 0.0;                                                                                                               
    int pid = std::abs(t.refTrkPdgId[itrk]);
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
    ROOT::Math::PtEtaPhiMVector v(t.refTrkPt [itrk], t.refTrkEta[itrk], t.refTrkPhi[itrk], mass);
    
    Int_t status = t.refTrkSta[itrk];                                                                                                
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


void MatchingTracksAggregation(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, tTree& t, Int_t ijet, Int_t ient){
  hadrons_4vec.clear();                                                                                                               
  hadrons_stat.clear();

  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {       //loop over trks
    // Track must belong to this jet
    if (t.trkJetId[itrk] != ijet) continue;
    // pT cut
    if (t.trkPt[itrk] < 1) continue;
    if (t.trkMatchSta[itrk] <100 ) continue;
    //if (t.trkSvtxId[itrk] <0 ) continue ;      //trk must be in SV
// Assign mass                                                                                                                    
    double mass = 0.0;                                                                                                                
    int pid = std::abs(t.trkMatchPdgId[itrk]);
    if      (pid == 211)  mass = 0.139570;                                                                                            
    else if (pid == 13)   mass = 0.105658;                                                                                            
    else if (pid == 11)   mass = 0.000510;                                                                                            
    else if (pid == 2212) mass = 0.938272;                                                                                            
    else if (pid == 321)  mass = 0.493677;                                                                                            
    else if (pid == 3112) mass = 1.19744;                                                                                             
    else if (pid == 3222) mass = 1.18937;                                                                                             
    else if (pid == 3312) mass = 1.32171;                                                                                             
    else if (pid == 3334) mass = 1.67245;  
    else if (pid == 1)   continue;
        //skip not matched                                                                                       
    else std::cout << "PDG:" << pid << std::endl;                                                                                     
                                                                                                                                      
    // 4 vector                                                                                                                       
    ROOT::Math::PtEtaPhiMVector v(t.trkPt [itrk], t.trkEta[itrk], t.trkPhi[itrk], mass);                                                                                              
    
    Int_t status = t.trkMatchSta[itrk];                                                                                                 
    
    if (status < 100) {   continue;                                                                                                            
      // I don't aggregate, just store each particle separately                                                                       
      //hadrons_stat.push_back(status);                                                                                                 
      //hadrons_4vec.push_back(v);                                                                                                      
    }                                                                                                                                 
    else {                                                                                                                            
      // I aggregate only if status >= 100                                                                                            
      auto it = std::find(hadrons_stat.begin(), hadrons_stat.end(), status);      //check if status alreday exists in hadrons_stat                                                                                                           
      if (it == hadrons_stat.end()) {                     //if not, create new hadrons                                                                             
        hadrons_stat.push_back(status);                                                                                               
        hadrons_4vec.push_back(v);                                                                                                    
      }                                                                                                                               
      else {                                                                                                                          
        size_t index = std::distance(hadrons_stat.begin(), it);                                                                       
        hadrons_4vec[index] += v;}                                                     //if yes, add to existing                               
    } 
  }
  
}

void groupVertexes1(std::vector<Vertex>& vertices)
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
                double dPhi = std::remainder(vertices[i].p4.Phi() - vertices[j].p4.Phi(), 2*M_PI);

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
	
	// merge match status lists
	auto &s1 = vertices[index1].trkMatchSta;
	auto &s2 = vertices[index2].trkMatchSta;
	s1.insert(s1.end(), s2.begin(), s2.end());

	//merge vertex id
	vertices[index1].trkSvtxId.insert(
	vertices[index1].trkSvtxId.end(),
        vertices[index2].trkSvtxId.begin(),
        vertices[index2].trkSvtxId.end());
	
	// remove merged vertex
	vertices.erase(vertices.begin() + index2);

    }
}


// ---- CHOOSE WETHER YOU WANT TO USE THE STATUS (makeSvtxs) OR THE BDT SCORE (makeSvtxs_withBDT) FOR THE TRACKS AGGREGATION ----

vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs(
  tTree& t,     //event
  Int_t& ijet,  //jet nb
  Long64_t& ient,
  double& agg_fail,     //dominant trkMatchSta fraction for vertex 0
  double& nb_sv ,    //   ..      ..           ..     ..    ..   1
  double& sv_fail,
  double& merge_fail){

  std::unordered_map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;
  std::unordered_map<Int_t, std::vector<int>> secVtxsMatchSta;
  vector<ROOT::Math::PtEtaPhiMVector> empty;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;
  no_sv_list.reserve(t.ntrk);
  no_sv_sta_list.reserve(t.ntrk);

    // Loop over tracks
    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //cuts for tracks
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;
	if (t.trkMatchSta[itrk] <100 ) continue;

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
	  no_sv_sta_list.push_back(t.trkMatchSta[itrk]); 
	} 
	//map vertex id to list of tracks
	else if (t.trkSvtxId[itrk] >= 0 ) {
	  secVtxs[t.trkSvtxId[itrk]].push_back(v1);
	  secVtxsMatchSta[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]); // store trkMatchSta of trk itrk as the itrk elem of the value associated to trkSvtxId
	} 
    }

    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      if (nb_sv < 15 ) {
      } 
      return empty;    }

    // Build Vertex objects
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      auto& v = secVtxsMatchSta[originalSvId];
      if (std::adjacent_find(v.begin(), v.end(),std::not_equal_to<int>()) != v.end()) {      // if not all tracks in the sv come from the same B 
	sv_fail += 1.0;
	//return empty;
      }
      Vertex vtx;
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      vtx.tracks.reserve(it.second.size());
      vtx.trkMatchSta.reserve(it.second.size());
      vtx.trkSvtxId.reserve(it.second.size());
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second)
	vtx.p4 += it.second[i];
	vtx.tracks.push_back(it.second[i]);
	vtx.trkMatchSta.push_back(v[i]);  // reuse already-looked-up reference
	vtx.trkSvtxId.push_back(originalSvId);
        }
      vertices.push_back(std::move(vtx));
    }
    // what we have created is tracks[i] ↔ trkMatchSta[i] ↔ trkSvtxId[i] correspondance :  originalSvId = vertex.trkSvtxId[i]

       
          

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices);
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl; 
       //return empty;
    }

    //vertices names
    auto &vertex0 = vertices[0];
    auto &vertex1 = vertices[1];

    //VERTEX 0 
    std::map<int,int> countMap0;
    for (size_t i = 0; i < vertex0.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex 
        int sta0 = vertex0.trkMatchSta[i];
        countMap0[sta0]++;
    }
    int totalTracks0 = vertex0.tracks.size();
    int maxCount0 = 0;
    for (const auto &kv : countMap0) {
      if (kv.second > maxCount0) maxCount0 = kv.second;
    }
    double prop0 = double(maxCount0) / double(totalTracks0);
    if (prop0 != 1) merge_fail += 1;

    //VERTEX 1
    std::map<int,int> countMap1;
    for (size_t i = 0; i < vertex1.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex
        int sta1 = vertex1.trkMatchSta[i];
        countMap1[sta1]++;
    }
    int totalTracks1 = vertex1.tracks.size();
    // std::cout << "Vertex 1 track proportions:\n";
    int maxCount1 = 0;
    for (const auto &kv : countMap1) {
      if (kv.second > maxCount1) maxCount1 = kv.second;
    }
    double prop1 = double(maxCount1) / double(totalTracks1);
    if (prop1 != 1) {
      merge_fail += 1;
      //printvtx( vertex1 , " wrong merging ", 1);
      //return empty ;
    }

    //Add in the particles that are not part of a sv
    int i = 0;
    for (const auto& v1 : no_sv_list) {
      double d0 = deltaR2(v1, vertex0.p4);
      double d1 = deltaR2(v1, vertex1.p4);
      if (d0 < d1) {
	vertex0.p4 += v1;
	if (no_sv_sta_list[i] != vertex0.trkMatchSta[0]) agg_fail += 1;
      }
      else if (d0 > d1) {
	vertex1.p4 += v1;
	if (no_sv_sta_list[i] != vertex1.trkMatchSta[0]) agg_fail += 1;
      }
      i += 1;
    }
    no_sv_list.clear();

    // Return summed 4-vectors of the 2 vertices
    vector<ROOT::Math::PtEtaPhiMVector> vecFinalSecVtxs;
    vecFinalSecVtxs.push_back(vertices[0].p4);
    vecFinalSecVtxs.push_back(vertices[1].p4);

    return vecFinalSecVtxs;
}


// Status is used for quality plots and BDT cut for the merging.
vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs_withBDT(
  tTree& t,     //event
  Int_t& ijet,  //jet nb
  Long64_t& ient, 
  double& agg_fail,     //dominant trkMatchSta fraction for vertex 0
  double& nb_sv ,    //   ..      ..           ..     ..    ..   1
  double& sv_fail,
  double& merge_fail,
  TH1D* h_score_bkg, TH1D* h_score_sg
){

  std::unordered_map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;
  std::unordered_map<Int_t, std::vector<int>> secVtxsMatchSta;
  vector<ROOT::Math::PtEtaPhiMVector> empty;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;
  no_sv_list.reserve(t.ntrk);
  no_sv_sta_list.reserve(t.ntrk);

    // Loop over tracks

    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //cuts for tracks 
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;
    //if (t.trkMatchSta[itrk] <100 ) continue;

    //if (t.trkMatchSta[itrk] <= 1) continue; //activate or not

    //control plots for the BDT score distribution in 0b and 1b jets
    //if (t.trkMatchSta[itrk] == 1)
      //{h_score_bkg->Fill(t.trkBdtScore[itrk], weight_tree);}
    //else if (t.trkMatchSta[itrk] >= 100)      
     // {h_score_sg->Fill(t.trkBdtScore[itrk], weight_tree);}
  


    // -- Add trkbdt score
    if (t.trkBdtScore[itrk]<=0.365) continue;
	
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
	  no_sv_sta_list.push_back(t.trkMatchSta[itrk]); 
	} 
	//map vertex id to list of tracks
	else if (t.trkSvtxId[itrk] >= 0 ) {
	  secVtxs[t.trkSvtxId[itrk]].push_back(v1);
	  secVtxsMatchSta[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]); // store trkMatchSta of trk itrk as the itrk elem of the value associated to trkSvtxId
	} 
    }

    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      if (nb_sv < 15 ) {
      } 
      return empty;    }

    // Build Vertex objects
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      auto& v = secVtxsMatchSta[originalSvId];
      if (std::adjacent_find(v.begin(), v.end(),std::not_equal_to<int>()) != v.end()) {      // if not all tracks in the sv come from the same B 
	sv_fail += 1.0;
	//return empty;
      }
      Vertex vtx;
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      vtx.tracks.reserve(it.second.size());
      vtx.trkMatchSta.reserve(it.second.size());
      vtx.trkSvtxId.reserve(it.second.size());
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second)
	vtx.p4 += it.second[i];
	vtx.tracks.push_back(it.second[i]);
	vtx.trkMatchSta.push_back(v[i]);  // reuse already-looked-up reference
	vtx.trkSvtxId.push_back(originalSvId);
        }
      vertices.push_back(std::move(vtx));
    }
    // what we have created is tracks[i] ↔ trkMatchSta[i] ↔ trkSvtxId[i] correspondance :  originalSvId = vertex.trkSvtxId[i]

       
          

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices);
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl; 
       //return empty;
    }

    //vertices names
    auto &vertex0 = vertices[0];
    auto &vertex1 = vertices[1];

    //VERTEX 0 
    std::map<int,int> countMap0;
    for (size_t i = 0; i < vertex0.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex 
        int sta0 = vertex0.trkMatchSta[i];
        countMap0[sta0]++;
    }
    int totalTracks0 = vertex0.tracks.size();
    int maxCount0 = 0;
    for (const auto &kv : countMap0) {
      if (kv.second > maxCount0) maxCount0 = kv.second;
    }
    double prop0 = double(maxCount0) / double(totalTracks0);
    if (prop0 != 1) merge_fail += 1;

    //VERTEX 1
    std::map<int,int> countMap1;
    for (size_t i = 0; i < vertex1.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex
        int sta1 = vertex1.trkMatchSta[i];
        countMap1[sta1]++;
    }
    int totalTracks1 = vertex1.tracks.size();
    // std::cout << "Vertex 1 track proportions:\n";
    int maxCount1 = 0;
    for (const auto &kv : countMap1) {
      if (kv.second > maxCount1) maxCount1 = kv.second;
    }
    double prop1 = double(maxCount1) / double(totalTracks1);
    if (prop1 != 1) {
      merge_fail += 1;
      //printvtx( vertex1 , " wrong merging ", 1);
      //return empty ;
    }

    //Add in the particles that are not part of a sv
    int i = 0;
    for (const auto& v1 : no_sv_list) {
      double d0 = deltaR2(v1, vertex0.p4);
      double d1 = deltaR2(v1, vertex1.p4);
      if (d0 < d1) {
	vertex0.p4 += v1;
	if (no_sv_sta_list[i] != vertex0.trkMatchSta[0]) agg_fail += 1;
      }
      else if (d0 > d1) {
	vertex1.p4 += v1;
	if (no_sv_sta_list[i] != vertex1.trkMatchSta[0]) agg_fail += 1;
      }
      i += 1;
    }
    no_sv_list.clear();

    // Return summed 4-vectors of the 2 vertices
    vector<ROOT::Math::PtEtaPhiMVector> vecFinalSecVtxs;
    vecFinalSecVtxs.push_back(vertices[0].p4);
    vecFinalSecVtxs.push_back(vertices[1].p4);

    return vecFinalSecVtxs;
}

  





// - - - - - - - - - - -  USING FUNCTIONS  - - - - - - - - - - -
void filter_b_bb(TString filename, TString output_folder, TString output_hist, TString domain, Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC, Int_t dataType) {
  
  // true level information aggregated to partial Bs
  tTree t;
  t.Init(filename, isMC);
  t.SetBranchStatus("*", 0);
    double lt2sv = 0 ;
    double totjt = 0;
    double beforeskip = 0 ;
    double vertskip = 0 ;
    double matchskip = 0;
    double merge_fail = 0 ;
    double sv_fail = 0 ;
    double agg_fail = 0;
    double nb_sv=0 ;
    int counter0 = 0 ;
    int counter1 = 0 ;
    int counter2 = 0 ;  

  std::vector<TString> active_branches = {
    "weight","jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "jtNtrk",
      "ntrk", "trkJetId", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "refTrkPdgId","refTrkSta", "refTrkMass",  "refmB", "refpt", "refeta", "refphi",
      "nrefTrk", "refTrkJetId", "refTrkPt", "refTrkEta", "refTrkPhi", "refTrkY", "refNtrk",
      "jtNsvtx", "trkSvtxId",  "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",  "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
      "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1"};
    t.SetBranchStatus(active_branches, 1);
  
   
  // Plots or histograms
  TH3D *h3D_2b = new TH3D("h3D_2b", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_1b = new TH3D("h3D_1b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  TH3D *h3D_0b = new TH3D("h3D_0b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  
  TH1D *h_0b_score_bkg = new TH1D("h_0b_score_bkg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_1b_score_bkg = new TH1D("h_1b_score_bkg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_2b_score_bkg = new TH1D("h_2b_score_bkg", "Distribution of Score", 1000, -1, 1);
  
  TH1D *h_0b_score_sg = new TH1D("h_0b_score_sg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_1b_score_sg = new TH1D("h_1b_score_sg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_2b_score_sg = new TH1D("h_2b_score_sg", "Distribution of Score", 1000, -1, 1);
  

  h3D_2b->Sumw2();
  h3D_1b->Sumw2();
  h3D_0b->Sumw2();

  Long64_t n_events = t.GetEntries();                                                                                                
  for (Long64_t ient = 0; ient < n_events; ient++) {
    //Progress
    if (ient % 50000 == 0) {                                                                                                         
      float percent = 100.0 * ient / n_events;
      std::cout << "\rProcessing: "  << percent << " %" << std::flush;
    }
    t.GetEntry(ient);
    //get the MC event weight
    double weight_tree = t.weight;

    //Add trigger stuff 
    //select HLT events in data depending if we are using the HighEG dataset or the LowEG dataset
    if(!isMC && dataType == 0){
      if(!(t.HLT_HIAK4PFJet100_v1 || t.HLT_HIAK4PFJet80_v1)) continue;
    }
	
	
    if(!isMC && dataType == -1){
      if(!((t.HLT_HIAK4PFJet60_v1 == 1 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) ||
	   (t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0))) continue;
    }

    //select HLT events with at least 40 GeV if MC
    //if(isMC){
      //if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;
    //}

    if(isMC){if(!(t.HLT_HIAK4PFJet80_v1 == 1)) continue;} //minitest

    // Loop over jets
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {  
      
      //some cuts
      if (std::abs(t.refeta[ijet]) > 1.9) continue;
      if ((isMC) && skipMC(t.refpt[ijet], t.pthat)) continue;
      if (t.refpt[ijet] < pT_low || t.refpt[ijet] > pT_high) continue;
      if (btag && t.discr_particleNet_BvsAll[ijet] <= 0.99) continue;
      
     

// ---------- CASE 0 B: ----------
if (t.jtNbHad[ijet] == 0){
  counter0 += 1 ;
  beforeskip += 1 ; 
        
  // new function to make vertices considering every track 
  vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec = makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail); //makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_0b_score_bkg, h_0b_score_sg);
        
  totjt += 1 ; 
        
  if (reco_hadrons_4vec.size() !=  2) continue;

  double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
  double pt1 = reco_hadrons_4vec[0].Pt();
  double pt2 = reco_hadrons_4vec[1].Pt();
  double eec = std::pow(pt1 * pt2, n);
  double jtpt = t.jtpt[ijet];
  double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

  h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree);
}

// ---------- CASE 1 B: ----------
      if (t.jtNbHad[ijet] == 1){
        counter1 += 1 ;
        beforeskip += 1 ; 
        
       // new function to make vertices considering every track 
       vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail); //makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_1b_score_bkg, h_1b_score_sg);
        
        totjt += 1 ; 
        
        
        if (reco_hadrons_4vec.size() !=  2) continue;

        double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
        double pt1 = reco_hadrons_4vec[0].Pt();
        double pt2 = reco_hadrons_4vec[1].Pt();
        double eec = std::pow(pt1 * pt2, n);
        double jtpt = t.jtpt[ijet];
        double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

      h3D_1b->Fill(mB, dr, jtpt, eec * weight_tree);
      }

      // ---------- CASE 2 B: ----------
      if (t.jtNbHad[ijet] == 2){
        counter2 += 1 ;
        beforeskip += 1 ; 
        vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =   makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail); //makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_2b_score_bkg, h_2b_score_sg);
        totjt += 1 ; 
      
      
        if (reco_hadrons_4vec.size() !=  2) continue;
        double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
        double pt1 = reco_hadrons_4vec[0].Pt();
        double pt2 = reco_hadrons_4vec[1].Pt();
        double eec = std::pow(pt1 * pt2, n);
        double jtpt = t.jtpt[ijet];
        double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

        h3D_2b->Fill(mB, dr, jtpt, eec * weight_tree);
    
      }

  } //close jet loop
  } //close event loop
  

  // -------- Save histograms -------
  TString label = "_btag";
  if(!btag) label = "_nobtag"; 

 label += "_truth"; //"BDT_cut_test_0365Cut_nonegstatus"; //if you want to add the BDT cut in the name of the output file

  TFile outFile( (output_folder + output_hist + label + domain).Data(), "RECREATE");                                                                                                                                             
  
  h3D_0b->Write();
  h3D_1b->Write();
  h3D_2b->Write();
  
  h_0b_score_bkg->Write();
  h_1b_score_bkg->Write();
  h_2b_score_bkg->Write();
  h_0b_score_sg->Write();
  h_1b_score_sg->Write();
  h_2b_score_sg->Write();


  outFile.Close();
}


void filter_b_bb_as_data_and_mc(TString filename_bjet, TString output_folder, TString output_hist, TString domain, Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC) {
  
  // true level information aggregated to partial Bs
  tTree t;
  t.Init(filename_bjet, isMC);
  t.SetBranchStatus("*", 0);
  double agg_fail = 0;
  double nb_sv=0 ;
  double lt2sv = 0 ;
  double totjt = 0;
  double beforeskip = 0 ;    
  double vertskip = 0 ;
  double matchskip = 0;
  double merge_fail = 0 ;
  double sv_fail = 0 ;
  double global_count = 0 ; 

  std::vector<TString> active_branches = {
      // common branches (data + MC)
      "jtpt", "jteta", "jtphi", "jtm", "nref", "jtNtrk", "jtNsvtx", "discr_particleNet_BvsAll",
      "ntrk", "trkJetId", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "trkSvtxId", "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt", "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
      "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1"};
  if (isMC) {
    std::vector<TString> mc_branches = {"weight", "pthat"};
    active_branches.insert(active_branches.end(), mc_branches.begin(), mc_branches.end());
  }
  t.SetBranchStatus(active_branches, 1);
  
   
  // Plots or histograms
  TH3D *h3D_bb_bkg = new TH3D("h3D_bb_bkg", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b_bkg = new TH3D("h3D_b_bkg", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  TH3D *h3D_bb_signal = new TH3D("h3D_bb_signal", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b_signal = new TH3D("h3D_b_signal", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    

  TH2F* global_histo = new TH2F("global_histo", "distribution of wrong reco with pt; pT; reco proportion", 40, 0, 180, 40, 0, 1.2);
  
  TH1D *h_0b_score = new TH1D("h_0b_score", "Distribution of Score", 50, -1.5, 1.5);
  TH1D *h_1b_score = new TH1D("h_1b_score", "Distribution of Score", 50, -1.5, 1.5);
 
  h3D_bb_bkg->Sumw2();
  h3D_b_bkg->Sumw2();
  h3D_bb_signal->Sumw2();
  h3D_b_signal->Sumw2();
    
  TRandom3 rnd(12345);
  Long64_t n_events = t.GetEntries();                                                                                                
  for (Long64_t ient = 0; ient < n_events; ient++) {
    //Progress
    if (ient % 50000 == 0) {                                                                                                         
      float percent = 100.0 * ient / n_events;
      std::cout << "\rProcessing: "  << percent << " %" << std::flush;
    }
    t.GetEntry(ient);
    double weight_tree = isMC ? t.weight : 1.0;

    bool isTemplate = (rnd.Uniform() < 0.5);

    // Loop over jets — identical logic for data and MC
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      if (std::abs(t.jteta[ijet]) > 1.9) continue;
      if (isMC && skipMC(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < pT_low || t.jtpt[ijet] > pT_high) continue;
      if (btag && t.discr_particleNet_BvsAll[ijet] <= 0.99) continue;

      // classify by number of reconstructed secondary vertices (same for data and MC)
      int nSV = t.jtNsvtx[ijet];
      if (nSV < 1) continue;

      beforeskip += 1;
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec = makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail);
      totjt += 1;
      if (reco_hadrons_4vec.empty()) lt2sv += 1;
      if (reco_hadrons_4vec.size() != 2) continue;

      double dr   = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
      double pt1  = reco_hadrons_4vec[0].Pt();
      double pt2  = reco_hadrons_4vec[1].Pt();
      double eec  = std::pow(pt1 * pt2, n);
      double jtpt = t.jtpt[ijet];
      double mB   = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

      if (nSV == 1) {
        if (!isTemplate) h3D_b_signal->Fill(mB, dr, jtpt, eec * weight_tree);
        else             h3D_b_bkg->Fill(mB, dr, jtpt, eec * weight_tree);
      } else {  // nSV >= 2
        if (!isTemplate) h3D_bb_signal->Fill(mB, dr, jtpt, eec * weight_tree);
        else             h3D_bb_bkg->Fill(mB, dr, jtpt, eec * weight_tree);
      }

  } //close jet loop
  } //close event loop


  // -------- Save histograms -------
  TString label = "_btag";
  if(!btag) label = "_nobtag"; 

  TFile outFile( (output_folder + output_hist + label + domain).Data(), "RECREATE");                                                                                                                                             
  h3D_b_signal->Write();
  h3D_b_bkg->Write();
  h3D_bb_signal->Write();
  h3D_bb_bkg->Write();


  h_0b_score->Write();
  h_1b_score->Write();


  outFile.Close();

}




// Build templates for the template fit.
// MC:   fills h3D_b (jtNbHad==1) and h3D_bb (jtNbHad==2) using reco-level cuts and reco SV reconstruction.
// Data: fills h3D_data with the same reco logic — no truth classification.
void make_templates(TString filename, TString output_folder, TString output_hist, TString domain,
                    Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC, Int_t dataType,
                    Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1) {

  tTree t;
  t.Init(filename, isMC);
  t.SetBranchStatus("*", 0);

  double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0;

  std::vector<TString> active_branches = {
    // reco branches — identical for data and MC
    "jtpt", "jteta", "nref", "jtNtrk", "jtNsvtx", "discr_particleNet_BvsAll",
    "ntrk", "trkJetId", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
    "trkSvtxId", "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",
    "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
    "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1"};
  if (isMC) {
    // MC-only: event weight, pthat for skipMC, and truth b-hadron count for template classification
    std::vector<TString> mc_branches = {"weight", "pthat", "jtNbHad"};
    active_branches.insert(active_branches.end(), mc_branches.begin(), mc_branches.end());
  }
  t.SetBranchStatus(active_branches, 1);

  // MC: separate 0b, b and bb templates
  TH3D *h3D_0b = new TH3D("h3D_0b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b  = new TH3D("h3D_b",  "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_bb = new TH3D("h3D_bb", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  // Data: single distribution to be fit
  TH3D *h3D_data = new TH3D("h3D_data", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);

  h3D_0b->Sumw2();   h3D_0b->SetCanExtend(TH1::kNoAxis);
  h3D_b->Sumw2();    h3D_b->SetCanExtend(TH1::kNoAxis);
  h3D_bb->Sumw2();   h3D_bb->SetCanExtend(TH1::kNoAxis);
  h3D_data->Sumw2(); h3D_data->SetCanExtend(TH1::kNoAxis);

  // Jet counts (no EEC weight): 3D (mB, dr, jtpt) — same axes as the EEC histograms
  TH3D *h_count_0b   = new TH3D("h_count_0b",   "jet counts 0b;m_{B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_b    = new TH3D("h_count_b",    "jet counts 1b;m_{B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_bb   = new TH3D("h_count_bb",   "jet counts 2b;m_{B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_data = new TH3D("h_count_data", "jet counts data;m_{B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h_count_0b->Sumw2();   h_count_0b->SetCanExtend(TH1::kNoAxis);
  h_count_b->Sumw2();    h_count_b->SetCanExtend(TH1::kNoAxis);
  h_count_bb->Sumw2();   h_count_bb->SetCanExtend(TH1::kNoAxis);
  h_count_data->Sumw2(); h_count_data->SetCanExtend(TH1::kNoAxis);

  Long64_t n_events = t.GetEntries();
  if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
  if (ev_first < 0) ev_first = 0;
  std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

  for (Long64_t ient = ev_first; ient < ev_last; ient++) {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
    t.GetEntry(ient);

    double weight_tree = isMC ? t.weight : 1.0;

    // trigger selection
    if (!isMC && dataType == 0) {
      if (!(t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1)) continue; }
    
      else if (!isMC && dataType == -1) {
      if (t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1) continue; 
      if (!(t.HLT_HIAK4PFJet40_v1 || t.HLT_HIAK4PFJet60_v1 )) continue;
    }
    else if (isMC) {if (!(t.HLT_HIAK4PFJet40_v1)) continue;}

    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // reco-level cuts — identical for data and MC
      if (std::abs(t.jteta[ijet]) > 1.9) continue;
      if (isMC && skipMC(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < pT_low || t.jtpt[ijet] > pT_high) continue;
      if (btag && t.discr_particleNet_BvsAll[ijet] <= 0.99) continue;

      // reco SV reconstruction — same for data and MC
      vector<ROOT::Math::PtEtaPhiMVector> reco_sv = makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr);
      if (reco_sv.size() < 2) continue;

      double dr   = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(), reco_sv[1].Eta(), reco_sv[1].Phi());
      double pt1  = reco_sv[0].Pt();
      double pt2  = reco_sv[1].Pt();
      double eec  = std::pow(pt1 * pt2, n);
      double jtpt = t.jtpt[ijet];
      double mB   = reco_sv[0].M() + reco_sv[1].M();
      if (mB > mb_max_fill) mB = mb_max_fill;  // fold overflow into last bin



      std::cout << "weight: " << weight_tree << std::endl;
      std::cout << "eec: " << eec << std::endl;

      //Fix the under/overflow
      if(dr < dr_min) dr = dr_min_fill;
      if(dr >= dr_max) dr = dr_max_fill;
      if(mB >= mb_max) mB = mb_max_fill;



      std::cout << "weight: " << weight_tree << std::endl;
      std::cout << "eec: " << eec << std::endl;

      if (isMC) {
        // use truth to classify: fill separate 0b, b and bb templates
        if      (t.jtNbHad[ijet] == 0) { h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree); h_count_0b->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 1) { h3D_b ->Fill(mB, dr, jtpt, eec * weight_tree); h_count_b ->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 2) { h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree); h_count_bb->Fill(mB, dr, jtpt, weight_tree); }
      } else {
        h3D_data->Fill(mB, dr, jtpt, eec * weight_tree);
        h_count_data->Fill(mB, dr, jtpt, weight_tree);
      }
    }
  }
  std::cout << std::endl;

  TString label = btag ? "_btag" : "_nobtag";
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
  TFile outFile((output_folder + output_hist + label + job_suffix + domain).Data(), "RECREATE");
  if (isMC) {
    h3D_0b->Write();
    h3D_b->Write();
    h3D_bb->Write();
    h_count_0b->Write();
    h_count_b->Write();
    h_count_bb->Write();
  } else {
    h3D_data->Write();
    h_count_data->Write();
  }
  outFile.Close();
}


//Step 1: filter bb from b. Only MC
//Step 2: filter bb from b, but split the sample in 2 and treat one as data and one as MC (to be used as template fit input)


void create_files_for_template_fit(Int_t dataType = 1, Float_t pT_low = 80, Float_t pT_high = 140, Int_t n = 1, bool btag = true, bool isMC = true){
 
TString filename;
TString output_hist;
TString output_folder = "/data_CMS/cms/zaidan/analysis_lise/";
TString domain = ".root";

//sanity check
if (isMC && dataType < 1) {
  std::cerr << "Invalid data type for MC sample" << std::endl;
  return;}

if (!isMC && dataType > 1) {
  std::cerr << "Invalid data type for data sample" << std::endl;
  return;}

 


if(dataType == -1){//________________________________data______________________________
  filename = "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root";
  output_hist = "template_for_fit_histos_3D_LowEG_f";
  isMC = false;
  cout<<"you chose data Low" <<endl;
  }

else if(dataType == 0) {
  filename = "/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
  output_hist = "template_for_fit_histos_3D_HighEG_f";
  isMC = false;
  cout<<"you chose data High" <<endl;       
  }      
                                                                                                                                                                                                                                                                        
else if(dataType == 1){//________________________________bjet______________________________
  filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
  output_hist = "template_for_fit_histos_3D_bjet_f";
  std::cout << "Creating files for template fit for bjet sample" << std::endl;
  cout<<"you chose bjet MC" <<endl;
  }

else if(dataType == 2){//________________________________dijet______________________________
  filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root"; 
  output_hist = "template_for_fit_histos_3D_qcd_f";
  std::cout << "Creating files for template fit for qcd sample" << std::endl;
  cout<<"you chose qcd MC" <<endl;
  }



else{
  cout<<" undefined data type"<<endl;
  return; 
  }




  make_templates(filename, output_folder, output_hist, domain, pT_low, pT_high, n, btag, isMC, dataType);
  //filter_b_bb(filename, output_folder, output_hist, domain, pT_low, pT_high, n, btag, isMC, dataType);
  //filter_b_bb_as_data_and_mc(filename, output_folder, output_hist, domain, pT_low, pT_high, n, btag, isMC);


  std::cout << "finished :)" << std::endl;
}

