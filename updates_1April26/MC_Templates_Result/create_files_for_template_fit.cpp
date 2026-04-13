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
#include "tTree.h"
#include "binning_histos_small.h"

// -- Global writing directory
TString smyData =  Form("/home/llr/cms/shatat/CMSAnalysis/EECs/MC_Templates_Result/");

struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    std::vector<int> trkMatchSta;
  vector<int> trkSvtxId;   
};


//given two vectors, computes deltaR squared 
auto deltaR2 = [](const ROOT::Math::PtEtaPhiMVector &a,
                  const ROOT::Math::PtEtaPhiMVector &b) {
    double dEta = a.Eta() - b.Eta();
    double dPhi = std::acos(std::cos(a.Phi() - b.Phi()));
    return dEta*dEta + dPhi*dPhi;
};


//printing function for a Vertex object 
void printvtx (const Vertex& vertex,
              const TString& txt,
	       int vtxnb,  Int_t& ijet, Long64_t& ient) {
  for (int itrk = 0; itrk < vertex.trkSvtxId.size(); itrk++   ) {
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
  std::vector<Int_t> vec ;
                                                                                                              
  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {       //loop over trks                                                                             
    // Track must belong to this jet
    if (t.trkJetId[itrk] != ijet) continue;
    
                                                                                
    // pT cut                                                                                                                         
    if (t.trkPt[itrk] < 1) continue;
    if (t.trkMatchSta[itrk] <100 ) continue;
    //if (t.trkSvtxId[itrk] <0 ) continue ;      //trk must be in SV 
    int value = t.trkSvtxId[itrk] ;
    if (std::find(vec.begin(), vec.end(), value) == vec.end()) {
      vec.push_back(value);
    }
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

  std::map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;    // map from vertex id to list of tracks :  secVtxs[svId] = list of track 4 vec in that vertex 
  std::map<Int_t, std::vector<int>> secVtxsMatchSta; // store trkMatchSta per vertex
  vector<ROOT::Math::PtEtaPhiMVector> empty;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;

   double weight_tree = t.weight ; 
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
        if      (pid == 211)  v1.SetM(0.139570); // pion 
        else if (pid == 13)   v1.SetM(0.105658); // muon 
        else if (pid == 11)   v1.SetM(0.000510); // electron 
        else if (pid == 2212) v1.SetM(0.938272); // proton 
        else if (pid == 321)  v1.SetM(0.493677); // kion 
        else if (pid == 3112) v1.SetM(1.19744);  // Sigma-                                                                                         
	      else if (pid == 3222) v1.SetM(1.18937);  // sigma +                                                                                  
	      else if (pid == 3312) v1.SetM(1.32171);  // xi -                                                                                           
	      else if (pid == 3334) v1.SetM( 1.67245); // Omega - 
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
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second) 
	vtx.p4 += it.second[i];                       //sums the tracks' 4 vector into vtx.p4   (it.second[i] = 4 vec of track i 
	vtx.tracks.push_back(it.second[i]);           //copies the track itself in vtx.tracks 
	vtx.trkMatchSta.push_back(secVtxsMatchSta[it.first][i]); //look at secVtxsMatchSta[vertex id][track i] and stores it in trkMatchSta (with the same index as the track saved above
	vtx.trkSvtxId.push_back(originalSvId);  //store provenance 
        }
      vertices.push_back(vtx);  //adds the Vertex vtx to the vector vertices 
    }
    // what we have created is tracks[i] ↔ trkMatchSta[i] ↔ trkSvtxId[i] correspondance :  originalSvId = vertex.trkSvtxId[i]

       
          

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices);
    if (vertices.size() != 2) {
       // std::cout<<" merging pb"<<std::endl; 
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
    // std::cout << "Vertex 0 track proportions:\n";
    int maxCount0 = 0;
    int maxSta0   = -10;
    for (const auto &kv : countMap0) {
      int sta0   = kv.first;
      int count0 = kv.second;
      if (count0 > maxCount0) {
        maxCount0 = count0;
        maxSta0   = sta0;
      }
    }
    double prop0 = double (maxCount0)/double(totalTracks0);
if (prop0 != 1 ) {
      merge_fail += 1 ;
      
    }
    
    //VERTEX 1
std::map<int,int> countMap1;
    for (size_t i = 0; i < vertex1.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex 
        int sta1 = vertex1.trkMatchSta[i];
        countMap1[sta1]++;
    }
    int totalTracks1 = vertex1.tracks.size();
    // std::cout << "Vertex 0 track proportions:\n";
    int maxCount1 = 0;
    int maxSta1   = -10;
    for (const auto &kv : countMap1) {
      int sta1   = kv.first;
      int count1 = kv.second;
      if (count1 > maxCount1) {
        maxCount1 = count1;
        maxSta1   = sta1;
      }
    }
    double prop1 = double (maxCount1)/double(totalTracks1);
    if (prop1 != 1 ) {
      merge_fail += 1 ;
      //printvtx( vertex1 , " wrong merging ", 1); 
      //return empty ;
    }

//Add in the particles that are not part of a sv
    int i = 0 ; 
    for (const auto& v1 : no_sv_list) {
      //std::cout<<"add tracks "<<no_sv_list.size()<<std::endl;
      double d0 = deltaR2 (v1, vertex0.p4);
      double d1 = deltaR2 (v1, vertex1.p4);
      if (d0 < d1 ) {
	//std::cout <<"adding coord "<< v1.Phi()<<v1.Eta() <<"  to v0 "<<std::endl; 
	vertex0.p4 += v1;

	if (no_sv_sta_list[i] != vertex0.trkMatchSta[0] ) {
	  agg_fail += 1 ;
	  //return empty;
	}
      }
      else if ( d0 > d1 ) {
	vertex1.p4 += v1 ;
	//std::cout <<"adding coord "<< v1.Phi()<<v1.Eta()<<"  to v1 "<<std::endl; 
	if (no_sv_sta_list[i] != vertex1.trkMatchSta[0] ) {
	  agg_fail += 1 ;
	  //return empty;
	}
      }
      i+=1 ; 
    }
    no_sv_list.clear();

std::set<int> originalVertices0(  vertices[0].trkSvtxId.begin(), vertices[0].trkSvtxId.end()    );   // keeps only one copy of each sv nb 
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

  std::map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;    // map from vertex id to list of tracks :  secVtxs[svId] = list of track 4 vec in that vertex 
  std::map<Int_t, std::vector<int>> secVtxsMatchSta; // store trkMatchSta per vertex
  vector<ROOT::Math::PtEtaPhiMVector> empty;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;

   double weight_tree = t.weight ; 
    // Loop over tracks

    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //cuts for tracks 
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;

    //control plots for the BDT score distribution in 0b and 1b jets
    if (t.trkMatchSta[itrk] == 1)
      {h_score_bkg->Fill(t.trkBdtScore[itrk], weight_tree);}
    else if (t.trkMatchSta[itrk] >= 100)      
      {h_score_sg->Fill(t.trkBdtScore[itrk], weight_tree);}
  


    // -- Add trkbdt score
    if (t.trkBdtScore[itrk]<=-0.3) continue;
	
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
	   return empty;
      }
      Vertex vtx;
  
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second) 
	vtx.p4 += it.second[i];                       //sums the tracks' 4 vector into vtx.p4   (it.second[i] = 4 vec of track i 
	vtx.tracks.push_back(it.second[i]);           //copies the track itself in vtx.tracks 
	vtx.trkMatchSta.push_back(secVtxsMatchSta[it.first][i]); //look at secVtxsMatchSta[vertex id][track i] and stores it in trkMatchSta (with the same index as the track saved above
	vtx.trkSvtxId.push_back(originalSvId);  //store provenance 
        }
      vertices.push_back(vtx);  //adds the Vertex vtx to the vector vertices 
    }
    // what we have created is tracks[i] ↔ trkMatchSta[i] ↔ trkSvtxId[i] correspondance :  originalSvId = vertex.trkSvtxId[i]

       
          

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices);
    if (vertices.size() != 2) {
       // std::cout<<" merging pb"<<std::endl; 
       return empty;
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
    // std::cout << "Vertex 0 track proportions:\n";
    int maxCount0 = 0;
    int maxSta0   = -10;
    for (const auto &kv : countMap0) {
      int sta0   = kv.first;
      int count0 = kv.second;
      if (count0 > maxCount0) {
        maxCount0 = count0;
        maxSta0   = sta0;
      }
    }
    double prop0 = double (maxCount0)/double(totalTracks0);
if (prop0 != 1 ) {
      merge_fail += 1 ;
      
    }
    
    //VERTEX 1
std::map<int,int> countMap1;
    for (size_t i = 0; i < vertex1.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex 
        int sta1 = vertex1.trkMatchSta[i];
        countMap1[sta1]++;
    }
    int totalTracks1 = vertex1.tracks.size();
    // std::cout << "Vertex 0 track proportions:\n";
    int maxCount1 = 0;
    int maxSta1   = -10;
    for (const auto &kv : countMap1) {
      int sta1   = kv.first;
      int count1 = kv.second;
      if (count1 > maxCount1) {
        maxCount1 = count1;
        maxSta1   = sta1;
      }
    }
    double prop1 = double (maxCount1)/double(totalTracks1);
    if (prop1 != 1 ) {
      merge_fail += 1 ;
      //printvtx( vertex1 , " wrong merging ", 1); 
      //return empty ;
    }

//Add in the particles that are not part of a sv
    int i = 0 ; 
    for (const auto& v1 : no_sv_list) {
      //std::cout<<"add tracks "<<no_sv_list.size()<<std::endl;
      double d0 = deltaR2 (v1, vertex0.p4);
      double d1 = deltaR2 (v1, vertex1.p4);
      if (d0 < d1 ) {
	//std::cout <<"adding coord "<< v1.Phi()<<v1.Eta() <<"  to v0 "<<std::endl; 
	vertex0.p4 += v1;

	if (no_sv_sta_list[i] != vertex0.trkMatchSta[0] ) {
	  agg_fail += 1 ;
	  //return empty;
	}
      }
      else if ( d0 > d1 ) {
	vertex1.p4 += v1 ;
	//std::cout <<"adding coord "<< v1.Phi()<<v1.Eta()<<"  to v1 "<<std::endl; 
	if (no_sv_sta_list[i] != vertex1.trkMatchSta[0] ) {
	  agg_fail += 1 ;
	  //return empty;
	}
      }
      i+=1 ; 
    }
    no_sv_list.clear();

    std::set<int> originalVertices0(  vertices[0].trkSvtxId.begin(), vertices[0].trkSvtxId.end()    );   // keeps only one copy of each sv nb 
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
    
  





// - - - - - - - - - - -  USING FUNCTIONS  - - - - - - - - - - -
void filter_b_bb(TString filename, TString output_folder, TString output_hist, TString domain, Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC) {
  
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
      "ntrk", "trkJetId", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "refTrkPdgId","refTrkSta", "refTrkMass",  "refmB", "refpt", "refeta", "refphi",                                                                                           
      "nrefTrk", "refTrkJetId", "refTrkPdgId", "refTrkPt", "refTrkEta", "refTrkPhi", "refTrkY", "refTrkSta", "refTrkMass", "refNtrk",
      "jtNsvtx", "trkSvtxId",  "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",  "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2"};
    t.SetBranchStatus(active_branches, 1);
  
   
  // Plots or histograms
  TH3D *h3D_2b = new TH3D("h3D_2b", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_1b = new TH3D("h3D_1b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  TH3D *h3D_0b = new TH3D("h3D_0b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  
  TH1D *h_0b_score_bkg = new TH1D("h_0b_score_bkg", "Distribution of Score", 50, -1, 1);
  TH1D *h_1b_score_bkg = new TH1D("h_1b_score_bkg", "Distribution of Score", 50, -1, 1);
  TH1D *h_2b_score_bkg = new TH1D("h_2b_score_bkg", "Distribution of Score", 50, -1, 1);
  
  TH1D *h_0b_score_sg = new TH1D("h_0b_score_sg", "Distribution of Score", 50, -1, 1);
  TH1D *h_1b_score_sg = new TH1D("h_1b_score_sg", "Distribution of Score", 50, -1, 1);
  TH1D *h_2b_score_sg = new TH1D("h_2b_score_sg", "Distribution of Score", 50, -1, 1);
  

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

    // Loop over jets
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {  
      
      //some cuts
      if (std::abs(t.refeta[ijet]) > 1.9) continue;
      if ((isMC) && skipMC(t.refpt[ijet], t.pthat)) continue;
      if (std::abs(t.refpt[ijet]) < pT_low || std::abs(t.refpt[ijet]) > pT_high) continue;
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;
      
     

// ---------- CASE 0 B: ----------
if (t.jtNbHad[ijet] == 0){
  counter0 += 1 ;
  beforeskip += 1 ; 
        
  // new function to make vertices considering every track 
  vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_0b_score_bkg, h_0b_score_sg);
        
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
       vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_1b_score_bkg, h_1b_score_sg);
        
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
        vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =   makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_2b_score_bkg, h_2b_score_sg);
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

 label += "_BDT_cut_test_LidasCut_nonegstatus"; //if you want to add the BDT cut in the name of the output file

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

  std::vector<TString> active_branches = {"weight","jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "jtNtrk",
      "ntrk", "trkJetId", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "refTrkPdgId","refTrkSta", "refTrkMass",  "refmB", "refpt", "refeta", "refphi",                                                                                           
      "nrefTrk", "refTrkJetId", "refTrkPdgId", "refTrkPt", "refTrkEta", "refTrkPhi", "refTrkY", "refTrkSta", "refTrkMass", "refNtrk",
      "jtNsvtx", "trkSvtxId",  "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",  "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",};
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
    //get the MC event weight
    double weight_tree = t.weight;


    bool isTemplate = (rnd.Uniform() < 0.5);

    // Loop over jets
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {  
      
      //some cuts
      if (std::abs(t.refeta[ijet]) > 1.9) continue;
      if ((isMC) && skipMC(t.refpt[ijet], t.pthat)) continue;
      if (std::abs(t.refpt[ijet]) < pT_low || std::abs(t.refpt[ijet]) > pT_high) continue;
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;
      
      std::vector<Int_t> hadrons_stat;
      std::vector<ROOT::Math::PtEtaPhiMVector> hadrons_4vec;
      PartialBsAggregation(hadrons_4vec, hadrons_stat, t, ijet);

      size_t nB_aggr = 0;
      for (size_t k = 0; k < hadrons_stat.size(); ++k) {
        if (hadrons_stat[k] >= 100)
          ++nB_aggr;
            }

      if (nB_aggr == 1){

      beforeskip += 1 ; 
    
      Float_t SVntrks1 = 0; 
      Float_t SVntrks2 = 0 ; 
      //tracks aggregation : takes (empty vector, t tree, jet nb) and fills the vector reco_hadrons_4vec with the2 reconstructed B
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail);
        
      //makeSvtxs( t,  ijet, SVntrks1, SVntrks2, global_histo, global_count );
      totjt += 1 ; 
      if (reco_hadrons_4vec.empty()) lt2sv += 1 ;
      
      if (reco_hadrons_4vec.size() !=  2) continue;
      double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
      double pt1 = reco_hadrons_4vec[0].Pt();
      double pt2 = reco_hadrons_4vec[1].Pt();
      double eec = std::pow(pt1 * pt2, n);
      double sumpt = pt1 + pt2;
      double jtpt = t.jtpt[ijet];
      double diffpt = std::abs(pt1-pt2);
      double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();


      if (!isTemplate) h3D_b_signal->Fill(mB, dr, jtpt, eec * weight_tree);
      else h3D_b_bkg->Fill(mB, dr, jtpt, eec * weight_tree);
  
      //if (dr < 0.1 ) continue ; 
      //h3D_b->Fill(mB, dr, jtpt, eec * weight_tree);                                                                                                                                                                                                                      
      }

      if (nB_aggr == 2){
       
         beforeskip += 1 ; 
      // -- now we start the reco step --
      Float_t SVntrks1 = 0; 
      Float_t SVntrks2 = 0 ; 
                //tracks aggregation : takes (empty vector, t tree, jet nb) and fills the vector reco_hadrons_4vec with the2 reconstructed B
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail);
         //makeSvtxs( t,  ijet, SVntrks1, SVntrks2, global_histo, global_count );
      totjt += 1 ; 
      if (reco_hadrons_4vec.empty()) {
	lt2sv += 1 ;
      } 

      
      if (reco_hadrons_4vec.size() !=  2) continue;
      double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(),
                            reco_hadrons_4vec[0].Phi(),
                            reco_hadrons_4vec[1].Eta(),
                            reco_hadrons_4vec[1].Phi()
                            );
      double pt1 = reco_hadrons_4vec[0].Pt();
      double pt2 = reco_hadrons_4vec[1].Pt();
      double eec = std::pow(pt1 * pt2, n);
      double sumpt = pt1 + pt2;
      double jtpt = t.jtpt[ijet];
      double diffpt = std::abs(pt1-pt2);
      double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

      //if (dr < 0.1 ) continue ; 

      if (!isTemplate) h3D_bb_signal->Fill(mB, dr, jtpt, eec * weight_tree);
      else h3D_bb_bkg->Fill(mB, dr, jtpt, eec * weight_tree);

      //h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree);                                                                                                                                                                                                                      
      
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
  //h_2b_score->Write();


  outFile.Close();

}

// -------------------------------------------------------------------
// ---------------- Version 2 : March 25th 2026 ---------------------------------------
// ---------------- makeSvtxs_BDTcut_v2() simialr to function used for data ---------- 
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
       // std::cout<<" merging pb"<<std::endl; 
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
   
void filter_b_bb_v2(TString filename, TString output_folder, TString output_hist, TString domain, Float_t pT_low, Float_t pT_high, double &BDTcutTh,  Int_t n, bool btag, bool isMC, bool apply_prescale, bool applyHLT40, bool applyHLT80) {
  
  /* Sample split histograms are also added, for cross check in mc as closure test.*/

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
      "ntrk", "trkJetId", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "refTrkPdgId","refTrkSta", "refTrkMass",  "refmB", "refpt", "refeta", "refphi",                                                                                           
      "nrefTrk", "refTrkJetId", "refTrkPdgId", "refTrkPt", "refTrkEta", "refTrkPhi", "refTrkY", "refTrkSta", "refTrkMass", "refNtrk",
      "jtNsvtx", "trkSvtxId",  "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",  "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
       "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
    };
    // Trigger info from HltTree (added as friend)

    t.SetBranchStatus(active_branches, 1);
  
   TH1::SetDefaultSumw2();  // enable Sumw2 globally, instead of one by one 
   
  // Plots or histograms
  TH3D *h3D_2b = new TH3D("h3D_2b", ";m_{2B} [GeV];#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_1b = new TH3D("h3D_1b", ";m_{2B} [GeV];#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  TH3D *h3D_0b = new TH3D("h3D_0b", ";m_{2B} [GeV];#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector); 
  TH3D *h3D_more2b = new TH3D("h3D_more2b", ";m_{2B} [GeV];#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector); 


  // -- For sampled templates: Templates and Pseudodata 
  TH3D *h3D_pseudodata_rns = new TH3D("h3D_pseudodata_rns", ";m_{2B} [GeV];#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_2b_rns = new TH3D("h3D_2b_rns", ";m_{2B} [GeV];#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_1b_rns = new TH3D("h3D_1b_rns", ";m_{2B} [GeV];#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  TH3D *h3D_0b_rns = new TH3D("h3D_0b_rns", ";m_{2B} [GeV];#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);  
  TH3D *h3D_more2b_rns = new TH3D("h3D_more2b_rns", ";m_{2B} [GeV];#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector); 
  TH3D *h3D_pseudodata_ALL_rns = new TH3D("h3D_pseudodata_ALL_rns", ";m_{2B} [GeV];#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);

  
  TH1D *h_0b_score_bkg = new TH1D("h_0b_score_bkg", "Distribution of Score", 50, -1, 1);
  TH1D *h_1b_score_bkg = new TH1D("h_1b_score_bkg", "Distribution of Score", 50, -1, 1);
  TH1D *h_2b_score_bkg = new TH1D("h_2b_score_bkg", "Distribution of Score", 50, -1, 1);
  
  TH1D *h_0b_score_sg = new TH1D("h_0b_score_sg", "Distribution of Score", 50, -1, 1);
  TH1D *h_1b_score_sg = new TH1D("h_1b_score_sg", "Distribution of Score", 50, -1, 1);
  TH1D *h_2b_score_sg = new TH1D("h_2b_score_sg", "Distribution of Score", 50, -1, 1);

  // -- Jet40 prescale 
  //Save the prescale factor (only for 40 GeV trigger)
    double prescale_pf40 = 33.917210;

  Long64_t n_events = t.GetEntries();   

  // --------------------
  // ------- TEST ------- 
  // n_events = 10e+03; // test

  for (Long64_t ient = 0; ient < n_events; ient++) {

    // Decide split sample: template or pseudodata 
    bool isTemplate; // true for templates, false for pseudodata 
    if (ient %2 == 0)  isTemplate = true;
    else  isTemplate = false;

    //Progress
    if (ient % 50000 == 0) {                                                                                                         
      float percent = 100.0 * ient / n_events;
      std::cout << "\rProcessing: "  << percent << " %" << std::flush;
    }
    t.GetEntry(ient);

    //get the MC event weight
    double weight_tree = t.weight;

    // if(!(applyHLT40 && t.HLT_HIAK4PFJet40_v1 == 1)) continue; // select ets of jet pt > 40 online  // WRONG!!! 
    
    //Select MC passing a trigger of at least 40 GeV
    if(applyHLT40 && t.HLT_HIAK4PFJet40_v1 != 1) { continue; }// select ets of jet pt > 40 online 

    // -- another trigger 
    if(applyHLT80 && t.HLT_HIAK4PFJet80_v1 != 1){ continue; } 


    // Loop over jets
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {  
      
      //some cuts
      if (std::abs(t.refeta[ijet]) > 1.9) continue;
      if ((isMC) && skipMC(t.refpt[ijet], t.pthat)) continue;
      if (std::abs(t.refpt[ijet]) < pT_low || std::abs(t.refpt[ijet]) > pT_high) continue;
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;
      
      // common counters 
      if (t.jtNbHad[ijet] == 0){counter0 += 1 ; beforeskip += 1 ;}
      else if (t.jtNbHad[ijet] == 1){ counter1 += 1 ; beforeskip += 1 ; }
      else if (t.jtNbHad[ijet] == 2){  counter2 += 1 ; beforeskip += 1 ; }

      // Common aggregation 
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec_v1 =  makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_0b_score_bkg, h_0b_score_sg); // Zoe version, with some useful histograms
      
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec = makeSvtxs_BDTcut_v2(t, ijet, ient, BDTcutTh);
      // common counters    
      totjt += 1; 
      // common condition: want 2B 
      if (reco_hadrons_4vec.size() !=  2) continue;

      // common calculations 
      double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
      double pt1 = reco_hadrons_4vec[0].Pt();
      double pt2 = reco_hadrons_4vec[1].Pt();
      double eec = std::pow(pt1 * pt2, n);
      double jtpt = t.jtpt[ijet];
      double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();


      eec*= weight_tree; // MC weight due to pthat 


      // and apply the prescale 40 for events of jet40 and not jet 60, 0, 100 
       if(applyHLT40 && 
           t.HLT_HIAK4PFJet40_v1 == 1 &&
           t.HLT_HIAK4PFJet60_v1 == 0 &&
            t.HLT_HIAK4PFJet80_v1 == 0 &&
             t.HLT_HIAK4PFJet100_v1 == 0)
              {
                  if (apply_prescale) eec*=prescale_pf40;
              }



      // treatement of under/over flow bins 
        //Fix the under/overflow
        if(dr < dr_min) dr = dr_min_fill;
        if(dr >= dr_max) dr = dr_max_fill;
        if(mB >= mb_max) mB = mb_max_fill;
        if(eec >= eec_max) eec = eec_max_fill;


    // -- Fill templates: 
    // ---------- CASE 1 B: ----------
    if (t.jtNbHad[ijet] == 1){

      h3D_1b->Fill(mB, dr, jtpt, eec );
      // split sample 
      if(isTemplate){ h3D_1b_rns->Fill(mB, dr, jtpt, eec ); }
      else{           h3D_pseudodata_rns->Fill(mB, dr, jtpt, eec );
                      h3D_pseudodata_ALL_rns->Fill(mB, dr, jtpt, eec );

          }

    }// case 1 


    // ---------- CASE 2 B: ----------
    else if (t.jtNbHad[ijet] == 2){
       h3D_2b->Fill(mB, dr, jtpt, eec );

      // split sample 
      if(isTemplate){ h3D_2b_rns->Fill(mB, dr, jtpt, eec ); }
      else{           h3D_pseudodata_rns->Fill(mB, dr, jtpt, eec );
                      h3D_pseudodata_ALL_rns->Fill(mB, dr, jtpt, eec );
          }

    }// case 2 
    
    // ---------- CASE 0 B: ----------
    else if (t.jtNbHad[ijet] == 0){
        // no sample split      
        h3D_0b->Fill(mB, dr, jtpt, eec);
        // split sample 
        if(isTemplate){ h3D_0b_rns->Fill(mB, dr, jtpt, eec ); }
        else{           h3D_pseudodata_rns->Fill(mB, dr, jtpt, eec );
                        h3D_pseudodata_ALL_rns->Fill(mB, dr, jtpt, eec );
            }
    } // case 0

    // --- Case: more than 2B ----------  
    else if (t.jtNbHad[ijet] > 2){
        h3D_more2b->Fill(mB, dr, jtpt, eec );
       if(isTemplate){ h3D_more2b_rns->Fill(mB, dr, jtpt, eec );}
      else{           h3D_pseudodata_ALL_rns->Fill(mB, dr, jtpt, eec );}
    }


    } //close jet loop
  } //close event loop
  

  // -------- Save histograms -------
  TString label = "_btag";
  if(!btag) label = "_nobtag"; 

 label += "_BDT_cut"; //if you want to add the BDT cut in the name of the output file

  TFile outFile( (output_folder + output_hist + label + domain).Data(), "RECREATE");                                                                                                                                             
  h3D_0b->Write();
  h3D_1b->Write();
  h3D_2b->Write();
  h3D_more2b->Write();

  
  h_0b_score_bkg->Write();
  h_1b_score_bkg->Write();
  h_2b_score_bkg->Write();
  h_0b_score_sg->Write();
  h_1b_score_sg->Write();
  h_2b_score_sg->Write();

  // Splitted samples histograms 
  h3D_pseudodata_rns ->Write();
  h3D_pseudodata_ALL_rns ->Write();

  h3D_0b_rns->Write();
  h3D_1b_rns->Write();
  h3D_2b_rns->Write();
  h3D_more2b_rns ->Write();



  outFile.Close();
}  


//Step 1: filter bb from b. Only MC
//Step 2: filter bb from b, but split the sample in 2 and treat one as data and one as MC (to be used as template fit input)


void create_files_for_template_fit( Float_t pT_low = 80, Float_t pT_high = 140, Int_t n = 1, bool btag = true, bool isMC = true, bool apply_prescale = false, bool applyHLT40 = false, bool applyHLT80 = false){
 //bjet = true => bjet sample,
 // bjet = false => sample is qcd generic sample
     // TString output_folder = gSystem->ExpandPathName("$mydata/analysis_lise/");
      TString domain = ".root";

      // make directory for aggregation using BDT cuts, with the value 
      double BDTcutTh = 0.365; // set first then run code 
      TString sDirname = Form("ALL_BDT_%03d", int(BDTcutTh * 1000));
      gSystem->mkdir(sDirname, kTRUE);
      TString output_folder = Form("%s%s/", smyData.Data(), sDirname.Data()); // personal directory + BDT value name 

      // bjet 
      TString filename_b = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      TString output_hist_b = "template_for_fit_histos_3D_step_1_bjet";
      std::cout << "Creating files for template fit for bjet sample" << std::endl;
     
      // dijet 
      TString filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      TString output_hist = "template_for_fit_histos_3D_step_1_qcd";
      std::cout << "Creating files for template fit for qcd sample" << std::endl;
/*
// -- code modified to have sample splitting in the same function ()
    // --- Run code 
    // test dijet 
    applyHLT40 = true; 
    apply_prescale = true;

    filter_b_bb_v2(filename, output_folder, output_hist, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40);
    // test bjet 
    filter_b_bb_v2(filename_b, output_folder, output_hist_b, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40);
    std::cout << "finished with prescale applied :)" << std::endl;


  // -- test turn off applying the prescales
  applyHLT40 = true; 
  apply_prescale = false;
    output_hist_b = "template_for_fit_histos_3D_step_1_bjet_HLT_noprescale";
    output_hist = "template_for_fit_histos_3D_step_1_qcd_HLT_noprescale";
    // test dijet 
    filter_b_bb_v2(filename, output_folder, output_hist, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40);
    // test bjet 
    filter_b_bb_v2(filename_b, output_folder, output_hist_b, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40);
    std::cout << "finished  without prescale applied:)" << std::endl;

  // -- test without the HLT40 
  applyHLT40 = false;
  apply_prescale = false;
    output_hist_b = "template_for_fit_histos_3D_step_1_bjet_noHLT";
    output_hist = "template_for_fit_histos_3D_step_1_qcd_noHLT";
    // test dijet 
    filter_b_bb_v2(filename, output_folder, output_hist, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40);
    // test bjet 
    filter_b_bb_v2(filename_b, output_folder, output_hist_b, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40);
    std::cout << "finished  without HLT applied:)" << std::endl;
*/

  // -- Re-run code with HLT 40 GeV
  applyHLT40 = true;
  apply_prescale = true;
  applyHLT80 = false;
  cout << "applyHLT40 " << applyHLT40 <<endl;
  cout << "apply_prescale " << apply_prescale <<endl;
  cout << "applyHLT80 " << applyHLT80 <<endl;
  output_hist_b = "template_for_fit_histos_3D_step_1_bjet_HLT40_andprescale40";
  output_hist = "template_for_fit_histos_3D_step_1_qcd_HLT40_andprescale40";
  filter_b_bb_v2(filename, output_folder, output_hist, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40, applyHLT80);
  filter_b_bb_v2(filename_b, output_folder, output_hist_b, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40, applyHLT80);

  // --additional testing: MC with Jet80 trigger, no prescale
  applyHLT40 = false;
  apply_prescale = false;
  applyHLT80 = true;
  output_hist_b = "template_for_fit_histos_3D_step_1_bjet_HLT80_noprescale";
  output_hist = "template_for_fit_histos_3D_step_1_qcd_HLT80_noprescale";
  filter_b_bb_v2(filename, output_folder, output_hist, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40, applyHLT80);
  filter_b_bb_v2(filename_b, output_folder, output_hist_b, domain, pT_low, pT_high, BDTcutTh, n, btag, isMC, apply_prescale, applyHLT40, applyHLT80);


}

// -- tests without prescale application 


