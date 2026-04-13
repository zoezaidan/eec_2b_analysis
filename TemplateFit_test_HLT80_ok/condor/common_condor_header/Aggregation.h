
// --- Zoe code -----
struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    // std::vector<int> trkMatchSta; // to test for data! 
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
    // std::cout<<"ient : "<<ient <<"  ijet : "<<ijet << txt << "vertex number: "<<vtxnb<<"  track original vertex : "<< vertex.trkSvtxId[itrk] <<"    and track original status : "<<vertex.trkMatchSta[itrk] << " pT = " << vertex.tracks[itrk].Pt() <<  std::endl;
  }
}
//Skip MC events that have a too large weight to stabilize the distribution                                                          
bool skipMC(double pt, double pthat) {
  if (pthat<0.35*pt) return true;                                                                                                  
  return false;
}                                              

// -- Zoe fnctions 
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
    // auto &s1 = vertices[index1].trkMatchSta;
    // auto &s2 = vertices[index2].trkMatchSta;
    // s1.insert(s1.end(), s2.begin(), s2.end());

    //merge vertex id
    vertices[index1].trkSvtxId.insert(
    vertices[index1].trkSvtxId.end(),
        vertices[index2].trkSvtxId.begin(),
        vertices[index2].trkSvtxId.end());
    
    // remove merged vertex
    vertices.erase(vertices.begin() + index2);

    }
}


// -- Afnan: this is test version of zoe function (removed match status info)
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
    // -- All track match status related lines --> not relevant for data --> comment it 

    // cout << "Hello, this is the aggreagtion fuction without matched sta info" << endl;
  std::unordered_map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;
  // std::unordered_map<Int_t, std::vector<int>> secVtxsMatchSta;
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
      // no_sv_sta_list.push_back(t.trkMatchSta[itrk]); 
    } 
    //map vertex id to list of tracks
    else if (t.trkSvtxId[itrk] >= 0 ) {
      secVtxs[t.trkSvtxId[itrk]].push_back(v1);
      // secVtxsMatchSta[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]); / not relevant to data
    } 
    }

    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      if (nb_sv < 15 ) {
      } 
        // cout << "secVtx size < 2!" << endl;
      // cout << "original step is to return empty vector " << endl;
      return empty; 
         }

    // Build Vertex objects
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      Vertex vtx;
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      vtx.tracks.reserve(it.second.size());
      // vtx.trkMatchSta.reserve(it.second.size());
      vtx.trkSvtxId.reserve(it.second.size());
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second)
    vtx.p4 += it.second[i];
    vtx.tracks.push_back(it.second[i]);
    // vtx.trkMatchSta.push_back(v[i]);  // reuse already-looked-up reference
    vtx.trkSvtxId.push_back(originalSvId);
        }
      vertices.push_back(std::move(vtx));
    }

       
          

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices);
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl; 
       //return empty;
    }

    //vertices names
    auto &vertex0 = vertices[0];
    auto &vertex1 = vertices[1];

    //Add in the particles that are not part of a sv
    int i = 0;
    for (const auto& v1 : no_sv_list) {
      double d0 = deltaR2(v1, vertex0.p4);
      double d1 = deltaR2(v1, vertex1.p4);
      if (d0 < d1) {
    vertex0.p4 += v1;
    // if (no_sv_sta_list[i] != vertex0.trkMatchSta[0]) agg_fail += 1;
      }
      else if (d0 > d1) {
    vertex1.p4 += v1;
    // if (no_sv_sta_list[i] != vertex1.trkMatchSta[0]) agg_fail += 1;
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