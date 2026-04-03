#include <TGraph.h>                                                                                                                   
#include "TString.h"                                                                                                                  
#include <TFile.h>                                                                                                                    
#include <TTree.h>                                                                                                                    
#include <TH1F.h>                                                                                                                     
#include <TCanvas.h>                                                                                                                  
#include <TLegend.h>                                                                                                                  
#include <iostream>                                                                                                                   
#include <Math/Vector4D.h> 
#include <map>
#include <string>                                                                                                           
#include "Math/VectorUtil.h"                                                                                                          
#include "tTree.h"                                                                                                                    
#include "binning_histos_small.h"
#include <vector>
#include "Math/Vector4D.h"

struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    std::vector<int> trkMatchSta;
  vector<int> trkSvtxId;
  double svtxnormchi2; 
};

int StatusToMarker(int status)
{
  if (status == -1 ) return 24;
    if (status < 100) return 20;          // circle
    if (status < 200) return 22;          // triangle
    if (status < 300) return 21;          // square
    if (status < 400) return 33;          // diamond
    return 36;                            // star for higher
}

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
void DrawJet ( tTree& t, Int_t ijet, Int_t ient ){

  // Canvas per  jet
  TCanvas *c = new TCanvas(Form("c_event_%d_jet_%d", ient, ijet),
			   Form("Event %d – Jet %d", ient, ijet),
			   900, 900
			   );
  TMultiGraph *mg = new TMultiGraph();   //create graph
  std::set<int> usedSV;
  std::set<int>usedStatus; 
  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {       //loop over trks                                                                             
    // Track must belong to this jet
    if (t.trkJetId[itrk] != ijet) continue;

    // std::cout<<"jet "<<ijet<< "Trk infos:  Status : "<<t.trkMatchSta[itrk]<<"  SV id :" <<t.trkSvtxId[itrk]
    //   <<"  4 vec : pt=" << t.trkPt [itrk]<<" eta ="<< t.trkEta[itrk]<<" phi = "<< t.trkPhi[itrk]<< " pid =" <<t.trkMatchPdgId[itrk]<<std::endl;
    //initialise characteristics 
    int    status = t.trkMatchSta[itrk];
    int    svId   = t.trkSvtxId[itrk];
    double pt     = t.trkPt[itrk];
    double eta    = t.trkEta[itrk];
    double phi    = t.trkPhi[itrk];

    usedSV.insert(svId);
    // one graph per track 
    TGraph *gr = new TGraph();
   
    gr->SetPoint(0, phi, eta);
    
    // shape from status
    gr->SetMarkerStyle(StatusToMarker(status)); 
    usedStatus.insert(status); 
    
    // color from SV id 
    
    int color =(svId == -1) ? kGray : 2 + (svId % 8);  
    gr->SetMarkerColorAlpha(color, 0.5);
    gr->SetLineColor(color);
    gr->SetLineWidth(1); 
    
    // size prop to pT
    double size = 0.7 + 0.2 * pt;
    gr->SetMarkerSize(size);
    
    mg->Add(gr);
  }
  
  // Draw
  gPad->SetBottomMargin(0.25); // plus d'espace sous l'axe X
  gPad->SetLeftMargin(0.15);   // garde un peu de marge à gauche
  gPad->SetRightMargin(0.1);  // marge droite
  
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("#phi");
  mg->GetYaxis()->SetTitle("#eta");



      // --- Legend pour Status (colonne gauche) ---
    TLegend* legStatus = new TLegend(0.3, 0.03, 0.51, 0.17);
    legStatus->SetBorderSize(1);
    legStatus->SetFillStyle(1001);
    
    std::vector<std::tuple<int,int,int,std::string>> statusMap = {
      {-100, 0, 24, "Status -1"}, 
      {0, 99, 20, "Status 1"},
      {100, 199, 22, "Status 100"},
      {200, 299, 21, "Status 200"},
      {300, 399, 33, "Status 300"},
      {400, 10000, 36, "Status > 300"}
    };
    
    for (auto [low, high, marker, label] : statusMap) {
      // Check if any usedStatus is in this range
      bool used = false;
      for (int s : usedStatus) {
	if (s >= low && s <= high) {
	  used = true;
	  break;
	}
      }
      if (used) {
	TGraph* g = new TGraph();
	g->SetMarkerStyle(marker);
	g->SetMarkerSize(2.0);
	legStatus->AddEntry(g, label.c_str(), "p");
      }
    }
    legStatus->Draw();
    
    // --- Legend pour SV (colonne droite) ---
    TLegend* legSV = new TLegend(0.52, 0.03, 0.75, 0.17);
    legSV->SetBorderSize(1);
    legSV->SetFillStyle(1001);

    for (int sv : usedSV) {
        TGraph* gsv = new TGraph();
        gsv->SetMarkerStyle(20);
        gsv->SetMarkerColor((sv == -1) ? kGray : 2 + (sv % 8));
	gsv->SetMarkerSize(2.0);
        legSV->AddEntry(gsv, Form("SV id %d", sv), "p");
    }

    legSV->Draw();

  
  c->Modified();
  c->Update();
}

void DrawJetwithB ( tTree& t, Int_t ijet, Int_t ient, vector<ROOT::Math::PtEtaPhiMVector> hadrons_4vec, vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec ){

  // Canvas per  jet
  TCanvas *c3 = new TCanvas(Form("Jet_with_B_event_%d_jet_%d", ient, ijet),
			   Form("Jet_with_B_Event %d – Jet %d", ient, ijet),
			   900, 900
			   );
  TMultiGraph *mg = new TMultiGraph();   //create graph
  std::set<int> usedSV;
  std::set<int>usedStatus; 
  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {       //loop over trks                                                                             
    // Track must belong to this jet
    if (t.trkJetId[itrk] != ijet) continue;

    // std::cout<<"jet "<<ijet<< "Trk infos:  Status : "<<t.trkMatchSta[itrk]<<"  SV id :" <<t.trkSvtxId[itrk]
    //   <<"  4 vec : pt=" << t.trkPt [itrk]<<" eta ="<< t.trkEta[itrk]<<" phi = "<< t.trkPhi[itrk]<< " pid =" <<t.trkMatchPdgId[itrk]<<std::endl;
    //initialise characteristics 
    int    status = t.trkMatchSta[itrk];
    int    svId   = t.trkSvtxId[itrk];
    double pt     = t.trkPt[itrk];
    double eta    = t.trkEta[itrk];
    double phi    = t.trkPhi[itrk];

    usedSV.insert(svId);
    // one graph per track 
    TGraph *gr = new TGraph();
   
    gr->SetPoint(0, phi, eta);
    
    // shape from status
    gr->SetMarkerStyle(StatusToMarker(status)); 
    usedStatus.insert(status); 
    
    // color from SV id 
    
    int color =(svId == -1) ? kGray : 2 + (svId % 8);  
    gr->SetMarkerColorAlpha(color, 0.5);
    gr->SetLineColor(color);
    gr->SetLineWidth(1); 
    
    // size prop to pT
    double size = 0.7 + 0.2 * pt;
    gr->SetMarkerSize(size);
    
    mg->Add(gr);
  }

  //add the Bs 
  for (const auto& p4 : reco_hadrons_4vec) {   
    double pt  = p4.Pt();
    double eta = p4.Eta();
    double phi = p4.Phi();
    
    TGraph* grReco = new TGraph();
    grReco->SetPoint(0, phi, eta);
    grReco->SetMarkerStyle(5);          // croix
    grReco->SetMarkerColor(kBlack);    // vert foncé
    grReco->SetLineColor(kBlack);
    grReco->SetLineWidth(3);
    
    double size = 0.7 + 0.2 * pt;        // même scaling que tracks
    grReco->SetMarkerSize(size);
    mg->Add(grReco);
  }
  
  
//  True hadrons (croix noires)
  for (const auto& p4 : hadrons_4vec) {
    double pt  = p4.Pt();
    double eta = p4.Eta();
    double phi = p4.Phi();
    
    TGraph* grTrue = new TGraph();
    grTrue->SetPoint(0, phi, eta);
    grTrue->SetMarkerStyle(2);      // croix
    grTrue->SetMarkerColor(kBlack);  // noir
    grTrue->SetLineColor(kBlack);
    grTrue->SetLineWidth(4);
    
    double size = 0.7 + 0.2 * pt;
    grTrue->SetMarkerSize(size);
    mg->Add(grTrue);
  }

  
  // Draw
  gPad->SetBottomMargin(0.25); // plus d'espace sous l'axe X
  gPad->SetLeftMargin(0.15);   // garde un peu de marge à gauche
  gPad->SetRightMargin(0.1);  // marge droite
  
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("#phi");
  mg->GetYaxis()->SetTitle("#eta");


  
  
  //  Legend Status (colonne gauche) 
  TLegend* legStatus = new TLegend(0.15, 0.03, 0.36, 0.17);
  legStatus->SetBorderSize(1);
  legStatus->SetFillStyle(1001);
  
  std::vector<std::tuple<int,int,int,std::string>> statusMap = {
    {-100, 0, 24, "Status -1"}, 
    {0, 99, 20, "Status 1"},
    {100, 199, 22, "Status 100"},
    {200, 299, 21, "Status 200"},
    {300, 399, 33, "Status 300"},
    {400, 10000, 36, "Status > 300"}
  };
  
  for (auto [low, high, marker, label] : statusMap) {
    // Check if any usedStatus is in this range
    bool used = false;
    for (int s : usedStatus) {
      if (s >= low && s <= high) {
	used = true;
	break;
      }
    }
    if (used) {
      TGraph* g = new TGraph();
      g->SetMarkerStyle(marker);
      g->SetMarkerSize(2.0);
      legStatus->AddEntry(g, label.c_str(), "p");
    }
  }
  legStatus->Draw();
  
  // Legend SV (colonne droite)
  TLegend* legSV = new TLegend(0.37, 0.03, 0.60, 0.17);
  legSV->SetBorderSize(1);
  legSV->SetFillStyle(1001);
  
  for (int sv : usedSV) {
    TGraph* gsv = new TGraph();
    gsv->SetMarkerStyle(20);
    gsv->SetMarkerColor((sv == -1) ? kGray : 2 + (sv % 8));
    gsv->SetMarkerSize(2.0);
    legSV->AddEntry(gsv, Form("SV id %d", sv), "p");
  }
  legSV->Draw();


  
  TLegend* legHad = new TLegend(0.61, 0.03, 0.80, 0.17);
  legHad->SetBorderSize(1);
  legHad->SetFillStyle(1001);
  
  TGraph* gRecoLeg = new TGraph();
  gRecoLeg->SetMarkerStyle(5);
  gRecoLeg->SetMarkerColor(kBlack);
  gRecoLeg->SetMarkerSize(2.0);
  legHad->AddEntry(gRecoLeg, "Reco B", "p");
  
  TGraph* gTrueLeg = new TGraph();
  gTrueLeg->SetMarkerStyle(2);
  gTrueLeg->SetMarkerColor(kBlack);
  gTrueLeg->SetMarkerSize(2.0);
  legHad->AddEntry(gTrueLeg, "Partial B", "p");
  
  legHad->Draw();
  
  
  c3->Modified();
  c3->Update();
}





void DrawGenEvent ( tTree& t, Int_t ijet, Int_t ient ){

  // Canvas per  jet
  TCanvas *c2 = new TCanvas(Form("c_Gen_event_%d_jet_%d", ient, ijet),
			   Form("c_Gen_Event %d – Jet %d", ient, ijet),
			   900, 900
			   );
  TMultiGraph *mg2 = new TMultiGraph();   //create graph
  std::set<int> usedSV;
  std::set<int>usedStatus;
  if (t.nrefTrk != t.ntrk){ std::cout<<"nref "<<t.nrefTrk<<" ntrk "<<t.ntrk<<std::endl;}
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {       //loop over trks                                                                             
  
    // std::cout<<"jet "<<ijet<< "Trk infos:  Status : "<<t.trkMatchSta[itrk]<<"  SV id :" <<t.trkSvtxId[itrk]
    //   <<"  4 vec : pt=" << t.trkPt [itrk]<<" eta ="<< t.trkEta[itrk]<<" phi = "<< t.trkPhi[itrk]<< " pid =" <<t.trkMatchPdgId[itrk]<<std::endl;
    //initialise characteristics 
    int    status = t.refTrkSta[itrk];
    int    svId   = t.refTrkJetId[itrk];
    double pt     = t.refTrkPt[itrk];
    double eta    = t.refTrkEta[itrk];
    double phi    = t.refTrkPhi[itrk];

    usedSV.insert(svId);
    // one graph per track 
    TGraph *gr2 = new TGraph();
   
    gr2->SetPoint(0, phi, eta);
    
    // shape from status
    gr2->SetMarkerStyle(StatusToMarker(status)); 
    usedStatus.insert(status); 
    
    // color from SV id 
    
    int color =(svId == -1) ? kGray : 2 + (svId % 8);  
    gr2->SetMarkerColorAlpha(color, 0.5);
    gr2->SetLineColor(color);
    gr2->SetLineWidth(1); 
    
    // size prop to pT
    double size = 0.7 + 0.2 * pt;
    gr2->SetMarkerSize(size);
    
    mg2->Add(gr2);
  }
  
  // Draw
  gPad->SetBottomMargin(0.25); // plus d'espace sous l'axe X
  gPad->SetLeftMargin(0.15);   // garde un peu de marge à gauche
  gPad->SetRightMargin(0.1);  // marge droite
  
  mg2->Draw("AP");
  mg2->GetXaxis()->SetTitle("#phi");
  mg2->GetYaxis()->SetTitle("#eta");



      // --- Legend pour Status (colonne gauche) ---
    TLegend* legStatus = new TLegend(0.3, 0.03, 0.51, 0.17);
    legStatus->SetBorderSize(1);
    legStatus->SetFillStyle(1001);
    
    std::vector<std::tuple<int,int,int,std::string>> statusMap = {
      {-100, 0, 24, "Status -1"}, 
      {0, 99, 20, "Status 1"},
      {100, 199, 22, "Status 100"},
      {200, 299, 21, "Status 200"},
      {300, 399, 33, "Status 300"},
      {400, 10000, 36, "Status > 300"}
    };
    
    for (auto [low, high, marker, label] : statusMap) {
      // Check if any usedStatus is in this range
      bool used = false;
      for (int s : usedStatus) {
	if (s >= low && s <= high) {
	  used = true;
	  break;
	}
      }
      if (used) {
	TGraph* g = new TGraph();
	g->SetMarkerStyle(marker);
	g->SetMarkerSize(2.0);
	legStatus->AddEntry(g, label.c_str(), "p");
      }
    }
    legStatus->Draw();
    
    // --- Legend pour SV (colonne droite) ---
    TLegend* legSV = new TLegend(0.52, 0.03, 0.75, 0.17);
    legSV->SetBorderSize(1);
    legSV->SetFillStyle(1001);

    for (int sv : usedSV) {
        TGraph* gsv = new TGraph();
        gsv->SetMarkerStyle(20);
        gsv->SetMarkerColor((sv == -1) ? kGray : 2 + (sv % 8));
	gsv->SetMarkerSize(2.0);
        legSV->AddEntry(gsv, Form("Jet id %d", sv), "p");
    }

    legSV->Draw();

  
  c2->Modified();
  c2->Update();
}



void DrawEvent ( tTree& t, Int_t ijet, Int_t ient ){
 std::cout << "DrawEvent called" << std::endl;

  // Canvas per  jet
  TCanvas *c1 = new TCanvas(Form("full_evt_c_event_%d_jet_%d", ient, ijet),
			   Form("full_evt_Event %d – Jet %d", ient, ijet),
			   900, 900
			   );
  c1->cd();
  TMultiGraph *mg1 = new TMultiGraph();   //create graph
  std::set<int> usedSV;
  std::set<int>usedStatus; 
  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {       //loop over trks                                                                             

    // std::cout<<"jet "<<ijet<< "Trk infos:  Status : "<<t.trkMatchSta[itrk]<<"  SV id :" <<t.trkSvtxId[itrk]
    //   <<"  4 vec : pt=" << t.trkPt [itrk]<<" eta ="<< t.trkEta[itrk]<<" phi = "<< t.trkPhi[itrk]<< " pid =" <<t.trkMatchPdgId[itrk]<<std::endl;
    //initialise characteristics 
    int    status = t.trkMatchSta[itrk];
    int    svId   = t.trkSvtxId[itrk];
    double pt     = t.trkPt[itrk];
    double eta    = t.trkEta[itrk];
    double phi    = t.trkPhi[itrk];

    usedSV.insert(svId);
    // one graph per track 
    TGraph *gr1 = new TGraph();
   
    gr1->SetPoint(0, phi, eta);
    
    // shape from status
    gr1->SetMarkerStyle(StatusToMarker(status)); 
    usedStatus.insert(status); 
    
    // color from SV id 
    
    int color =(svId == -1) ? kGray : 2 + (svId % 8);  
    gr1->SetMarkerColorAlpha(color, 0.5);
    gr1->SetLineColor(color);
    gr1->SetLineWidth(1); 
    
    // size prop to pT
    double size = 0.7 + 0.2 * pt;
    gr1->SetMarkerSize(size);
    
    mg1->Add(gr1);
  }
  
  // Draw
  gPad->SetBottomMargin(0.25); // plus d'espace sous l'axe X
  gPad->SetLeftMargin(0.15);   // garde un peu de marge à gauche
  gPad->SetRightMargin(0.1);  // marge droite
  
  mg1->Draw("AP");
  mg1->GetXaxis()->SetTitle("#phi");
  mg1->GetYaxis()->SetTitle("#eta");
  mg1->SetTitle("FULL EVENT");




      // --- Legend pour Status (colonne gauche) ---
    TLegend* legStatus = new TLegend(0.3, 0.03, 0.51, 0.17);
    legStatus->SetBorderSize(1);
    legStatus->SetFillStyle(1001);
    
    std::vector<std::tuple<int,int,int,std::string>> statusMap = {
      {-100, 0, 24, "Status -1"}, 
      {0, 99, 20, "Status 1"},
      {100, 199, 22, "Status 100"},
      {200, 299, 21, "Status 200"},
      {300, 399, 33, "Status 300"},
      {400, 10000, 36, "Status > 300"}
    };
    
    for (auto [low, high, marker, label] : statusMap) {
      // Check if any usedStatus is in this range
      bool used = false;
      for (int s : usedStatus) {
	if (s >= low && s <= high) {
	  used = true;
	  break;
	}
      }
      if (used) {
	TGraph* g = new TGraph();
	g->SetMarkerStyle(marker);
	g->SetMarkerSize(2.0);
	legStatus->AddEntry(g, label.c_str(), "p");
      }
    }
    legStatus->Draw();
    
    // --- Legend pour SV (colonne droite) ---
    TLegend* legSV = new TLegend(0.52, 0.03, 0.75, 0.17);
    legSV->SetBorderSize(1);
    legSV->SetFillStyle(1001);

    for (int sv : usedSV) {
        TGraph* gsv = new TGraph();
        gsv->SetMarkerStyle(20);
        gsv->SetMarkerColor((sv == -1) ? kGray : 2 + (sv % 8));
	gsv->SetMarkerSize(2.0);
        legSV->AddEntry(gsv, Form("SV id %d", sv), "p");
    }

    legSV->Draw();

  
  c1->Modified();
  c1->Update();
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
    
    //  std::cout<<"jet "<<ijet<< "Trk infos:  Status : "<<t.trkMatchSta[itrk]<<"  SV id :" <<t.trkSvtxId[itrk]
    //    <<"  4 vec : pt=" << t.trkPt [itrk]<<" eta ="<< t.trkEta[itrk]<<" phi = "<< t.trkPhi[itrk]<< " pid =" <<t.trkMatchPdgId[itrk]<<std::endl;
                                                                                           
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
	      
                double dist = deltaR2 (vertices[i].p4, vertices[j].p4) ; 

                if (dist < min_distance) {
                    min_distance = dist;
                    index1 = i;
                    index2 = j;
                }
            }
        }

	// merge vertex index2 into index1 
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


vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs(
					      tTree& t,     //event
					      Int_t& ijet,  //jet nb
					      Long64_t& ient, 
					      double& agg_fail,     //dominant trkMatchSta fraction for vertex 0
					      double& nb_sv ,    //   ..      ..           ..     ..    ..   1
					      TH2F* global_histo,   // histogram with nb of seconcary vertex per jet 
					      double& sv_fail,
					      double& merge_fail,
					      TH1D* chi_single_sta,
					      TH1D* chi_mixed_sta ){
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
	if      (pid == 211)  v1.SetM(0.139570);                                                                                            
	else if (pid == 13)   v1.SetM(0.105658);                                                                                            
	else if (pid == 11)   v1.SetM(0.000510);                                                                                            
	else if (pid == 2212) v1.SetM(0.938272);                                                                                            
	else if (pid == 321)  v1.SetM(0.493677);                                                                                            
	else if (pid == 3112) v1.SetM(1.19744);                                                                                             
	else if (pid == 3222) v1.SetM(1.18937);                                                                                             
	else if (pid == 3312) v1.SetM(1.32171);                                                                                             
	else if (pid == 3334) v1.SetM( 1.67245);
	else if (pid == 1 ) continue ;  //bug
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
    //if (secVtxs.size() < 2 and  t.nsvtx >= 2 ) {
    // nb_sv += 1 ; 
    //} 
    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      if (nb_sv < 15 ) {
	//DrawGenEvent(t, ijet, ient);
	//DrawEvent(t, ijet, ient) ;
	//DrawJet(t, ijet, ient);
      } 
      return empty;
    }

    // Build Vertex objects
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      auto& v = secVtxsMatchSta[originalSvId];
      if (std::adjacent_find(v.begin(), v.end(),std::not_equal_to<int>()) != v.end()) {      // if not all tracks in the sv come from the same B 
	sv_fail += 1.0;
	chi_mixed_sta->Fill(t.svtxnormchi2[originalSvId], weight_tree); 
	//return empty;
      }
      else {chi_single_sta->Fill(t.svtxnormchi2[originalSvId], weight_tree);}  //all tracks come from the same B 
      Vertex vtx;
      vtx.svtxnormchi2 = t.svtxnormchi2[originalSvId];
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
      //DrawJet (t, ijet, ient); 
      // printvtx( vertex0  , " wrong merging ", 0); 
      //return empty ;
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


    //printvtx(vertex0, "sv_fail ", 0, ijet, ient );
    //printvtx(vertex1, "sv_fail ", 1, ijet, ient );
    
// counts the number of original vertex merged into one final vertex / one final B 
    std::set<int> originalVertices0(  vertices[0].trkSvtxId.begin(), vertices[0].trkSvtxId.end()    );   // keeps only one copy of each sv nb 
    int nOriginalVertices0 = originalVertices0.size();
    // vertex 1
    std::set<int> originalVertices1(  vertices[1].trkSvtxId.begin(), vertices[1].trkSvtxId.end() );
    int nOriginalVertices1 = originalVertices1.size();

    // if (nOriginalVertices1 +  nOriginalVertices0 == 3 ) {
    //   group3vertices(t, ijet, vertices, global_count) ;        
    //    } 
    
    //  if (nOriginalVertices0 +  nOriginalVertices1 == 4  )  {
    //  global_count +=1 ; 
    // group4vertices( t,  ijet, vertices) ;
    // vertex0 = vertices[0];
    // vertex1 = vertices[1];
    // }  
    
      
      // Return summed 4-vectors of the 2 vertices
    vector<ROOT::Math::PtEtaPhiMVector> vecFinalSecVtxs;
    vecFinalSecVtxs.push_back(vertices[0].p4);
    vecFinalSecVtxs.push_back(vertices[1].p4);
    /* SVntrks1 = prop ;
    SVntrks2 = prop1;
    double weight_tree = t.weight;
    global_histo ->Fill(  t.jtpt[ijet], (prop+prop1)/2, weight_tree); 
    if (prop != 1 or prop1 != 1){
      if (prop != 1) {global_histo ->Fill(  vertices[0].p4.Pt(), prop, weight_tree); }
      if (prop1 != 1) {global_histo ->Fill(  vertices[1].p4.Pt(), prop1, weight_tree); }
    } 
    if (prop == 1 and prop1 == 1) {
	global_histo ->Fill(  vertices[0].p4.Pt(), prop, weight_tree);
	global_histo ->Fill(  vertices[1].p4.Pt(), prop1, weight_tree);
	} */

    return vecFinalSecVtxs;
}


void reco_tracks_eec_sec_ver(TString filename_bjet, TString output_folder, TString output_hist, TString domain, Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC) {
    
    tTree t;
    t.Init(filename_bjet, isMC);
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
    int counter = 0 ;
    int counter1 = 0 ;
    int counter2 = 0 ; 

    //activates branches 
    std::vector<TString> active_branches = {"weight","jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "jtNtrk",
      "ntrk", "trkJetId", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "refTrkPdgId","refTrkSta", "refTrkMass",  "refmB", "refpt", "refeta", "refphi",                                                                                           
      "nrefTrk", "refTrkJetId", "refTrkPdgId", "refTrkPt", "refTrkEta", "refTrkPhi", "refTrkY", "refTrkSta", "refTrkMass", "refNtrk",
      "jtNsvtx", "trkSvtxId",  "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",  "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",};
    t.SetBranchStatus(active_branches, 1);

    // Plots or histograms

    TH1D *partial_h_tracks_2b = new TH1D("partial_h_tracks_2b", "#DeltaR;EEC", bins_dr, dr_binsVector);
    partial_h_tracks_2b->Sumw2();
    TH1D *partial_M = new TH1D( "partial_M",";Reconstructed mass; ",   200,    0.0,  50);
    partial_M->Sumw2();    
    TH1D *partial_pT = new TH1D( "partial_pT",";Reconstructed p_{T}; ",   100,    0.0, 100 );
    partial_pT->Sumw2();    
    TH1D *partial_h_dR = new TH1D( "partial_h_dR",";#DeltaR between B hadrons; ",   100,    0.0,  1 );
    partial_h_dR->Sumw2();    
    TH1D *partial_d_pT = new TH1D( "partial_d_pT",";#Delta p_{T} between B hadrons; ",   200,    0.0,  200 );
    partial_d_pT->Sumw2();    
    TH1D *partial_np = new TH1D( "partial_np",";Number of particles in the jet's SVs; ",   50,    0.0, 50 );
    partial_np->Sumw2();  
    
    //all_jets 
    TH1D *alljets_h_tracks_2b = new TH1D("alljets_h_tracks_2b", "#DeltaR;EEC", bins_dr, dr_binsVector);
    alljets_h_tracks_2b->Sumw2();
    TH1D *alljets_M = new TH1D( "alljets_M",";Reconstructed mass; ",   200,    0.0,  50);
    alljets_M->Sumw2();    
    TH1D *alljets_pT = new TH1D( "alljets_pT",";Reconstructed p_{T}; ",   100,    0.0, 100 );
    alljets_pT->Sumw2();    
    TH1D *alljets_h_dR = new TH1D( "alljets_h_dR",";#DeltaR between B hadrons; ",   100,    0.0,  1 );
    alljets_h_dR->Sumw2();    
    TH1D *alljets_d_pT = new TH1D( "alljets_d_pT",";#Delta p_{T} between B hadrons; ",   200,    0.0,  200 );
    alljets_d_pT->Sumw2();    
    TH1D *alljets_np = new TH1D( "alljets_np",";Number of particles in the jet's SVs; ",   50,    0.0, 50 );
    alljets_np->Sumw2();    
    
    //matched2goodSV = mat2gdSV
    TH1D *mat2gdSV_h_tracks_2b = new TH1D("mat2gdSV_h_tracks_2b", "#DeltaR;EEC", bins_dr, dr_binsVector);
    mat2gdSV_h_tracks_2b->Sumw2();    
    TH1D *mat2gdSV_M = new TH1D( "mat2gdSV_M",";Reconstructed mass; ",   200,    0.0,  50 );
    mat2gdSV_M->Sumw2();    
    TH1D *mat2gdSV_pT = new TH1D( "mat2gdSV_pT",";Reconstructed p_{T}; ",   100,    0.0, 100);
    mat2gdSV_pT->Sumw2();    
    TH1D *mat2gdSV_h_dR = new TH1D( "mat2gdSV_h_dR",";#DeltaR between B hadrons; ",   100,    0.0,  1 );
    mat2gdSV_h_dR->Sumw2();    
    TH1D *mat2gdSV_d_pT = new TH1D( "mat2gdSV_d_pT",";#Delta p_{T} between B hadrons; ",   200,    0.0,  200 );
    mat2gdSV_d_pT->Sumw2();    
    TH1D *mat2gdSV_np = new TH1D( "mat2gdSV_np",";Number of particles in the jet's SVs; ",   50,    0.0, 50 );
    mat2gdSV_np->Sumw2();
    
        //reco 2 good SV = rec2gdSV
    TH1D *rec2gdSV_h_tracks_2b = new TH1D("rec2gdSV_h_tracks_2b", "#DeltaR;EEC", bins_dr, dr_binsVector);
    rec2gdSV_h_tracks_2b->Sumw2();    
    TH1D *rec2gdSV_M = new TH1D( "rec2gdSV_M",";Reconstructed mass; ",   200,    0.0,  50 );
    rec2gdSV_M->Sumw2();    
    TH1D *rec2gdSV_pT = new TH1D( "rec2gdSV_pT",";Reconstructed p_{T}; ",   100,    0.0, 100);
    rec2gdSV_pT->Sumw2();    
    TH1D *rec2gdSV_h_dR = new TH1D( "rec2gdSV_h_dR",";#DeltaR between B hadrons; ",   100,    0.0,  1 );
    rec2gdSV_h_dR->Sumw2();    
    TH1D *rec2gdSV_d_pT = new TH1D( "rec2gdSV_d_pT",";#Delta p_{T} between B hadrons; ",   200,    0.0,  200 );
    rec2gdSV_d_pT->Sumw2();    
    TH1D *rec2gdSV_np = new TH1D( "rec2gdSV_np",";Number of particles in the jet's SVs; ",   50,    0.0, 50 );
    rec2gdSV_np->Sumw2();    
    
    //wrong jets = wrgjets 
    TH1D *wrgjets_h_tracks_2b = new TH1D("wrgjets_h_tracks_2b", "#DeltaR;EEC", bins_dr, dr_binsVector);
    wrgjets_h_tracks_2b->Sumw2();    
    TH1D *wrgjets_M = new TH1D( "wrgjets_M",";Reconstructed mass; ",  200,    0.0, 50 );
    wrgjets_M->Sumw2();    
    TH1D *wrgjets_pT = new TH1D( "wrgjets_pT",";Reconstructed p_{T}; ",   100,    0.0, 100 );
    wrgjets_pT->Sumw2();    
    TH1D *wrgjets_h_dR = new TH1D( "wrgjets_h_dR",";#DeltaR between B hadrons; ",   100,    0.0,  1 );
    wrgjets_h_dR->Sumw2();    
    TH1D *wrgjets_d_pT = new TH1D( "wrgjets_d_pT",";#Delta p_{T} between B hadrons; ",   200,    0.0,  200 );
    wrgjets_d_pT->Sumw2();    
    TH1D *wrgjets_np = new TH1D( "wrgjets_np",";Number of particles in the jet's SVs; ",   50,    0.0, 50 );
    wrgjets_np->Sumw2();
 


    //old histos 
    TH1D *bpt_frac = new TH1D(  "bpt_frac", ";(p_{T1}+p_{T2})/jet p_{T}; ",    100,     0,     1.1  );
    bpt_frac->Sumw2(); 
    TH2F* global_histo = new TH2F("global_histo", "distribution of wrong reco with pt; pT; reco proportion", 40, 0, 180, 40, 0, 1.2);
    global_histo->Sumw2(); 
    TH2F *pt_comp = new TH2F(  "pt_comp",  "pt comparison; #Delta p_{T}_{reco} ;#Delat p_{T}_{ref}", 
    500, 0, 100,                // X bins, min, max
    500, 0, 100                // Y bins, min, max
);
    TH2F *asym_vs_fracpt = new TH2F(  "asym_vs_fracpt",  "asym vs frac pt ;pt1+pt2/jtpt ; asym", 
    50, 0, 1,                // X bins, min, max
    20, 0, 1                 // Y bins, min, max
);
    TH1D *d_matched_frac = new TH1D( "d_matched_frac",  "; ; ",   100,  0.0,  5.0  );
    d_matched_frac->Sumw2();
    TH1D *chi_single_sta = new TH1D( "chi_single_sta",  "; ; ",   1000,  0.0,  100  );
    chi_single_sta->Sumw2();
    TH1D *chi_mixed_sta = new TH1D( "chi_mixed_sta",  "; ; ",   1000,  0.0,  100  );
    chi_mixed_sta->Sumw2(); 

    int i = 0 ; 
      

    Long64_t n_bjet = t.GetEntries();
    for (Long64_t ient = 0; ient < n_bjet; ient++) {        //event loop (ient)
      //Progress    
    if (ient % 50000 == 0) {                                                                   
      float percent = 100.0 * ient / n_bjet;
      std::cout << "\rProcessing reco_sec_ver_eec: " << percent << " %" << std::flush;                                                                                                                                         
    }
    //if (ient % 10 != 0) continue;
    t.GetEntry(ient);
    //get the MC event weight                                                                                                                                  
    double weight_tree = t.weight;

    for (Int_t ijet = 0; ijet < t.nref; ijet++) {            //jet loop (ijet)
      //some cuts
      if (std::abs(t.jteta[ijet]) > 1.9) continue;
      if (std::abs(t.jtpt[ijet]) < pT_low || std::abs(t.jtpt[ijet]) > pT_high) continue;
      if ((isMC) && skipMC(t.jtpt[ijet], t.pthat)) continue;
      
      //Select jets passing the b-jet tagging (if needed)                                                                                                        
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;
      if (t.jtNbHad[ijet] !=  2 ) {
	counter += 1 ;
	//DrawEvent (t, ijet, ient) ; 
	continue;  
      }

    

      //checks that the jets contains two Bs (in our definition, contsining 2 Bs = we reconstruct 2 from gen lvl 
      std::vector<Int_t> hadrons_stat;
      std::vector<ROOT::Math::PtEtaPhiMVector> hadrons_4vec;
    
      PartialBsAggregation(hadrons_4vec, hadrons_stat, t, ijet);
      if (hadrons_4vec.size() < 2) continue;  //skip if single B 
    
      
      //calculate the nb of B had (if >= 2) 
      size_t nB_aggr = 0;
      for (size_t k = 0; k < hadrons_stat.size(); ++k) {
        if (hadrons_stat[k] >= 100) ++nB_aggr;
      }

      if (nB_aggr != 2  ) {
	continue;                    //cuts jets that don t have 2 b
      }
      
      
      counter2 += 1 ; 

      double partial_dr = t.calc_dr(hadrons_4vec[0].Eta(),
				    hadrons_4vec[0].Phi(),
				    hadrons_4vec[1].Eta(),
				    hadrons_4vec[1].Phi()
				    );
      double partial_pt1 = hadrons_4vec[0].Pt();
      double partial_pt2 = hadrons_4vec[1].Pt();
      double partial_eec = std::pow(partial_pt1 * partial_pt2, n);
      double partial_sumpt = partial_pt1 + partial_pt2;
      double partial_diffpt = std::abs(partial_pt1 - partial_pt2);
      double partialM0 = hadrons_4vec[0].M();
      double partialM1 = hadrons_4vec[1].M();


      
      partial_h_tracks_2b->Fill(partial_dr, partial_eec * weight_tree);
      partial_M->Fill(partialM0, weight_tree);
      partial_M->Fill(partialM1, weight_tree);
      partial_pT->Fill(partial_pt1, weight_tree);
      partial_pT->Fill(partial_pt2, weight_tree);
      partial_h_dR->Fill(partial_dr, weight_tree);
      partial_d_pT->Fill(partial_diffpt, weight_tree);
      
      
      // -- now we start the reco step --
      Float_t SVntrks2 = 0 ;
      beforeskip += 1 ;
      // if (beforeskip < 1 ) DrawJet(t, ijet, ient) ; 
      double merge_fail1;
      merge_fail1 =  merge_fail ; 
                //tracks aggregation : takes (empty vector, t tree, jet nb) and fills the vector reco_hadrons_4vec with the2 reconstructed B
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs( t,  ijet, ient, agg_fail , nb_sv, global_histo, sv_fail, merge_fail, chi_single_sta, chi_mixed_sta);
    
      if (reco_hadrons_4vec.empty()) {
	lt2sv += 1 ;
      }
      if (reco_hadrons_4vec.size() !=  2) continue;
      totjt += 1 ; 

      
      if (   matchskip < 0 ) { 
	//DrawJet (t, ijet, ient);
	DrawJetwithB (t, ijet, ient, hadrons_4vec, reco_hadrons_4vec) ;
	DrawEvent (t, ijet, ient);
	DrawGenEvent (t, ijet, ient) ;
	matchskip += 1 ; 
      }
     
      double alljets_dr = t.calc_dr(reco_hadrons_4vec[0].Eta(),
                reco_hadrons_4vec[0].Phi(),
                reco_hadrons_4vec[1].Eta(),
                reco_hadrons_4vec[1].Phi()
                );
      double alljets_pt1 = reco_hadrons_4vec[0].Pt();
      double alljets_pt2 = reco_hadrons_4vec[1].Pt();
      double alljets_eec = std::pow(alljets_pt1 * alljets_pt2, n);
      double alljets_sumpt = alljets_pt1 +alljets_pt2;
      double jtpt = t.jtpt[ijet];
      double alljets_diffpt = std::abs(alljets_pt1-alljets_pt2);
      double alljetsM0 = reco_hadrons_4vec[0].M();
      double alljetsM1 = reco_hadrons_4vec[1].M();

      
      alljets_h_tracks_2b -> Fill (alljets_dr, alljets_eec * weight_tree);
      alljets_M -> Fill (alljetsM0, weight_tree) ;
      alljets_M -> Fill (alljetsM1, weight_tree) ;
      alljets_pT -> Fill (alljets_pt1 , weight_tree) ;
      alljets_pT -> Fill (alljets_pt2, weight_tree) ;
      alljets_h_dR -> Fill (alljets_dr, weight_tree) ;
      alljets_d_pT -> Fill (alljets_diffpt , weight_tree) ;
      
      //if (dr < 0.1  ) continue ;

      //histogram filling 
      
      //bpt_frac->Fill( sumpt/jtpt  , weight_tree);
      //h_dR->Fill(dr, weight_tree);
      //d_matched_frac->Fill(SVntrks1, weight_tree);
      //d_matched_frac->Fill(SVntrks2, weight_tree); 
      //h_tracks_2b ->Fill(dr, eec*weight_tree); 
 
      
      //compare with results from matched tracks aggregation 
      vector<ROOT::Math::PtEtaPhiMVector> mat_reco_hadrons_4vec ;
      std::vector<Int_t> mat_reco_hadrons_stat;
      MatchingTracksAggregation(mat_reco_hadrons_4vec, mat_reco_hadrons_stat, t, ijet, ient);
      // std::cout << "size = " << mat_reco_hadrons_4vec.size() << std::endl;

      //print
      for (size_t i = 0; i < mat_reco_hadrons_4vec.size(); ++i) {
	const auto& v = mat_reco_hadrons_4vec[i];
	//	std::cout << "  [" << i << "] " << "pt="  << v.Pt()  << ", eta=" << v.Eta() << ", phi=" << v.Phi() << ", m="   << v.M()  << std::endl;
      }

      
      if (mat_reco_hadrons_4vec.size() !=  2  ){
	//std::cout<<"drawing from main..."<<mat_reco_hadrons_4vec.size()<<std::endl;
        if (t.jtNbHad[ijet] == 1) {
	  vertskip += 1 ; 
	  //std::cout<<"for event " << ient << "   jet "<<ijet<<"  JtNbHad = " << t.jtNbHad[ijet] <<std::endl;
	}
	//	DrawJet(t, ijet, ient ) ;
	/*	if (i < 10 ) {
	  std::cout<<"entering loop with i " <<i <<" and nB aggr : "<<nB_aggr <<std::endl;
	  for (size_t i = 0; i < hadrons_stat.size(); ++i) {
	    std::cout << "hadrons_stat[" << i << "] = "
		      << hadrons_stat[i] << std::endl;
	  }
	  DrawGenEvent (t, ijet, ient); 
	  DrawEvent (t, ijet, ient) ;
	  DrawJet(t, ijet, ient) ; 
	  i+= 1 ; 
	  } */
	wrgjets_h_tracks_2b -> Fill (alljets_dr, alljets_eec * weight_tree);
	wrgjets_M -> Fill (alljetsM0, weight_tree) ;
	wrgjets_M -> Fill (alljetsM1, weight_tree) ;
	wrgjets_pT -> Fill (alljets_pt1 , weight_tree) ;
	wrgjets_pT -> Fill (alljets_pt2, weight_tree) ;
	wrgjets_h_dR -> Fill (alljets_dr, weight_tree) ;
	wrgjets_d_pT -> Fill (alljets_diffpt , weight_tree) ;
	continue ;
      } //close mat reco had 4 vec size != 2 
    

      

      double mat2gdSV_dr = t.calc_dr(mat_reco_hadrons_4vec[0].Eta(),
                mat_reco_hadrons_4vec[0].Phi(),
                mat_reco_hadrons_4vec[1].Eta(),
                mat_reco_hadrons_4vec[1].Phi()
                );
      double mat2gdSV_pt1 = mat_reco_hadrons_4vec[0].Pt();
      double mat2gdSV_pt2 = mat_reco_hadrons_4vec[1].Pt();
      double mat2gdSV_eec = std::pow(mat2gdSV_pt1 * mat2gdSV_pt2, n);
      double mat2gdSV_sumpt = mat2gdSV_pt1 +mat2gdSV_pt2;
      double mat2gdSV_diffpt = std::abs(mat2gdSV_pt1-mat2gdSV_pt2);
      double mat2gdSV_M0 = mat_reco_hadrons_4vec[0].M();
      double mat2gdSV_M1 = mat_reco_hadrons_4vec[1].M();

      mat2gdSV_h_tracks_2b -> Fill (mat2gdSV_dr, mat2gdSV_eec * weight_tree);
      mat2gdSV_M -> Fill (mat2gdSV_M0, weight_tree) ;
      mat2gdSV_M -> Fill (mat2gdSV_M1, weight_tree) ;
      mat2gdSV_pT -> Fill (mat2gdSV_pt1 , weight_tree) ;
      mat2gdSV_pT -> Fill (mat2gdSV_pt2, weight_tree) ;
      mat2gdSV_h_dR -> Fill (mat2gdSV_dr, weight_tree) ;
      mat2gdSV_d_pT -> Fill (mat2gdSV_diffpt , weight_tree) ;


      rec2gdSV_h_tracks_2b -> Fill (alljets_dr, alljets_eec * weight_tree);
      rec2gdSV_M -> Fill (alljetsM0, weight_tree) ;
      rec2gdSV_M -> Fill (alljetsM1, weight_tree) ;
      rec2gdSV_pT -> Fill (alljets_pt1 , weight_tree) ;
      rec2gdSV_pT -> Fill (alljets_pt2, weight_tree) ;
      rec2gdSV_h_dR -> Fill (alljets_dr, weight_tree) ;
      rec2gdSV_d_pT -> Fill (alljets_diffpt , weight_tree) ;
      

      
      //if (std::abs(diffpt-_1diffpt) == 0  ) {
	// }
      // if (std::abs(diffpt - _1diffpt) != 0 ) {
	//matchskip += 1;
	//pt_comp->Fill(_1diffpt , diffpt, weight_tree);
      //  }
      
      //pt_comp->Fill(_1diffpt , diffpt, weight_tree);     
      //asym_vs_fracpt->Fill(_1dr, alpha - _1alpha , weight_tree); 
      //std::cout << "filling with eec : " << eec<<"  and   dr : "<< dr << std::endl;
    }
    }//close event loop

    std::cout<<"nb of sv fail : "<< sv_fail <<"  nb of merge fail  "<< merge_fail <<"  agg fail : " << agg_fail<<"  nb_sv  "<<nb_sv <<"  beforeskip "<< beforeskip <<" lt2sv "<<lt2sv<< "   verskip =  " << vertskip<<"  totjt " << totjt<< " matchskip = "<<matchskip<< " counter "<<counter<<  std::endl;
     std::cout<<"counter1 "<<counter1 << "  counter2 : " <<counter2 <<std::endl; 
    

    //save histos 
    TString label = "_btag";
    if(!btag) label = "_nobtag";
    TFile outFile( (output_folder + output_hist + label + domain).Data(), "RECREATE");
    
    //all_jets
    alljets_h_tracks_2b->Write();
    alljets_M->Write();
    alljets_pT->Write();
    alljets_h_dR->Write();
    alljets_d_pT->Write();
    //alljets_np->Write();
    
    //matched2goodSV = mat2gdSV
    mat2gdSV_h_tracks_2b->Write();
    mat2gdSV_M->Write();
    mat2gdSV_pT->Write();
    mat2gdSV_h_dR->Write();
    mat2gdSV_d_pT->Write();
    //mat2gdSV_np->Write();
    
    //reco 2 good SV = rec2gdSV
    rec2gdSV_h_tracks_2b->Write();
    rec2gdSV_M->Write();
    rec2gdSV_pT->Write();
    rec2gdSV_h_dR->Write();
    rec2gdSV_d_pT->Write();
    //rec2gdSV_np->Write();
    
    //wrong jets = wrgjets
    wrgjets_h_tracks_2b->Write();
    wrgjets_M->Write();
    wrgjets_pT->Write();
    wrgjets_h_dR->Write();
    wrgjets_d_pT->Write();
    //wrgjets_np->Write();

    //partila
    partial_h_tracks_2b->Write();
    partial_M->Write();
    partial_pT->Write();
    partial_h_dR->Write();
    partial_d_pT->Write();
   
    /*
      bpt_frac->Write();
      pt_comp->Write(); 
      global_histo ->Write();     
      asym_vs_fracpt->Write();
      chi_mixed_sta -> Write();
      chi_single_sta -> Write(); 
    */
    
    outFile.Close();                                                                                                                                                                                                                                                                                                             
}   























void create_histos_EEC(Float_t pT_low = 80, Float_t pT_high = 140, Int_t n = 1, bool isMC = true, bool btag = true) {

  TString output_folder = gSystem->ExpandPathName("$mydata/");
  TString output_hist_had = "EEC_at_had_level"; 
  //TString output_hist_gen = "EEC_at_gen_level";
  TString output_hist_reco_tracks = "EEC_at_reco_tracks_level_sv_only";
  TString output_hist_reco_tracks_2svcut = "EEC_at_reco_tracks_level_2svcut";
  TString output_hist_reco_tracks_closest = "EEC_at_reco_tracks_closest";
  TString output_hist_reco_sec_ver = "reco_vs_matched_allBtrks_jtnbhaddef";

  TString filename_bjet = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
  TString domain = ".root";
  
  reco_tracks_eec_sec_ver(filename_bjet, output_folder, output_hist_reco_sec_ver, domain, pT_low, pT_high, n, btag, isMC);
  
}
